#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(testthat))

file_argument <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][[1L]]
test_file <- sub("^--file=", "", file_argument)
root <- normalizePath(file.path(dirname(test_file), "..", "..", ".."), mustWork = TRUE)
ci_dir <- file.path(root, ".github", "ci")
rscript <- file.path(R.home("bin"), "Rscript")

run_script <- function(script, args, env = character()) {
  stdout <- tempfile("stdout-")
  stderr <- tempfile("stderr-")
  on.exit(unlink(c(stdout, stderr)), add = TRUE)
  status <- system2(
    rscript,
    c("--vanilla", shQuote(file.path(ci_dir, script)), vapply(args, shQuote, "")),
    stdout = stdout,
    stderr = stderr,
    env = env
  )
  list(
    status = if (is.null(status)) 0L else as.integer(status),
    stdout = paste(readLines(stdout, warn = FALSE), collapse = "\n"),
    stderr = paste(readLines(stderr, warn = FALSE), collapse = "\n")
  )
}

make_source <- function(directory) {
  package <- file.path(directory, "colocboost")
  dir.create(package, recursive = TRUE)
  writeLines(c(
    "Package: colocboost",
    "Version: 1.0.9",
    "Depends: R (>= 4.0.0)",
    "Imports: Rfast, matrixStats",
    "Suggests: testthat (>= 3.0.0), knitr, rmarkdown, ashr, MASS, susieR",
    "VignetteBuilder: knitr"
  ), file.path(package, "DESCRIPTION"))
  tarball <- file.path(directory, "colocboost_1.0.9.tar.gz")
  old <- setwd(directory)
  on.exit(setwd(old), add = TRUE)
  utils::tar(tarball, files = "colocboost", compression = "gzip", tar = "internal")
  tarball
}

test_that("dependency plans distinguish full MKL from exact noSuggests policy", {
  directory <- tempfile("dependency-plan-")
  dir.create(directory)
  tarball <- make_source(directory)

  mkl_evidence <- file.path(directory, "mkl.txt")
  mkl <- run_script("prepare-rhub-dependencies.R", c(
    "--environment-id=mkl",
    paste0("--tarball=", tarball),
    paste0("--library=", file.path(directory, "mkl-lib")),
    paste0("--evidence=", mkl_evidence),
    "--install-package=false",
    "--plan-only=true"
  ), env = "CI_DEPENDENCY_PREP_TESTING=1")
  expect_identical(mkl$status, 0L, info = mkl$stderr)
  mkl_doc <- jsonlite::read_json(mkl_evidence, simplifyVector = TRUE)
  expect_identical(mkl_doc$dependency_policy, "all")
  expect_identical(mkl_doc$ci_tooling, c("devtools", "jsonlite", "yaml"))
  expect_false(mkl_doc$install_package)
  expect_true(mkl_doc$plan_only)

  no_evidence <- file.path(directory, "nosuggests.txt")
  no <- run_script("prepare-rhub-dependencies.R", c(
    "--environment-id=nosuggests",
    paste0("--tarball=", tarball),
    paste0("--library=", file.path(directory, "no-lib")),
    paste0("--evidence=", no_evidence),
    "--install-package=true",
    "--plan-only=true"
  ), env = "CI_DEPENDENCY_PREP_TESTING=1")
  expect_identical(no$status, 0L, info = no$stderr)
  no_doc <- jsonlite::read_json(no_evidence, simplifyVector = TRUE)
  expect_identical(no_doc$dependency_policy, "hard-plus-testing-and-vignette-builder")
  expect_identical(no_doc$hard_dependencies, c("R", "Rfast", "matrixStats"))
  expect_identical(no_doc$selected_soft_dependencies, c("testthat", "knitr"))
  expect_identical(no_doc$excluded_suggests, c("rmarkdown", "ashr", "MASS", "susieR"))
  expect_true(no_doc$install_package)
})

test_that("dependency prep rejects wrong install mode and reused libraries", {
  directory <- tempfile("dependency-reject-")
  dir.create(directory)
  tarball <- make_source(directory)
  library <- file.path(directory, "library")
  dir.create(library)
  writeLines("owned", file.path(library, "stale"))
  base <- c(
    paste0("--tarball=", tarball),
    paste0("--library=", library),
    paste0("--evidence=", file.path(directory, "evidence.txt")),
    "--plan-only=true"
  )
  reused <- run_script("prepare-rhub-dependencies.R", c(
    "--environment-id=mkl", base, "--install-package=false"
  ), env = "CI_DEPENDENCY_PREP_TESTING=1")
  expect_failure(expect_identical(reused$status, 0L))

  fresh_base <- c(
    paste0("--tarball=", tarball),
    paste0("--library=", file.path(directory, "fresh")),
    paste0("--evidence=", file.path(directory, "other-evidence.txt")),
    "--plan-only=true"
  )
  wrong <- run_script("prepare-rhub-dependencies.R", c(
    "--environment-id=nosuggests",
    fresh_base,
    "--install-package=false"
  ), env = "CI_DEPENDENCY_PREP_TESTING=1")
  expect_failure(expect_identical(wrong$status, 0L))
})

make_maps <- function(directory, extra = character(), omit = character()) {
  required <- c(
    "/opt/intel/lib/libmkl_intel_lp64.so",
    "/opt/intel/lib/libmkl_core.so.2",
    "/opt/intel/lib/libmkl_sequential.so"
  )
  required <- setdiff(required, omit)
  path <- file.path(directory, "maps")
  writeLines(sprintf("0000-1111 r-xp 0000 00:00 0 %s", c(required, extra)), path)
  path
}

run_verify <- function(directory, maps, threads = c("1", "1"), verbose = "MKL_VERBOSE oneMKL 2025 GEMM") {
  verbose_file <- file.path(directory, "verbose.log")
  writeLines(verbose, verbose_file)
  evidence <- file.path(directory, paste0("evidence-", length(list.files(directory)), ".json"))
  result <- run_script("verify-mkl.R", c(
    paste0("--policy=", file.path(ci_dir, "check-policy.yml")),
    paste0("--evidence=", evidence),
    paste0("--maps-file=", maps),
    paste0("--verbose-file=", verbose_file)
  ), env = c(
    "CI_VERIFY_MKL_TESTING=1",
    paste0("MKL_NUM_THREADS=", threads[[1]]),
    paste0("OMP_NUM_THREADS=", threads[[2]]),
    "MKL_VERBOSE=1"
  ))
  list(result = result, evidence = evidence)
}

test_that("MKL verifier accepts exact serial required mappings in fixture mode", {
  directory <- tempfile("mkl-pass-")
  dir.create(directory)
  verified <- run_verify(directory, make_maps(directory))
  expect_identical(verified$result$status, 0L, info = verified$result$stderr)
  evidence <- jsonlite::read_json(verified$evidence, simplifyVector = TRUE)
  expect_true(evidence$status == "pass")
  expect_identical(evidence$threads$MKL_NUM_THREADS, "1")
  expect_identical(evidence$threads$OMP_NUM_THREADS, "1")
  expect_match(evidence$mkl_verbose, "MKL_VERBOSE")
  expect_true(isTRUE(evidence$blas_operation$completed))
})

test_that("MKL verifier rejects missing, threaded, and alternative BLAS mappings", {
  forbidden <- c(
    "/opt/intel/lib/libmkl_intel_thread.so",
    "/usr/lib/libopenblas.so.0",
    "/usr/lib/libatlas.so",
    "/usr/lib/libblis.so"
  )
  for (library in forbidden) {
    directory <- tempfile("mkl-forbidden-")
    dir.create(directory)
    verified <- run_verify(directory, make_maps(directory, extra = library))
    expect_failure(expect_identical(verified$result$status, 0L), info = library)
  }
  directory <- tempfile("mkl-missing-")
  dir.create(directory)
  missing <- run_verify(
    directory,
    make_maps(directory, omit = "/opt/intel/lib/libmkl_core.so.2")
  )
  expect_failure(expect_identical(missing$result$status, 0L))

  directory <- tempfile("mkl-thread-")
  dir.create(directory)
  threaded <- run_verify(directory, make_maps(directory), threads = c("2", "1"))
  expect_failure(expect_identical(threaded$result$status, 0L))
})

test_that("the local non-MKL runtime fails the production verifier", {
  directory <- tempfile("not-mkl-")
  dir.create(directory)
  result <- run_script("verify-mkl.R", c(
    paste0("--policy=", file.path(ci_dir, "check-policy.yml")),
    paste0("--evidence=", file.path(directory, "evidence.json"))
  ), env = c("MKL_NUM_THREADS=1", "OMP_NUM_THREADS=1", "MKL_VERBOSE=1"))
  expect_failure(expect_identical(result$status, 0L))
})
