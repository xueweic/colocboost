#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
required <- c("environment-id", "tarball", "library", "evidence", "install-package")
optional <- "plan-only"

parse_cli <- function(values) {
  if (any(!grepl("^--[a-z][a-z0-9-]*=.+$", values))) {
    stop("Every argument must use --name=value syntax.")
  }
  keys <- sub("^--([^=]+)=.*$", "\\1", values)
  vals <- sub("^--[^=]+=", "", values)
  if (anyDuplicated(keys)) stop("Duplicate argument: --", keys[duplicated(keys)][[1L]])
  unknown <- setdiff(keys, c(required, optional))
  if (length(unknown)) stop("Unknown argument: --", unknown[[1L]])
  missing <- setdiff(required, keys)
  if (length(missing)) stop("Missing argument: --", missing[[1L]])
  result <- as.list(vals)
  names(result) <- keys
  result
}

scalar_bool <- function(value, label) {
  if (!value %in% c("true", "false")) stop(label, " must be true or false.")
  identical(value, "true")
}

absolute <- function(value, label, must_work = FALSE) {
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", value)) stop(label, " must be absolute.")
  normalizePath(value, winslash = "/", mustWork = must_work)
}

split_dependencies <- function(value) {
  if (is.na(value) || !nzchar(trimws(value))) return(character())
  pieces <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  pieces <- trimws(gsub("\\s*\\([^)]*\\)\\s*$", "", pieces))
  pieces[nzchar(pieces)]
}

write_evidence <- function(path, fields) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to write dependency evidence.")
  }
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
    stop("Cannot create evidence directory.")
  }
  temporary <- tempfile("dependency-evidence-", tmpdir = parent)
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(
    fields, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  if (!file.rename(temporary, path)) stop("Cannot finalize dependency evidence.")
}

description_from_tarball <- function(tarball) {
  members <- utils::untar(tarball, list = TRUE)
  if (!length(members) || any(grepl("(^/|(^|/)\\.\\.(/|$)|\\\\)", members))) {
    stop("Tarball member paths are unsafe.")
  }
  roots <- unique(sub("/.*$", "", members))
  if (length(roots) != 1L) stop("Tarball must have exactly one top-level directory.")
  description_member <- paste0(roots[[1L]], "/DESCRIPTION")
  if (sum(members == description_member) != 1L) {
    stop("Tarball must contain exactly one top-level DESCRIPTION.")
  }
  temporary <- tempfile("description-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  utils::untar(tarball, files = description_member, exdir = temporary)
  description <- file.path(temporary, description_member)
  if (!file.exists(description) || nzchar(Sys.readlink(description))) {
    stop("DESCRIPTION must extract as a regular non-symlink file.")
  }
  read.dcf(description, all = TRUE)[1L, , drop = FALSE]
}

status <- tryCatch({
  cli <- parse_cli(args)
  environment_id <- cli[["environment-id"]]
  if (!environment_id %in% c("mkl", "nosuggests")) {
    stop("environment-id must be mkl or nosuggests.")
  }
  install_package <- scalar_bool(cli[["install-package"]], "install-package")
  plan_only_value <- if (is.null(cli[["plan-only"]])) "false" else cli[["plan-only"]]
  plan_only <- scalar_bool(plan_only_value, "plan-only")
  if (plan_only && !identical(Sys.getenv("CI_DEPENDENCY_PREP_TESTING"), "1")) {
    stop("plan-only is test-only.")
  }
  if (environment_id == "mkl" && install_package) {
    stop("MKL source-mode prep must not install the target package.")
  }
  if (environment_id == "nosuggests" && !install_package) {
    stop("noSuggests prep must install the target with tests.")
  }

  tarball <- absolute(cli$tarball, "tarball", must_work = TRUE)
  tarball_info <- file.info(tarball)
  if (nzchar(Sys.readlink(tarball)) || is.na(tarball_info$isdir) || tarball_info$isdir) {
    stop("tarball must be a regular non-symlink file.")
  }
  library <- absolute(cli$library, "library")
  evidence <- absolute(cli$evidence, "evidence")
  if (dir.exists(library)) {
    if (length(list.files(library, all.files = TRUE, no.. = TRUE))) {
      stop("library must be initially empty.")
    }
  } else if (!dir.create(library, recursive = TRUE)) {
    stop("Cannot create library.")
  }

  description <- description_from_tarball(tarball)
  field <- function(name) if (name %in% colnames(description)) description[1L, name] else NA_character_
  hard <- unique(c(
    split_dependencies(field("Depends")),
    split_dependencies(field("Imports")),
    split_dependencies(field("LinkingTo"))
  ))
  suggests <- split_dependencies(field("Suggests"))
  vignette_builders <- split_dependencies(field("VignetteBuilder"))
  testing <- c("testthat", "RUnit", "tinytest")
  selected <- suggests[suggests %in% unique(c(testing, vignette_builders))]
  excluded <- setdiff(suggests, selected)
  dependency_policy <- if (environment_id == "mkl") {
    "all"
  } else {
    "hard-plus-testing-and-vignette-builder"
  }
  ci_tooling <- if (environment_id == "mkl") {
    c("devtools", "jsonlite", "yaml")
  } else {
    character()
  }

  installed_here <- character()
  selected_refs <- selected
  availability <- NULL
  target_installed <- FALSE
  installed_tests <- FALSE
  if (!plan_only) {
    if (!requireNamespace("pak", quietly = TRUE)) stop("Package 'pak' is required.")
    old_paths <- .libPaths()
    on.exit(.libPaths(old_paths), add = TRUE)
    .libPaths(unique(c(library, old_paths)))
    ref <- paste0("deps::", tarball)
    if (environment_id == "mkl") {
      pak::pkg_install(ref, dependencies = TRUE, lib = library)
      pak::pkg_install(ci_tooling, lib = library)
    } else {
      suggested <- pak::pkg_deps(
        ref,
        dependencies = list(direct = "Suggests", indirect = character())
      )
      suggested <- suggested[-1L, , drop = FALSE]
      selected_refs <- suggested[["ref"]][suggested[["package"]] %in% selected]
      if (length(selected_refs) != length(selected)) {
        stop("Could not resolve every allowed direct Suggests reference.")
      }
      pak::pkg_install(c(ref, selected_refs), lib = library)
      r_binary <- file.path(R.home("bin"), "R")
      install_status <- system2(
        r_binary,
        c(
          "CMD", "INSTALL", "--install-tests",
          paste0("--library=", shQuote(library)), shQuote(tarball)
        )
      )
      if (!identical(as.integer(install_status), 0L)) {
        stop("R CMD INSTALL --install-tests failed.")
      }
      installed_here <- sort(rownames(installed.packages(lib.loc = library)))
      if (!all(selected %in% installed_here)) {
        stop("Allowed testing or vignette packages were not installed.")
      }
      if (any(excluded %in% installed_here)) {
        stop("An excluded direct Suggests package was installed in the lane library.")
      }
      if (any(c("devtools", "ashr", "susieR") %in% installed_here)) {
        stop("A forbidden noSuggests package was installed in the lane library.")
      }
      if (any(vapply(c("devtools", "ashr", "susieR"), requireNamespace, logical(1), quietly = TRUE))) {
        stop("A forbidden noSuggests package is available on the effective library path.")
      }
      if (!"colocboost" %in% installed_here) {
        stop("The target package was not installed with its tests.")
      }
      target_installed <- TRUE
      installed_tests <- dir.exists(system.file("tests", package = "colocboost"))
      if (!installed_tests) stop("Installed package tests are missing.")
    }
    installed_here <- sort(rownames(installed.packages(lib.loc = library)))
    availability_names <- if (environment_id == "mkl") {
      unique(c(setdiff(hard, "R"), suggests, ci_tooling))
    } else {
      unique(c(setdiff(hard, "R"), selected, "devtools", "ashr", "susieR"))
    }
    availability <- setNames(
      lapply(availability_names, requireNamespace, quietly = TRUE),
      availability_names
    )
  }

  fields <- list(
    schema_version = 1L,
    kind = "rhub-dependency-proof",
    status = if (plan_only) "plan" else "pass",
    environment_id = environment_id,
    dependency_policy = dependency_policy,
    ci_tooling = ci_tooling,
    hard_dependencies = hard,
    suggested_dependencies = suggests,
    recognized_testing_frameworks = testing,
    vignette_builders = vignette_builders,
    selected_soft_dependencies = if (environment_id == "mkl") suggests else selected,
    selected_refs = selected_refs,
    excluded_suggests = if (environment_id == "mkl") character() else excluded,
    install_package = install_package,
    plan_only = plan_only,
    library = library,
    tarball = tarball,
    installed_packages = installed_here,
    availability = availability,
    target_installed = target_installed,
    installed_tests = installed_tests
  )
  write_evidence(evidence, fields)
  0L
}, error = function(error) {
  message("dependency prep error: ", conditionMessage(error))
  2L
})

quit(save = "no", status = status, runLast = FALSE)
