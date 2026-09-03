#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
required <- c("policy", "evidence")
optional <- c("maps-file", "verbose-file")

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

absolute <- function(value, label, must_work = FALSE) {
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", value)) stop(label, " must be absolute.")
  normalizePath(value, winslash = "/", mustWork = must_work)
}

write_json_atomic <- function(document, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for MKL evidence.")
  }
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
    stop("Cannot create MKL evidence directory.")
  }
  temporary <- tempfile("mkl-evidence-", tmpdir = parent, fileext = ".json")
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(document, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null")
  if (!file.rename(temporary, path)) stop("Cannot finalize MKL evidence.")
}

read_regular <- function(path, label) {
  path <- absolute(path, label, must_work = TRUE)
  info <- file.info(path)
  if (nzchar(Sys.readlink(path)) || is.na(info$isdir) || info$isdir) {
    stop(label, " must be a regular non-symlink file.")
  }
  readLines(path, warn = FALSE)
}

run_verbose_probe <- function() {
  rscript <- file.path(R.home("bin"), "Rscript")
  expression <- paste(
    "x <- matrix(as.double(seq_len(4096)), 64L, 64L);",
    "y <- crossprod(x);",
    "stopifnot(all(is.finite(y)));",
    "cat('probe-checksum=', format(sum(y), scientific=FALSE), '\\n', sep='')"
  )
  output <- system2(
    rscript,
    c("--vanilla", "-e", shQuote(expression)),
    stdout = TRUE,
    stderr = TRUE,
    env = c("MKL_VERBOSE=1", "MKL_NUM_THREADS=1", "OMP_NUM_THREADS=1")
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("MKL verbose child probe failed.")
  if (!any(grepl("MKL_VERBOSE", output, fixed = TRUE))) {
    stop("MKL verbose child probe did not emit MKL_VERBOSE evidence.")
  }
  output
}

evidence_path <- NULL
result <- tryCatch({
  cli <- parse_cli(args)
  evidence_path <- absolute(cli$evidence, "evidence")
  policy_path <- absolute(cli$policy, "policy", must_work = TRUE)
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Package 'yaml' is required.")
  policy <- yaml::read_yaml(policy_path)
  backend <- policy$numerical_backends$mkl
  required_patterns <- backend$required_library_patterns
  forbidden_patterns <- backend$forbidden_library_patterns
  if (
    !is.character(required_patterns) || length(required_patterns) != 3L ||
    !is.character(forbidden_patterns) || length(forbidden_patterns) != 4L
  ) stop("MKL policy patterns are malformed.")

  threads <- list(
    MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS", unset = ""),
    OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", unset = "")
  )
  if (!identical(threads$MKL_NUM_THREADS, "1") || !identical(threads$OMP_NUM_THREADS, "1")) {
    stop("MKL_NUM_THREADS and OMP_NUM_THREADS must both equal 1.")
  }
  if (!identical(Sys.getenv("MKL_VERBOSE", unset = ""), "1")) {
    stop("MKL_VERBOSE must equal 1.")
  }

  fixture_args <- !is.null(cli[["maps-file"]]) || !is.null(cli[["verbose-file"]])
  fixture_mode <- identical(Sys.getenv("CI_VERIFY_MKL_TESTING"), "1")
  if (fixture_args && !fixture_mode) stop("MKL fixture inputs are test-only.")
  if (fixture_args && (is.null(cli[["maps-file"]]) || is.null(cli[["verbose-file"]]))) {
    stop("Both MKL fixture inputs are required together.")
  }

  operation <- crossprod(matrix(as.double(seq_len(4096)), 64L, 64L))
  if (!all(is.finite(operation))) stop("BLAS matrix operation returned non-finite values.")
  maps_path <- if (fixture_args) cli[["maps-file"]] else "/proc/self/maps"
  maps <- read_regular(maps_path, "process maps")
  maps_text <- paste(maps, collapse = "\n")
  missing <- required_patterns[!vapply(
    required_patterns, grepl, logical(1), x = maps_text, perl = TRUE
  )]
  forbidden <- forbidden_patterns[vapply(
    forbidden_patterns, grepl, logical(1), x = maps_text, perl = TRUE
  )]
  if (length(missing)) stop("Required serial MKL process mapping is missing: ", missing[[1L]])
  if (length(forbidden)) stop("Forbidden numerical library is loaded: ", forbidden[[1L]])

  verbose <- if (fixture_args) {
    read_regular(cli[["verbose-file"]], "MKL verbose fixture")
  } else {
    run_verbose_probe()
  }
  if (!any(grepl("MKL_VERBOSE", verbose, fixed = TRUE))) {
    stop("MKL_VERBOSE execution evidence is missing.")
  }

  setvars_path <- "/opt/intel/oneapi/setvars.sh"
  setvars <- list(
    path = setvars_path,
    exists = file.exists(setvars_path) && !nzchar(Sys.readlink(setvars_path)),
    SETVARS_ARGS = Sys.getenv("SETVARS_ARGS", unset = ""),
    MKLROOT = Sys.getenv("MKLROOT", unset = ""),
    LD_LIBRARY_PATH = Sys.getenv("LD_LIBRARY_PATH", unset = "")
  )
  if (!fixture_mode && (!isTRUE(setvars$exists) || !nzchar(setvars$MKLROOT))) {
    stop("Intel oneAPI setvars evidence is missing.")
  }

  r_binary <- file.path(R.home("bin"), "R")
  blas_libs <- system2(r_binary, c("CMD", "config", "BLAS_LIBS"), stdout = TRUE, stderr = TRUE)
  document <- list(
    schema_version = 1L,
    kind = "mkl-runtime-proof",
    status = "pass",
    r_executable = normalizePath(r_binary, mustWork = TRUE),
    r_version = as.character(getRversion()),
    session_info = capture.output(sessionInfo()),
    ext_soft_version = as.list(extSoftVersion()),
    la_library = tryCatch(La_library(), error = function(error) NA_character_),
    la_version = tryCatch(as.character(La_version()), error = function(error) NA_character_),
    blas_libs = blas_libs,
    threads = threads,
    setvars = setvars,
    maps_path = normalizePath(maps_path, mustWork = TRUE),
    required_library_patterns = required_patterns,
    forbidden_library_patterns = forbidden_patterns,
    loaded_library_lines = maps[grepl("\\.so", maps)],
    mkl_verbose = paste(verbose, collapse = "\n"),
    blas_operation = list(completed = TRUE, checksum = sum(operation)),
    fixture_mode = fixture_mode
  )
  write_json_atomic(document, evidence_path)
  list(status = 0L, document = document)
}, error = function(error) {
  message("MKL verification error: ", conditionMessage(error))
  if (!is.null(evidence_path) && requireNamespace("jsonlite", quietly = TRUE)) {
    try(write_json_atomic(list(
      schema_version = 1L,
      kind = "mkl-runtime-proof",
      status = "fail",
      message = conditionMessage(error)
    ), evidence_path), silent = TRUE)
  }
  list(status = 2L)
})

quit(save = "no", status = result$status, runLast = FALSE)
