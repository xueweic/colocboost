#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
required <- c(
  "environment-id", "runtime-profile", "r-executable", "source-sha",
  "event-sha", "tarball-sha256", "evidence"
)

parse_cli <- function(values) {
  if (any(!grepl("^--[a-z][a-z0-9-]*=.+$", values))) {
    stop("Every argument must use --name=value syntax.")
  }
  keys <- sub("^--([^=]+)=.*$", "\\1", values)
  vals <- sub("^--[^=]+=", "", values)
  if (anyDuplicated(keys)) stop("Duplicate argument: --", keys[duplicated(keys)][[1L]])
  unknown <- setdiff(keys, required)
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

config_value <- function(r, name) {
  output <- system2(r, c("CMD", "config", name), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("R CMD config failed for ", name, ".")
  paste(output, collapse = " ")
}

compiler_evidence <- function(configuration, label) {
  command <- strsplit(trimws(configuration), "[[:space:]]+")[[1L]][[1L]]
  path <- unname(Sys.which(command))
  if (!nzchar(path)) stop(label, " compiler executable is unavailable.")
  path <- absolute(path, paste(label, "compiler"), must_work = TRUE)
  output <- system2(path, "--version", stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || !length(output) || !nzchar(output[[1L]])) {
    stop(label, " compiler version probe failed.")
  }
  list(path = path, version = output[[1L]])
}

os_release_evidence <- function() {
  path <- "/etc/os-release"
  if (!file.exists(path)) stop("/etc/os-release is required for Linux identity.")
  lines <- readLines(path, warn = TRUE)
  value <- function(key) {
    match <- lines[startsWith(lines, paste0(key, "="))]
    if (length(match) != 1L) stop("Missing or duplicate ", key, " in /etc/os-release.")
    result <- sub("^[^=]+=", "", match)
    result <- sub('^"(.*)"$', "\\1", result)
    if (!nzchar(result) || grepl("[\\r\\n]", result)) stop("Invalid ", key, " in /etc/os-release.")
    result
  }
  list(id = tolower(value("ID")), version = value("VERSION_ID"))
}

write_evidence <- function(path, fields) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.")
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
    stop("Cannot create evidence directory.")
  }
  temporary <- tempfile("runtime-evidence-", tmpdir = parent)
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(fields, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null")
  if (!file.rename(temporary, path)) stop("Cannot finalize runtime evidence.")
}

status <- tryCatch({
  cli <- parse_cli(args)
  selected_r <- absolute(cli[["r-executable"]], "r-executable", must_work = TRUE)
  actual_r <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
  actual_r <- absolute(actual_r, "active R", must_work = TRUE)
  if (!identical(selected_r, actual_r)) stop("Active R does not match selected r-executable.")
  evidence <- absolute(cli$evidence, "evidence")

  matrix_dimension <- 3L
  lhs <- matrix(seq_len(matrix_dimension^2L), nrow = matrix_dimension)
  rhs <- matrix(rev(seq_len(matrix_dimension^2L)), nrow = matrix_dimension)
  product <- lhs %*% rhs
  matrix_checksum <- sum(product)
  if (!is.finite(matrix_checksum) || matrix_checksum == 0) stop("BLAS operation failed.")

  maps_path <- "/proc/self/maps"
  if (!file.exists(maps_path)) stop("/proc/self/maps is required for in-process proof.")
  maps <- readLines(maps_path, warn = TRUE)
  paths <- sub("^.*[[:space:]](/[^[:space:]]+)$", "\\1", maps)
  paths <- sort(unique(paths[grepl("^/", paths)]))
  if (!length(paths)) stop("No loaded libraries were observed in this R process.")

  sysname <- tolower(unname(Sys.info()[["sysname"]]))
  os_name <- if (identical(sysname, "darwin")) "darwin" else if (grepl("windows", sysname)) "windows" else sysname
  architecture <- R.version$arch
  if (is.null(architecture) || !nzchar(architecture)) architecture <- .Platform$r_arch
  architecture <- sub("^amd64$", "x86_64", architecture)
  cc <- config_value(selected_r, "CC")
  cxx <- config_value(selected_r, "CXX")
  fc <- config_value(selected_r, "FC")
  cc_evidence <- compiler_evidence(cc, "C")
  cxx_evidence <- compiler_evidence(cxx, "C++")
  fc_evidence <- compiler_evidence(fc, "Fortran")
  distribution <- os_release_evidence()
  compiler <- cc_evidence$version
  la_library <- if (exists("La_library", mode = "function")) La_library() else "unavailable"
  la_version <- if (exists("La_version", mode = "function")) La_version() else "unavailable"
  environment_names <- c(
    "PATH", "R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE",
    "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "CC", "CXX", "FC", "F77",
    "LANG", "LC_ALL", "LC_CTYPE", "ASAN_OPTIONS", "UBSAN_OPTIONS",
    "LD_PRELOAD", "_R_CHECK_DONTTEST_EXAMPLES_"
  )
  environment <- as.list(Sys.getenv(environment_names, unset = NA_character_))
  names(environment) <- environment_names

  fields <- list(
    schema_version = 1L,
    kind = "r-runtime-proof",
    status = "pass",
    environment_id = cli[["environment-id"]],
    runtime_profile = cli[["runtime-profile"]],
    source_sha = cli[["source-sha"]],
    event_sha = cli[["event-sha"]],
    tarball_sha256 = cli[["tarball-sha256"]],
    r_executable = selected_r,
    r_resolved = actual_r,
    r_version = R.version.string,
    r_status = R.version$status,
    r_platform = R.version$platform,
    os = os_name,
    os_release = unname(Sys.info()[["release"]]),
    distribution_id = distribution$id,
    distribution_version = distribution$version,
    architecture = architecture,
    compiler = compiler,
    cc = cc,
    cc_path = cc_evidence$path,
    cc_version = cc_evidence$version,
    cxx = cxx,
    cxx_path = cxx_evidence$path,
    cxx_version = cxx_evidence$version,
    fc = fc,
    fc_path = fc_evidence$path,
    fc_version = fc_evidence$version,
    locale = Sys.getlocale(),
    session_info = paste(capture.output(sessionInfo()), collapse = " | "),
    ext_soft_version = as.list(extSoftVersion()),
    la_library = as.character(la_library),
    la_version = as.character(la_version),
    blas_libs = config_value(selected_r, "BLAS_LIBS"),
    matrix_dimension = matrix_dimension,
    matrix_checksum = matrix_checksum,
    maps_method = "proc-self-maps",
    loaded_libraries = paths,
    long_double = .Machine$sizeof.longdouble > .Machine$sizeof.double,
    environment = environment
  )
  write_evidence(evidence, fields)
  0L
}, error = function(error) {
  message("runtime probe error: ", conditionMessage(error))
  2L
})

quit(save = "no", status = status, runLast = FALSE)
