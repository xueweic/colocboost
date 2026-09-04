#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
required <- c(
  "environment-id", "purpose", "tarball", "library", "evidence",
  "source-sha", "event-sha", "tarball-sha256", "r-executable"
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

absolute_lexical <- function(value, label) {
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", value)) stop(label, " must be absolute.")
  if (grepl("[\r\n]", value)) stop(label, " must be single-line.")
  gsub("\\\\", "/", value)
}

split_dependencies <- function(value) {
  if (is.na(value) || !nzchar(trimws(value))) return(character())
  pieces <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  pieces <- trimws(gsub("\\s*\\([^)]*\\)\\s*$", "", pieces))
  pieces[nzchar(pieces)]
}

description_from_tarball <- function(tarball) {
  members <- utils::untar(tarball, list = TRUE)
  if (!length(members) || any(grepl("(^/|(^|/)\\.\\.(/|$)|\\\\)", members))) {
    stop("Tarball member paths are unsafe.")
  }
  roots <- unique(sub("/.*$", "", members))
  if (length(roots) != 1L) stop("Tarball must have exactly one top-level directory.")
  member <- paste0(roots[[1L]], "/DESCRIPTION")
  if (sum(members == member) != 1L) stop("Tarball must contain one DESCRIPTION.")
  temporary <- tempfile("description-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  utils::untar(tarball, files = member, exdir = temporary)
  path <- file.path(temporary, member)
  info <- file.info(path)
  if (nzchar(Sys.readlink(path)) || is.na(info$isdir) || info$isdir) {
    stop("DESCRIPTION must be a regular non-symlink file.")
  }
  read.dcf(path, all = TRUE)[1L, , drop = FALSE]
}

write_evidence <- function(path, fields) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.")
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
    stop("Cannot create evidence directory.")
  }
  temporary <- tempfile("full-dependency-evidence-", tmpdir = parent)
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(fields, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null")
  if (!file.rename(temporary, path)) stop("Cannot finalize dependency evidence.")
}

status <- tryCatch({
  cli <- parse_cli(args)
  if (!cli$purpose %in% c("unit", "check")) stop("purpose must be unit or check.")
  tarball <- absolute(cli$tarball, "tarball", must_work = TRUE)
  tarball_info <- file.info(tarball)
  if (nzchar(Sys.readlink(tarball)) || is.na(tarball_info$isdir) || tarball_info$isdir) {
    stop("tarball must be a regular non-symlink file.")
  }
  library <- absolute(cli$library, "library")
  evidence <- absolute(cli$evidence, "evidence")
  selected_r <- absolute_lexical(cli[["r-executable"]], "r-executable")
  selected_r_resolved <- absolute(selected_r, "r-executable", must_work = TRUE)
  selected_r_home <- system2(
    selected_r, "RHOME", stdout = TRUE, stderr = FALSE
  )
  rhome_status <- attr(selected_r_home, "status")
  if (!is.null(rhome_status) && rhome_status != 0L) {
    stop("Selected R failed its RHOME identity probe.")
  }
  selected_r_home <- selected_r_home[nzchar(trimws(selected_r_home))]
  if (length(selected_r_home) != 1L) stop("Selected R returned ambiguous RHOME identity.")
  selected_r_home <- absolute(
    trimws(selected_r_home[[1L]]), "selected R home", must_work = TRUE
  )
  active_r_home <- absolute(R.home(), "active R home", must_work = TRUE)
  if (!identical(selected_r_home, active_r_home)) {
    stop("Active R home does not match selected r-executable.")
  }
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
    split_dependencies(field("Depends")), split_dependencies(field("Imports")),
    split_dependencies(field("LinkingTo"))
  ))
  suggests <- split_dependencies(field("Suggests"))
  builders <- split_dependencies(field("VignetteBuilder"))
  tooling <- c("jsonlite", "yaml")

  old_paths <- .libPaths()
  on.exit(.libPaths(old_paths), add = TRUE)
  .libPaths(unique(c(library, .Library)))
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  if (!requireNamespace("pak", quietly = TRUE)) {
    utils::install.packages("pak", lib = library, quiet = TRUE)
  }
  if (!requireNamespace("pak", quietly = TRUE)) stop("Package 'pak' is required.")
  required_here <- unique(c(
    setdiff(hard, "R"), setdiff(suggests, "MASS"), tooling
  ))
  pak::pkg_install(paste0("deps::", tarball), dependencies = TRUE, lib = library)
  pak::pkg_install(required_here, lib = library)

  installed <- sort(rownames(installed.packages(lib.loc = library)))
  availability_names <- unique(c(setdiff(hard, "R"), suggests, tooling))
  availability <- setNames(
    lapply(availability_names, requireNamespace, quietly = TRUE),
    availability_names
  )
  package_origins <- setNames(
    lapply(availability_names, function(name) {
      normalizePath(find.package(name, quiet = FALSE), winslash = "/", mustWork = TRUE)
    }),
    availability_names
  )
  if (!all(unlist(availability, use.names = FALSE))) {
    stop("At least one full dependency is unavailable after installation.")
  }
  if (!all(required_here %in% installed)) {
    stop("Fresh lane library lacks a required explicitly installed package.")
  }

  fields <- list(
    schema_version = 1L,
    kind = "full-dependency-proof",
    status = "pass",
    environment_id = cli[["environment-id"]],
    purpose = cli$purpose,
    dependency_policy = "all",
    source_sha = cli[["source-sha"]],
    event_sha = cli[["event-sha"]],
    tarball_sha256 = cli[["tarball-sha256"]],
    r_executable = selected_r,
    r_resolved = selected_r_resolved,
    r_version = R.version.string,
    r_platform = R.version$platform,
    library = library,
    tarball = tarball,
    hard_dependencies = hard,
    suggested_dependencies = suggests,
    vignette_builders = I(builders),
    tooling = tooling,
    installed_packages = installed,
    availability = availability,
    package_origins = package_origins,
    plan_only = FALSE
  )
  write_evidence(evidence, fields)
  0L
}, error = function(error) {
  message("full dependency prep error: ", conditionMessage(error))
  2L
})

quit(save = "no", status = status, runLast = FALSE)
