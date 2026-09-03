#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
initial_working_directory <- getwd()

required_keys <- c(
  "log", "reported-log-path", "policy", "output", "environment-id",
  "source-sha", "event-sha", "tarball-sha256", "check-exit-code"
)

absolute_path <- function(path, must_work = FALSE) {
  expanded <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", expanded)) {
    expanded <- file.path(initial_working_directory, expanded)
  }
  normalizePath(expanded, winslash = "/", mustWork = must_work)
}

extract_recoverable <- function(key) {
  prefix <- paste0("--", key, "=")
  values <- args[startsWith(args, prefix)]
  if (length(values) != 1L) return(NULL)
  value <- substring(values, nchar(prefix) + 1L)
  if (nzchar(value)) value else NULL
}

recovered <- list(
  output = extract_recoverable("output"),
  environment_id = extract_recoverable("environment-id"),
  source_sha = extract_recoverable("source-sha"),
  event_sha = extract_recoverable("event-sha"),
  tarball_sha256 = extract_recoverable("tarball-sha256")
)

is_identity_recoverable <- function(identity) {
  is.character(identity$output) && length(identity$output) == 1L &&
    nzchar(identity$output) &&
    is.character(identity$environment_id) &&
    length(identity$environment_id) == 1L && nzchar(identity$environment_id) &&
    is.character(identity$source_sha) &&
    grepl("^[0-9a-f]{40}$", identity$source_sha) &&
    is.character(identity$event_sha) &&
    grepl("^[0-9a-f]{40}$", identity$event_sha) &&
    is.character(identity$tarball_sha256) &&
    grepl("^[0-9a-f]{64}$", identity$tarball_sha256)
}

write_json_atomic <- function(document, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to write the check result.")
  }
  output_path <- absolute_path(path)
  parent <- dirname(output_path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
    stop(sprintf("Cannot create result directory: %s", parent))
  }
  temporary <- tempfile("check-result-", tmpdir = parent, fileext = ".json")
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(
    document,
    temporary,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null"
  )
  if (!file.rename(temporary, output_path)) {
    stop(sprintf("Cannot atomically move check result into place: %s", output_path))
  }
}

infrastructure_result <- function(message_text) {
  list(
    schema_version = 1L,
    result_kind = "infrastructure",
    environment_id = recovered$environment_id,
    source_sha = recovered$source_sha,
    event_sha = recovered$event_sha,
    tarball_sha256 = recovered$tarball_sha256,
    status = "fail",
    classification = "check-log-classification",
    message = message_text
  )
}

fail_infrastructure <- function(message_text) {
  if (is_identity_recoverable(recovered)) {
    tryCatch(
      write_json_atomic(infrastructure_result(message_text), recovered$output),
      error = function(error) {
        message("Unable to write infrastructure result: ", conditionMessage(error))
      }
    )
  }
  message(message_text)
  quit(save = "no", status = 2L, runLast = FALSE)
}

parse_cli <- function(values) {
  malformed <- values[!grepl("^--[A-Za-z][A-Za-z0-9-]*=.+$", values)]
  if (length(malformed) > 0L) {
    stop(sprintf("Malformed argument: %s", malformed[[1L]]))
  }
  keys <- sub("^--([^=]+)=.*$", "\\1", values)
  parsed_values <- sub("^--[^=]+=", "", values)
  unknown <- setdiff(unique(keys), required_keys)
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown argument: --%s", unknown[[1L]]))
  }
  duplicated <- unique(keys[duplicated(keys)])
  if (length(duplicated) > 0L) {
    stop(sprintf("Duplicate argument: --%s", duplicated[[1L]]))
  }
  missing <- setdiff(required_keys, keys)
  if (length(missing) > 0L) {
    stop(sprintf("Missing required argument: --%s", missing[[1L]]))
  }
  parsed <- as.list(parsed_values)
  names(parsed) <- keys
  if (!grepl("^[0-9]+$", parsed[["check-exit-code"]])) {
    stop("--check-exit-code must be a nonnegative integer.")
  }
  parsed[["check-exit-code"]] <- as.integer(parsed[["check-exit-code"]])
  if (is.na(parsed[["check-exit-code"]])) {
    stop("--check-exit-code is too large.")
  }
  parsed
}

cli <- tryCatch(parse_cli(args), error = identity)
if (inherits(cli, "error")) {
  fail_infrastructure(conditionMessage(cli))
}

recovered <- list(
  output = cli$output,
  environment_id = cli[["environment-id"]],
  source_sha = cli[["source-sha"]],
  event_sha = cli[["event-sha"]],
  tarball_sha256 = cli[["tarball-sha256"]]
)
if (!is_identity_recoverable(recovered)) {
  fail_infrastructure("The result identity is not valid.")
}

log_path <- tryCatch(absolute_path(cli$log, must_work = TRUE), error = identity)
if (inherits(log_path, "error") || !file.exists(log_path) || dir.exists(log_path)) {
  fail_infrastructure(sprintf("Check log is not a regular file: %s", cli$log))
}
policy_path <- tryCatch(absolute_path(cli$policy, must_work = TRUE), error = identity)
if (inherits(policy_path, "error") || !file.exists(policy_path) || dir.exists(policy_path)) {
  fail_infrastructure(sprintf("Check policy is not a regular file: %s", cli$policy))
}
if (!nzchar(cli[["reported-log-path"]])) {
  fail_infrastructure("--reported-log-path must be nonempty.")
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  fail_infrastructure("Package 'yaml' is required to read the check policy.")
}

is_nonempty_string <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

parse_expiry <- function(value) {
  if (!is_nonempty_string(value) ||
      !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", value)) {
    return(as.Date(NA))
  }
  parsed <- suppressWarnings(as.Date(value, format = "%Y-%m-%d"))
  if (is.na(parsed) || format(parsed, "%Y-%m-%d") != value) as.Date(NA) else parsed
}

validate_note_policy <- function(policy) {
  if (!is.list(policy) || is.null(names(policy))) {
    stop("Policy root must be a named mapping.")
  }
  if (!is.numeric(policy$version) || length(policy$version) != 1L ||
      is.na(policy$version) || policy$version != 1) {
    stop("Policy version must be exactly 1.")
  }
  section <- policy$r_cmd_check
  if (!is.list(section) || is.null(names(section)) ||
      !identical(sort(names(section)), "allowed_notes")) {
    stop("Policy r_cmd_check must contain only allowed_notes.")
  }
  rules <- section$allowed_notes
  if (!is.list(rules)) {
    stop("Policy r_cmd_check.allowed_notes must be a list.")
  }
  required <- c("id", "check", "output", "rationale", "expires")
  for (index in seq_along(rules)) {
    rule <- rules[[index]]
    if (!is.list(rule) || is.null(names(rule)) ||
        !identical(sort(names(rule)), sort(required))) {
      stop(sprintf("NOTE allowance %d has unexpected or missing fields.", index))
    }
    invalid <- required[!vapply(rule[required], is_nonempty_string, logical(1))]
    if (length(invalid) > 0L) {
      stop(sprintf(
        "NOTE allowance %d field '%s' must be a nonempty string.",
        index,
        invalid[[1L]]
      ))
    }
    expiry <- parse_expiry(rule$expires)
    if (is.na(expiry)) {
      stop(sprintf("NOTE allowance %d has an invalid expiry date.", index))
    }
    if (expiry < Sys.Date()) {
      stop(sprintf("NOTE allowance '%s' expired on %s.", rule$id, rule$expires))
    }
  }
  if (length(rules) > 0L) {
    ids <- vapply(rules, `[[`, character(1), "id")
    if (anyDuplicated(ids)) {
      stop(sprintf("Duplicate NOTE allowance id: %s", ids[duplicated(ids)][[1L]]))
    }
    keys <- vapply(rules, function(rule) {
      paste(rule$check, rule$output, sep = "\r")
    }, character(1))
    if (anyDuplicated(keys)) {
      stop("Multiple NOTE allowances match the same exact check chunk.")
    }
  }
  rules
}

policy <- tryCatch(yaml::read_yaml(policy_path), error = identity)
if (inherits(policy, "error")) {
  fail_infrastructure(paste("Unable to parse check policy:", conditionMessage(policy)))
}
allowed_notes <- tryCatch(validate_note_policy(policy), error = identity)
if (inherits(allowed_notes, "error")) {
  fail_infrastructure(conditionMessage(allowed_notes))
}

read_log_lines <- function(path) {
  size <- file.info(path)$size
  if (is.na(size)) stop("Unable to determine check-log size.")
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  text <- readChar(connection, nchars = size, useBytes = TRUE)
  text <- gsub("\r\n", "\n", text, fixed = TRUE)
  if (grepl("\r", text, fixed = TRUE)) {
    stop("Check log contains unsupported bare carriage returns.")
  }
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  while (length(lines) > 0L && !nzchar(lines[[length(lines)]])) {
    lines <- lines[-length(lines)]
  }
  lines
}

parse_footer <- function(lines) {
  if (length(lines) == 0L) stop("Check log is empty.")
  done <- which(lines == "* DONE")
  status_positions <- which(startsWith(lines, "Status:"))
  if (length(done) != 1L || length(status_positions) != 1L) {
    stop("Check log must contain one unique '* DONE' and terminal Status footer.")
  }
  status_position <- status_positions[[1L]]
  if (status_position != done[[1L]] + 1L) {
    stop("Status footer must immediately follow '* DONE'.")
  }
  tail <- if (status_position < length(lines)) {
    lines[seq.int(status_position + 1L, length(lines))]
  } else {
    character()
  }
  if (length(tail) > 1L ||
      (length(tail) == 1L && !grepl("^\\* elapsed time .+$", tail))) {
    stop("Check log contains an unrecognized footer postlude.")
  }

  status_text <- substring(lines[[status_position]], nchar("Status: ") + 1L)
  if (identical(lines[[status_position]], "Status: OK")) {
    counts <- c(errors = 0L, warnings = 0L, notes = 0L)
  } else {
    pieces <- strsplit(status_text, ", ", fixed = TRUE)[[1L]]
    pattern <- "^([1-9][0-9]*) (ERROR|WARNING|NOTE)(s?)$"
    matched <- regexec(pattern, pieces, perl = TRUE)
    captures <- regmatches(pieces, matched)
    if (any(lengths(captures) != 4L)) {
      stop("Status footer has malformed counts.")
    }
    values <- vapply(captures, function(value) as.integer(value[[2L]]), integer(1))
    kinds <- vapply(captures, `[[`, character(1), 3L)
    plurals <- vapply(captures, `[[`, character(1), 4L)
    expected_plural <- ifelse(values == 1L, "", "s")
    if (any(plurals != expected_plural) || any(is.na(values))) {
      stop("Status footer has invalid singular or plural counts.")
    }
    expected_order <- c("ERROR", "WARNING", "NOTE")
    order_indices <- match(kinds, expected_order)
    if (anyDuplicated(kinds) || any(is.na(order_indices)) ||
        !identical(order_indices, sort(order_indices))) {
      stop("Status footer kinds must be unique and ordered ERROR, WARNING, NOTE.")
    }
    counts <- c(errors = 0L, warnings = 0L, notes = 0L)
    names_by_kind <- c(ERROR = "errors", WARNING = "warnings", NOTE = "notes")
    counts[names_by_kind[kinds]] <- values
  }
  list(body = lines[seq_len(done[[1L]] - 1L)], counts = counts)
}

parse_chunks <- function(lines) {
  if (length(lines) == 0L) stop("Check log contains no check headings.")
  drop_parent <- rep(FALSE, length(lines))
  if (length(lines) > 1L) {
    for (index in seq_len(length(lines) - 1L)) {
      drop_parent[[index]] <-
        (startsWith(lines[[index]], "* loading checks for arch") &&
           startsWith(lines[[index + 1L]], "** checking")) ||
        (startsWith(lines[[index]], "* checking examples") &&
           startsWith(lines[[index + 1L]], "** running examples for arch")) ||
        (startsWith(lines[[index]], "* checking tests") &&
           startsWith(lines[[index + 1L]], "** running tests for arch"))
    }
  }
  lines <- lines[!drop_parent]

  heading_pattern <- paste0(
    "^\\*\\*? ((checking|creating|running examples for arch|",
    "running tests for arch) .*) \\.\\.\\.( (\\[[^ ]*\\]))?( (.*)|)$"
  )
  matches <- regexec(heading_pattern, lines, perl = TRUE)
  captured <- regmatches(lines, matches)
  heading_positions <- which(lengths(captured) > 0L)
  suspicious <- which(
    grepl(
      "^\\*\\*? (checking|creating|running examples for arch|running tests for arch)( |$)",
      lines,
      perl = TRUE
    ) & lengths(captured) == 0L
  )
  if (length(suspicious) > 0L) {
    stop(sprintf("Malformed check heading at line %d.", suspicious[[1L]]))
  }
  if (length(heading_positions) == 0L) {
    stop("Check log contains no check headings.")
  }

  chunks <- vector("list", length(heading_positions))
  allowed_statuses <- c("OK", "NONE", "SKIPPED", "INFO", "ERROR", "WARNING", "NOTE")
  for (chunk_index in seq_along(heading_positions)) {
    position <- heading_positions[[chunk_index]]
    values <- captured[[position]]
    status <- values[[7L]]
    if (!nzchar(status)) {
      stop(sprintf("Check heading at line %d has no terminal status.", position))
    }
    is_extension_type <- identical(values[[2L]], "checking extension type") &&
      identical(status, "Package")
    if (!status %in% allowed_statuses && !is_extension_type) {
      stop(sprintf("Check heading at line %d has unknown status '%s'.", position, status))
    }
    next_position <- if (chunk_index < length(heading_positions)) {
      heading_positions[[chunk_index + 1L]]
    } else {
      length(lines) + 1L
    }
    output_lines <- if (position + 1L < next_position) {
      lines[seq.int(position + 1L, next_position - 1L)]
    } else {
      character()
    }
    chunks[[chunk_index]] <- list(
      check = values[[2L]],
      status = status,
      output = paste(output_lines, collapse = "\n")
    )
  }
  chunks
}

parsed <- tryCatch({
  footer <- parse_footer(read_log_lines(log_path))
  chunks <- parse_chunks(footer$body)
  observed_counts <- c(
    errors = sum(vapply(chunks, `[[`, character(1), "status") == "ERROR"),
    warnings = sum(vapply(chunks, `[[`, character(1), "status") == "WARNING"),
    notes = sum(vapply(chunks, `[[`, character(1), "status") == "NOTE")
  )
  observed_counts <- as.integer(observed_counts)
  names(observed_counts) <- c("errors", "warnings", "notes")
  if (!identical(observed_counts, footer$counts)) {
    stop(sprintf(
      "Status footer counts (%s) do not match heading counts (%s).",
      paste(footer$counts, collapse = "/"),
      paste(observed_counts, collapse = "/")
    ))
  }
  list(chunks = chunks, counts = observed_counts)
}, error = identity)
if (inherits(parsed, "error")) {
  fail_infrastructure(conditionMessage(parsed))
}

note_chunks <- parsed$chunks[vapply(
  parsed$chunks,
  function(chunk) identical(chunk$status, "NOTE"),
  logical(1)
)]
note_matches <- if (length(note_chunks) == 0L) {
  integer()
} else {
  vapply(note_chunks, function(chunk) {
    sum(vapply(allowed_notes, function(rule) {
      identical(chunk$check, rule$check) && identical(chunk$output, rule$output)
    }, logical(1)))
  }, integer(1))
}
rule_matches <- if (length(allowed_notes) == 0L) {
  integer()
} else {
  vapply(allowed_notes, function(rule) {
    sum(vapply(note_chunks, function(chunk) {
      identical(chunk$check, rule$check) && identical(chunk$output, rule$output)
    }, logical(1)))
  }, integer(1))
}
notes_allowed <- all(note_matches == 1L) && all(rule_matches == 1L)

package_pass <- parsed$counts[["errors"]] == 0L &&
  parsed$counts[["warnings"]] == 0L && notes_allowed &&
  cli[["check-exit-code"]] == 0L
document <- list(
  schema_version = 1L,
  result_kind = "package-check",
  environment_id = cli[["environment-id"]],
  source_sha = cli[["source-sha"]],
  event_sha = cli[["event-sha"]],
  tarball_sha256 = cli[["tarball-sha256"]],
  status = if (package_pass) "pass" else "fail",
  check = list(
    errors = parsed$counts[["errors"]],
    warnings = parsed$counts[["warnings"]],
    notes = parsed$counts[["notes"]],
    log_path = cli[["reported-log-path"]]
  )
)

write_error <- tryCatch({
  write_json_atomic(document, cli$output)
  NULL
}, error = identity)
if (inherits(write_error, "error")) {
  message(conditionMessage(write_error))
  quit(save = "no", status = 2L, runLast = FALSE)
}

quit(save = "no", status = if (package_pass) 0L else 1L, runLast = FALSE)
