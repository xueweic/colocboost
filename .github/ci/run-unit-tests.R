#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
initial_working_directory <- getwd()

required_keys <- c("package", "load-package", "context", "policy", "output")
optional_keys <- "filter"
allowed_keys <- c(required_keys, optional_keys)
forbidden_skip_sentinel <- paste0(
  "One run returned no colocalization; ",
  "entropy-vs-uniform comparison skipped"
)

absolute_path <- function(path, must_work = FALSE) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path <- file.path(initial_working_directory, path)
  }
  normalizePath(path, winslash = "/", mustWork = must_work)
}

extract_recoverable <- function(key) {
  prefix <- paste0("--", key, "=")
  values <- args[startsWith(args, prefix)]
  if (length(values) != 1L) {
    return(NULL)
  }
  value <- substring(values, nchar(prefix) + 1L)
  if (!nzchar(value)) NULL else value
}

recovered_output <- extract_recoverable("output")
if (!is.null(recovered_output)) {
  recovered_output <- absolute_path(recovered_output)
}
recovered_context <- extract_recoverable("context")
recovered_load_package <- extract_recoverable("load-package")
recovered_package <- extract_recoverable("package")
recovered_filter <- extract_recoverable("filter")

empty_summary <- function() {
  list(
    total = 0L,
    pass = 0L,
    failure = 0L,
    error = 0L,
    warning = 0L,
    empty = 0L,
    skip_like_success = 0L,
    skip = 0L,
    observations = list(
      total = 0L,
      error = 0L,
      failure = 0L,
      warning = 0L,
      empty = 0L,
      skip_like_success = 0L,
      skip = 0L,
      pass = 0L
    )
  )
}

new_report <- function(
  status,
  context,
  load_package,
  package,
  summary = empty_summary(),
  cases = list(),
  matched_allowances = list(),
  unused_allowances = list(),
  violations = list(),
  filter = recovered_filter
) {
  list(
    schema_version = 1L,
    kind = "unit-tests",
    context = context,
    status = status,
    environment = list(
      id = context,
      load_package = load_package,
      package = package,
      r_version = as.character(getRversion()),
      suite = if (is.null(filter)) "full" else "filtered",
      filter = filter
    ),
    summary = summary,
    cases = unname(cases),
    matched_allowances = unname(matched_allowances),
    unused_allowances = unname(unused_allowances),
    violations = unname(violations)
  )
}

write_report <- function(report, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to write the unit-test report.")
  }
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
    stop(sprintf("Cannot create report directory: %s", parent))
  }
  temporary <- tempfile("unit-report-", tmpdir = parent, fileext = ".json")
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(
    report,
    temporary,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null"
  )
  if (!file.rename(temporary, path)) {
    stop(sprintf("Cannot move unit-test report into place: %s", path))
  }
}

exit_with_diagnostic <- function(type, message_text, status = 2L) {
  report <- new_report(
    status = "error",
    context = recovered_context,
    load_package = recovered_load_package,
    package = recovered_package,
    violations = list(list(type = type, message = message_text))
  )
  if (!is.null(recovered_output)) {
    tryCatch(
      write_report(report, recovered_output),
      error = function(error) {
        message("Unable to write diagnostic JSON: ", conditionMessage(error))
      }
    )
  }
  message(message_text)
  quit(save = "no", status = status, runLast = FALSE)
}

parse_cli <- function(values) {
  malformed <- values[!grepl("^--[A-Za-z][A-Za-z0-9-]*=.+$", values)]
  if (length(malformed) > 0L) {
    stop(sprintf("Malformed argument: %s", malformed[[1L]]))
  }
  keys <- sub("^--([^=]+)=.*$", "\\1", values)
  parsed_values <- sub("^--[^=]+=", "", values)
  unknown <- setdiff(unique(keys), allowed_keys)
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown argument: --%s", unknown[[1L]]))
  }
  duplicated_keys <- unique(keys[duplicated(keys)])
  if (length(duplicated_keys) > 0L) {
    stop(sprintf("Duplicate argument: --%s", duplicated_keys[[1L]]))
  }
  missing <- setdiff(required_keys, keys)
  if (length(missing) > 0L) {
    stop(sprintf("Missing required argument: --%s", missing[[1L]]))
  }
  result <- as.list(parsed_values)
  names(result) <- keys
  if (!result[["load-package"]] %in% c("source", "installed")) {
    stop("--load-package must be 'source' or 'installed'.")
  }
  result
}

cli <- tryCatch(parse_cli(args), error = identity)
if (inherits(cli, "error")) {
  exit_with_diagnostic("cli-error", conditionMessage(cli))
}

package_path <- tryCatch(
  absolute_path(cli$package, must_work = TRUE),
  error = identity
)
if (inherits(package_path, "error") || !dir.exists(package_path)) {
  exit_with_diagnostic("filesystem-error", sprintf(
    "Package directory does not exist: %s",
    cli$package
  ))
}

policy_path <- tryCatch(
  absolute_path(cli$policy, must_work = TRUE),
  error = identity
)
if (inherits(policy_path, "error") || !file.exists(policy_path)) {
  exit_with_diagnostic("filesystem-error", sprintf(
    "Policy file does not exist: %s",
    cli$policy
  ))
}

output_path <- absolute_path(cli$output)
recovered_output <- output_path
recovered_context <- cli$context
recovered_load_package <- cli[["load-package"]]
recovered_package <- package_path
recovered_filter <- cli$filter

if (!requireNamespace("yaml", quietly = TRUE)) {
  exit_with_diagnostic("infrastructure-error", "Package 'yaml' is required.")
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  exit_with_diagnostic("infrastructure-error", "Package 'testthat' is required.")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  exit_with_diagnostic("infrastructure-error", "Package 'jsonlite' is required.")
}
suppressWarnings(suppressPackageStartupMessages(library(testthat)))

is_nonempty_string <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

parse_expiry <- function(value) {
  if (!is_nonempty_string(value) || !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", value)) {
    return(as.Date(NA))
  }
  parsed <- tryCatch(as.Date(value, format = "%Y-%m-%d"), error = function(error) as.Date(NA))
  if (is.na(parsed) || format(parsed, "%Y-%m-%d") != value) as.Date(NA) else parsed
}

validate_policy <- function(policy) {
  if (!is.list(policy) || is.null(names(policy))) {
    stop("Policy root must be a named mapping.")
  }
  if (!is.numeric(policy$version) || length(policy$version) != 1L ||
      is.na(policy$version) || policy$version != 1) {
    stop("Policy version must be exactly 1.")
  }
  if (!is.list(policy$unit_tests) || is.null(policy$unit_tests$allowed_skips)) {
    stop("Policy unit_tests.allowed_skips must be a list.")
  }
  allowances <- policy$unit_tests$allowed_skips
  if (!is.list(allowances)) {
    stop("Policy unit_tests.allowed_skips must be a list.")
  }
  required <- c(
    "id", "context", "file", "test_title", "reason", "expected_count",
    "rationale", "expires"
  )
  optional <- character()
  for (index in seq_along(allowances)) {
    allowance <- allowances[[index]]
    if (!is.list(allowance) || is.null(names(allowance))) {
      stop(sprintf("Skip allowance %d must be a named mapping.", index))
    }
    missing <- setdiff(required, names(allowance))
    unknown <- setdiff(names(allowance), c(required, optional))
    if (length(missing) > 0L) {
      stop(sprintf(
        "Skip allowance %d is missing field '%s'.",
        index,
        missing[[1L]]
      ))
    }
    if (length(unknown) > 0L) {
      stop(sprintf(
        "Skip allowance %d has unknown field '%s'.",
        index,
        unknown[[1L]]
      ))
    }
    string_fields <- setdiff(required, "expected_count")
    invalid_string <- string_fields[!vapply(
      allowance[string_fields],
      is_nonempty_string,
      logical(1)
    )]
    if (length(invalid_string) > 0L) {
      stop(sprintf(
        "Skip allowance %d field '%s' must be a nonempty string.",
        index,
        invalid_string[[1L]]
      ))
    }
    if (!is.numeric(allowance$expected_count) ||
        length(allowance$expected_count) != 1L ||
        is.na(allowance$expected_count) ||
        allowance$expected_count < 1 ||
        allowance$expected_count %% 1 != 0) {
      stop(sprintf(
        "Skip allowance %d expected_count must be a positive integer.",
        index
      ))
    }
    if (basename(allowance$file) != allowance$file) {
      stop(sprintf("Skip allowance %d file must be a basename.", index))
    }
    expiry <- parse_expiry(allowance$expires)
    if (is.na(expiry)) {
      stop(sprintf("Skip allowance %d expires must be a valid YYYY-MM-DD date.", index))
    }
    allowances[[index]]$expected_count <- as.integer(allowance$expected_count)
    allowances[[index]]$.expiry <- expiry
  }
  if (length(allowances) > 0L) {
    ids <- vapply(allowances, `[[`, character(1), "id")
    if (anyDuplicated(ids)) {
      stop(sprintf("Duplicate skip allowance id: %s", ids[duplicated(ids)][[1L]]))
    }
    keys <- vapply(allowances, function(allowance) paste(
      allowance$context,
      allowance$file,
      allowance$test_title,
      allowance$reason,
      sep = "\r"
    ), character(1))
    if (anyDuplicated(keys)) {
      stop("Duplicate skip allowance matching key.")
    }
  }
  allowances
}

policy <- tryCatch(yaml::read_yaml(policy_path), error = identity)
if (inherits(policy, "error")) {
  exit_with_diagnostic("policy-error", paste(
    "Unable to parse policy:", conditionMessage(policy)
  ))
}
allowances <- tryCatch(validate_policy(policy), error = identity)
if (inherits(allowances, "error")) {
  exit_with_diagnostic("policy-error", conditionMessage(allowances))
}

condition_message <- function(condition) {
  message_text <- condition$message
  if (is.null(message_text)) {
    message_text <- tryCatch(conditionMessage(condition), error = function(error) "")
  }
  paste(as.character(message_text), collapse = "\n")
}

condition_call <- function(condition) {
  call <- condition$call
  if (is.null(call)) {
    call <- tryCatch(conditionCall(condition), error = function(error) NULL)
  }
  if (is.null(call)) NULL else paste(deparse(call), collapse = " ")
}

condition_line <- function(condition) {
  reference <- condition$srcref
  if (is.null(reference)) {
    return(NA_integer_)
  }
  values <- suppressWarnings(as.integer(reference))
  if (length(values) == 0L || is.na(values[[1L]])) NA_integer_ else values[[1L]]
}

condition_file <- function(condition) {
  reference <- condition$srcref
  if (is.null(reference)) {
    return(NULL)
  }
  filename <- tryCatch(
    getSrcFilename(reference, full.names = TRUE),
    error = function(error) NULL
  )
  if (is.null(filename) || length(filename) == 0L || !nzchar(filename[[1L]])) {
    NULL
  } else {
    basename(filename[[1L]])
  }
}

active_source_evidence <- function() {
  calls <- sys.calls()
  for (index in rev(seq_along(calls))) {
    reference <- attr(calls[[index]], "srcref")
    if (!is.null(reference)) {
      filename <- tryCatch(
        getSrcFilename(reference, full.names = TRUE),
        error = function(error) NULL
      )
      return(list(
        call = paste(deparse(calls[[index]]), collapse = " "),
        file = if (is.null(filename) || length(filename) == 0L ||
                   !nzchar(filename[[1L]])) NULL else basename(filename[[1L]]),
        line = {
          values <- suppressWarnings(as.integer(reference))
          if (length(values) == 0L || is.na(values[[1L]])) {
            NA_integer_
          } else {
            values[[1L]]
          }
        }
      ))
    }
  }
  list(call = NULL, file = NULL, line = NA_integer_)
}

contains_sentinel <- function(condition) {
  text <- paste(condition_message(condition), condition_call(condition), sep = "\n")
  grepl(forbidden_skip_sentinel, text, fixed = TRUE)
}

classify_observation <- function(condition) {
  message_text <- condition_message(condition)
  if (inherits(condition, "expectation_error")) {
    "error"
  } else if (inherits(condition, "expectation_failure")) {
    "failure"
  } else if (inherits(condition, "expectation_warning")) {
    "warning"
  } else if (inherits(condition, "expectation_skip") &&
             identical(message_text, "Reason: empty test")) {
    "empty"
  } else if (inherits(condition, "expectation_success") &&
             contains_sentinel(condition)) {
    "skip-like-success"
  } else if (inherits(condition, "expectation_skip")) {
    "skip"
  } else {
    "pass"
  }
}

count_observations <- function(classes) {
  list(
    total = as.integer(length(classes)),
    error = as.integer(sum(classes == "error")),
    failure = as.integer(sum(classes == "failure")),
    warning = as.integer(sum(classes == "warning")),
    empty = as.integer(sum(classes == "empty")),
    skip_like_success = as.integer(sum(classes == "skip-like-success")),
    skip = as.integer(sum(classes == "skip")),
    pass = as.integer(sum(classes == "pass"))
  )
}

classify_case <- function(result_list, row) {
  classes <- vapply(result_list, classify_observation, character(1))
  error_indices <- which(classes == "error")
  failure_indices <- which(classes == "failure")
  warning_indices <- which(classes == "warning")
  empty_indices <- which(classes == "empty")
  skip_like_indices <- which(classes == "skip-like-success")
  skip_indices <- which(classes == "skip")

  if (length(error_indices) > 0L || isTRUE(row$error)) {
    list(class = "error", index = if (length(error_indices)) error_indices[[1L]] else 1L)
  } else if (length(failure_indices) > 0L || row$failed > 0L) {
    list(class = "failure", index = if (length(failure_indices)) failure_indices[[1L]] else 1L)
  } else if (length(warning_indices) > 0L || row$warning > 0L) {
    list(class = "warning", index = if (length(warning_indices)) warning_indices[[1L]] else 1L)
  } else if (length(empty_indices) > 0L) {
    list(class = "empty", index = empty_indices[[1L]])
  } else if (length(skip_like_indices) > 0L) {
    list(class = "skip-like-success", index = skip_like_indices[[1L]])
  } else if (length(skip_indices) > 0L || isTRUE(row$skipped)) {
    list(class = "skip", index = if (length(skip_indices)) skip_indices[[1L]] else 1L)
  } else {
    list(class = "pass", index = if (length(result_list)) 1L else NA_integer_)
  }
}

run_tests <- function() {
  old <- setwd(package_path)
  on.exit(setwd(old), add = TRUE)
  test_path <- if (cli[["load-package"]] == "source") {
    "tests/testthat"
  } else {
    system.file("tests/testthat", package = "colocboost")
  }
  if (!nzchar(test_path) || !dir.exists(test_path)) {
    stop(sprintf(
      "Test directory not found for load mode '%s'.",
      cli[["load-package"]]
    ))
  }
  arguments <- list(
    path = test_path,
    package = "colocboost",
    load_package = cli[["load-package"]],
    reporter = testthat::ListReporter$new(),
    stop_on_failure = FALSE,
    stop_on_warning = FALSE
  )
  if (!is.null(cli$filter)) {
    arguments$filter <- cli$filter
  }
  do.call(testthat::test_dir, arguments)
}

captured_run_warnings <- list()
capture_run_warning <- function(condition) {
  source_evidence <- active_source_evidence()
  call <- condition_call(condition)
  file <- condition_file(condition)
  line <- condition_line(condition)
  captured_run_warnings[[length(captured_run_warnings) + 1L]] <<- list(
    message = condition_message(condition),
    call = if (is.null(call)) source_evidence$call else call,
    file = if (is.null(file)) source_evidence$file else file,
    line = if (is.na(line)) source_evidence$line else line
  )
  restart <- findRestart("muffleWarning")
  if (!is.null(restart)) {
    invokeRestart(restart)
  }
}

test_results <- tryCatch(
  withCallingHandlers(run_tests(), warning = capture_run_warning),
  error = identity
)
if (inherits(test_results, "error")) {
  exit_with_diagnostic("execution-error", paste(
    "Unable to execute unit tests:", conditionMessage(test_results)
  ))
}

test_data <- tryCatch(as.data.frame(test_results), error = identity)
if (inherits(test_data, "error")) {
  exit_with_diagnostic("execution-error", paste(
    "Unable to convert unit-test results:", conditionMessage(test_data)
  ))
}

cases <- vector("list", nrow(test_data) + length(captured_run_warnings))
for (index in seq_len(nrow(test_data))) {
  row <- test_data[index, , drop = FALSE]
  result_list <- test_results[[index]]$results
  classification <- classify_case(result_list, row)
  representative <- if (is.na(classification$index) || !length(result_list)) {
    NULL
  } else {
    result_list[[classification$index]]
  }
  message_text <- if (is.null(representative)) "" else condition_message(representative)
  reason <- if (classification$class == "skip") {
    sub("^Reason: ", "", message_text)
  } else {
    NULL
  }
  cases[[index]] <- list(
    class = classification$class,
    message = message_text,
    call = if (is.null(representative)) NULL else condition_call(representative),
    file = basename(row$file[[1L]]),
    line = if (is.null(representative)) NA_integer_ else condition_line(representative),
    test_title = row$test[[1L]],
    test_context = if (is.null(row$context[[1L]]) || !nzchar(row$context[[1L]])) {
      NULL
    } else {
      row$context[[1L]]
    },
    reason = reason,
    observations = count_observations(vapply(
      result_list,
      classify_observation,
      character(1)
    ))
  )
}

if (length(captured_run_warnings) > 0L) {
  for (warning_index in seq_along(captured_run_warnings)) {
    warning <- captured_run_warnings[[warning_index]]
    cases[[nrow(test_data) + warning_index]] <- list(
      class = "warning",
      message = warning$message,
      call = warning$call,
      file = warning$file,
      line = warning$line,
      test_title = "<testthat setup/load>",
      test_context = NULL,
      reason = NULL,
      observations = count_observations("warning")
    )
  }
}

case_classes <- vapply(cases, `[[`, character(1), "class")
summary <- list(
  total = length(cases),
  pass = sum(case_classes == "pass"),
  failure = sum(case_classes == "failure"),
  error = sum(case_classes == "error"),
  warning = sum(case_classes == "warning"),
  empty = sum(case_classes == "empty"),
  skip_like_success = sum(case_classes == "skip-like-success"),
  skip = sum(case_classes == "skip")
)
summary <- lapply(summary, as.integer)
observation_names <- names(count_observations(character()))
summary$observations <- setNames(lapply(observation_names, function(name) {
  as.integer(sum(vapply(cases, function(case) case$observations[[name]], integer(1))))
}), observation_names)

violations <- list()
matched_allowances <- list()
unused_allowances <- list()

append_violation <- function(value) {
  violations[[length(violations) + 1L]] <<- value
}

strict_classes <- c("error", "failure", "warning", "empty", "skip-like-success")
for (case in cases[case_classes %in% strict_classes]) {
  append_violation(list(
    type = case$class,
    file = case$file,
    line = case$line,
    test_title = case$test_title,
    message = case$message
  ))
}

skip_indices <- which(case_classes == "skip")
skip_accounted_for <- rep(FALSE, length(skip_indices))

public_allowance <- function(allowance, observed_count) {
  allowance$.expiry <- NULL
  allowance$observed_count <- as.integer(observed_count)
  allowance
}

if (cli[["load-package"]] == "installed") {
  selected_allowances <- allowances[vapply(
    allowances,
    function(allowance) identical(allowance$context, cli$context),
    logical(1)
  )]
  for (allowance in selected_allowances) {
    matching <- vapply(skip_indices, function(case_index) {
      case <- cases[[case_index]]
      identical(case$file, allowance$file) &&
        identical(case$test_title, allowance$test_title) &&
        identical(case$reason, allowance$reason)
    }, logical(1))
    observed_count <- sum(matching)
    expired <- allowance$.expiry < Sys.Date()
    if (expired) {
      unused_allowances[[length(unused_allowances) + 1L]] <- public_allowance(
        allowance,
        observed_count
      )
      append_violation(list(
        type = "expired-allowance",
        id = allowance$id,
        expires = allowance$expires,
        observed_count = as.integer(observed_count),
        message = sprintf("Skip allowance '%s' expired on %s.", allowance$id, allowance$expires)
      ))
      next
    }
    if (observed_count == 0L) {
      unused_allowances[[length(unused_allowances) + 1L]] <- public_allowance(
        allowance,
        observed_count
      )
      append_violation(list(
        type = "unused-allowance",
        id = allowance$id,
        message = sprintf("Skip allowance '%s' was not observed.", allowance$id)
      ))
    } else if (observed_count != allowance$expected_count) {
      skip_accounted_for[matching] <- TRUE
      append_violation(list(
        type = "skip-count-mismatch",
        id = allowance$id,
        expected_count = allowance$expected_count,
        observed_count = as.integer(observed_count),
        message = sprintf(
          "Skip allowance '%s' expected %d observation(s), but saw %d.",
          allowance$id,
          allowance$expected_count,
          observed_count
        )
      ))
    } else {
      skip_accounted_for[matching] <- TRUE
      matched_allowances[[length(matched_allowances) + 1L]] <- public_allowance(
        allowance,
        observed_count
      )
    }
  }
}

for (position in which(!skip_accounted_for)) {
  case <- cases[[skip_indices[[position]]]]
  append_violation(list(
    type = "unexpected-skip",
    file = case$file,
    line = case$line,
    test_title = case$test_title,
    reason = case$reason,
    message = sprintf(
      "Unexpected skip in %s: %s [%s]",
      case$file,
      case$test_title,
      case$reason
    )
  ))
}

status <- if (length(violations) == 0L) "pass" else "fail"
report <- new_report(
  status = status,
  context = cli$context,
  load_package = cli[["load-package"]],
  package = package_path,
  summary = summary,
  cases = cases,
  matched_allowances = matched_allowances,
  unused_allowances = unused_allowances,
  violations = violations
)

tryCatch(
  write_report(report, output_path),
  error = function(error) {
    message(conditionMessage(error))
    quit(save = "no", status = 2L, runLast = FALSE)
  }
)

quit(
  save = "no",
  status = if (identical(status, "pass")) 0L else 1L,
  runLast = FALSE
)
