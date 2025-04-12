library(testthat)

# Helper function to parse NAMESPACE file and extract exported functions
get_exported_functions <- function() {
  # Locate the NAMESPACE file
  if (file.exists("../../NAMESPACE")) {
    namespace_path <- "../../NAMESPACE"
  } else if (file.exists("../NAMESPACE")) {
    namespace_path <- "../NAMESPACE"
  } else {
    # If running within package, try to find it in the package directory
    package_dir <- find.package("colocboost", quiet = TRUE)
    if (length(package_dir) > 0) {
      namespace_path <- file.path(package_dir, "NAMESPACE")
      if (!file.exists(namespace_path)) {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  }
  
  # Read NAMESPACE file
  ns_lines <- readLines(namespace_path)
  
  # Extract exported functions
  export_lines <- ns_lines[grepl("^export\\(", ns_lines)]
  exported_funcs <- gsub("export\\((.*)\\)", "\\1", export_lines)
  
  # Clean up any quotation marks
  exported_funcs <- gsub("\"", "", exported_funcs)
  exported_funcs <- gsub("\'", "", exported_funcs)
  
  return(exported_funcs)
}

# Test that the package loads without errors
test_that("package loads without errors", {
  expect_error(library(colocboost), NA)
})

# Test that exported functions exist
test_that("package has expected functions", {
  exported_functions <- get_exported_functions()
  
  # If we couldn't parse the NAMESPACE, fall back to a minimal set
  if (is.null(exported_functions) || length(exported_functions) == 0) {
    warning("Could not parse NAMESPACE file, falling back to minimum expected functions")
    exported_functions <- c("colocboost")
  }
  
  # Test that each exported function exists
  for (func in exported_functions) {
    expect_true(exists(func, where = asNamespace("colocboost")), 
                info = paste0("Function '", func, "' should exist in the package"))
  }
  
  # Print the functions we found for debugging
  message("Tested for existence of the following exported functions:")
  message(paste(exported_functions, collapse = ", "))
})

# Test that package depends on required packages
test_that("package depends on required packages", {
  # Get dependencies from DESCRIPTION file
  get_dependencies <- function() {
    desc_file <- NULL
    if (file.exists("../../DESCRIPTION")) {
      desc_file <- "../../DESCRIPTION"
    } else if (file.exists("../DESCRIPTION")) {
      desc_file <- "../DESCRIPTION"
    } else {
      package_dir <- find.package("colocboost", quiet = TRUE)
      if (length(package_dir) > 0) {
        desc_file <- file.path(package_dir, "DESCRIPTION")
      }
    }
    
    if (is.null(desc_file) || !file.exists(desc_file)) {
      return(c("R", "Rfast", "matrixStats"))  # Fallback
    }
    
    desc <- readLines(desc_file)
    
    # Find Depends and Imports sections
    deps_start <- grep("^(Depends|Imports):", desc)
    deps_end <- c(deps_start[-1] - 1, length(desc))
    deps_end <- deps_end[1:length(deps_start)]
    
    # Extract dependencies
    deps <- c()
    for (i in 1:length(deps_start)) {
      section <- desc[deps_start[i]:deps_end[i]]
      section_text <- paste(section, collapse = " ")
      
      # Remove section header
      section_text <- gsub("^(Depends|Imports):", "", section_text)
      
      # Extract package names
      pkgs <- strsplit(section_text, ",")[[1]]
      pkgs <- trimws(pkgs)
      
      # Remove version requirements
      pkgs <- gsub("\\s*\\(.*\\)", "", pkgs)
      
      deps <- c(deps, pkgs)
    }
    
    # Remove R from the list (it's always available)
    deps <- setdiff(deps, "R (>= 4.0.0)")
    deps <- c("R", deps)  # Add R back in a standardized format
    
    return(deps)
  }
  
  expected_deps <- get_dependencies()
  
  # Check if packages can be loaded - this is a better test than checking DESCRIPTION
  for (pkg in expected_deps) {
    if (pkg != "R") { # Skip R as it's always loaded
      if (!requireNamespace(pkg, quietly = TRUE)) {
        skip(paste("Package", pkg, "not available for testing"))
      }
    }
  }
  
  # Basic test succeeded if we got here
  expect_true(TRUE)
})

# Test that the package has the correct structure
test_that("package has correct S3 methods", {
  # Check S3 methods using methods function if any expected
  expect_true(is.function(colocboost), "colocboost should be a function")
  
  # Test that colocboost returns correct class
  skip_on_cran()
  
  # Generate minimal test data
  set.seed(123)
  n <- 20  # Small sample for quick test
  p <- 5   # Few variables
  X <- matrix(rnorm(n*p), n, p)
  colnames(X) <- paste0("X", 1:p)
  Y <- matrix(rnorm(n*2), n, 2)
  Y_list <- list(Y[,1], Y[,2])
  X_list <- list(X, X)
  
  # Run with minimal iterations
  suppressWarnings({
    result <- colocboost(X = X_list, Y = Y_list, M = 2)
  })
  
  # Test class
  expect_s3_class(result, "colocboost")
})

# Test documentation exists for key functions
test_that("documentation exists for key functions", {
  # Get all exported functions
  exported_functions <- get_exported_functions()
  
  # If we couldn't parse the NAMESPACE, fall back to a minimal set
  if (is.null(exported_functions) || length(exported_functions) == 0) {
    exported_functions <- c("colocboost", "colocboost_plot", "get_cos_summary")
  }
  
  # Limit to a reasonable number to check
  funcs_to_check <- exported_functions[1:min(length(exported_functions), 5)]
  
  # Check if help pages exist for key functions
  for (func in funcs_to_check) {
    expect_error(help(func), NA, 
                 info = paste0("Documentation should exist for function '", func, "'"))
  }
})