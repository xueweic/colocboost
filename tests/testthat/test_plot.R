library(testthat)

# Utility function to generate a simple colocboost results 
generate_test_result <- function(n = 500, p = 60, L = 4, seed = 42, output_level = 3) {
  set.seed(seed)
  
  # Generate X with LD structure
  sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] <- 0.9^abs(i - j)
    }
  }
  X <- MASS::mvrnorm(n, rep(0, p), sigma)
  colnames(X) <- paste0("SNP", 1:p)
  
  # Generate true effects based on the number of traits
  true_beta <- matrix(0, p, L)
  
  if (L == 1) {
    # Single trait case
    true_beta[5, 1] <- 0.5  # SNP5 affects the trait
    true_beta[30, 1] <- 0.3 # SNP10 also affects the trait
  } else if (L == 2) {
    # Simple multi-trait case
    true_beta[5, 1] <- 0.5  # SNP5 affects trait 1
    true_beta[5, 2] <- 0.4  # SNP5 also affects trait 2 (colocalized)
    true_beta[30, 2] <- 0.3 # SNP10 only affects trait 2
  } else if (L == 4) {
    # Complex multi-trait case with multiple colocalization patterns
    # SNP5 affects traits 1, 2, and 3 (colocalized across 3 traits)
    true_beta[5, 1] <- 0.5  
    true_beta[5, 2] <- 0.5
    true_beta[5, 3] <- 0.5  
    
    # SNP10 affects traits 2 and 4 (colocalized across 2 traits)
    true_beta[30, 2] <- 0.5
    true_beta[30, 4] <- 0.5
    
    # SNP15 only affects trait 3 (trait-specific effect)
    true_beta[40, 3] <- 0.6
    
    # SNP18 only affects trait 4 (trait-specific effect)
    true_beta[55, 4] <- 0.5
  }
  
  # Generate Y with some noise
  Y <- matrix(0, n, L)
  for (l in 1:L) {
    Y[, l] <- X %*% true_beta[, l] + rnorm(n, 0, 1)
  }
  
  # Prepare input for colocboost
  if (L == 1) {
    # For single trait, Y should be a vector
    Y_input <- Y[,1]
    X_input <- X
  } else {
    # For multiple traits, convert to list format
    Y_input <- lapply(1:L, function(l) Y[,l])
    X_input <- replicate(L, X, simplify = FALSE)
  }
    
  
  # Run colocboost with minimal parameters to get a model object
  suppressWarnings({
    result <- colocboost(
      X = X_input, 
      Y = Y_input,
      focal_outcome_idx = L,
      output_level = output_level
    )
  })
  
  return(result)
}

# Generate a test colocboost result
cb_res <- generate_test_result()

# Test colocboost_plot function with basic options
test_that("colocboost_plot basic functionality works", {
  
  # Basic plot call
  expect_error(suppressWarnings(colocboost_plot(cb_res)), NA)
  
  # Test with non-colocboost object
  expect_error(colocboost_plot("not_a_colocboost_object"), 
               "Input of colocboost_plot must be a 'colocboost' object!")
})

# Test colocboost_plot with different y-axis options
test_that("colocboost_plot handles different y-axis options", {
  
  # Test with different y-axis values
  expect_error(suppressWarnings(colocboost_plot(cb_res, y = "log10p")), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, y = "z_original")), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, y = "cos_vcp")), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, y = "vcp")), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, y = "coef")), NA)
  
  # Test with invalid y-axis value
  expect_error(suppressWarnings(colocboost_plot(cb_res, y = "invalid")), 
               "Invalid y value! Choose from 'log10p', 'z_original', 'vcp', 'coef', or 'cos_vcp'!")
})

# Test colocboost_plot with plot filtering options
test_that("colocboost_plot handles filtering options", {
  
  # Test with outcome index filtering
  expect_error(suppressWarnings(colocboost_plot(cb_res, outcome_idx = 1)), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, outcome_idx = 2)), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, outcome_idx = 1:2)), NA)
  
  # Test with plot_cos_idx filtering
  cos_count <- length(cb_res$cos_details$cos$cos_variables)
  if (cos_count > 0) {
    expect_error(suppressWarnings(colocboost_plot(cb_res, plot_cos_idx = 1)), NA)
    
    # Test with invalid plot_cos_idx
    if (cos_count < 10) {
      expect_error(suppressWarnings(colocboost_plot(cb_res, plot_cos_idx = 10)),
                   "Please check plot_cos_idx!")
    }
  }
  
  # Test with plot_focal_only
  expect_error(suppressWarnings(colocboost_plot(cb_res, plot_focal_only = TRUE)), NA)
  
  # Test with plot_focal_cos_outcome_only
  expect_error(suppressWarnings(colocboost_plot(cb_res, plot_focal_cos_outcome_only = TRUE)), NA)
})

# Test colocboost_plot with visual customization options
test_that("colocboost_plot handles visual customization options", {
  
  # Test with custom colors
  expect_error(suppressWarnings(colocboost_plot(cb_res, points_color = "red")), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, cos_color = c("blue", "green", "orange"))), NA)
  
  # Test with custom styling
  expect_error(suppressWarnings(colocboost_plot(cb_res, lab_style = c(3, 2))), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, axis_style = c(2.5, 2))), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, title_style = c(3, 3))), NA)
  
  # Test with legend position options
  expect_error(suppressWarnings(colocboost_plot(cb_res, outcome_legend_pos = "top")), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, outcome_legend_pos = "bottom")), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, outcome_legend_pos = "left")), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, outcome_legend_pos = "right")), NA)
  
  # Test with custom legend sizes
  expect_error(suppressWarnings(colocboost_plot(cb_res, outcome_legend_size = 2.0)), NA)
})

# Test colocboost_plot with layout options
test_that("colocboost_plot handles layout options", {
  
  # Test with different plot_cols values
  expect_error(suppressWarnings(colocboost_plot(cb_res, plot_cols = 1)), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, plot_cols = 3)), NA)
  
  # Test with ylim_each option
  expect_error(suppressWarnings(colocboost_plot(cb_res, ylim_each = TRUE)), NA)
  # When ylim_each is FALSE, we need to provide a ylim parameter
  expect_error(suppressWarnings(colocboost_plot(cb_res, ylim_each = FALSE, ylim = c(0, 10))), NA)
  
  # Test with title_specific option
  expect_error(suppressWarnings(colocboost_plot(cb_res, title_specific = "BRCA1")), NA)
  
  # Test with variant_coord option
  expect_error(suppressWarnings(colocboost_plot(cb_res, variant_coord = FALSE)), NA)
})

# Test colocboost_plot with additional visualization options
test_that("colocboost_plot handles additional visualization options", {
  
  # Test with vertical line options
  expect_error(suppressWarnings(colocboost_plot(cb_res, add_vertical = TRUE, add_vertical_idx = c(5, 10))), NA)
  
  # Test with show_top_variables option
  expect_error(suppressWarnings(colocboost_plot(cb_res, show_top_variables = TRUE)), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res, show_top_variables = FALSE)), NA)
  
  # Test with show_variable option
  expect_error(suppressWarnings(colocboost_plot(cb_res, show_variable = TRUE)), NA)
})

# Test colocboost_plot with custom outcome names
test_that("colocboost_plot handles custom outcome names", {
  
  # Test with custom outcome names
  expect_error(suppressWarnings(colocboost_plot(cb_res, outcome_names = c("Trait1", "Trait2"))), NA)
})

# Test colocboost_plot with a specific range
test_that("colocboost_plot handles zoom-in with grange", {
  
  # Test with grange option to zoom in
  expect_error(suppressWarnings(colocboost_plot(cb_res, grange = 5:15)), NA)
})

# Test colocboost_plot with focal outcome in L=4 case
test_that("colocboost_plot handles focal outcome in complex cases", {
  
  # Generate a test colocboost result with 4 traits and focal outcome set
  cb_res_focal <- generate_test_result(L = 4, output_level = 3)
  
  # Basic plot call with focal outcome
  expect_error(suppressWarnings(colocboost_plot(cb_res_focal)), NA)
  
  # Test plot_focal_only option
  expect_error(suppressWarnings(colocboost_plot(cb_res_focal, plot_focal_only = TRUE)), NA)
  
  # Test plot_focal_cos_outcome_only option
  expect_error(suppressWarnings(colocboost_plot(cb_res_focal, plot_focal_cos_outcome_only = TRUE)), NA)
  
  # Combine focal outcome filtering with other options
  expect_error(suppressWarnings(colocboost_plot(cb_res_focal, 
                                               plot_focal_only = TRUE,
                                               y = "cos_vcp")), NA)
  
  expect_error(suppressWarnings(colocboost_plot(cb_res_focal, 
                                               plot_focal_cos_outcome_only = TRUE,
                                               plot_ucos = TRUE)), NA)
  
  # Test focusing only on outcomes colocalized with focal outcome
  expect_error(suppressWarnings(colocboost_plot(cb_res_focal, 
                                               plot_focal_cos_outcome_only = TRUE,
                                               outcome_idx = 1:3)), NA)
})


# Test colocboost_plot with single trait (finemapping) results
test_that("colocboost_plot handles single trait results", {
  
  # Generate a single trait colocboost result
  cb_res_single <- generate_test_result(L = 1)
  
  # Basic plot call for single trait
  expect_error(suppressWarnings(colocboost_plot(cb_res_single)), NA)
  
  # Test custom options with single trait
  expect_error(suppressWarnings(colocboost_plot(cb_res_single, y = "vcp")), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res_single, plot_cols = 1)), NA)
})

# Test colocboost_plot with L=4 case for complex colocalization and trait-specific effects
test_that("colocboost_plot handles L=4 case with complex colocalization patterns", {
  
  # Generate a test colocboost result with 4 traits and high output level
  cb_res_complex <- generate_test_result(L = 4, output_level = 3)
  
  # Basic plot call for complex case
  expect_error(suppressWarnings(colocboost_plot(cb_res_complex)), NA)
  
  # Test y-axis options for complex case
  expect_error(suppressWarnings(colocboost_plot(cb_res_complex, y = "log10p")), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res_complex, y = "cos_vcp")), NA)
  
  # Test filtering for specific outcomes
  expect_error(suppressWarnings(colocboost_plot(cb_res_complex, outcome_idx = 1:2)), NA)
  expect_error(suppressWarnings(colocboost_plot(cb_res_complex, outcome_idx = c(1,3))), NA)
  
  # Test plot_ucos for visualizing trait-specific effects
  expect_error(suppressWarnings(colocboost_plot(cb_res_complex, plot_ucos = TRUE)), NA)
  
  # Test with specific plot_ucos_idx if available
  # This is conditional because it depends on the actual number of ucos detected
  if (!is.null(cb_res_complex$ucos_details) && 
      !is.null(cb_res_complex$ucos_details$ucos) && 
      length(cb_res_complex$ucos_details$ucos$ucos_index) > 0) {
    n_ucos <- length(cb_res_complex$ucos_details$ucos$ucos_index)
    if (n_ucos > 0) {
      expect_error(suppressWarnings(colocboost_plot(cb_res_complex, 
                                                   plot_ucos = TRUE, 
                                                   plot_ucos_idx = 1:min(n_ucos, 2))), NA)
    }
  }
  
  # Test visualization of both colocalization and trait-specific effects together
  expect_error(suppressWarnings(colocboost_plot(cb_res_complex, 
                                               plot_ucos = TRUE, 
                                               show_cos_to_uncoloc = TRUE)), NA)
})


# Generate a test colocboost result with high output level to include ucos details
cb_res <- generate_test_result(output_level = 3)

# Test colocboost_plot with uncolocalized visualization options
test_that("colocboost_plot handles uncolocalized visualization options", {
  
  
  # Test with plot_ucos options
  expect_error(suppressWarnings(colocboost_plot(cb_res, plot_ucos = TRUE)), NA)
  
  # Test with show_cos_to_uncoloc options
  expect_error(suppressWarnings(colocboost_plot(cb_res, show_cos_to_uncoloc = TRUE)), NA)
  
  # Generate a different colocboost result to test the warning for plot_ucos
  cb_res_low <- generate_test_result(output_level = 1)
  # This should give a warning but not an error
  expect_warning(colocboost_plot(cb_res_low, plot_ucos = TRUE),
                 "Since you want to plot trait-specific \\(uncolocalized\\) sets with plot_ucos = TRUE")
})

# Test colocboost_plot with varying cutoff settings from get_robust_colocalization
test_that("colocboost_plot handles varying cutoff settings", {
  
  # Generate a test colocboost result
  cb_res <- generate_test_result(output_level = 3)
  
  # Test with different cutoff settings
  cutoff_settings <- list(
    list(cos_npc_cutoff = 0.5, npc_outcome_cutoff = 0.2),
    list(cos_npc_cutoff = 0.7, npc_outcome_cutoff = 0.3),
    list(cos_npc_cutoff = 0.9, npc_outcome_cutoff = 0.5),
    list(cos_npc_cutoff = 1.0, npc_outcome_cutoff = 0.5),  # Corner case: cos_npc_cutoff = 1
    list(cos_npc_cutoff = 0.5, npc_outcome_cutoff = 1.0),  # Corner case: npc_outcome_cutoff = 1
    list(cos_npc_cutoff = 1.0, npc_outcome_cutoff = 1.0)   # Corner case: both cutoffs = 1
  )
  
  for (cutoff in cutoff_settings) {
    # Apply robust colocalization filtering
    filter_res <- get_robust_colocalization(
      cb_res, 
      cos_npc_cutoff = cutoff$cos_npc_cutoff, 
      npc_outcome_cutoff = cutoff$npc_outcome_cutoff
    )
    
    # Test basic plot call with filtered results
    expect_error(suppressWarnings(colocboost_plot(filter_res)), NA)
    
    # Test y-axis options with filtered results
    expect_error(suppressWarnings(colocboost_plot(filter_res, y = "log10p")), NA)
    expect_error(suppressWarnings(colocboost_plot(filter_res, y = "cos_vcp")), NA)
    
    # Test plot_focal_only option
    expect_error(suppressWarnings(colocboost_plot(filter_res, plot_focal_only = TRUE)), NA)
    
    # Test plot_ucos option
    expect_error(suppressWarnings(colocboost_plot(filter_res, plot_ucos = TRUE)), NA)
    
    # Test show_cos_to_uncoloc option
    expect_error(suppressWarnings(colocboost_plot(filter_res, show_cos_to_uncoloc = TRUE)), NA)
    
    # Test combined options
    expect_error(suppressWarnings(colocboost_plot(filter_res, 
                                                  plot_focal_only = TRUE, 
                                                  plot_ucos = TRUE, 
                                                  y = "cos_vcp")), NA)
  }
})


