make_coordinate_ucos_result <- function() {
  variables <- c(
    "chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G", "chr1:400:A:G"
  )
  ucos_index <- list("ucos1:y2" = c(2L, 3L))

  result <- list(
    data_info = list(
      n_outcomes = 2L,
      n_variables = length(variables),
      outcome_info = list(
        outcome_names = c("BP", "MDD"),
        is_focal = c(FALSE, FALSE)
      ),
      variables = variables,
      z = list(BP = c(0.1, 0.2, 0.3, 0.4), MDD = c(0.4, 0.3, 0.2, 0.1)),
      coef = list(BP = rep(0, length(variables)), MDD = rep(0, length(variables)))
    ),
    vcp = rep(0, length(variables)),
    cos_summary = data.frame(cos_id = character(), focal_outcome = logical()),
    cos_details = list(
      cos = list(
        cos_variables = NULL,
        cos_index = NULL,
        cos_top_variables = data.frame(top_index = integer(), top_variables = character()),
        cos_outcomes = list(outcome_index = list()),
        cos_vcp = list()
      )
    ),
    ucos_details = list(
      ucos = list(
        ucos_index = ucos_index,
        ucos_variables = list("ucos1:y2" = variables[c(2L, 3L)])
      ),
      ucos_outcomes = list(
        outcome_index = list("ucos1:y2" = 2L),
        outcome_name = list("ucos1:y2" = "MDD")
      ),
      ucos_weight = list("ucos1:y2" = c(0, 0.6, 0.4, 0)),
      ucos_top_variables = data.frame(
        top_index = 2L,
        top_variables = variables[2L],
        row.names = "ucos1:y2"
      )
    )
  )
  class(result) <- "colocboost"
  result
}

test_that("coordinate plots map appended uCoS indices to genomic positions", {
  plot_input <- suppressWarnings(colocboost:::get_input_plot(
    make_coordinate_ucos_result(),
    plot_ucos = TRUE,
    variant_coord = TRUE
  ))

  expect_equal(plot_input$cos[["ucos1:y2"]], c(200, 300))
  expect_equal(plot_input$cos_hits[["ucos1:y2"]], 200)
  expect_true(all(plot_input$cos[["ucos1:y2"]] %in% plot_input$x$pos))
})
