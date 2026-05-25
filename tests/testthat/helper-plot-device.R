.colocboost_plot_device_env <- new.env(parent = emptyenv())
.colocboost_plot_device_env$file <- tempfile("colocboost-test-plots-", fileext = ".pdf")
grDevices::pdf(.colocboost_plot_device_env$file)

reg.finalizer(.colocboost_plot_device_env, function(e) {
  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  unlink(e$file)
  unlink("Rplots.pdf")
  unlink(file.path("tests", "testthat", "Rplots.pdf"))
}, onexit = TRUE)
