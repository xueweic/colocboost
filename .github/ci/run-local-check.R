#!/usr/bin/env Rscript

tex_root <- tryCatch(tinytex::tinytex_root(), error = function(error) "")
if (nzchar(tex_root)) {
  tex_bins <- list.dirs(file.path(tex_root, "bin"), recursive = FALSE)
  complete_bins <- Filter(function(tex_bin) {
    file.exists(file.path(tex_bin, "pdflatex")) &&
      file.exists(file.path(tex_bin, "makeindex"))
  }, tex_bins)
  if (length(complete_bins)) {
    Sys.setenv(PATH = paste(complete_bins[[1L]], Sys.getenv("PATH"), sep = .Platform$path.sep))
  }
}
missing_tools <- names(which(!nzchar(Sys.which(c(pdflatex = "pdflatex", makeindex = "makeindex")))))
if (length(missing_tools)) {
  stop("A complete TeX installation is required for the CRAN manual check: ",
       paste(missing_tools, collapse = ", "))
}

devtools::check(
  pkg = ".",
  args = "--as-cran",
  build_args = character(),
  manual = TRUE
)
