run_make_report <- function(){ if(requireNamespace("quarto", quietly=TRUE)) quarto::quarto_render("reports/preeclampsia_report.qmd") else message("quarto package unavailable; skipping report") }
