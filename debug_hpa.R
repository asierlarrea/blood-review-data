sink("debug_output.txt")
setwd("C:/Users/HP zBook 15v/Documents/ebi-work/blood-review-data")

cat("=== DEBUGGING HPA LOADING ===\n\n")

source("scripts/config/analysis_config.R")
source("scripts/utilities/data_loader.R")

cat("Available sources:", paste(get_all_data_sources(), collapse=", "), "\n\n")

cat("Testing HPA MS...\n")
result_ms <- tryCatch({
  hpa_ms <- load_data_source("hpa_ms", force_mapping = FALSE)
  cat("SUCCESS: HPA MS loaded with", nrow(hpa_ms), "rows\n")
  cat("Sample abundance values:", paste(head(hpa_ms$abundance, 5), collapse=", "), "\n")
  hpa_ms
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  NULL
})

cat("\nTesting HPA PEA...\n")
result_pea <- tryCatch({
  hpa_pea <- load_data_source("hpa_pea", force_mapping = FALSE)
  cat("SUCCESS: HPA PEA loaded with", nrow(hpa_pea), "rows\n")
  cat("Sample abundance values:", paste(head(hpa_pea$abundance, 5), collapse=", "), "\n")
  hpa_pea
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  NULL
})

sink()



























