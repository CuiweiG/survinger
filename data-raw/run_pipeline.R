.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all()

ecdc <- readRDS("data-raw/ecdc_surveillance.rds")
cat("=== Running full pipeline on ECDC real-world data ===\n\n")

design <- surv_design(
  data = ecdc$sequences,
  strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population,
  source_type = "source_type"
)
cat("--- Design ---\n")
print(design)

cat("\n--- System Report ---\n")
report <- surv_report(design)

lineages <- sort(table(ecdc$sequences$lineage), decreasing = TRUE)
target_lin <- names(lineages)[!grepl("Other", names(lineages))][1]
cat("\n--- Prevalence for", target_lin, "---\n")
weighted <- surv_lineage_prevalence(design, target_lin, method = "hajek")
naive <- surv_naive_prevalence(design, target_lin)
cat("Weighted mean:", round(mean(weighted$estimates$prevalence, na.rm = TRUE), 4), "\n")
cat("Naive mean:", round(mean(naive$estimates$prevalence, na.rm = TRUE), 4), "\n")
cat("Mean |diff|:", round(mean(abs(
  weighted$estimates$prevalence - naive$estimates$prevalence
), na.rm = TRUE), 4), "\n")

cat("\n--- Allocation ---\n")
alloc <- surv_optimize_allocation(design, "min_mse", total_capacity = 1000)
print(alloc)

cat("\n--- Delay ---\n")
delay <- surv_estimate_delay(design)
print(delay)

cat("\n--- Nowcast ---\n")
nc <- surv_nowcast_lineage(design, delay, target_lin)
print(nc)

cat("\n--- Adjusted ---\n")
adj <- surv_adjusted_prevalence(design, delay, target_lin)
print(adj)

cat("\n--- Detection ---\n")
det <- surv_detection_probability(design, 0.01)
cat("Weekly detection at 1%:", round(det$overall, 3), "\n")

# Generate plots
png_dir <- "vignettes/figures"
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

grDevices::png(file.path(png_dir, "ecdc_compare.png"), width = 800, height = 500)
print(surv_compare_estimates(weighted, naive,
  title = paste("ECDC Real Data:", target_lin, "- Weighted vs Naive")))
grDevices::dev.off()

grDevices::png(file.path(png_dir, "ecdc_rates.png"), width = 800, height = 400)
print(surv_plot_sequencing_rates(design))
grDevices::dev.off()

grDevices::png(file.path(png_dir, "ecdc_allocation.png"), width = 800, height = 400)
print(plot(alloc))
grDevices::dev.off()

grDevices::png(file.path(png_dir, "ecdc_delay.png"), width = 800, height = 400)
print(plot(delay))
grDevices::dev.off()

grDevices::png(file.path(png_dir, "ecdc_nowcast.png"), width = 800, height = 500)
print(plot(nc))
grDevices::dev.off()

grDevices::png(file.path(png_dir, "ecdc_adjusted.png"), width = 800, height = 500)
print(plot(adj))
grDevices::dev.off()

cat("\nPlots saved:", paste(list.files(png_dir), collapse = ", "), "\n")
cat("\n=== FULL PIPELINE ON REAL DATA: SUCCESS ===\n")
