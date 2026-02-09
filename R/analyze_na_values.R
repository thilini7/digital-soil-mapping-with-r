# Analyze NA values in raster files to find best reference raster
library(terra)

cov_path <- "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r/Data/data_in/soil_covariates"
cov_files <- list.files(cov_path, pattern = ".tif$", full.names = TRUE)

results <- data.frame(
  file = character(),
  total_cells = numeric(),
  na_cells = numeric(),
  valid_cells = numeric(),
  na_percent = numeric(),
  stringsAsFactors = FALSE
)

cat("Analyzing", length(cov_files), "raster files for NA values...\n\n")

for (f in cov_files) {
  r <- rast(f)
  total <- ncell(r)
  na_count <- global(is.na(r), "sum", na.rm = FALSE)[[1]]
  valid <- total - na_count
  na_pct <- (na_count / total) * 100
  
  results <- rbind(results, data.frame(
    file = basename(f),
    total_cells = total,
    na_cells = na_count,
    valid_cells = valid,
    na_percent = round(na_pct, 2)
  ))
}

# Sort by NA percentage (ascending)
results <- results[order(results$na_percent), ]

cat("=== RASTERS SORTED BY NA PERCENTAGE (LOWEST FIRST) ===\n\n")
print(results, row.names = FALSE)

cat("\n=== BEST REFERENCE RASTER (MINIMUM NA) ===\n")
cat("File:", results$file[1], "\n")
cat("NA Percentage:", results$na_percent[1], "%\n")
cat("Valid Cells:", format(results$valid_cells[1], big.mark = ","), "\n")
cat("\nFull path:\n")
cat(file.path(cov_path, results$file[1]), "\n")
