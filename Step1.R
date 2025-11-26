if (!require(archive)) install.packages("archive")
if (!require(data.table)) install.packages("data.table")
if (!require(dplyr)) install.packages("dplyr")

library(archive)
library(data.table)
library(dplyr)

pqtl_path <- "/BiO/hae/000001_UKB_pQTL/"
saved_path <- "/BiO/hae/000005_MRSPI/00_pqtl_sig/"

if (!dir.exists(saved_path)) dir.create(saved_path)

pqtl_files <- list.files(pqtl_path, pattern = "\\.tar$", full.names = TRUE)

pval_threshold <- 1.100015e-12

for (tar_file in pqtl_files) {
  file_name <- tools::file_path_sans_ext(basename(tar_file))
  archive_contents <- archive::archive(tar_file)
  chromosome_files <- archive_contents$path[grepl("chr[0-9]+|chr1[0-9]|chr2[0-2]", archive_contents$path)]
  
  chromosome_data <- lapply(chromosome_files, function(file_path) {
    archive_read <- archive::archive_read(tar_file, file = file_path)
    con <- gzcon(archive_read)
    raw_text <- readLines(con, warn = FALSE)
    text_input <- paste(raw_text, collapse = "\n")
    data <- fread(text_input, na.strings = c("", "NA"))
    close(archive_read)
    return(data)
  })
  
  combined_data <- rbindlist(chromosome_data, fill = TRUE)
  combined_data$MAF <- pmin(combined_data$A1FREQ, 1 - combined_data$A1FREQ)
  combined_data <- combined_data[combined_data$MAF > 0.01, ]
  combined_data$P <- 10^(-combined_data$LOG10P)
  filtered_data <- combined_data[combined_data$P < pval_threshold, ]
  
  if (nrow(filtered_data) > 0) {
    save_file <- paste0(saved_path, file_name, "_sig.tsv")
    fwrite(filtered_data, save_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Saved:", save_file, "(", nrow(filtered_data), "rows )\n")
  } else {
    cat("No significant data:", file_name, "\n")
  }
}
