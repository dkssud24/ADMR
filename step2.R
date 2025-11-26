#Step2
#Clumping
library(data.table)
library(dplyr)

# Load 1000G reference BIM file once
root_bim <- fread("/BiO/hae/000006_ref_1000G/all_bim.txt.v2")

# Set paths
pqtl_path <- "/BiO/hae/000005_MRSPI/00_pqtl_sig/"
output_path <- "/BiO/hae/000005_MRSPI/11_pqtl_sig_clump/"
plink_path <- "plink"
ref_data <- "/BiO/hae/000006_ref_1000G/ref"

# Create output directory if not exists
if (!dir.exists(output_path)) dir.create(output_path)

# List all pqtl files
pqtl_files <- list.files(pqtl_path, pattern = "\\.tsv$", full.names = TRUE)

# Loop through each pqtl file
for (pqtl_file in pqtl_files) {
  pqtl_sig <- fread(pqtl_file)
  pqtl_sig$original <- sub("(:[^:]+){2}$", "", pqtl_sig$ID)
  
  merged_data <- left_join(pqtl_sig, root_bim, by = "original")
  clumping_input <- merged_data[, c("CHROM", "V2", "GENPOS", "ALLELE0", "ALLELE1", "P", "BETA", "SE")]
  colnames(clumping_input) <- c("CHR", "SNP", "BP", "A1", "A2", "P", "BETA", "SE")
  clumping_input <- na.omit(clumping_input)
  
  file_base <- tools::file_path_sans_ext(basename(pqtl_file))
  plink_input_path <- file.path(output_path, paste0(file_base, "_sig.txt"))
  result_output_path <- file.path(output_path, paste0("result_", file_base, "_sig.txt"))

  fwrite(clumping_input, plink_input_path, sep = "\t", quote = FALSE, row.names = FALSE)
  
  plink_command <- paste(plink_path,
                         "--bfile", ref_data,
                         "--clump", plink_input_path,
                         "--clump-kb 1000",
                         "--clump-r2 0.01",
                         "--clump-p1 1",
                         "--clump-p2 1",
                         "--out", result_output_path)

  system(plink_command)

  clumped_file <- paste0(result_output_path, ".clumped")
  
  if (file.exists(clumped_file)) {
    clumped_result <- fread(clumped_file)
    fwrite(clumped_result, paste0(result_output_path, "_clumped.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Clumping completed for:", pqtl_file, "\n")
  } else {
    cat("No significant clumping results for:", pqtl_file, "\n")
  }
}

