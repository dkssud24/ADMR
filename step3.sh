#step3_single
#pQTL_clumping_QC_for_input_MR_tmp
# Load necessary libraries
library(data.table)
library(dplyr)

# Load pQTL clumped data
pqtl <- fread("/BiO//hae/000005_MRSPI/11_pqtl_sig_clump/ZP3_P21754_OID30265_v1_Cardiometabolic_II_sig_sig.txt")
pqtl <- data.frame(pqtl)

# Load clumped SNPs
pqtl_clump <- fread("/BiO/hae/000005_MRSPI/11_pqtl_sig_clump/result_ZP3_P21754_OID30265_v1_Cardiometabolic_II_sig_sig.txt_clumped.tsv")
pqtl_clump_snp <- data.frame(SNP = pqtl_clump$SNP)

# Merge clumped SNPs with pQTL data
m <- left_join(pqtl_clump_snp, pqtl, by = 'SNP')

# Load original pQTL data
pqtl_ori <- fread('/BiO/hae/000005_MRSPI/00_pqtl_sig/ZP3_P21754_OID30265_v1_Cardiometabolic_II_sig.tsv')

# Create SNP_ID for matching
pqtl_ori$SNP_ID <- paste0(pqtl_ori$CHROM, ":", pqtl_ori$GENPOS, ":", pqtl_ori$ALLELE0, ":", pqtl_ori$ALLELE1)
m$SNP_ID <- paste0(m$CHR, ":", m$BP, ":", m$A1, ":", m$A2)  # Ensure matching alleles correctly

# Match datasets using SNP_ID
m_matched <- left_join(m, pqtl_ori, by = "SNP_ID")

# Prepare MR data with necessary columns
mr_data <- m_matched %>%
  select(SNP, CHR, BP, A1 = ALLELE1, A2 = ALLELE0, BETA = BETA.x, SE = SE.x, P = P.x, A1FREQ)

# Check final data
head(mr_data)
dim(mr_data)

#step3_bulk
library(data.table)
library(dplyr)

# Define directories
pqtl_sig_dir <- "/BiO/hae/000005_MRSPI/00_pqtl_sig/"
pqtl_clump_dir <- "/BiO/hae/000005_MRSPI/11_pqtl_sig_clump/"
output_dir <- "/BiO/hae/000005_MRSPI/22_pqtl_sig_clump_QC/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir)

# List all pQTL sig and clump files
pqtl_sig_files <- list.files(pqtl_sig_dir, pattern = "_sig.tsv$", full.names = TRUE)
pqtl_clump_files <- list.files(pqtl_clump_dir, pattern = "_sig_sig.txt$", full.names = TRUE)

# Process each file
for (pqtl_sig_file in pqtl_sig_files) {
  file_base <- tools::file_path_sans_ext(basename(pqtl_sig_file))
  clump_file <- pqtl_clump_files[grepl(file_base, pqtl_clump_files)]
  
  if (length(clump_file) == 0) next
  
  # Load data
  pqtl <- fread(clump_file)
  if (nrow(pqtl) == 0) {
    cat("Skipping empty pqtl file:", clump_file, "\n")
    next
  }
  
  pqtl_clump <- fread(clump_file)
  if (nrow(pqtl_clump) == 0) {
    cat("Skipping empty clump file:", clump_file, "\n")
    next
  }
  
  pqtl_clump_snp <- data.frame(SNP = pqtl_clump$SNP)
  m <- left_join(pqtl_clump_snp, pqtl, by = 'SNP')
  if (nrow(m) == 0) {
    cat("No matching SNPs in:", file_base, "\n")
    next
  }
  
  pqtl_ori <- fread(pqtl_sig_file)
  if (nrow(pqtl_ori) == 0) {
    cat("Skipping empty original pqtl file:", pqtl_sig_file, "\n")
    next
  }

  # Create SNP_ID for matching
  pqtl_ori$SNP_ID <- paste0(pqtl_ori$CHROM, ":", pqtl_ori$GENPOS, ":", pqtl_ori$ALLELE0, ":", pqtl_ori$ALLELE1)
  m$SNP_ID <- paste0(m$CHR, ":", m$BP, ":", m$A1, ":", m$A2)

  # Match datasets
  m_matched <- left_join(m, pqtl_ori, by = "SNP_ID")
  if (nrow(m_matched) == 0) {
    cat("No matched SNPs in:", file_base, "\n")
    next
  }

  # Prepare MR data
  mr_data <- m_matched %>% select(SNP, CHR, BP, A1 = ALLELE1, A2 = ALLELE0, BETA = BETA.x, SE = SE.x, P = P.x, A1FREQ)

  # Save result
  output_file <- paste0(output_dir, file_base, "_MR_ready.tsv")
  fwrite(mr_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Processed:", file_base, "->", output_file, "\n")
}

