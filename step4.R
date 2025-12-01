# Load required libraries
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(MR.SPI)

# Define file paths
pqtl_dir <- "/BiO/hae/000005_MRSPI/33_pqtl_sig_clump_QC_list/GCST002245/"
gwas_files <- "/BiO/hae/000005_MRSPI/99_System_AD/00_GWAS_data/GCST002245"
output_base_dir <- "/BiO/hae/000005_MRSPI/99_System_AD/04_SNP_Pro_AD/GCST002245/"

if (!dir.exists(output_base_dir)) dir.create(output_base_dir)

# Read list of GWAS files
#gwas_files <- readLines(gwas_list_file)

# List all pQTL files
pqtl_files <- list.files(pqtl_dir, pattern = "\\.tsv$", full.names = TRUE)

# Iterate over each GWAS file
for (gwas_file in gwas_files) {
  
  # Load GWAS outcome data
  cat("Loading GWAS file", gwas_file, "\n")
  outcome_data <- fread(gwas_file)
#  colnames(outcome_data) <- c("CHR", "SNP", "BP", "other_allele.outcome", "effect_allele.outcome",
#                              "beta.outcome", "se.outcome", "neg_log10_p_value", "odds_ratio", "ci_lower",
#                              "ci_upper", "eaf.outcome", "pval.outcome")
  outcome_data$outcome <- "Brain MRI Volume"
  outcome_data$id.outcome <- outcome_data$SNP
  outcome_data <- data.frame(outcome_data)
  
  # Extract GWAS ID for folder naming (e.g., GCST90002621)
  gwas_id <- gsub(".*(GCST[0-9]+).*", "\\1", gwas_file)
  output_dir <- file.path(output_base_dir)
  
  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  cat("Created directory", output_dir, "for results\n")
  
  # Iterate over each pQTL file for current GWAS
  for (pqtl_file in pqtl_files) {
    
    cat("Processing pQTL file", pqtl_file, "with GWAS", gwas_file, "\n")
    
    # Load pQTL data
    pqtl_data <- fread(pqtl_file)
    colnames(pqtl_data) <- c("SNP", "CHR", "BP", "effect_allele.exposure", "other_allele.exposure",
                             "beta.exposure", "se.exposure", "pval.exposure", "eaf.exposure")
    pqtl_data$exposure <- "Protein"
    pqtl_data$id.exposure <- pqtl_data$SNP
    pqtl_data <- data.frame(pqtl_data)
    
    # Harmonize data
    dat <- harmonise_data(pqtl_data, outcome_data, action = 1)
    
    # Prepare MR-SPI input and count SNPs
    mr_spi_input <- dat %>%
      dplyr::select(SNP, beta.exposure, se.exposure, beta.outcome, se.outcome, eaf.exposure) %>%
      na.omit()
    
    # Count harmonized SNPs
    snp_count <- nrow(mr_spi_input)
    
    # Skip if less than 3 SNPs
    if (snp_count < 3) {
      cat("Insufficient SNPs less than 3 for MR-SPI in", pqtl_file, "with GWAS", gwas_file, "- Skipping.\n")
      next
    }
    
    # Perform MR-SPI analysis
    n1 <- 33995  # GWAS sample size
    n2 <- 46668  # pQTL sample size
    
    mr_spi_result <- MR.SPI(
      gamma = mr_spi_input$beta.exposure,
      Gamma = mr_spi_input$beta.outcome,
      se_gamma = mr_spi_input$se.exposure,
      se_Gamma = mr_spi_input$se.outcome,
      n1 = n1,
      n2 = n2,
      freq = mr_spi_input$eaf.exposure,
      max_clique = TRUE,
      verbose = TRUE
    )
    
    # Extract MR-SPI results
    beta <- mr_spi_result$betaHat[[1]]
    se <- mr_spi_result$beta.sdHat[[1]]
    ci_lower <- mr_spi_result$ci[[1]][1]
    ci_upper <- mr_spi_result$ci[[1]][2]
    z_score <- beta / se
    p_value <- 2 * pnorm(-abs(z_score))
    
    # Create result dataframe with SNP count
    result_df <- data.frame(
      GWAS_File = basename(gwas_file),
      PQTL_File = basename(pqtl_file),
      SNP_Count = snp_count,
      Beta = beta,
      SE = se,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      P_Value = p_value
    )
    
    # Define output file name
    gwas_name <- gsub(".tsv", "", basename(gwas_file))
    pqtl_name <- gsub(".tsv", "", basename(pqtl_file))
    output_file <- file.path(output_dir, paste0("MR_SPI_", gwas_name, "_", pqtl_name, ".csv"))
    
    # Save results to CSV
    fwrite(result_df, output_file, sep = ",", quote = FALSE, row.names = FALSE)
    
    cat("Saved result to", output_file, "\n")
  }
}

cat("All analyses completed.\n")
