############################################
# ABF colocalisation (coloc) between:
#  - deCODE plasma pQTLs (protein)
#  - MS GWAS (OpenGWAS: ieu-b-18)
#
# Notes:
# - Requires pQTL summary stats that include eaf (effect allele frequency)
# - Uses harmonise_data(action=2) to align effects prior to coloc
############################################

# ---- Libraries ----
library(TwoSampleMR)
library(coloc)     # ABF coloc: coloc.abf()
library(dplyr)
library(stringr)
library(readr)

# ---- User parameters ----

# 1) deCODE pQTL regional summary stats
# Must include:
# phenotype, SNP, chr_name, chrom_start, effect_allele, other_allele, beta, se, pval, eaf
# pqtl_sumstats file is already cis-filtered per protein
pqtl_sumstats_filepath <- "/path/to/decode_pqtl_cis_sumstats.txt"

# 2) List of proteins to test for colocalisation.
# provide a simple text file with one id.exposure (phenotype string) per line
proteins_list_filepath <- "/path/to/proteins_to_coloc.txt"

# 3) MS OpenGWAS outcome ID
ms_outcome_id <- "ieu-b-18"

# 4) MS GWAS sample sizes (needed for coloc case-control dataset)
ms_n_cases    <- 47429
ms_n_controls <- 68374

# 5) Coloc threshold for "strong evidence"
pp_h4_strong_threshold <- 0.8

# 6) Output directory
output_dir <- "/path/to/output/coloc_results"

# ---- Load pQTL summary stats ----
# Expecting a tab-delimited file
pqtl <- read.delim(pqtl_sumstats_filepath, stringsAsFactors = FALSE)

# Standardise/rename to TwoSampleMR-like columns
pqtl <- pqtl %>%
 rename(
  phenotype = phenotype,
  SNP = SNP,
  chr_name = chr_name,
  chrom_start = chrom_start,
  effect_allele.exposure = effect_allele,
  other_allele.exposure  = other_allele,
  beta.exposure = beta,
  se.exposure   = se,
  pval.exposure = pval,
  eaf.exposure  = eaf
 )

pqtl$id.exposure <- pqtl$phenotype
pqtl$exposure <- pqtl$phenotype

# QC: remove missing essentials
pqtl <- pqtl %>%
 filter(
  !is.na(SNP),
  !is.na(beta.exposure),
  !is.na(se.exposure),
  !is.na(pval.exposure),
  !is.na(effect_allele.exposure),
  !is.na(other_allele.exposure),
  !is.na(eaf.exposure)
 )

# ---- Determine proteins to coloc ----
proteins_to_coloc <- unique(mr_output_IVW_df_sig$id.exposure)

message("Number of proteins to coloc: ", length(proteins_to_coloc))

# ---- Coloc runner per protein ----
run_coloc_for_protein <- function(protein_id) {
 
 # Subset pQTL region for this protein.
 exp_dat <- pqtl %>% filter(id.exposure == protein_id)
 
 # Pull MS outcome stats for these SNPs from OpenGWAS
 out_dat <- extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = ms_outcome_id
 )

 # Harmonise alleles so exposure/outcome betas correspond to the same effect allele
 dat_h <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat  = out_dat,
  action = 2
 )
 
 # coloc needs: beta, varbeta, snp, position, MAF, N, type (+ s for cc)
 # For MAF use min(eaf, 1-eaf)
 dat_h <- dat_h %>%
  mutate(
   maf.exposure = pmin(eaf.exposure, 1 - eaf.exposure),
   maf.outcome  = if ("eaf.outcome" %in% colnames(dat_h)) pmin(eaf.outcome, 1 - eaf.outcome) else NA_real_
  ) %>%
  filter(
   !is.na(maf.exposure),
   !is.na(maf.outcome),
   maf.exposure > 0, maf.outcome > 0
  )

 # Assemble coloc datasets
 d1 <- list(
  beta    = dat_h$beta.exposure,
  varbeta = (dat_h$se.exposure)^2,
  snp     = dat_h$SNP,
  position= dat_h$chrom_start.exposure,   
  MAF     = dat_h$maf.exposure,
  N       = nrow(dat_h),                 
  type    = "quant"
 )
 
 totalN <- ms_n_cases + ms_n_controls
 s <- ms_n_cases / totalN
 
 d2 <- list(
  beta    = dat_h$beta.outcome,
  varbeta = (dat_h$se.outcome)^2,
  snp     = dat_h$SNP,
  position= dat_h$chrom_start.outcome,    
  MAF     = dat_h$maf.outcome,
  N       = totalN,
  s       = s,
  type    = "cc"
 )
 
 # Run coloc ABF
 coloc_res <- coloc.abf(dataset1 = d1, dataset2 = d2)
 
 # Extract summary posteriors
 summ <- coloc_res$summary
 # Identify lead SNP by max SNP.PP.H4 (posterior for being shared causal)
 # (coloc returns per-SNP posteriors in $results)
 lead <- coloc_res$results %>%
  arrange(desc(SNP.PP.H4)) %>%
  slice(1)
 
 out <- list(
  protein_id = protein_id,
  nsnps = nrow(dat_h),
  PP.H0 = unname(summ["PP.H0.abf"]),
  PP.H1 = unname(summ["PP.H1.abf"]),
  PP.H2 = unname(summ["PP.H2.abf"]),
  PP.H3 = unname(summ["PP.H3.abf"]),
  PP.H4 = unname(summ["PP.H4.abf"]),
  lead_snp = lead$snp,
  lead_snp_pp_h4 = lead$SNP.PP.H4,
  coloc_object = coloc_res
 )
 
 return(out)
}

# ---- Run coloc across proteins ----
results_list <- vector("list", length(proteins_to_coloc))
names(results_list) <- proteins_to_coloc

for (i in seq_along(proteins_to_coloc)) {
 prot <- proteins_to_coloc[i]
 message("Coloc: ", prot, " (", i, "/", length(proteins_to_coloc), ")")
 results_list[[prot]] <- run_coloc_for_protein(prot)
}

# ---- Build results table ----
coloc_summary_df <- bind_rows(lapply(results_list, function(x) {
 data.frame(
  id.exposure = x$protein_id,
  nsnps = x$nsnps,
  PP.H0 = x$PP.H0,
  PP.H1 = x$PP.H1,
  PP.H2 = x$PP.H2,
  PP.H3 = x$PP.H3,
  PP.H4 = x$PP.H4,
  lead_snp = x$lead_snp,
  lead_snp_pp_h4 = x$lead_snp_pp_h4,
  strong_coloc = (x$PP.H4 >= pp_h4_strong_threshold),
  stringsAsFactors = FALSE
 )
})) %>% arrange(desc(PP.H4))

# ---- Save outputs ----
save(results_list, coloc_summary_df,
     file = file.path(output_dir, "coloc_ms_decode_results.Robj"))

write.table(coloc_summary_df,
            file = file.path(output_dir, "coloc_ms_decode_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)