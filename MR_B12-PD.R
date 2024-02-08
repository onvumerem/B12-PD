# install.packages("dplyr")
# install.packages("tibble")
# install.packages("psych")
# install.packages("plyr")
# install.packages("markdown")

install.packages("devtools")
devtools::install_github("rondolab/MR-PRESSO")

install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR") 
# remotes::install_github("MRCIEU/TwoSampleMR@0.4.26")

library(ggplot2)
library(dplyr)
library(MRPRESSO)
library(data.table)
library(googledrive)
library(TwoSampleMR)


PA_full_table_wo_casecontrol <- fread("PA_B12.txt")
PD_full_table_yesUKBB <- fread("META_no23_yesUKBB.txt")
PD_full_table_noUKBB <- fread("IPDGC_sum_stats_no_UKB.txt")
B12_full_table <- fread("Vi.B12.fastGWA")

PA_full_table <- PA_full_table_wo_casecontrol %>%
  mutate(n_cases = case_when(
    n_studies == 1 & n_samples == 390780 ~ 754,
    n_studies == 2 & n_samples == 529365 ~ 1132,
    n_studies == 2 & n_samples == 523097 ~ 1788,
    n_studies == 3 & n_samples == 661682 ~ 2166,
    n_studies == 1 & n_samples == 138585 ~ 378,
    n_studies == 2 & n_samples == 270902 ~ 1412,
    n_studies == 1 & n_samples == 132317 ~ 1034,
    TRUE ~ NA_real_
  ),
  n_controls = case_when(
    n_studies == 1 & n_samples == 390780 ~ 390026,
    n_studies == 2 & n_samples == 529365 ~ 528233,
    n_studies == 2 & n_samples == 523097 ~ 521309,
    n_studies == 3 & n_samples == 661682 ~ 659516,
    n_studies == 1 & n_samples == 138585 ~ 138207,
    n_studies == 2 & n_samples == 270902 ~ 269490,
    n_studies == 1 & n_samples == 132317 ~ 131283,
    TRUE ~ NA_real_))

write.table(PA_full_table, file="PA_full_table.txt", row.names=FALSE, quote=FALSE, sep = "\t")


### PART 1: B12 deficiency or Pernicious anemia as exposure, Parkinson's disease risk as outcome ###

PA_full_table_SNP <- PA_full_table %>% filter(PA_full_table[["p-value"]] < 5e-8 & PA_full_table[["p-value"]] != 0)
PA_full_table_SNP_1 <- format_data(PA_full_table_SNP_1, type = "exposure", snps = NULL, header = TRUE, snp_col = "rs_number", beta_col = "beta", se_col = "se", eaf_col = "eaf", effect_allele_col = "reference_allele", other_allele_col = "other_allele", pval_col = "p-value", z_col = "z", samplesize_col = "n_samples", ncase_col = "n_cases",  ncontrol_col = "n_controls", log_pval = FALSE)
PA_full_table_SNP_1_clumped <- clump_data(PA_full_table_SNP_1, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1) #, pop = "EUR")

B12_full_table_SNP_1 <- B12_full_table %>% filter(B12_full_table[["p_value"]] < 5e-8 & B12_full_table[["p_value"]] != 0)
B12_full_table_SNP_1 <- format_data(B12_full_table_SNP_1, type = "exposure", snps = NULL, header = TRUE, snp_col = "variant_id", beta_col = "BETA", se_col = "SE",  eaf_col = "AF1", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "p_value", log_pval = FALSE, samplesize_col = "N", chr_col = "chromosome", pos_col = "base_pair_location")
B12_full_table_SNP_1_clumped <- clump_data(B12_full_table_SNP_1, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1) #, pop = "EUR")

PD_fullsumstats_noUKBB_for_PA <- read_outcome_data(
  snps = PA_full_table_SNP_1_clumped$SNP,
  sep = "\t",
  filename = "IPDGC_sum_stats_no_UKB.txt",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "FreqSE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P-value",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol")

PD_fullsumstats_noUKBB_for_B12 <- read_outcome_data(
  snps = B12_full_table_SNP_1_clumped$SNP,
  sep = "\t",
  filename = "IPDGC_sum_stats_no_UKB.txt",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "FreqSE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P-value",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol")

PD_fullsumstats_yesUKBB_for_B12 <- read_outcome_data(
  snps = B12_full_table_SNP_1_clumped$SNP,
  sep = "\t",
  filename = "META_no23_yesUKBB.txt",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "Freq",
  pval_col = "p",
  ncase_col = "N_cases",
  ncontrol_col = "N_controls")

exp <- list(B12_full_table_SNP_1_clumped, B12_full_table_SNP_1_clumped, PA_full_table_SNP_1_clumped)
exp_names <- c("B12", "B12", "PA")
out <- list(PD_fullsumstats_yesUKBB_for_B12, PD_fullsumstats_noUKBB_for_B12, PD_fullsumstats_noUKBB_for_PA)
out_names <- c("PD yes UKBB", "PD no UKBB", "PD no UKB")

for(i in 1:length(exp)){
  folder <- paste(exp_names[i], "VS", out_names[i])
  unlink(paste("/kaggle/working/", folder), recursive = TRUE)
  unlink(paste("/kaggle/working/", folder, ".zip", sep=""), recursive = TRUE)
  dir.create(folder)
  exp_data <- exp[i]
  out_data <-out[i]
  dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
  dat$units.outcome <-"log odds"
  dat$units.exposure <-"log odds"
  dat1 <-subset(dat, dat$eaf.exposure!="NA")
  if (i == 3){
    dat1$prevalence <- round(dat1$ncase.exposure / dat1$ncontrol.exposure, 1)
    dat1$rsq.exposure <- get_r_from_lor(
      lor = dat1$beta.exposure,
      af = dat1$eaf.exposure,
      ncase = dat1$ncase.exposure,
      ncontrol = dat1$ncontrol.exposure,
      prevalence = dat1$prevalence,
      model = "logit",
      correction = FALSE)} 
  else {
    dat1$rsq.exposure <- get_r_from_pn(dat1$pval.exposure, dat1$samplesize.exposure)}
  dat1$rsq.outcome <- get_r_from_pn(dat1$pval.outcome, dat1$samplesize.outcome)
  steiger <- steiger_filtering(dat1)
  sig <- subset(steiger, steiger$steiger_dir==TRUE)
  presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig, NbDistribution = 1000,  SignifThreshold = 0.05)
  capture.output(print(presso), file = paste(folder, "presso.txt", sep = "/"))
  mr_res <- mr(sig)
  mr_res_with_OR <- generate_odds_ratios(mr_res)
  write.table(mr_res_with_OR, file=file.path(folder, "mr_res_with_OR.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
  R2 <- mean(sig$rsq.exposure)
  capture.output(print(R2), file = paste(folder, "r2.txt", sep = "/"))
  n <- mean(sig$samplesize.exposure)
  k <- nrow(subset(sig, sig$ambiguous == FALSE))
  F <- (R2*(n-1-k))/((1-R2)*k)
  capture.output(print(F), file = paste(folder, "f.txt", sep = "/"))
  mr_report(sig, study = folder,output_path = folder)
  res_single <- mr_singlesnp(sig)
  p5 <- mr_forest_plot(res_single)
  p5[[1]]
  ggsave(p5[[1]], file= "plot.jpg", path=folder, width=7, height=12)
  zip(zipfile = paste(folder, ".zip", sep=""), files = folder)
  gc()
}


### PART 2: B12 deficiency or Pernicious anemia as exposure, Parkinson's disease age at onsent or progression as outcome ###
### 2.1. Using rsID as SNP column 


PD_age <- fread("IPDGC_AAO_GWAS_sumstats_april_2018.txt")

# link changes, find it here: https://pdgenetics.shinyapps.io/pdprogmetagwasbrowser/ --> Download all results --> Download allele reference (131.5 MB)
url <- "https://pdgenetics.shinyapps.io/pdprogmetagwasbrowser/_w_1ac0e7ce/session/5d19d1791e6288d9b9986efa5fc113af/download/DLref?w=1ac0e7ce"
download.file(url, destfile = "DLref.csv", mode = "wb")
ref <- read.table("DLref.csv", header = TRUE, sep = ",", na.strings = "NA")

out_p <- c("PD_Hoehn and Yahr scale.txt", "cont_UPDRS3_scaled.txt", "PD_Cognitive_impairment.txt", "cont_MMSE.txt")
for(j in 1:length(out_p)){
  unlink(paste("/kaggle/working/", out_p[j], "_merged_with_ref.txt"), recursive = TRUE)
  df1 <- fread(out_p[j])
  merged_df <- merge(df1, ref, by="SNP")
  write.table(merged_df, file = paste(out_p[j], "_merged_with_ref.txt", sep=""), sep="\t", row.names=FALSE)
  gc()}


exp <- list(B12_full_table_SNP_1_clumped, PA_full_table_SNP_1_clumped)
exp_names <- c('B12', 'PA')
out_phenotypes <- c("PD_Hoehngend Yahr scale.txt", "cont_UPDRS3_scaled.txt", "PD_Cognitive_impairment.txt", "cont_MMSE.txt")
out_phenotypes_names <- c("PD_Hoehn_and_Yahr_scale", "PD_UPDRS3", "PD_Cognitive_impairment", "PD_MMSE")

for(i in 1:length(exp)){
  for(j in 2:length(out_phenotypes)){
    folder <- paste(exp_names[i], " VS ", out_phenotypes_names[j], sep="")
    unlink(paste("/kaggle/working/", folder, sep=""), recursive = TRUE)
    unlink(paste("/kaggle/working/", folder, ".zip", sep=""), recursive = TRUE)
    dir.create(folder)
    exp_data <- exp[i]
    if (j == 1){
      out_data <- read_outcome_data(snps = exp_data$SNP, filename = paste(out_phenotypes[j], "_merged_with_ref.txt", sep=""), sep="\t",  snp_col = "RSID.x", beta_col = "BETA.x", se_col = "SE.x", eaf_col = "MAF.x",effect_allele_col = "ALT.x", other_allele_col = "REF.x", pval_col = "P.x", samplesize_col = "N.x") #, chr_col = "CHR.x", pos_col = "START.x")
    }
    else {
      out_data <- read_outcome_data(snps = exp_data$SNP, filename = paste(out_phenotypes[j], "_merged_with_ref.txt", sep=""), sep="\t",  snp_col = "RSID", beta_col = "BETA", se_col = "SE", eaf_col = "MAF",effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "P", samplesize_col = "N")} #, chr_col = "CHR.x", pos_col = "START.x")
    dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
    gc()
    dat$units.outcome <-"log odds"
    dat$units.exposure <-"log odds"
    dat1 <-subset(dat, dat$eaf.exposure!="NA")
    if (i == 2){
      dat1$rsq.exposure <- get_r_from_lor(
        lor = dat1$beta.exposure,
        af = dat1$eaf.exposure,
        ncase = dat1$ncase.exposure,
        ncontrol = dat1$ncontrol.exposure,
        prevalence = dat1$prevalence,
        model = "logit",
        correction = FALSE)}
    else {
      dat1$rsq.exposure <- get_r_from_pn(dat1$pval.exposure, dat1$samplesize.exposure)}
    dat1$rsq.outcome <- get_r_from_pn(dat1$pval.outcome, dat1$samplesize.outcome)
    steiger <- steiger_filtering(dat1)
    sig <- subset(steiger, steiger$steiger_dir==TRUE)
    presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig, NbDistribution = 1000,  SignifThreshold = 0.05)
    capture.output(print(presso), file = paste(folder, "presso.txt", sep = "/"))
    mr_res <- mr(sig)
    mr_res_with_OR <- generate_odds_ratios(mr_res)
    write.table(mr_res_with_OR, file=file.path(folder, "mr_res_with_OR.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
    R2 <- mean(sig$rsq.exposure)
    capture.output(print(R2), file = paste(folder, "r2.txt", sep = "/"))
    n <- mean(sig$samplesize.exposure)
    k <- nrow(subset(sig, sig$ambiguous == FALSE))
    F <- (R2*(n-1-k))/((1-R2)*k)
    capture.output(print(F), file = paste(folder, "f.txt", sep = "/"))
    mr_report(sig, study = folder,output_path = folder)
    res_single <- mr_singlesnp(sig)
    p5 <- mr_forest_plot(res_single)
    p5[[1]]
    ggsave(p5[[1]], file= "plot.jpg", path=folder , width=7, height=12)
    zip(zipfile = paste(folder, ".zip", sep=""), files = folder)
    gc()
  }
}


### 2.1. Using chr:pos as SNP column 

# chromosome = c("chr1", "chr1", "chr2", "chr6", "chr6", "chr10", "chr21", "chr6", "chr6", "chr6", "chr6")
# base_pair_location = c(113761186, 113893860, 55581879, 32393072, 32658674, 6053445, 44294411, 32559257, 32655666, 32737819, 32633199)
# chrpos <- c(
#   "chr1:113761186",
#   "chr1:113893860",
#   "chr2:55581879",
#   "chr6:32393072",
#   "chr6:32658674",
#   "chr10:6053445",
#   "chr21:44294411",
#   "chr6:32559257",
#   "chr6:32655666",
#   "chr6:32737819",
#   "chr6:32633199"
# )

# PA_full_table_SNP_1_clumped$chrpos <- chrpos
# PA_full_table_SNP_1_clumped$chromosome <- chromosome
# PA_full_table_SNP_1_clumped$base_pair_location <- base_pair_location
# PA_full_table_SNP_1_clumped

# PA_chrpos_asID <- format_data(PA_full_table_SNP_1_clumped, type = "exposure", snps = NULL, header = TRUE, snp_col = "chrpos", beta_col = "beta.exposure", se_col = "se.exposure",  eaf_col = "eaf.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", pval_col = "pval.exposure", log_pval = FALSE, samplesize_col = "samplesize.exposure", chr_col = "chromosome", pos_col = "base_pair_location", ncase_col = "ncase.exposure", ncontrol_col = "ncontrol.exposure", z_col = "z.exposure")
# PD_age_for_PA <- read_outcome_data(
#   snps = PA_chrpos_asID$SNP,
#   sep = "\t",
#   filename = "IPDGC_AAO_GWAS_sumstats_april_2018.txt",
#   snp_col = "MarkerName",
#   beta_col = "Effect",
#   se_col = "FreqSE",
#   effect_allele_col = "Allele1",
#   other_allele_col = "Allele2",
#   eaf_col = "Freq1",
#   pval_col = "P-value")

B12_chrpos_asID <- format_data(B12_full_table_SNP_1_clumped, type = "exposure", snps = NULL, header = TRUE, snp_col = "chrpos", beta_col = "beta.exposure", se_col = "se.exposure",  eaf_col = "eaf.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", pval_col = "pval.exposure", log_pval = FALSE, samplesize_col = "samplesize.exposure", chr_col = "chromosome", pos_col = "pos.exposure")
PD_age_for_B12 <- read_outcome_data(
  snps = B12_chrpos_asID$SNP,
  sep = "\t",
  filename = "IPDGC_AAO_GWAS_sumstats_april_2018.txt",
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "FreqSE",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P-value")

exp <- list(B12_chrpos_asID)
exp_names <- c("B12 (chr:pos as ID)")
out <- list(PD_age_for_B12)
out_names <- c("PD age")

for(i in 1:length(exp)){
  folder <- paste(exp_names[i], "VS", out_names[i])
  unlink(paste("/kaggle/working/", folder), recursive = TRUE)
  unlink(paste("/kaggle/working/", folder, ".zip", sep=""), recursive = TRUE)
  dir.create(folder)
  exp_data <- exp[i]
  out_data <-out[i]
  dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
  dat$units.outcome <-"log odds"
  dat$units.exposure <-"log odds"
  dat1 <-subset(dat, dat$eaf.exposure!="NA")
  dat1$rsq.exposure <- get_r_from_pn(dat1$pval.exposure, dat1$samplesize.exposure)
  dat1$samplesize.outcome <- 28568
  dat1$rsq.outcome <- get_r_from_pn(dat1$pval.outcome, dat1$samplesize.outcome)
  steiger <- steiger_filtering(dat1)
  sig <- subset(steiger, steiger$steiger_dir==TRUE)
  presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig, NbDistribution = 1000,  SignifThreshold = 0.05)
  capture.output(print(presso), file = paste(folder, "presso.txt", sep = "/"))
  mr_res <- mr(sig)
  mr_res_with_OR <- generate_odds_ratios(mr_res)
  write.table(mr_res_with_OR, file=file.path(folder, "mr_res_with_OR.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
  R2 <- mean(sig$rsq.exposure)
  capture.output(print(R2), file = paste(folder, "r2.txt", sep = "/"))
  n <- mean(sig$samplesize.exposure)
  k <- nrow(subset(sig, sig$ambiguous == FALSE))
  F <- (R2*(n-1-k))/((1-R2)*k)
  capture.output(print(F), file = paste(folder, "f.txt", sep = "/"))
  mr_report(sig, study = folder,output_path = folder)
  res_single <- mr_singlesnp(sig)
  p5 <- mr_forest_plot(res_single)
  p5[[1]]
  ggsave(p5[[1]], file= "plot.jpg", path=folder, width=7, height=12)
  zip(zipfile = paste(folder, ".zip", sep=""), files = folder)
  gc()
}


exp <- list(B12_chrpos_asID)
exp_names <- c("B12 (chr:pos as ID)")
out_phenotypes <- c("PD_Hoehngend Yahr scale.txt", "cont_UPDRS3_scaled.txt", "PD_Cognitive_impairment.txt", "cont_MMSE.txt")
out_phenotypes_names <- c("PD_Hoehn_and_Yahr_scale", "PD_UPDRS3", "PD_Cognitive_impairment", "PD_MMSE")


for(i in 1:length(exp)){
  for(j in 2:length(out_phenotypes)){
    folder <- paste(exp_names[i], " VS ", out_phenotypes_names[j], sep="")
    unlink(paste("/kaggle/working/", folder, sep=""), recursive = TRUE)
    unlink(paste("/kaggle/working/", folder, ".zip", sep=""), recursive = TRUE)
    dir.create(folder)
    exp_data <- exp[i]
    if (j == 1){
      out_data <- read_outcome_data(snps = exp_data$SNP, filename = paste(out_phenotypes[j], "_merged_with_ref.txt", sep=""), sep="\t",  snp_col = "RSID.x", beta_col = "BETA.x", se_col = "SE.x", eaf_col = "MAF.x",effect_allele_col = "ALT.x", other_allele_col = "REF.x", pval_col = "P.x", samplesize_col = "N.x") #, chr_col = "CHR.x", pos_col = "START.x")
    }
    else {
      out_data <- read_outcome_data(snps = exp_data$SNP, filename = paste(out_phenotypes[j], "_merged_with_ref.txt", sep=""), sep="\t",  snp_col = "RSID", beta_col = "BETA", se_col = "SE", eaf_col = "MAF",effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "P", samplesize_col = "N")} #, chr_col = "CHR.x", pos_col = "START.x")
    dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
    gc()
    dat$units.outcome <-"log odds"
    dat$units.exposure <-"log odds"
    dat1 <-subset(dat, dat$eaf.exposure!="NA")
    if (i == 2){
      dat1$rsq.exposure <- get_r_from_lor(
        lor = dat1$beta.exposure,
        af = dat1$eaf.exposure,
        ncase = dat1$ncase.exposure,
        ncontrol = dat1$ncontrol.exposure,
        prevalence = dat1$prevalence,
        model = "logit",
        correction = FALSE)}
    else {
      dat1$rsq.exposure <- get_r_from_pn(dat1$pval.exposure, dat1$samplesize.exposure)}
    dat1$rsq.outcome <- get_r_from_pn(dat1$pval.outcome, dat1$samplesize.outcome)
    steiger <- steiger_filtering(dat1)
    sig <- subset(steiger, steiger$steiger_dir==TRUE)
    presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig, NbDistribution = 1000,  SignifThreshold = 0.05)
    capture.output(print(presso), file = paste(folder, "presso.txt", sep = "/"))
    mr_res <- mr(sig)
    mr_res_with_OR <- generate_odds_ratios(mr_res)
    write.table(mr_res_with_OR, file=file.path(folder, "mr_res_with_OR.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
    R2 <- mean(sig$rsq.exposure)
    capture.output(print(R2), file = paste(folder, "r2.txt", sep = "/"))
    n <- mean(sig$samplesize.exposure)
    k <- nrow(subset(sig, sig$ambiguous == FALSE))
    F <- (R2*(n-1-k))/((1-R2)*k)
    capture.output(print(F), file = paste(folder, "f.txt", sep = "/"))
    mr_report(sig, study = folder,output_path = folder)
    res_single <- mr_singlesnp(sig)
    p5 <- mr_forest_plot(res_single)
    p5[[1]]
    ggsave(p5[[1]], file= "plot.jpg", path=folder , width=7, height=12)
    zip(zipfile = paste(folder, ".zip", sep=""), files = folder)
    gc()
  }
}


### PART 3: REVERSE MR ###

### 3.1 Parkinson's disease risk as exposure, B12 deficiency or Pernicious anemia as outcome

PD_full_table_noUKBB_SNP_1 <- PD_full_table_noUKBB %>% filter(PD_full_table_noUKBB[["P-value"]] < 5e-8 & PD_full_table_noUKBB[["P-value"]] != 0)
PD_full_table_noUKBB_SNP_1 <- format_data(PD_full_table_noUKBB_SNP_1, type = "exposure", snps = NULL, header = TRUE, snp_col = "SNP", beta_col = "Effect", se_col = "FreqSE", effect_allele_col = "Allele1", other_allele_col = "Allele2", eaf_col = "Freq1", pval_col = "P-value", ncase_col = "ncase", ncontrol_col = "ncontrol", log_pval = FALSE)
PD_full_table_noUKBB_SNP_1_clumped <- clump_data(PD_full_table_noUKBB_SNP_1, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1) #, pop = "EUR")

PD_full_table_yesUKBB_SNP_1 <- PD_full_table_yesUKBB %>% filter(PD_full_table_yesUKBB[["p"]] < 5e-8 & PD_full_table_yesUKBB[["p"]] != 0)
PD_full_table_yesUKBB_SNP_1 <- format_data(PD_full_table_yesUKBB_SNP_1, type = "exposure", snps = NULL, header = TRUE, snp_col = "SNP", beta_col = "b", se_col = "se", effect_allele_col = "A1", other_allele_col = "A2", eaf_col = "Freq", pval_col = "p", ncase_col = "N_cases", ncontrol_col = "N_controls", log_pval = FALSE)
PD_full_table_yesUKBB_SNP_1_clumped <- clump_data(PD_full_table_yesUKBB_SNP_1, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1) #, pop = "EUR")

PA_fullsumstats_for_PD <- read_outcome_data(snps = PD_full_table_noUKBB_SNP_1_clumped$SNP, sep = "\t", filename = "PA_B12.txt", snp_col = "rs_number", beta_col = "beta", se_col = "se", eaf_col = "eaf", effect_allele_col = "reference_allele", other_allele_col = "other_allele", pval_col = "p-value", samplesize_col = "n_samples",ncase_col = "n_cases", ncontrol_col = "n_controls", log_pval = FALSE)#  z_col = "z",
B12_fullsumstats_for_PD <- read_outcome_data(snps = PD_full_table_noUKBB_SNP_1_clumped$SNP, sep = "\t", filename = "Vi.B12.fastGWA", snp_col = "variant_id", beta_col = "BETA", se_col = "SE",  eaf_col = "AF1", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "p_value", log_pval = FALSE, samplesize_col = "N") #chr_col = "chromosome", pos_col = "base_pair_location")

exp <- list(PD_full_table_noUKBB_SNP_1_clumped, PD_full_table_noUKBB_SNP_1_clumped)
exp_names <- c("PD no UKBB", "PD no UKB")
out <- list(PA_fullsumstats_for_PD, B12_fullsumstats_for_PD)
out_names <- c("PA", "B12")


for(i in 1:length(out)){
  folder <- paste(exp_names[i], "VS", out_names[i])
  unlink(paste("/kaggle/working/", folder), recursive = TRUE)
  unlink(paste("/kaggle/working/", folder, ".zip", sep=""), recursive = TRUE)
  dir.create(folder)
  exp_data <- exp[i]
  out_data <-out[i]
  dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
  dat$units.outcome <-"log odds"
  dat$units.exposure <-"log odds"
  dat1 <-subset(dat, dat$eaf.exposure!="NA")
  dat1$rsq.exposure <- get_r_from_pn(dat1$pval.exposure, dat1$samplesize.exposure)
  if (i == 1){
    dat1$prevalence <- round(dat1$ncase.outcome / dat1$ncontrol.outcome, 1)
    dat1$rsq.outcome <- get_r_from_lor(
      lor = dat1$beta.outcome,
      af = dat1$eaf.outcome,
      ncase = dat1$ncase.outcome,
      ncontrol = dat1$ncontrol.outcome,
      prevalence = dat1$prevalence,
      model = "logit",
      correction = FALSE)} 
  else {
    dat1$rsq.outcome <- get_r_from_pn(dat1$pval.outcome, dat1$samplesize.outcome)}
  steiger <- steiger_filtering(dat1)
  sig <- subset(steiger, steiger$steiger_dir==TRUE)
  presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig, NbDistribution = 1000,  SignifThreshold = 0.05)
  capture.output(print(presso), file = paste(folder, "presso.txt", sep = "/"))
  mr_res <- mr(sig)
  mr_res_with_OR <- generate_odds_ratios(mr_res)
  write.table(mr_res_with_OR, file=file.path(folder, "mr_res_with_OR.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
  R2 <- mean(sig$rsq.exposure)
  capture.output(print(R2), file = paste(folder, "r2.txt", sep = "/"))
  n <- mean(sig$samplesize.exposure)
  k <- nrow(subset(sig, sig$ambiguous == FALSE))
  F <- (R2*(n-1-k))/((1-R2)*k)
  capture.output(print(F), file = paste(folder, "f.txt", sep = "/"))
  mr_report(sig, study = folder,output_path = folder)
  res_single <- mr_singlesnp(sig)
  p5 <- mr_forest_plot(res_single)
  p5[[1]]
  ggsave(p5[[1]], file= "plot.jpg", path=folder, width=7, height=12)
  zip(zipfile = paste(folder, ".zip", sep=""), files = folder)
  gc()
}


### 3.2. Parkinson's disease penetrance or age of onsent (LRRK2-associated) as exposure, B12 deficiency or Pernicious anemia as outcome

PD_penetrance_LRRK2 <- fread("Penetrance model.csv")
PD_age_LRRK2 <- fread("Age-at-onset model.csv")
PD_penetrance_LRRK2$Ncase <- 776
PD_penetrance_LRRK2$Ncontrol <- 1103
PD_age_LRRK2$Ncase <- 776
PD_age_LRRK2$Ncontrol <- 1103

PD_age_LRRK2_out <- format_data(PD_age_LRRK2, type = "exposure", snps = NULL, header = TRUE,  snp_col = "rsid", beta_col = "BETA",
                                se_col = "SE",
                                eaf_col = "MAF",
                                effect_allele_col = "Allele 1",
                                other_allele_col = "Allele 2",
                                pval_col = "P-Value",
                                ncase_col = "Ncase",
                                ncontrol_col = "Ncontrol",
                                gene_col = "Gene",
                                min_pval = 1e-200,
                                chr_col = "CHR",
                                pos_col = "BP",
                                log_pval = FALSE
)

PD_age_LRRK2_out_clumped <- clump_data(PD_age_LRRK2_out, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1) #, pop = "EUR")

PD_penetrance_LRRK2_out <- format_data(PD_penetrance_LRRK2, type = "exposure", snps = NULL, header = TRUE,  snp_col = "rsid", beta_col = "BETA",
                                       se_col = "SE",
                                       eaf_col = "MAF",
                                       effect_allele_col = "Allele 1",
                                       other_allele_col = "Allele 2",
                                       pval_col = "P-Value",
                                       ncase_col = "Ncase",
                                       ncontrol_col = "Ncontrol",
                                       gene_col = "Gene",
                                       min_pval = 1e-200,
                                       chr_col = "CHR",
                                       pos_col = "BP",
                                       log_pval = FALSE
)

PD_penetrance_LRRK2_out_clumped <- clump_data(PD_penetrance_LRRK2_out, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1) #, pop = "EUR")


exp <- list(PD_penetrance_LRRK2_out, PD_age_LRRK2_out, PD_penetrance_LRRK2_out_clumped, PD_age_LRRK2_out_clumped)
exp_names <- c('PD_penetrance_LRRK2', 'PD_age_LRRK2','PD_penetrance_LRRK2_out_clumped', 'PD_age_LRRK2_out_clumped')
out_phenotypes <- c("PA_full_table.txt", "Vi.B12.fastGWA")
out_phenotypes_names <- c("PA", "B12")

for(i in 1:length(exp)){
  for(i in list(1, 2)){
    for(j in 1:length(out_phenotypes)){
      folder <- paste(exp_names[i], " VS ", out_phenotypes_names[j], sep="")
      unlink(paste("/kaggle/working/", folder, sep=""), recursive = TRUE)
      unlink(paste("/kaggle/working/", folder, ".zip", sep=""), recursive = TRUE)
      dir.create(folder)
      exp_data <- exp[i]
      if (j == 1){
        out_data <- read_outcome_data(filename = out_phenotypes[j], snps = exp_data$SNP, sep="\t", snp_col = "rs_number", beta_col = "beta", se_col = "se", eaf_col = "eaf", effect_allele_col = "reference_allele", other_allele_col = "other_allele", pval_col = "p-value", samplesize_col = "n_samples", ncase_col = "n_cases",  ncontrol_col = "n_controls")
      }
      else {
        out_data <- read_outcome_data(filename = out_phenotypes[j], snps = exp_data$SNP, sep = "\t", snp_col = "variant_id", beta_col = "BETA", se_col = "SE",  eaf_col = "AF1", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "p_value", samplesize_col = "N")} #chr_col = "chromosome", pos_col = "base_pair_location")
      dat <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
      gc()
      dat$units.outcome <-"log odds"
      dat$units.exposure <-"log odds"
      dat1 <-subset(dat, dat$eaf.exposure!="NA")
      
      if (i == 1 || i == 3){
        dat1$prevalenceexp <- round(dat1$ncase.exposure / dat1$ncontrol.exposure, 1)
        dat1$rsq.exposure <- get_r_from_lor(
          lor = dat1$beta.exposure,
          af = dat1$eaf.exposure,
          ncase = dat1$ncase.exposure,
          ncontrol = dat1$ncontrol.exposure,
          prevalence = dat1$prevalenceexp,
          model = "logit",
          correction = FALSE)}
      else {
        dat1$rsq.exposure <- get_r_from_pn(dat1$pval.exposure, dat1$samplesize.exposure)}
      
      if (j == 1){
        dat1$prevalenceout <- round(dat1$ncase.outcome / dat1$ncontrol.outcome, 1)
        dat1$rsq.outcome <- get_r_from_lor(
          lor = dat1$beta.outcome,
          af = dat1$eaf.outcome,
          ncase = dat1$ncase.outcome,
          ncontrol = dat1$ncontrol.outcome,
          prevalence = dat1$prevalenceout,
          model = "logit",
          correction = FALSE)}
      else {
        dat1$rsq.outcome <- get_r_from_pn(dat1$pval.outcome, dat1$samplesize.outcome)}
      
      steiger <- steiger_filtering(dat1)
      sig <- subset(steiger, steiger$steiger_dir==TRUE)
      presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig, NbDistribution = 1000,  SignifThreshold = 0.05)
      capture.output(print(presso), file = paste(folder, "presso.txt", sep = "/"))
      mr_res <- mr(sig)
      mr_res_with_OR <- generate_odds_ratios(mr_res)
      write.table(mr_res_with_OR, file=file.path(folder, "mr_res_with_OR.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
      R2 <- mean(sig$rsq.exposure)
      capture.output(print(R2), file = paste(folder, "r2.txt", sep = "/"))
      n <- mean(sig$samplesize.exposure)
      k <- nrow(subset(sig, sig$ambiguous == FALSE))
      F <- (R2*(n-1-k))/((1-R2)*k)
      capture.output(print(F), file = paste(folder, "f.txt", sep = "/"))
      mr_report(sig, study = folder,output_path = folder)
      res_single <- mr_singlesnp(sig)
      p5 <- mr_forest_plot(res_single)
      p5[[1]]
      ggsave(p5[[1]], file= "plot.jpg", path=folder , width=7, height=12)
      zip(zipfile = paste(folder, ".zip", sep=""), files = folder)
      gc()
    }
  } 