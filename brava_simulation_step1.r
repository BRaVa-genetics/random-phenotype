#!/usr/bin/env Rscript

include_fixed_effect_covariates_in_fit <- FALSE

if (include_fixed_effect_covariates_in_fit) {
    covariates <- "age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
    categorical_covariates <- "sex"
    destination <- "/Duncan/simulation_study/outputs/step1/"
} else {
    covariates <- ""
    categorical_covariates <- ""
    destination <- "/Duncan/simulation_study/outputs/step1_no_covariates/"
}

create_step1_sim_cmd <- function(
    phenotype, trait_type, covariates, categorical_covariates,
    sim_phenos_file="/Duncan/simulation_study/random_phenos_EUR_including_covariates.tsv.gz",
    instance_type="mem3_ssd1_v2_x4",
    priority="low",
    destination="/Duncan/simulation_study/outputs/step1/")
{
    name <- paste0("step1_", phenotype, "_EUR")
    output_prefix <- paste0("out/", phenotype, "_EUR")
    cmd <- paste0("dx run saige-universal-step-1",
        " -i output_prefix=", output_prefix,
        " -i sample_ids=/brava/inputs/ancestry_sample_ids/qced_EUR_sample_IDs.txt",
        " -i genotype_bed=/brava/outputs/step0/brava_EUR.plink_for_var_ratio.bed",
        " -i genotype_bim=/brava/outputs/step0/brava_EUR.plink_for_var_ratio.bim",
        " -i genotype_fam=/brava/outputs/step0/brava_EUR.plink_for_var_ratio.fam",
        " -i pheno_list=", sim_phenos_file, 
        " -i pheno=",phenotype, 
        " -i GRM=/brava/inputs/GRM/brava_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx",
        " -i GRM_samples=/brava/inputs/GRM/brava_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt",
        " -i covariates=", covariates,
        " -i categorical_covariates=", categorical_covariates,
        " -i trait_type=", trait_type,
        " --instance-type ", instance_type,
        " --priority ", priority,
        " --destination ", destination,
        " -y",
        " --name ", name)
}

# Required packages
packages <- c('data.table', 'dplyr')

for (p in packages) {
    if (!require(p, character.only = TRUE)) {
        install.packages(p, repos = "http://cran.us.r-project.org")
    }
}

# The covariates need to be added to the phenotype file
local_sim_pheno_list <- "random_phenos_EUR.tsv.gz"
covariates_vec <- c("age", "age2", "age_sex", "age2_sex", "sex", paste0("PC", seq(1,10)))

dt_pheno <- fread(local_sim_pheno_list)
dt_pheno[, IID := userID]
dt_pheno[, userID := NULL]

system("dx download brava/inputs/phenotypes/BRaVa_phenotypes_with_superpopulation_labels_updated_combined.tsv")
dt_covar <- fread("BRaVa_phenotypes_with_superpopulation_labels_updated_combined.tsv")
dt_covar <- dt_covar %>% select(all_of(c("IID", covariates_vec)))
dt_covar <- data.table(dt_covar)
setkey(dt_pheno, "IID")
setkey(dt_covar, "IID")

dt <- merge(dt_pheno, dt_covar)
cts_phenos <- grep("continuous", names(dt), value=TRUE)
binary_phenos <- setdiff(names(dt), c("IID", covariates_vec, cts_phenos))
print(binary_phenos)
fwrite(dt, file="random_phenos_EUR_including_covariates.tsv.gz", sep="\t")

# Ensure that the phenotype file is present on the RAP
system("dx upload random_phenos_EUR_including_covariates.tsv.gz --path /Duncan/simulation_study/")
system("dx build -f ~/Repositories/universal-saige-dnanexus/saige-universal-step-1")

for (pheno in binary_phenos) {
    cmd <- create_step1_sim_cmd(pheno, trait_type="binary",
        covariates=covariates,
        categorical_covariates=categorical_covariates,
        destination=destination)
    print(cmd)
    system(cmd)
}

for (pheno in cts_phenos) {
    cmd <- create_step1_sim_cmd(pheno, trait_type="quantitative",
        covariates=covariates,
        categorical_covariates=categorical_covariates,
        destination=destination)
    print(cmd)
    system(cmd)
}
