#!/usr/bin/env Rscript

# Required packages
packages <- c('data.table', 'dplyr')

for (p in packages) {
    if (!require(p, character.only = TRUE)) {
        install.packages(p, repos = "http://cran.us.r-project.org")
    }
}

system("dx build -f ~/Repositories/universal-saige-dnanexus/saige-universal-step-1")

create_step1_sim_cmd <- function(
    phenotype, trait_type, covariates, categorical_covariates,
    population_label, sim_phenos_file,
    instance_type="mem3_ssd1_v2_x4", priority="low",
    destination="/Duncan/simulation_study/outputs/step1/")
{
    pop <- population_label
    name <- paste0("step1_", phenotype, "_", pop)
    output_prefix <- paste0("out/", phenotype, "_", pop)
    cmd <- paste0("dx run saige-universal-step-1",
        " -i output_prefix=", output_prefix,
        " -i sample_ids=/brava/inputs/ancestry_sample_ids/qced_", pop, "_sample_IDs.txt",
        " -i genotype_bed=/brava/outputs/step0/brava_", pop ,".plink_for_var_ratio.bed",
        " -i genotype_bim=/brava/outputs/step0/brava_", pop, ".plink_for_var_ratio.bim",
        " -i genotype_fam=/brava/outputs/step0/brava_", pop, ".plink_for_var_ratio.fam",
        " -i pheno_list=", sim_phenos_file, 
        " -i pheno=",phenotype, 
        " -i GRM=/brava/outputs/step0/brava_", pop, "_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx",
        " -i GRM_samples=/brava/outputs/step0/brava_", pop, "_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt",
        " -i covariates=", covariates,
        " -i categorical_covariates=", categorical_covariates,
        " -i trait_type=", trait_type,
        " --instance-type ", instance_type,
        " --priority ", priority,
        " --destination ", destination,
        " -y",
        " --name ", name)
}

merge_in_covariates_for_population_label_run <- function(population_label, upload=FALSE)
{
    pop <- population_label
    # The covariates need to be added to the phenotype file
    local_sim_pheno_list <- paste0("random_phenos_", pop, ".tsv.gz")
    covariates_vec <- c("age", "age2", "age_sex", "age2_sex", "sex", paste0("PC", seq(1,10)))

    dt_pheno <- fread(local_sim_pheno_list)
    dt_pheno[, IID := userID]
    dt_pheno[, userID := NULL]

    if (!file.exists("BRaVa_phenotypes_with_superpopulation_labels_updated_combined.tsv")) {
        system("dx download brava/inputs/phenotypes/BRaVa_phenotypes_with_superpopulation_labels_updated_combined.tsv")
    }

    dt_covar <- fread("BRaVa_phenotypes_with_superpopulation_labels_updated_combined.tsv")
    dt_covar <- dt_covar %>% select(all_of(c("IID", covariates_vec)))
    dt_covar <- data.table(dt_covar)
    setkey(dt_pheno, "IID")
    setkey(dt_covar, "IID")

    dt <- merge(dt_pheno, dt_covar)
    cts_phenos <- grep("continuous", names(dt), value=TRUE)
    binary_phenos <- setdiff(names(dt), c("IID", covariates_vec, cts_phenos))
    print(binary_phenos)
    fwrite(dt, file=paste0("random_phenos_", pop, "_including_covariates.tsv.gz"), sep="\t")

    # Ensure that the phenotype file is present on the RAP
    if (upload) {
        system(paste0("dx upload random_phenos_", pop, "_including_covariates.tsv.gz --path /Duncan/simulation_study/"))
    }
    return(list(binary=binary_phenos, cts=cts_phenos))
}

include_fixed_effect_covariates_in_fit <- TRUE

if (include_fixed_effect_covariates_in_fit) {
    covariates <- "age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
    categorical_covariates <- "sex"
    destination <- "/Duncan/simulation_study/outputs/step1/"
} else {
    covariates <- ""
    categorical_covariates <- ""
    destination <- "/Duncan/simulation_study/outputs/step1_no_covariates/"
}

for (pop in c("AFR", "AMR", "EAS", "SAS"))
{
    phenos <- merge_in_covariates_for_population_label_run(pop)
    binary_phenos <- phenos$binary
    cts_phenos <- phenos$cts
    sim_phe_file <- paste0("/Duncan/simulation_study/random_phenos_", pop, "_including_covariates.tsv.gz")
    for (pheno in binary_phenos) {
        cmd <- create_step1_sim_cmd(pheno, trait_type="binary",
            covariates=covariates,
            categorical_covariates=categorical_covariates,
            population_label=pop, sim_phenos_file=sim_phe_file,
            destination=destination)
        print(cmd)
        system(cmd)
    }

    for (pheno in cts_phenos) {
        print(pheno)
        cmd <- create_step1_sim_cmd(pheno, trait_type="quantitative",
            covariates=covariates,
            categorical_covariates=categorical_covariates,
            population_label=pop, sim_phenos_file=sim_phe_file,
            destination=destination)
        print(cmd)
        system(cmd)
    }
}

