#!/usr/bin/env Rscript

# Required packages
packages <- c('data.table')

for (p in packages) {
    if (!require(p, character.only = TRUE)) {
        install.packages(p, repos = "http://cran.us.r-project.org")
    }
}

system("dx build -f ~/Repositories/universal-saige-dnanexus/saige-universal-step-2")

create_step2_sim_cmd <- function(
    chr, phenotype, pop,
    GRM, GRM_samples,
    test_type="group",
    model_file_location="/Duncan/simulation_study/outputs/step1/",
    group_file="/brava/inputs/annotations/v7/ukb_wes_450k.july.qced.brava_common_rare.v7.chr@.saige.txt.gz",
    exome_file="/Barney/wes/sample_filtered/ukb_wes_450k.qced.chr@",
    instance_type="mem3_ssd1_v2_x8", priority="low",
    destination="/Duncan/simulation_study/outputs/step2/",
    annotations="pLoF,damaging_missense_or_protein_altering,other_missense_or_protein_altering,synonymous,pLoF:damaging_missense_or_protein_altering,pLoF:damaging_missense_or_protein_altering:other_missense_or_protein_altering:synonymous",
    split=NULL
    ) {

    if (is.null(split)) {
        output_prefix <- paste0("out/chr", chr, "_", phenotype, "_", pop)
        name <- paste(chr, phenotype, pop, sep="_")
    } else {
        output_prefix <- paste0("out/chr", chr, "_split_", split, "_", phenotype, "_", pop)
        name <- paste(chr, "split", split, phenotype, pop, sep="_")
    }

    # Step 1 outputs
    model_file <- paste0(model_file_location, phenotype, "_", pop, ".rda")
    variance_ratio <- paste0(model_file_location, phenotype, "_", pop, ".varianceRatio.txt")
    
    group_file <- gsub("@", chr, group_file)

    # Exome plink files
    exome_bed <- paste0(gsub("@", chr, exome_file), ".bed")
    exome_bim <- paste0(gsub("@", chr, exome_file), ".bim")
    exome_fam <- paste0(gsub("@", chr, exome_file), ".fam")

    cmd <- paste0("dx run saige-universal-step-2",
        " -i output_prefix=", output_prefix,
        " -i model_file=", model_file, 
        " -i variance_ratio=", variance_ratio,
        " -i chrom=", chr,
        " -i group_file=", group_file,
        " -i annotations=", annotations,
        " -i test_type=", test_type,
        " -i exome_bed=", exome_bed,
        " -i exome_bim=", exome_bim,
        " -i exome_fam=", exome_fam,
        " -i GRM=", GRM,
        " -i GRM_samples=", GRM_samples,
        " --instance-type ", instance_type,
        " --priority ", priority,
        " --destination ", destination,
        " -y",
        " --name ", name
    )
    return(cmd)
}

include_fixed_effect_covariates_in_fit <- TRUE
if (include_fixed_effect_covariates_in_fit) {
    model_file_location <-"/Duncan/simulation_study/outputs/step1/"
    destination <- "/Duncan/simulation_study/outputs/step2/"
} else {
    model_file_location <- "/Duncan/simulation_study/outputs/step1_no_covariates/"
    destination <- "/Duncan/simulation_study/outputs/step2_no_covariates/"
}

# Ensure that the phenotype files locally and on the RAP match
# system("dx upload random_phenos_EUR_including_covariates.tsv.gz --path /Duncan/simulation_study/")
for (pop in c("AFR", "AMR", "EAS", "SAS"))
{
    dt <- fread(paste0("random_phenos_", pop, "_including_covariates.tsv.gz"))
    covariates <- c("age", "age2", "age_sex", "age2_sex", "sex", paste0("PC", seq(1,10)))
    phenos <- setdiff(names(dt), c("IID", covariates))
    outputs <- system(paste("dx ls", destination), intern=TRUE)
    outputs <- outputs[-grep("singleAssoc.txt.gz", outputs)]

    grm <- paste0("/brava/outputs/step0/brava_", pop,
        "_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx")
    grm_samples <- paste0("/brava/outputs/step0/brava_", pop,
        "_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")

    for (pheno in phenos) {
        for (chr in c(seq(2,22), "X")) {
            if (!any(grepl(paste0("chr", chr, "_", pheno, "_", pop), outputs))) {
                cmd <- create_step2_sim_cmd(
                    chr=chr, phenotype=pheno,
                    pop=pop, GRM=grm, GRM_samples=grm_samples,
                    destination=destination,
                    model_file_location=model_file_location)
                print(cmd)
                system(cmd)
            }
        }
    }
}
