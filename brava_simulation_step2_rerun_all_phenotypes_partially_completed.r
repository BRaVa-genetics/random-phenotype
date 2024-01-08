#!/usr/bin/env Rscript

system("dx build -f ~/Repositories/universal-saige-dnanexus/saige-universal-step-2")

include_fixed_effect_covariates_in_fit <- FALSE
if (include_fixed_effect_covariates_in_fit) {
    model_file_location <-"/Duncan/simulation_study/outputs/step1/"
    destination <- "/Duncan/simulation_study/outputs/step2/"
} else {
    model_file_location <- "/Duncan/simulation_study/outputs/step1_no_covariates/"
    destination <- "/Duncan/simulation_study/outputs/step2_no_covariates/"
}

create_step2_sim_cmd <- function(
    chr, phenotype,
    test_type="group",
    model_file_location="/Duncan/simulation_study/outputs/step1/",
    group_file="/brava/inputs/annotations/v7/ukb_wes_450k.july.qced.brava_common_rare.v7.chr@.saige.txt.gz",
    exome_file="/Barney/wes/sample_filtered/ukb_wes_450k.qced.chr@",
    GRM="/brava/inputs/GRM/brava_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx",
    GRM_samples="/brava/inputs/GRM/brava_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt",
    instance_type="mem3_ssd1_v2_x8", priority="low",
    destination="/Duncan/simulation_study/outputs/step2/",
    annotations="pLoF,damaging_missense_or_protein_altering,other_missense_or_protein_altering,synonymous,pLoF:damaging_missense_or_protein_altering,pLoF:damaging_missense_or_protein_altering:other_missense_or_protein_altering:synonymous",
    split=NULL
    ) {

    pop <- "EUR"
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

# Required packages
packages <- c('data.table', 'dplyr')

for (p in packages) {
    if (!require(p, character.only = TRUE)) {
        install.packages(p, repos = "http://cran.us.r-project.org")
    }
}

RAP_outputs_folder <- destination
outputs <- system(paste("dx ls", RAP_outputs_folder), intern=TRUE)
outputs <- outputs[-grep("singleAssoc.txt.gz", outputs)]
# Ensure reruns are not counted as distinct
outputs <- unique(gsub(" : file-.*", "", outputs))
rerun_high_priority <- FALSE

dt <- data.table(
	chr=gsub("^(chr[0-9,X]+)_.*", "\\1", outputs),
	sex=ifelse(grepl("_F.txt", outputs), "F", ifelse(grepl("_M.txt", outputs), "M", "both")),
	pop=gsub(".*_([A-Z]{3})_*[A-Z]*.txt.*", "\\1", outputs),
	phenotype=gsub("^chr[0-9,X]+_(.*)_[0-9]+_([A-Z]{3})_*[A-Z]*.txt.*", "\\1", outputs),
	replicate=as.integer(gsub("^chr[0-9,X]+_.*_([0-9]+)_.*", "\\1", outputs))
)

priority_rerun <- "low"

# Create the equivalent data.table that would have been created if they had all completed
# Define the set of things to rerun as the set difference between the two.
dt_full <- unique(dt %>% select(-chr))
dt_check <- dt_full %>% mutate(chr="chrX")
for (chr in seq(2,22)) {
	dt_check <- rbind(dt_check, dt_full %>% mutate(chr=paste0("chr", chr)))
}

dt_rerun <- setdiff(dt_check, dt)

if (rerun_high_priority) {
	priority_rerun <- "high"
} else {
	priority_rerun <- "low"
}

# Run them all by extracting from the relevant entries of this data.table
for (row in 1:nrow(dt_rerun)) {
	cmd <- create_step2_sim_cmd(
		chr=gsub("chr", "", dt_rerun$chr[row]),
		phenotype=paste(dt_rerun$phenotype[row], dt_rerun$replicate[row], sep="_"),
        destination=destination,
        model_file_location=model_file_location,
		priority=priority_rerun
	)
	print(cmd)
	system(cmd)
}
