#!/usr/bin/env Rscript
library(data.table)
library(argparse)
library(dplyr)

main <- function(args)
{
    create_gene_annotation_CAFs <- function(
        dt, annotation_group, gene_group)
    {
        annotation_group <- enquo(annotation_group)
        gene_group <- enquo(gene_group)

        dt_CAF <- dt %>% group_by(!!annotation_group, !!gene_group) %>% 
            summarise(CAF = sum(MAF), CAC = sum(MAC), n_variants=n())
        dt_CAF <- data.table(dt_CAF)
        return(dt_CAF)
    }

    AC_path <- args$AC_path
    vep_spliceAI_path <- args$vep_spliceAI_processed
    
    dt_brava_annot <- fread(vep_spliceAI_path, key=c("SNP_ID"))
    if (!all(dt_brava_annot$annotation %in% c(
        "pLoF",
        "damaging_missense_or_protein_altering",
        "other_missense_or_protein_altering",
        "synonymous",
        "non_coding"))) {
        stop("An annotation is present in the file which is not in the following:
            - pLoF
            - damaging_missense_or_protein_altering
            - other_missense_or_protein_altering
            - synonymous
            - non_coding")
    }
    if (all(is.na(dt_brava_annot$max_DS))) {
        stop("splice AI information has not been incorporated.
            All entries of the max_DS column are empty. 
            Check variant ID merging between spliceAI and VEP information. 
            Perhaps chrCHR:POS:REF:ALT vs CHR:POS:REF:ALT.")
    }
    dt_brava_annot[, SNP_ID := gsub("chr", "", SNP_ID)]
    dt_AC <- fread(AC_path)

    dt_AC[, SNP_ID := gsub("chr", "", SNP)]
    setkey(dt_AC, "SNP_ID")

    # Check to ensure that all SNP_IDs are off the form CHR:POS:REF:ALT
    AC_correct_format <- all(
        grepl("^[1-9,X]{1}[0-9]*:[0-9]+:[A,C,G,T]+:[A,C,G,T]+", dt_AC$SNP_ID))
    brava_annot_correct_format <- all(
        grepl("^[1-9,X]{1}[0-9]*:[0-9]+:[A,C,G,T]+:[A,C,G,T]+", dt_brava_annot$SNP_ID)
        )
    if (!AC_correct_format | !brava_annot_correct_format) {
        stop(paste("One or both SNP/SNP_ID columns in", AC_path,
            "and", vep_spliceAI_path, "is not formatted to CHR:POS:REF:ALT"))
    }

    dt_AC[, AC_A1 := 2*`C(HOM A1)` + `C(HET)`]
    dt_AC[, AC_A2 := 2*`C(HOM A2)` + `C(HET)`]
    dt_AC[, MAC := pmin(AC_A1, AC_A2)]
    dt_AC[, check := `C(HOM A1)` + `C(HOM A2)` + `C(HET)` + `C(MISSING)`]
    dt_AC <- dt_AC[MAC > 0]
    n_samples <- dt_AC$check[1]
    dt_AC[, MAF := MAC/(AC_A1 + AC_A2)]

    # Filter based on MAF (given that gene-based tests will be MAF dependent)
    dt_CAF_list <- list()
    for (MAF_cutoff in c(1e-4, 1e-3, 0.01)) {
        dt_AC_tmp <- dt_AC[MAF < MAF_cutoff,]
        dt_variant_gene_AC_tmp <- merge(dt_AC_tmp, dt_brava_annot)
        dt_CAF_list[[as.character(MAF_cutoff)]] <- create_gene_annotation_CAFs(dt_variant_gene_AC_tmp, annotation, GENE)
        dt_CAF_list[[as.character(MAF_cutoff)]][, MAF_upper := MAF_cutoff]
    }
    dt_CAF <- rbindlist(dt_CAF_list)

    # Ensure the output file is gzipped .tsv
    out <- gsub(".tsv.gz$", "", args$out)
    out <- gsub(".gz$", "", args$out)
    out <- paste0(out, ".tsv.gz")

    fwrite(dt_CAF, file=out, sep="\t", quote=FALSE)
    cat("created gene CAF split by annotation group\n")
    return(dt_CAF)
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--AC_path", default=NULL, required=TRUE,
    help="Path to allele count information, output from plink --freqx")
parser$add_argument("--vep_spliceAI_processed", default=NULL, required=TRUE,
    help=paste0("Path to the processed VEP file with spliceAI information ",
        "included. This is the 'long' output of brava_create_annot.py"))
parser$add_argument("--out", default=NULL, required=TRUE,
    help="Output filepath")
args <- parser$parse_args()

main(args)
