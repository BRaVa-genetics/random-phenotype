#!/usr/bin/env Rscript
library(data.table)
library(dplyr)

# Combine all files together to create a single file comprising all CAFs across 
# all superpopulation labellings.
directory <- "gene_CAFs"

dt_list <- list()
for (file in dir(directory, full.names=TRUE)) {
	# Extract population label
	dt_list[[file]] <- fread(file)
	dt_list[[file]][, super_population := gsub(".*\\.([A-Z]+)\\.chr[1-9,X]{1}[0-9]*.tsv.gz", "\\1", file)]
}

dt <- rbindlist(dt_list)
fwrite(dt, file="gene_CAFs/combined_gene_CAFs.tsv.gz", sep='\t')
