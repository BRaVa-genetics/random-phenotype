#!/usr/bin/env Rscript

packages <- c('data.table', 'dplyr', 'ggplot2', 'tidyr', 'latex2exp', 'argparse', 'ggrastr')

for (p in packages) {
    if (!require(p, character.only = TRUE)) {
        install.packages(p, repos = "http://cran.us.r-project.org")
    }
}

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(argparse)
library(ggrastr)

# mkdir data
# dx download "/Duncan/simulation_study/outputs/step2/*AMR*" -o data/AMR/
annotation_names <- list(
    "damaging missense or protein altering" = "damaging_missense_or_protein_altering",
    "non-coding" = "non_coding", 
    "other missense or protein altering" = "other_missense_or_protein_altering",
    "pLoF" = "pLoF",
    "synonymous" = "synonymous",
    "pLoF; damaging missense or protein altering" = "pLoF;damaging_missense_or_protein_altering",
    "pLoF; damaging missense or protein altering;\nother missense or protein altering; synonymous" = "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous"
    )

extract_phenotype_name <- function(filename)
{
    if (grepl("singleAssoc", filename)) {
        gsub("^.*chr[0-9,X]+_(p_.*).txt.singleAssoc.txt.gz", "\\1", filename)
    } else {
        gsub("^.*chr[0-9,X]+_(p_.*).txt.gz", "\\1", filename)
    }
}

extract_phenotype_type <- function(filename)
{
    if (grepl("continuous", filename)) {
        return("continuous")
    } else {
        return("binary")
    }
}

extract_prevalence <- function(filename)
{
    if (grepl("continuous", filename)) {
        suppressWarnings(NA)
    } else {
        if (grepl("singleAssoc", filename)) {
            as.numeric(gsub("^.*chr[0-9,X]+_p_[0-9,\\.]+_([0-9,\\.]+).*.txt.singleAssoc.txt.gz", "\\1", filename))
        } else {
            as.numeric(gsub("^.*chr[0-9,X]+_p_[0-9,\\.]+_([0-9,\\.]+).*.txt.gz", "\\1", filename))
        }
    }	
}

extract_test_type <- function(filename)
{
    if (grepl("singleAssoc", filename)) {
        return("variant")
    } else {
        return("group")
    }
}

extract_population <- function(filename)
{
    if (grepl("singleAssoc", filename)) {
        gsub("^.*([A-Z]{3}).txt.singleAssoc.*", "\\1", filename)
    } else {
        gsub("^.*([A-Z]{3}).txt.gz", "\\1", filename)
    }
}

extract_run <- function(filename)
{
    if (grepl("singleAssoc", filename)) {
        as.integer(gsub("^.*_([0-9]+)+_[A-Z]{3}.txt.singleAssoc.*", "\\1",
            filename))
    } else {
        as.integer(gsub("^.*_([0-9]+)+_[A-Z]{3}.*", "\\1", filename))
    }
}

extract_file_information <- function(filename) {
    data.table(
        file = filename,
        phenotype = extract_phenotype_name(filename),
        type = extract_phenotype_type(filename),
        test = extract_test_type(filename),
        population = extract_population(filename),
        run = extract_run(filename),
        prevalence = extract_prevalence(filename)
    )
}

fread_and_prevalence <- function(filename)
{
    if (grepl("singleAssoc", filename)) {
        dt <- fread(filename, select=c("p.value", "Allele1", "AF_Allele2"),
            col.names=c("p.value", "Allele1", "MAF"))
        dt <- dt[Allele1 != "UR", ]
        dt[, MAF:=pmin(MAF, 1-MAF)]
        dt[, Allele1:=NULL]
    } else {
        dt <- fread(filename,
            select=c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "Region",
                "Group", "max_MAF"),
            col.names=c("SKAT-O", "Burden", "SKAT", "GENE",
                "annotation", "MAF_upper"))
        dt <- data.table(dt %>% pivot_longer(
            c("SKAT-O", "SKAT", "Burden"),
            names_to = "test", values_to = "p.value"))
    }
    dt[, prevalence:=extract_prevalence(filename)]
    dt[, prevalence:=ifelse(is.na(prevalence), "continuous", prevalence)]
    dt[, run:=extract_run(filename)]
    return(dt)
}

main <- function(args)
{
    dt_results <- rbindlist(
        lapply(dir(args$directory, full.name=TRUE),
            extract_file_information))
    if (!is.null(args$single_run)) {
        dt_results <- dt_results %>% filter(run == as.integer(args$single_run))
    }

    # Single variant tests
    dt_variant <- rbindlist(
        lapply(dt_results[test=="variant",]$file, fread_and_prevalence))
    dt_variant[, p.value := as.numeric(p.value)]
    dt_variant[, Interval := cut(
        MAF,
        breaks=c(0, 0.0001, 0.001, 0.01, 0.1, Inf),
        labels=c("<0.0001", "[0.0001, 0.001)", "[0.001, 0.01)", "[0.01, 0.1)",
            ">0.1"),
        include.lowest=TRUE)
    ]

    if (nrow(dt_variant) > 0) {
        if (is.null(args$single_run)) {
            # Find the minimum number of variants analysed across the runs
            dt_variant <- dt_variant %>% 
                group_by(Interval, prevalence, run) %>% 
                mutate(count=n()) %>% 
                group_by(Interval, prevalence) %>% 
                mutate(min_count = min(count), .groups="drop_last")
            dt_variant <- dt_variant %>% 
                group_by(Interval, prevalence, run) %>%
                reframe(quantile_p = quantile(p.value,
                            probs=seq(0,1,length.out=(min_count[1]+1)[-(min_count[1]+1)])),
                ) %>% group_by(Interval, prevalence, run) %>%
                mutate(order = order(quantile_p))
            dt_variant <- dt_variant %>% 
                group_by(Interval, prevalence, order) %>% 
                summarise(p.value = median(quantile_p), .groups="drop_last") %>% 
                mutate(p.value.expected = seq(1,n())/n())
        } else {
            dt_variant <- dt_variant[order(p.value)] %>% 
                group_by(Interval, prevalence) %>% 
                mutate(p.value.expected = seq(1,n())/n(), .groups="drop_last")
        }

        cat("Creating variant plots...\n")
        pdf(file=gsub(".pdf", "_variant.pdf", args$out))
        dt_variant <- dt_variant %>% mutate(p.value = as.numeric(p.value)+1e-300)
        p <- ggplot(dt_variant,
            aes(x=-log10(p.value.expected), y=-log10(p.value), col=Interval)) +
        rasterise(geom_point(), dpi=500) + geom_abline(intercept=0, slope=1, color='black') + 
        facet_wrap(~prevalence, scales="free") + theme_classic() + 
        xlab(TeX("$-\\log_{10}(P_{observed})$")) +
        ylab(TeX("$-\\log_{10}(P_{expected})$")) + 
        labs(title = "Variant based testing") #+
        # coord_fixed(ratio = 1)
        print(p)
        dev.off()
        cat("Variant plots created.\n")
    }

    # Group tests
    CAFs_master <- args$CAF_file
    dt_CAF <- fread(CAFs_master)[super_population == args$pop,]
    setkeyv(dt_CAF, c("GENE", "annotation", "MAF_upper"))
    dt_group <- rbindlist(
        lapply(dt_results[test=="group",]$file, fread_and_prevalence))[
            annotation != "Cauchy"
        ]

    # Determine the collection of combinations of annotations
    dt_group[, annotation:=factor(annotation)]
    combined_annotations <- grep(";", levels(dt_group$annotation), value=TRUE)

    dt_CAF_to_add <- list()
    for (c in combined_annotations)
    {
        dt_CAF_to_add[[c]] <- dt_CAF %>% 
            filter(annotation %in% strsplit(c, split=";")[[1]]) %>%
            group_by(GENE, MAF_upper) %>% 
            summarise(
                CAF = sum(CAF),
                CAC=sum(CAC),
                n_variants=sum(n_variants), .groups="drop_last"
                ) %>% mutate(annotation=c, super_population=args$pop)
    }
    dt_CAF_to_add <- rbindlist(dt_CAF_to_add)
    dt_CAF <- rbind(dt_CAF, dt_CAF_to_add)

    setkeyv(dt_group, c("GENE", "annotation", "MAF_upper"))
    dt_group <- merge(dt_group, dt_CAF)[,
        c("CAC", "n_variants", "super_population") := NULL]
    dt_group[, Interval := cut(
        CAF,
        breaks=c(0, 0.0001, 0.001, 0.01, 0.1, Inf),
        labels=c("<0.0001", "[0.0001, 0.001)", "[0.001, 0.01)", "[0.01, 0.1)",
            ">0.1"),
        include.lowest=TRUE)
    ]
    dt_group[, annotation := factor(annotation)]
    levels(dt_group$annotation) <- annotation_names

    if (is.null(args$single_run)) {
        # Find the minimum number of variants analysed across the runs
        dt_group_split <- dt_group %>% 
            group_by(Interval, prevalence, MAF_upper, annotation, test, run) %>% 
            mutate(count=n()) %>% 
            group_by(Interval, prevalence, MAF_upper, annotation, test) %>% 
            mutate(min_count = min(count))
        dt_group_split <- dt_group_split %>% 
            group_by(Interval, prevalence, MAF_upper, annotation, test, run) %>%
            reframe(quantile_p = quantile(as.numeric(p.value),
                probs=seq(0,1,length.out=(min_count[1]+1)[-(min_count[1]+1)]))) %>%
            group_by(Interval, prevalence, MAF_upper, annotation, test, run) %>%
            mutate(order = order(quantile_p))
        dt_group_split <- dt_group_split %>% 
            group_by(Interval, prevalence, MAF_upper, annotation, test, order) %>% 
            summarise(p.value = median(quantile_p), .groups="drop_last") %>% 
            mutate(p.value.expected = seq(1,n())/n())
        dt_group_split <- data.table(dt_group_split)
    } else {
        dt_group_split <- dt_group[order(as.numeric(p.value))] %>% 
        group_by(Interval, prevalence, MAF_upper, annotation, test) %>% 
        mutate(p.value.expected = seq(1,n())/n(), .groups="drop_last")
        dt_group_split <- data.table(dt_group_split)
    }

    if (is.null(args$single_run)) {
        # Find the minimum number of variants analysed across the runs
        dt_group <- dt_group %>% 
            group_by(Interval, prevalence, MAF_upper, test, run) %>% 
            mutate(count=n()) %>% 
            group_by(Interval, prevalence, MAF_upper, test) %>% 
            mutate(min_count = min(count))
        dt_group <- dt_group %>% 
            group_by(Interval, prevalence, MAF_upper, test, run) %>%
            reframe(quantile_p = quantile(as.numeric(p.value),
                probs=seq(0,1,length.out=(min_count[1]+1)[-(min_count[1]+1)]))) %>%
            group_by(Interval, prevalence, MAF_upper, test, run) %>%
            mutate(order = order(quantile_p))
        dt_group <- dt_group %>% 
            group_by(Interval, prevalence, MAF_upper, test, order) %>% 
            summarise(p.value = median(quantile_p), .groups="drop_last") %>% 
            mutate(p.value.expected = seq(1,n())/n())
        dt_group <- data.table(dt_group)
    } else {
        dt_group <- dt_group[order(as.numeric(p.value))] %>% 
            group_by(Interval, prevalence, MAF_upper, test) %>% 
            mutate(p.value.expected = seq(1,n())/n(), .groups="drop_last")
        dt_group <- data.table(dt_group)
    }

    cat("Creating gene plots...\n")
    pdf(file=gsub(".pdf", "_group.pdf", args$out), width=5, height=10)
    # At the level of tests and annotations
    for (t in unique(dt_group_split$test)) {
        cat(paste0("test: ",t, "\n"))
        for (ann in unique(dt_group_split$annotation)) {
            cat(paste0("annotation: ", ann, "\n"))
            p <- ggplot(dt_group_split[(annotation==ann & test==t),],
                aes(x=-log10(p.value.expected),
                    y=-log10(p.value), col=Interval)
                ) +
            rasterise(geom_point(), dpi=500) + 
            geom_abline(intercept=0, slope=1, color='black') + 
            facet_grid(prevalence~MAF_upper, scales="free") + 
            geom_abline(intercept=0, slope=1, color='black') + 
            theme_classic() +
            xlab(TeX("$-\\log_{10}(P_{observed})$")) +
            ylab(TeX("$-\\log_{10}(P_{expected})$")) + 
            labs(title = "Gene based testing",
                subtitle = paste0(t, ", ", ann)) #+
            # coord_fixed(ratio = 1)
            print(p)
        }
    }

    # At the level of tests (combine over annotations)
    for (max_MAF in unique(dt_group$MAF_upper)) {
        cat(paste0("max MAF: ", max_MAF, "\n"))
        p <- ggplot(dt_group[MAF_upper==max_MAF],
            aes(x=-log10(p.value.expected), y=-log10(p.value), col=Interval)) +
        rasterise(geom_point(), dpi=500) + geom_abline(intercept=0, slope=1, color='black') + 
        facet_grid(prevalence~test, scales="free") + 
        geom_abline(intercept=0, slope=1, color='black') + 
        theme_classic() +
        xlab(TeX("$-\\log_{10}(P_{observed})$")) +
        ylab(TeX("$-\\log_{10}(P_{expected})$")) +
        labs(title = "Gene based testing",
            subtitle = paste("max MAF:", max_MAF)) #+
        # coord_fixed(ratio = 1)
        print(p)
    }
    dev.off()
    cat("Gene plots created.\n")
    print(warnings())
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--directory", default="data", required=FALSE,
    help=paste0("Path to directory containing step 2 output split by",
        " chromosome to combine"))
parser$add_argument("--CAF_file",
    default="gene_CAFs/combined_gene_CAFs.tsv.gz", required=FALSE,
    help=paste0("Path to file containing the combined allele frequency split ",
    " by annotation"))
parser$add_argument("--out", default="simulation_plots.pdf", required=FALSE,
    help="Output filepath")
parser$add_argument("--single_run", default=NULL, required=FALSE,
    help=paste0("Pass an integer for the required run to plot"))
parser$add_argument("--pop", default="EUR", required=FALSE,
    help=paste0("Population label used for combined allele frequency binning"))
args <- parser$parse_args()

main(args)

# To generate all of the plots, run
# Rscript plot_simulation_study.r
# Rscript plot_simulation_study.r --directory data/AFR --out simulation_plots_no_covariates_AFR.pdf --pop AFR --single_run 1
Rscript plot_simulation_study.r --directory data/AMR --out simulation_plots_no_covariates_AMR.pdf --pop AMR --single_run 1
# Rscript plot_simulation_study.r --directory data/EAS --out simulation_plots_no_covariates_EAS.pdf --pop EAS --single_run 1
# Rscript plot_simulation_study.r --directory data/EUR --out simulation_plots_no_covariates_EUR.pdf --pop EUR --single_run 1
# Rscript plot_simulation_study.r --directory data/SAS --out simulation_plots_no_covariates_SAS.pdf --pop SAS --single_run 1
