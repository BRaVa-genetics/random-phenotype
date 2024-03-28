# Random phenotype simulation study
This repository contains functions and scripts required to generate random phenotypes under the null using the GRM as structured noise, run SAIGE, and plot the results.

First, construct the sparse GRM. This should have already been done along the path to running SAIGE for on the pilot phenotypes.

Next, generate random phenotypes under the null using `random_phenos.R`.

`Rscript random_phenos.R --n-reps 1 --grm ${EUR_sparseGRM} --samples ${EUR_sampleIDs} --out random_phenos_EUR`

Once this is done, run step 1 and step 2 using SAIGE on the resultant phenotypes.

The R scripts:
- `brava_simulation_step1.r`
- `brava_simulation_step2.r`
- `brava_simulation_step2_rerun_all_partially_completed.r`

provide examples of how to do this on the RAP using their submission system for UK Biobank data.

In the case of your cohort/biobank, you will need to create an equivalent submission script to loop over all of the phenotype data.

These are wrappers around calls to step 1 and 2 of SAIGE in [universal-saige](https://github.com/BRaVa-genetics/universal-saige).
Templates for these steps are in the [templates folder](https://github.com/BRaVa-genetics/universal-saige/tree/main/templates) of that repository.

Next, we wish to split the genes based on 'combined allele count', so that downstream we can see how lambdaGC varies as a function of expected combined allele count in cases.

`extract_CAF.r` provides a function to do this which you can call from the command line.

You need to provide the path to the allele count information (output from `plink --freqx`) and the vep and spliceAI annotated VCF file (the file created as output during variant annotation ending with `long.csv.gz` (Note that both of these files will have been generated when summarising the counts of classes of variation in the annotations which you can generate [here](https://github.com/BRaVa-genetics/BRaVa_curation/tree/main/QC/annotation_summary).

Finally, combine the resultant CAF files across ancestry labels and chromosomes. `combine_extracted_CAFs.r` does exactly this, combining across ancestry labellings for the files contained in a single folder named `gene_CAFs`. Note that this is a very naive script which assumes that only the required files are present in the `gene_CAFs` folder in your current working directory.

Once the results have been created, we can plot the results.

This is acheived by running 
`Rscript plot_simulation_study.r`

Note that the function here used to determine various parameters assumes an explicit filename structure which I used in the `brava_simulation_step2.r` above. Namely:

`chr${chr}_${phenotype}_${label}`

where ${phenotype} is the column name in the output of the call to `random_phenos.R`, and `${label}` is one of AFR, AMR, EAS, EUR, SAS.
...

