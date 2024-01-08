# Restrict sample IDs to the Europeans
# Run the simulate phenotypes code
# Upload the phenotypes file to the RAP

# Important - the files being downloaded have been restricted
# to the collection of qced Europeans - the naming of the file 
# is misleading
RAP_location="/brava/inputs/GRM/"
EUR_sampleIDs="brava_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
EUR_sparseGRM="brava_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"

# Step 0 is already complete and generates the required input for the phenotypes to be simulated
# Download to interactive node and generate random phenotypes
dx download ${RAP_location}${EUR_sampleIDs}
dx download ${RAP_location}${EUR_sparseGRM}

# These inputs are used to generate null phenotypes
Rscript random_phenos.r --grm ${EUR_sparseGRM} --samples ${EUR_sampleIDs} --out random_phenos_EUR

# Submit SAIGE steps as jobs
Rscript brava_simulation_step1.r
Rscript brava_simulation_step2.r
Rscript brava_simulation_step2_rerun_all_phenotypes_partially_completed.r
