# Restrict sample IDs to the Europeans
# Run the simulate phenotypes code
# Upload the phenotypes file to the RAP

RAP_location="/brava/outputs/step0/"
pops=("AFR" "AMR" "EAS" "EUR" "SAS")
pops=("AFR" "AMR" "EAS" "SAS")
filter_samples=true

for pop in ${pops[@]}; do
	# Sample IDs and sparse GRM filenames
	sampleIDs="brava_${pop}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
	sparseGRM="brava_${pop}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
	# Step 0 is already complete and generates the required input for the phenotypes to be simulated
	# Download to interactive node and generate random phenotypes
	# Note that before this step, you need to ensure that dx has been installed (pip3 install dxpy)
	# and that you have run dx login, linked to the appropriate project which contains that sparse GRMs
	dx download -f ${RAP_location}${sampleIDs}
	dx download -f ${RAP_location}${sparseGRM}
done

if $filter_samples
then
	Rscript filter_to_cleaned_samples.r
fi

for pop in ${pops[@]}; do
	# These inputs are used to generate null phenotypes
	if $filter_samples
	then
		sampleIDs="brava_${pop}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.filtered.mtx.sampleIDs.txt"
		sparseGRM="brava_${pop}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.filtered.mtx"
		echo "Rscript random_phenos.R --grm ${sparseGRM} --samples ${sampleIDs} --out random_phenos_${pop} --n_reps 1"
		Rscript random_phenos.R --grm ${sparseGRM} --samples ${sampleIDs} --out random_phenos_${pop} --n_reps 1
	else
		sampleIDs="brava_${pop}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
		sparseGRM="brava_${pop}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
		echo "Rscript random_phenos.R --grm ${sparseGRM} --samples ${sampleIDs} --out random_phenos_${pop} --n_reps 1"
		Rscript random_phenos.R --grm ${sparseGRM} --samples ${sampleIDs} --out random_phenos_${pop} --n_reps 1
	fi

	# Upload the result
	dx upload random_phenos_${pop}.tsv.gz --path /Duncan/simulation_study/

	# Submit SAIGE steps as jobs
	# Rscript brava_simulation_step1.r
	# Rscript brava_simulation_step2.r
	# Rscript brava_simulation_step2_rerun_all_phenotypes_partially_completed.r
done

