#!/bin/bash

# dx upload extract_CAF.r --path Duncan/counts/scripts/
rscript_remote="Duncan/counts/scripts/extract_CAF.r"

in_dir="/mnt/project/Duncan/long_annotations"
out_dir="/Duncan/gene_CAFs"
dx mkdir -p ${out_dir}

# Including spliceAI information in the annotation
for anc in AFR AMR EAS EUR SAS; do
   for CHR in {{1..22},X}; do
      out_prefix="ukb_wes_450k.${anc}.chr${CHR}"
      AC_path=/mnt/project/Duncan/counts/plink.frqx.chr${CHR}.${anc}.gz
      vep_processed_long_path="${in_dir}/ukb_wes_450k.july.qced.brava_common_rare.v7.chr${CHR}.saige.txt.long.csv.gz"
      dx run app-swiss-army-knife \
        -iimage_file="/docker/rsuite.tar.gz"\
        -icmd="
           Rscript /mnt/project/${rscript_remote} \
            --AC_path ${AC_path} \
            --vep_spliceAI_processed ${vep_processed_long_path} \
            --out ${out_prefix} \
          "\
        --instance-type mem1_ssd1_v2_x8 \
        --folder=".${out_dir}" \
        --priority normal \
        --name gene_CAFs -y
  done
done

# Download the required summary data, looping over chromosomes.
# Sum them up and create barplots of the results
# mkdir -p gene_CAFs 
# dx download /Duncan/gene_CAFs/* -f -o gene_CAFs/
# Combine the resultant files
# Rscript combine_extracted_CAFs.r