#!/bin/bash

for chr in `seq 1 22`;
do
	echo "processing chromosome: ${chr}"
	#create the snp list to keep based on the freeze 1 sample NWD1000295 to save time extracting unnecessary snps in the first step
	bcftools query -f '%CHROM\t%POS\n' ../output/LAI/NWD100295/NWD100295_filtered_phased_chr${chr}.vcf.gz > topmed_chr${chr}_1000g_common.txt

done
