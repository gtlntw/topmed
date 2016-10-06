#!/bin/bash

for chr in `seq 1 22`;
do
	echo "processing chromosome: ${chr}"
	#create the snp list to keep based on the freeze 1 sample NWD1000295 to save time extracting unnecessary snps in the first step
	cat ../output/LAI/NWD100295/NWD100295_HGDP_chr${chr}.pos > topmed_freeze1_chr${chr}_HGDP_common.txt
done

for chr in `seq 1 22`;
do
	echo "processing chromosome: ${chr}"
	#create the snp list to keep based on the freeze 1 sample NWD1000295 to save time extracting unnecessary snps in the first step
	cat ../output/LAI/NWD999941/NWD999941_HGDP_chr${chr}.pos > topmed_freeze2_chr${chr}_HGDP_common.txt
done

python ./2016-07-11\ intersection_freeze12.py

