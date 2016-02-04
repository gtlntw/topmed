#!/bin/sh

##extract the necessary snps
# for chr in `seq 1 22`;
# do
# 	bcftools view /net/topmed3/working/hmkang/freeze2/10597.v2/final/topmed_freeze2_10597.chr${chr}.overlap_removed.svm_pass.genotypes.vcf.gz -T /net/topmed2/working/khlin/common_site/topmed_chr${chr}_1000g_common.txt --types snps -M2 --exclude-uncalled -f PASS --output-type z --output-file /net/topmed2/working/khlin/topmed.freeze2.subset/topmed_freeze2_10597.chr${chr}.subset.filtered.vcf.gz &
# done


#tabix
for chr in `seq 1 22`;
do
	bcftools index -t -f /net/topmed2/working/khlin/topmed.freeze2.subset/topmed_freeze2_10597.chr${chr}.subset.filtered.vcf.gz &
done
