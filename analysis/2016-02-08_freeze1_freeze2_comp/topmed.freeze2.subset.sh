#!/bin/sh

#extract the necessary snps
for chr in `seq 1 22`;
do
	#take out '-f PASS' so that snps that are overlap_indel and overlap_vntr are kept
	bcftools view /net/topmed3/working/hmkang/freeze2/10597.v2/final/topmed_freeze2_10597.chr${chr}.overlap_removed.svm_pass.genotypes.vcf.gz --types snps -M2 --exclude-uncalled -f PASS -T common_site/freeze2_1000g_chr${chr}_common.txt --output-type z --output-file topmed.freeze2.subset/topmed_freeze2_10597.chr${chr}.subset.filtered.vcf.gz &
done


# ## tabix
# for chr in `seq 1 22`;
# do
# 	bcftools index -t -f /net/topmed2/working/khlin/topmed.freeze2.subset/topmed_freeze2_10597.chr${chr}.subset.filtered.vcf.gz &
# done
