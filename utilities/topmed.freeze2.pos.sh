#!/bin/sh

## extract position
for chr in `seq 1 22`;
do
	bcftools query -f '%CHROM\t%POS\n' /net/topmed3/working/hmkang/freeze2/10597.v2/final/topmed_freeze2_10597.chr${chr}.overlap_removed.svm_pass.sites.vcf.gz > chr${chr}.freeze2.new.pos &
done
