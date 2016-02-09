#!/bin/sh

for chr in {1..22}
do
	echo "processing chromosome ${chr}"
	zcat /net/topmed2/working/khlin/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz | awk '{print $2}' | sort > 1000g_chr${chr}.pos
	bcftools query /net/topmed3/working/hmkang/freeze2/10597.v2/final/topmed_freeze2_10597.chr${chr}.overlap_removed.svm_pass.sites.vcf.gz -f '%POS\n' | sort > freeze2_chr${chr}.pos
	#common in both
	comm -12 1000g_chr${chr}.pos freeze2_chr${chr}.pos | sort -n |awk '{print '${chr}'"\t"$1}' > common_site/freeze2_1000g_chr${chr}_common.txt
done

echo "total number of snps used"
for chr in {1..22}
do
	echo "processing chromosome ${chr}"
	cut -f 2 common_site/freeze2_1000g_chr${chr}_common.txt | wc -l | awk '{print '${chr}'"\t"$1}' >> freeze2_1000g_common.txt
done

awk 'BEGIN {SUM=0} {SUM=SUM+$2} END {print 0"\t"SUM}' freeze2_1000g_common.txt >> freeze2_1000g_common.txt


rm *.pos
rm *.snp