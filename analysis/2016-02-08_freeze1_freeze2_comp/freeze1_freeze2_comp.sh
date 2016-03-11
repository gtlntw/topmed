#!/bin/sh
chr=20
rm freeze1_freeze2_diff.txt
rm hgdp_freeze2_diff.txt
rm *.pos
rm *.snp


#generate data
for chr in {1..22}
do
	echo "processing chromosome ${chr}"
	bcftools query /net/topmed2/working/khlin/output/LAI/NWD100295/NWD100295_filtered_phased_chr${chr}.vcf.gz -f '%POS\n' |sort > freeze1_chr${chr}.pos
	bcftools query /net/topmed3/working/hmkang/freeze2/10597.v2/final/topmed_freeze2_10597.chr${chr}.overlap_removed.svm_pass.sites.vcf.gz -f '%POS\n' |sort > freeze2_chr${chr}.pos
	cut -f 2 /net/topmed2/working/khlin/common_site/topmed_chr${chr}_HGDP_common.txt | sort > hgdp_overlap_chr${chr}.sorted.pos
	#unique in freeze 1
	comm -23 freeze1_chr${chr}.pos freeze2_chr${chr}.pos | wc -l | awk '{print '${chr}'"\t"$1}'>> freeze1_freeze2_diff.txt
	#unique in hgdP
	comm -23 hgdp_overlap_chr${chr}.sorted.pos freeze2_chr${chr}.pos | wc -l | awk '{print '${chr}'"\t"$1}'>> hgdp_freeze2_diff.txt
done

echo "total number of snps used"
for chr in {1..22}
do
	echo "processing chromosome ${chr}"
	bcftools query /net/topmed2/working/khlin/output/LAI/NWD100295/NWD100295_filtered_phased_chr${chr}.vcf.gz -f '%POS\n' | wc -l | awk '{print $1}' >> freeze1.snp
	cut -f 2 /net/topmed2/working/khlin/common_site/topmed_chr${chr}_HGDP_common.txt | wc -l | awk '{print $1}' >> hgdp.snp
done

paste freeze1_freeze2_diff.txt freeze1.snp > freeze1_freeze2_diff_final.txt
paste hgdp_freeze2_diff.txt hgdp.snp > hgdp_freeze2_diff_final.txt

#freeze 1 SNPs used in the phasing step
awk 'BEGIN {SUM=0;SUM1=0} {SUM=SUM+$2;SUM1=SUM1+$3} END {print 0"\t"SUM"\t"SUM1}' freeze1_freeze2_diff_final.txt >> freeze1_freeze2_diff_final.txt
#freeze 1 SNPs used in the local ancestry calling part
awk 'BEGIN {SUM=0;SUM1=0} {SUM=SUM+$2;SUM1=SUM1+$3} END {print 0"\t"SUM"\t"SUM1}' hgdp_freeze2_diff_final.txt >> hgdp_freeze2_diff_final.txt

rm *.pos
rm *.snp
