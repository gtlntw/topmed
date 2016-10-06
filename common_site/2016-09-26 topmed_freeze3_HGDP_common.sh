#!/bin/bash

##extract snp information only and then remove multiallelic snps
# for chr in `seq 1 22`;
# do
# 	echo "processing chromosome: ${chr}"

# 	bcftools view /net/fantasia/home/sayantan/DATABASE/TOPMED/REFERENCE_DATA/VCF/ALL.chr${chr}.freeze3a.pass.gtonly.genotypes.remove.missing.remove.monomorphic.Eagle.Phased.vcf.gz -G --types snps -M2 --exclude-uncalled | bcftools norm --multiallelics +snps | bcftools view --types snps -M2 --exclude-uncalled --output-type z --output-file /net/topmed2/working/khlin/topmed.freeze3.phased.hgdp.snp/topmed.freeze3.chr${chr}.sites.snp.biallelic.vcf.gz &
# done


##compare with HGDP and extract the overlapped snps
for chr in `seq 1 22`;
do
	echo "processing chromosome: ${chr}"

	/net/snowwhite/home/khlin/bin/vcftools --gzvcf /net/topmed2/working/khlin/topmed.freeze3.phased.hgdp.snp/topmed.freeze3.chr${chr}.sites.snp.biallelic.vcf.gz --gzdiff /net/topmed2/working/khlin/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz --diff-site --stdout | awk '{if($4 == "B")  print $1 "\t" $2}' > /net/topmed2/working/khlin/common_site/topmed_freeze3_chr${chr}_HGDP_common.txt &
done