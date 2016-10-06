#!/bin/sh
toolDir="/net/snowwhite/home/khlin/tools/bcftools-1.3"

##extract the necessary snps
for chr in `seq 1 22`;
do
	#take out '-f PASS' so that snps that are overlap_indel and overlap_vntr are kept
	${toolDir}/bcftools view /net/fantasia/home/sayantan/DATABASE/TOPMED/REFERENCE_DATA/VCF/ALL.chr${chr}.freeze3a.pass.gtonly.genotypes.remove.missing.remove.monomorphic.Eagle.Phased.vcf.gz -T /net/topmed2/working/khlin/common_site/topmed_freeze3_chr${chr}_HGDP_common.txt --types snps -M2 --exclude-uncalled --output-type z --output-file /net/topmed2/working/khlin/topmed.freeze3.phased.hgdp.snp/topmed.freeze3.chr${chr}.phased.hgdp.snp.filtered.vcf.gz &
done


## tabix
# for chr in `seq 1 22`;
# do
# 	bcftools index -t -f /net/topmed2/working/khlin/topmed.freeze2.subset/topmed_freeze2_10597.chr${chr}.subset.filtered.vcf.gz &
# done

