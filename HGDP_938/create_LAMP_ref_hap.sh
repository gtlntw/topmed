#!/bin/bash
while getopts s: opt; do
  case $opt in
  s)
      subset=$OPTARG
      ;;
  esac
done

topmed_dr="/net/topmed2/working/khlin"
target_dr="$topmed_dr/HGDP_938/LAMPLD"


for chr in `seq 1 22`;
do
	#europe
	vcftools --gzvcf $topmed_dr/HGDP_938/HGDP_938_chr${chr}_phased.vcf.gz --chr ${chr} --positions $topmed_dr/common_site/topmed_chr${chr}_HGDP_subset_common${subset}.txt --keep $topmed_dr/HGDP_938/HGDP_europe.txt --IMPUTE --out ${target_dr}/HGDP_europe_chr${chr} &&  \
	rm -f ${target_dr}/HGDP_europe_chr${chr}.impute.hap.indv ${target_dr}/HGDP_europe_chr${chr}.impute.legend ${target_dr}/HGDP_europe_chr${chr}.log
	awk -f /net/snowwhite/home/khlin/tools/transpose.awk < ${target_dr}/HGDP_europe_chr${chr}.impute.hap > ${target_dr}/HGDP_europe_chr${chr}.hap.tmp && mv ${target_dr}/HGDP_europe_chr${chr}.hap.tmp ${target_dr}/HGDP_europe_chr${chr}.impute.hap
	#africa
	vcftools --gzvcf $topmed_dr/HGDP_938/HGDP_938_chr${chr}_phased.vcf.gz --chr ${chr} --positions $topmed_dr/common_site/topmed_chr${chr}_HGDP_subset_common${subset}.txt --keep $topmed_dr/HGDP_938/HGDP_africa.txt --IMPUTE --out ${target_dr}/HGDP_africa_chr${chr} &&  \
	rm -f ${target_dr}/HGDP_africa_chr${chr}.impute.hap.indv ${target_dr}/HGDP_africa_chr${chr}.impute.legend ${target_dr}/HGDP_africa_chr${chr}.log
	awk -f /net/snowwhite/home/khlin/tools/transpose.awk < ${target_dr}/HGDP_africa_chr${chr}.impute.hap > ${target_dr}/HGDP_africa_chr${chr}.hap.tmp && mv ${target_dr}/HGDP_africa_chr${chr}.hap.tmp ${target_dr}/HGDP_africa_chr${chr}.impute.hap

	#native america
	vcftools --gzvcf $topmed_dr/HGDP_938/HGDP_938_chr${chr}_phased.vcf.gz --chr ${chr} --positions $topmed_dr/common_site/topmed_chr${chr}_HGDP_subset_common${subset}.txt --keep $topmed_dr/HGDP_938/HGDP_native_america.txt --IMPUTE --out ${target_dr}/HGDP_native_america_chr${chr} &&  \
	rm -f ${target_dr}/HGDP_native_america_chr${chr}.impute.hap.indv ${target_dr}/HGDP_native_america_chr${chr}.impute.legend ${target_dr}/HGDP_native_america_chr${chr}.log
	awk -f /net/snowwhite/home/khlin/tools/transpose.awk < ${target_dr}/HGDP_native_america_chr${chr}.impute.hap > ${target_dr}/HGDP_native_america_chr${chr}.hap.tmp && mv ${target_dr}/HGDP_native_america_chr${chr}.hap.tmp ${target_dr}/HGDP_native_america_chr${chr}.impute.hap

	#pos file
	awk '{print $2}' $topmed_dr/common_site/topmed_chr${chr}_HGDP_subset_common${subset}.txt > ${target_dr}/chr${chr}.pos

done
