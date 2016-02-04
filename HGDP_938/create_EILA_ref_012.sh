#!/bin/bash
while getopts s: opt; do
  case $opt in
  s)
      subset=$OPTARG
      ;;
  esac
done

topmedDir="/net/snowwhite/home/khlin/topmed"
targetDir="$topmedDir/HGDP_938/EILA"


for chr in `seq 20 22`;
do
	#europe
	vcftools --gzvcf $topmedDir/HGDP_938/HGDP_938_chr${chr}_phased.vcf.gz --chr ${chr} --positions $topmedDir/common_site/1000g_chr${chr}_HGDP_common.txt --keep $topmedDir/HGDP_938/HGDP_europe.txt --012 --out ${targetDir}/HGDP_europe_chr${chr} &&  \
	awk -f /net/snowwhite/home/khlin/tools/transpose_withSpace.awk < ${targetDir}/HGDP_europe_chr${chr}.012 | tail -n +2 > ${targetDir}/HGDP_europe_chr${chr}.012.tmp && \
    mv ${targetDir}/HGDP_europe_chr${chr}.012.tmp ${targetDir}/HGDP_europe_chr${chr}.012
	
	#africa
	vcftools --gzvcf $topmedDir/HGDP_938/HGDP_938_chr${chr}_phased.vcf.gz --chr ${chr} --positions $topmedDir/common_site/1000g_chr${chr}_HGDP_common.txt --keep $topmedDir/HGDP_938/HGDP_africa.txt --012 --out ${targetDir}/HGDP_africa_chr${chr} &&  \
	awk -f /net/snowwhite/home/khlin/tools/transpose_withSpace.awk < ${targetDir}/HGDP_africa_chr${chr}.012 | tail -n +2 > ${targetDir}/HGDP_africa_chr${chr}.012.tmp && \
    mv ${targetDir}/HGDP_africa_chr${chr}.012.tmp ${targetDir}/HGDP_africa_chr${chr}.012

	#native_america
	vcftools --gzvcf $topmedDir/HGDP_938/HGDP_938_chr${chr}_phased.vcf.gz --chr ${chr} --positions $topmedDir/common_site/1000g_chr${chr}_HGDP_common.txt --keep $topmedDir/HGDP_938/HGDP_native_america.txt --012 --out ${targetDir}/HGDP_native_america_chr${chr} &&  \
	awk -f /net/snowwhite/home/khlin/tools/transpose_withSpace.awk < ${targetDir}/HGDP_native_america_chr${chr}.012 | tail -n +2 > ${targetDir}/HGDP_native_america_chr${chr}.012.tmp && \
    mv ${targetDir}/HGDP_native_america_chr${chr}.012.tmp ${targetDir}/HGDP_native_america_chr${chr}.012
done

rm -rf ${targetDir}/*.log ${targetDir}/*.pos ${targetDir}/*.indv
