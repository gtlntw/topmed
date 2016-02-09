#!/bin/bash
set pipefail; /net/snowwhite/home/khlin/bin/vcftools --remove-indels --remove-filtered-all --max-alleles 2 --gzvcf /net/topmed2/working/hmkang/snps/4318/filtered/chr6.filtered.rehdr.gt2.vcf.gz --gzdiff /net/topmed2/working/khlin/HGDP_938/HGDP_938_chr6_filtered_phased.vcf.gz --diff-site --stdout | awk '{if($4 == "B")  print $1 "\t" $2}' > /net/topmed2/working/khlin/common_site/topmed_chr6_HGDP_common_global.txt