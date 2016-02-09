#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr17_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD373188/NWD373188_filtered_phased_chr17.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD373188/HGDP_938_chr17_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD373188/NWD373188_HGDP_subset_chr17_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD373188/NWD373188_HGDP_subset_chr17_filtered_phased.vcf.gz