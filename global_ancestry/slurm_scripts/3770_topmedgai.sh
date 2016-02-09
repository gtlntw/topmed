#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr8_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD260830/NWD260830_filtered_phased_chr8.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD260830/HGDP_938_chr8_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD260830/NWD260830_HGDP_subset_chr8_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD260830/NWD260830_HGDP_subset_chr8_filtered_phased.vcf.gz