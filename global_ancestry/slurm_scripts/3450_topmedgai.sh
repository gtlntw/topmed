#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr18_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD516954/NWD516954_filtered_phased_chr18.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD516954/HGDP_938_chr18_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD516954/NWD516954_HGDP_subset_chr18_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD516954/NWD516954_HGDP_subset_chr18_filtered_phased.vcf.gz