#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr3_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD242497/NWD242497_filtered_phased_chr3.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/HGDP_938_chr3_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_subset_chr3_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_subset_chr3_filtered_phased.vcf.gz