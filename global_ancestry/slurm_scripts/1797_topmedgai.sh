#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr15_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD231341/NWD231341_filtered_phased_chr15.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD231341/HGDP_938_chr15_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD231341/NWD231341_HGDP_subset_chr15_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD231341/NWD231341_HGDP_subset_chr15_filtered_phased.vcf.gz