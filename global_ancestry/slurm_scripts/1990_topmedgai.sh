#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr10_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD241576/NWD241576_filtered_phased_chr10.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD241576/HGDP_938_chr10_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD241576/NWD241576_HGDP_subset_chr10_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD241576/NWD241576_HGDP_subset_chr10_filtered_phased.vcf.gz