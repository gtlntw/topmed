#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr6_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD607532/NWD607532_filtered_phased_chr6.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD607532/HGDP_938_chr6_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD607532/NWD607532_HGDP_subset_chr6_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD607532/NWD607532_HGDP_subset_chr6_filtered_phased.vcf.gz