#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr16_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD985300/NWD985300_filtered_phased_chr16.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD985300/HGDP_938_chr16_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD985300/NWD985300_HGDP_subset_chr16_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD985300/NWD985300_HGDP_subset_chr16_filtered_phased.vcf.gz