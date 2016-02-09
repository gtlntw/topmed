#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr14_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD245359/NWD245359_filtered_phased_chr14.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD245359/HGDP_938_chr14_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD245359/NWD245359_HGDP_subset_chr14_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD245359/NWD245359_HGDP_subset_chr14_filtered_phased.vcf.gz