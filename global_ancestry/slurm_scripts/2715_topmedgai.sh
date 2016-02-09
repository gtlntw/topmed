#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr9_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD248269/NWD248269_filtered_phased_chr9.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD248269/HGDP_938_chr9_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD248269/NWD248269_HGDP_subset_chr9_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD248269/NWD248269_HGDP_subset_chr9_filtered_phased.vcf.gz