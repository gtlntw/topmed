#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr4_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD260039/NWD260039_filtered_phased_chr4.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD260039/HGDP_938_chr4_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD260039/NWD260039_HGDP_subset_chr4_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD260039/NWD260039_HGDP_subset_chr4_filtered_phased.vcf.gz