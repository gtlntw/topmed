#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr13_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD913138/NWD913138_filtered_phased_chr13.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/HGDP_938_chr13_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_subset_chr13_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_subset_chr13_filtered_phased.vcf.gz