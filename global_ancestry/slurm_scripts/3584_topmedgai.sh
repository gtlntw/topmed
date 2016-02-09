#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr20_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD259800/NWD259800_filtered_phased_chr20.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD259800/HGDP_938_chr20_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD259800/NWD259800_HGDP_subset_chr20_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD259800/NWD259800_HGDP_subset_chr20_filtered_phased.vcf.gz