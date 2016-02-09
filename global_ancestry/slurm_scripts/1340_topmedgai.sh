#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr20_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD228196/NWD228196_filtered_phased_chr20.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD228196/HGDP_938_chr20_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD228196/NWD228196_HGDP_subset_chr20_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD228196/NWD228196_HGDP_subset_chr20_filtered_phased.vcf.gz