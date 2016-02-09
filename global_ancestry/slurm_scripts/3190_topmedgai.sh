#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr22_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD256045/NWD256045_filtered_phased_chr22.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD256045/HGDP_938_chr22_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD256045/NWD256045_HGDP_subset_chr22_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD256045/NWD256045_HGDP_subset_chr22_filtered_phased.vcf.gz