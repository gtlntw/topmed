#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr1_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD360786/NWD360786_filtered_phased_chr1.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD360786/HGDP_938_chr1_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD360786/NWD360786_HGDP_subset_chr1_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD360786/NWD360786_HGDP_subset_chr1_filtered_phased.vcf.gz