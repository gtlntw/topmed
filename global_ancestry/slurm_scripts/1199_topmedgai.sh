#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr11_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD151923/NWD151923_filtered_phased_chr11.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD151923/HGDP_938_chr11_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD151923/NWD151923_HGDP_subset_chr11_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD151923/NWD151923_HGDP_subset_chr11_filtered_phased.vcf.gz