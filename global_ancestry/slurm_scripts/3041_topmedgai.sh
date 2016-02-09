#!/bin/bash
set pipefail; bcftools merge -O u -R /net/topmed2/working/khlin/common_site/topmed_chr5_HGDP_common_global.txt /net/topmed2/working/khlin/output/LAI/NWD368534/NWD368534_filtered_phased_chr5.vcf.gz /net/topmed2/working/khlin/global_ancestry/temp/NWD368534/HGDP_938_chr5_subset_filtered_phased.vcf.gz | bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o /net/topmed2/working/khlin/global_ancestry/temp/NWD368534/NWD368534_HGDP_subset_chr5_filtered_phased.vcf.gz && bcftools index -t -f /net/topmed2/working/khlin/global_ancestry/temp/NWD368534/NWD368534_HGDP_subset_chr5_filtered_phased.vcf.gz