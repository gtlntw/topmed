#!/bin/bash
set pipefail; bcftools view /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_subset_chr4_filtered_phased.vcf.gz | bcftools query -f '[%GT]\n' - | sed 's/|//g' > /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_chr4.alleles && bcftools view /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_subset_chr4_filtered_phased.vcf.gz | bcftools query -f '%POS\n' - | Rscript /net/topmed2/working/khlin/utilities/generate_markerLocations_file.R stdin /net/topmed2/working/khlin/genetic_map_GRCh37/genetic_map_chr4_combined_b37.txt /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_chr4.locations && bcftools query -l /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_subset_chr4_filtered_phased.vcf.gz | Rscript /net/topmed2/working/khlin/utilities/make_classes_file.R stdin /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_chr4.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_ref.txt && bcftools query -f '%CHROM	%POS\n' /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_subset_chr4_filtered_phased.vcf.gz > /net/topmed2/working/khlin/global_ancestry/temp/NWD242497/NWD242497_HGDP_chr4.pos