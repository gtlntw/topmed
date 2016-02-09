#!/bin/bash
set pipefail; bcftools view /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_subset_chr17_filtered_phased.vcf.gz | bcftools query -f '[%GT]\n' - | sed 's/|//g' > /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_chr17.alleles && bcftools view /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_subset_chr17_filtered_phased.vcf.gz | bcftools query -f '%POS\n' - | Rscript /net/topmed2/working/khlin/utilities/generate_markerLocations_file.R stdin /net/topmed2/working/khlin/genetic_map_GRCh37/genetic_map_chr17_combined_b37.txt /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_chr17.locations && bcftools query -l /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_subset_chr17_filtered_phased.vcf.gz | Rscript /net/topmed2/working/khlin/utilities/make_classes_file.R stdin /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_chr17.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_ref.txt && bcftools query -f '%CHROM	%POS\n' /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_subset_chr17_filtered_phased.vcf.gz > /net/topmed2/working/khlin/global_ancestry/temp/NWD913138/NWD913138_HGDP_chr17.pos