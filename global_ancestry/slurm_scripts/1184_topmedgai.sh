#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD357867/NWD357867_HGDP_chr18.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD357867/NWD357867_HGDP_chr18.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD357867/NWD357867_HGDP_chr18.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD357867/NWD357867_HGDP_chr18 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD357867 chr 18