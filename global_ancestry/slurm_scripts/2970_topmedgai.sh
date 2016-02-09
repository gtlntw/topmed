#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD249716/NWD249716_HGDP_chr22.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD249716/NWD249716_HGDP_chr22.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD249716/NWD249716_HGDP_chr22.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD249716/NWD249716_HGDP_chr22 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD249716 chr 22