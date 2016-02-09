#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD318475/NWD318475_HGDP_chr13.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD318475/NWD318475_HGDP_chr13.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD318475/NWD318475_HGDP_chr13.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD318475/NWD318475_HGDP_chr13 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD318475 chr 13