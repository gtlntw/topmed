#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD878008/NWD878008_HGDP_chr3.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD878008/NWD878008_HGDP_chr3.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD878008/NWD878008_HGDP_chr3.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD878008/NWD878008_HGDP_chr3 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD878008 chr 3