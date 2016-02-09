#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD182356/NWD182356_HGDP_chr19.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD182356/NWD182356_HGDP_chr19.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD182356/NWD182356_HGDP_chr19.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD182356/NWD182356_HGDP_chr19 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD182356 chr 19