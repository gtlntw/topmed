#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD248269/NWD248269_HGDP_chr2.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD248269/NWD248269_HGDP_chr2.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD248269/NWD248269_HGDP_chr2.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD248269/NWD248269_HGDP_chr2 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD248269 chr 2