#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD282166/NWD282166_HGDP_chr15.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD282166/NWD282166_HGDP_chr15.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD282166/NWD282166_HGDP_chr15.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD282166/NWD282166_HGDP_chr15 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD282166 chr 15