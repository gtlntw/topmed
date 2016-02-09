#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD606831/NWD606831_HGDP_chr7.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD606831/NWD606831_HGDP_chr7.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD606831/NWD606831_HGDP_chr7.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD606831/NWD606831_HGDP_chr7 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD606831 chr 7