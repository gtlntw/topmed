#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD229074/NWD229074_HGDP_chr12.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD229074/NWD229074_HGDP_chr12.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD229074/NWD229074_HGDP_chr12.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD229074/NWD229074_HGDP_chr12 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD229074 chr 12