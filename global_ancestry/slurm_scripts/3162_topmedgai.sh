#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD256485/NWD256485_HGDP_chr16.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD256485/NWD256485_HGDP_chr16.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD256485/NWD256485_HGDP_chr16.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD256485/NWD256485_HGDP_chr16 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD256485 chr 16