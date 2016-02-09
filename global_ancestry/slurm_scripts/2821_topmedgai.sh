#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD985300/NWD985300_HGDP_chr5.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD985300/NWD985300_HGDP_chr5.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD985300/NWD985300_HGDP_chr5.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD985300/NWD985300_HGDP_chr5 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD985300 chr 5