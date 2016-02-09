#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD259800/NWD259800_HGDP_chr21.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD259800/NWD259800_HGDP_chr21.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD259800/NWD259800_HGDP_chr21.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD259800/NWD259800_HGDP_chr21 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD259800 chr 21