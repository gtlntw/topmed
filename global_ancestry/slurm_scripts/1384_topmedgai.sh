#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD228196/NWD228196_HGDP_chr20.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD228196/NWD228196_HGDP_chr20.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD228196/NWD228196_HGDP_chr20.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD228196/NWD228196_HGDP_chr20 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD228196 chr 20