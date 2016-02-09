#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD972722/NWD972722_HGDP_chr22.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD972722/NWD972722_HGDP_chr22.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD972722/NWD972722_HGDP_chr22.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD972722/NWD972722_HGDP_chr22 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD972722 chr 22