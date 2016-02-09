#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD490037/NWD490037_HGDP_chr1.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD490037/NWD490037_HGDP_chr1.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD490037/NWD490037_HGDP_chr1.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD490037/NWD490037_HGDP_chr1 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD490037 chr 1