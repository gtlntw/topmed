#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD243371/NWD243371_HGDP_chr4.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD243371/NWD243371_HGDP_chr4.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD243371/NWD243371_HGDP_chr4.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD243371/NWD243371_HGDP_chr4 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD243371 chr 4