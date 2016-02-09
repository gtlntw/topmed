#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD360786/NWD360786_HGDP_chr9.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD360786/NWD360786_HGDP_chr9.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD360786/NWD360786_HGDP_chr9.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD360786/NWD360786_HGDP_chr9 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD360786 chr 9