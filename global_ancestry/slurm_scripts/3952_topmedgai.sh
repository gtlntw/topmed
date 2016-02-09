#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD260039/NWD260039_HGDP_chr14.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD260039/NWD260039_HGDP_chr14.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD260039/NWD260039_HGDP_chr14.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD260039/NWD260039_HGDP_chr14 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD260039 chr 14