#!/bin/bash
set pipefail; cd /net/snowwhite/home/khlin/tools/RFMix_v1.5.4 && python ./RunRFMix.py PopPhased /net/topmed2/working/khlin/global_ancestry/temp/NWD499673/NWD499673_HGDP_chr8.alleles /net/topmed2/working/khlin/global_ancestry/temp/NWD499673/NWD499673_HGDP_chr8.classes /net/topmed2/working/khlin/global_ancestry/temp/NWD499673/NWD499673_HGDP_chr8.locations -o /net/topmed2/working/khlin/global_ancestry/temp/NWD499673/NWD499673_HGDP_chr8 --forward-backward --num-threads 1 && Rscript /net/topmed2/working/khlin/global_ancestry/global_convert_class2plot.R id NWD499673 chr 8