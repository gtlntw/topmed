#!/bin/bash
set pipefail; cat /net/topmed2/working/khlin/global_ancestry/temp/NWD318475/NWD318475_filtered_chr10.bim | awk '{print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > /net/topmed2/working/khlin/global_ancestry/temp/NWD318475/NWD318475_filtered_chr10.bim.tmp && mv /net/topmed2/working/khlin/global_ancestry/temp/NWD318475/NWD318475_filtered_chr10.bim.tmp /net/topmed2/working/khlin/global_ancestry/temp/NWD318475/NWD318475_filtered_chr10.bim