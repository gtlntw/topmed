#!/bin/bash
set pipefail; cat /net/topmed2/working/khlin/global_ancestry/temp/NWD373188/NWD373188_filtered_chr6.bim | awk '{print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > /net/topmed2/working/khlin/global_ancestry/temp/NWD373188/NWD373188_filtered_chr6.bim.tmp && mv /net/topmed2/working/khlin/global_ancestry/temp/NWD373188/NWD373188_filtered_chr6.bim.tmp /net/topmed2/working/khlin/global_ancestry/temp/NWD373188/NWD373188_filtered_chr6.bim