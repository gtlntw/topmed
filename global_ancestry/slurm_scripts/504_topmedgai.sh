#!/bin/bash
set pipefail; cat /net/topmed2/working/khlin/global_ancestry/temp/NWD937291/NWD937291_filtered_chr20.bim | awk '{print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > /net/topmed2/working/khlin/global_ancestry/temp/NWD937291/NWD937291_filtered_chr20.bim.tmp && mv /net/topmed2/working/khlin/global_ancestry/temp/NWD937291/NWD937291_filtered_chr20.bim.tmp /net/topmed2/working/khlin/global_ancestry/temp/NWD937291/NWD937291_filtered_chr20.bim