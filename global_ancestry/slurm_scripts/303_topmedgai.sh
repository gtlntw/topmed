#!/bin/bash
set pipefail; cat /net/topmed2/working/khlin/global_ancestry/temp/NWD241576/NWD241576_filtered_chr17.bim | awk '{print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > /net/topmed2/working/khlin/global_ancestry/temp/NWD241576/NWD241576_filtered_chr17.bim.tmp && mv /net/topmed2/working/khlin/global_ancestry/temp/NWD241576/NWD241576_filtered_chr17.bim.tmp /net/topmed2/working/khlin/global_ancestry/temp/NWD241576/NWD241576_filtered_chr17.bim