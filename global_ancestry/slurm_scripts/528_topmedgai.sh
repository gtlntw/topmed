#!/bin/bash
set pipefail; cat /net/topmed2/working/khlin/global_ancestry/temp/NWD579612/NWD579612_filtered_chr22.bim | awk '{print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' > /net/topmed2/working/khlin/global_ancestry/temp/NWD579612/NWD579612_filtered_chr22.bim.tmp && mv /net/topmed2/working/khlin/global_ancestry/temp/NWD579612/NWD579612_filtered_chr22.bim.tmp /net/topmed2/working/khlin/global_ancestry/temp/NWD579612/NWD579612_filtered_chr22.bim