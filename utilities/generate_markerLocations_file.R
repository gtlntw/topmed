#courtesy of Gregory Zajac
#example: vcftools --gzvcf /net/wonderland/home/gzajac/ChromoPainter_test/GFG_Interim_20150608/HGDP_938_imputed_GFG_Interim_20150608_chr22_variant_QC_phased.vcf.gz --max-missing 1 --remove-filtered-all --phased --recode -c | Rscript /net/wonderland/home/gzajac/RFmix_test/generate_markerLocations_file.R stdin /net/wonderland/home/gzajac/genetic_map_b37/genetic_map_chr22_combined_b37.txt /net/wonderland/home/gzajac/RFmix_test/GFG_Interim_20150608/HGDP_938_imputed_GFG_Interim_20150608_chr22_variant_QC_phased_markerLocations.txt
#example: Rscript /net/wonderland/home/gzajac/RFmix_test/generate_markerLocations_file.R /net/wonderland/home/gzajac/RFmix_test/GFG_Interim_20150608/HGDP_938_imputed_GFG_Interim_20150608_chr19_variant_QC_phased_positions.txt /net/wonderland/home/gzajac/genetic_map_b37/genetic_map_chr19_combined_b37.txt /net/wonderland/home/gzajac/RFmix_test/GFG_Interim_20150608/HGDP_938_imputed_GFG_Interim_20150608_chr19_variant_QC_phased_markerLocations.txt
#argv = c("/net/wonderland/home/gzajac/RFmix_test/GFG_Interim_20150608/HGDP_938_imputed_GFG_Interim_20150608_chr19_variant_QC_phased_positions.txt", "/net/wonderland/home/gzajac/genetic_map_b37/genetic_map_chr19_combined_b37.txt", "/net/wonderland/home/gzajac/RFmix_test/GFG_Interim_20150608/HGDP_938_imputed_GFG_Interim_20150608_chr19_variant_QC_phased_markerLocations.txt")

argv <- commandArgs(TRUE)

vcf_positions_file = argv[1]
genetic_map_file = argv[2]
output_file = argv[3]

vcf_marker_positions = read.table(vcf_positions_file, header=F, stringsAsFactors=F, col.names = "position")
vcf_marker_positions$vcf = T
genetic_map = read.table(genetic_map_file, header=T, stringsAsFactors=F, check.names=F)

vcf_map = merge(vcf_marker_positions, genetic_map, by = "position", all=T)
#write.table(file = paste0(output_file, ".vcf_map.txt"), vcf_map, row.names=F, col.names=F)
for(i in 1:nrow(vcf_map)){
#cat(file = paste0(output_file, ".debugging.txt"), i, unlist(vcf_map[i,]), "\n")
	if(is.na(vcf_map[i,"Genetic_Map(cM)"])){
		j = 1
		k = 1
		while(i-j>0 && is.na(vcf_map[i-j, "Genetic_Map(cM)"])){j = j + 1}
		while(i+k <= nrow(vcf_map) && is.na(vcf_map[i+k, "Genetic_Map(cM)"])){k = k + 1}
		if(i - j <= 0){
			vcf_map[i,"interpolated_map_distance"] = 0
		}
		else if (i + k > nrow(vcf_map)){
			vcf_map[i,"interpolated_map_distance"] = vcf_map[i-j,"Genetic_Map(cM)"]
		}
		else{
			#formula y = y2 + (y2 - y1) * (x - x2) / (x2 - x1)
			vcf_map[i,"interpolated_map_distance"] = vcf_map[i+k,"Genetic_Map(cM)"] + (vcf_map[i+k,"Genetic_Map(cM)"] - vcf_map[i-j,"Genetic_Map(cM)"]) * (vcf_map[i,"position"] - vcf_map[i+k,"position"]) / (vcf_map[i+k,"position"] - vcf_map[i-j,"position"])
		}
	}
	else {
		vcf_map[i,"interpolated_map_distance"] = vcf_map[i,"Genetic_Map(cM)"]
	}

}

write.table(file = output_file, vcf_map[which(vcf_map$vcf == T),"interpolated_map_distance"], row.names=F, col.names=F)

#for debugging purposes
#write.table(file = paste0(output_file, ".full_file.txt"), vcf_map, row.names=F, col.names=F)