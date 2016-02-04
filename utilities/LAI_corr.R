id <- c("NWD104907", "NWD536211", "NWD128368", "NWD430673", "NWD457384", "NWD114636", "NWD408313", "NWD598963", "NWD786333", "NWD919036")

# result.list <- list()

combinedData <- data.frame(x=NULL, y=NULL)
for(chr in 1:22) {
	for(i in 1:10) {
	  chr_RFMix <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id[i],"/",id[i],"_HGDP_subset_chr",chr,".0.Viterbi.txt", sep=""))
	  chr_LAMPLD <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id[i],"/",id[i],"_chr",chr,"_lampped_plot.out", sep=""))
	  
	  pos_RFMix <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id[i],"/",id[i],"_HGDP_subset_chr",chr,".pos", sep=""))
	  pos_LAMPLD <- read.table(paste("/net/topmed2/working/khlin/common_site/topmed_chr",chr,"_HGDP_common.txt", sep=""))
	  
	  pos_overlap <- which(pos_LAMPLD[,2] %in% pos_RFMix[,2])
	   
	  chr_RFMix_dosage <- apply(chr_RFMix==1, 1, sum)
	  chr_LAMPLD_dosage <- apply(chr_LAMPLD[pos_overlap,]==1, 1, sum)
	  combinedData <- rbind(combinedData, data.frame(x=chr_RFMix_dosage, y=chr_LAMPLD_dosage))
	  
	  print(dim(combinedData))
	#   cor_RFMix_LAMPLD <- cor(chr_RFMix_dosage, chr_LAMPLD_dosage)
	#   
	#   result.list[[chr]] <- c(cor_RFMix_LAMPLD)
	}
}

# result.list
cor(combinedData[,1], combinedData[,2])

