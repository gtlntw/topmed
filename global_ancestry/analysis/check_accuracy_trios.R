#####################################################
##full ref.panel vs adpative subset vs subject subset
#determine the switch point
check_switch <- function(lai_tmp) {
	hap1 <- diff(lai_tmp[,1])
	hap2 <- diff(lai_tmp[,2])
	switch_point <- which((hap1!=0 | hap2!=0)==T) + 1 #switching marker id
	return(switch_point)
}

#count the switch point in trios
count_switch <- function(switch_point_data_family, tolerence=0, family_member_id) {
	n_family_member <- length(switch_point_data_family)
	result <- list()

	for (member_idx in 3:n_family_member) {
		f1_n_switch <- length(unlist(switch_point_data_family[[1]])) #number of switch in founder1
		f2_n_siwtch <- length(unlist(switch_point_data_family[[2]])) #number of switch in founder2
		f1_n_switch_found <- 0
		f2_n_switch_found <- 0
		for (chr in 1:22) {
			switch_point_child <- sapply(switch_point_data_family[[member_idx]][[chr]], function(x) x + -tolerence:tolerence)
			#check the number of switch in founder1 also found in child
			f1_n_switch_found <- f1_n_switch_found + sum(switch_point_data_family[[1]][[chr]] %in% switch_point_child)
			#check the number of switch in founder2 also found in child
			f2_n_switch_found <- f2_n_switch_found + sum(switch_point_data_family[[2]][[chr]] %in% switch_point_child)
		}
		total_switch <- f1_n_switch+f2_n_siwtch
		total_switch_found <- f1_n_switch_found+f2_n_switch_found
		result[[member_idx-2]] <- data.frame(id=paste(family_member_id[1],"_",family_member_id[member_idx+1],sep=""), f1_n_switch, f1_n_switch_found, f2_n_siwtch, f2_n_switch_found, total_switch, total_switch_found, prop_found=total_switch_found/total_switch)
	}
	return(do.call(rbind, result))
}

#read in trio ID
trios <- readLines("/net/topmed2/working/khlin/infereed_trios.barnes.txt")

switch_analysis <- function(tolerence=0, n_trios=30) {
  tolerence <- tolerence #search the nearby switch point within # markers
  result_rfmix <- list()
  result_rfmix_subset <- list()
  result_rfmix_gai_subset <- list()
  for (family_id in 1:n_trios) { #
  	family_member_id <- unlist(strsplit(trios[family_id], split="\t"))
  	#read in trios file
  	n_family_member <- length(family_member_id) - 1
  	switch_point_data_family_rfmix <- list()
  	switch_point_data_family_rfmix_subset <- list()
  	switch_point_data_family_rfmix_gai_subset <- list()
  	for (family_member.idx in 1:(n_family_member)) {
  		id <- family_member_id[family_member.idx + 1] #the id starts at the second element
  		switch_point_chr_rfmix <- list()
  		switch_point_chr_rfmix_subset <- list()
  		switch_point_chr_rfmix_gai_subset <- list()
  		for (chr in 1:22) {
  			#RFMix needs to determine the switchpoint
  			lai_tmp_rfmix <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
  			switch_point_chr_rfmix[[chr]] <- check_switch(lai_tmp_rfmix)
  			#RFMix subset needs to determine the switchpoint
  			lai_tmp_rfmix_subset <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id,"/",id,"_HGDP_subset_chr",chr,".0.Viterbi.txt", sep=""))
  			switch_point_chr_rfmix_subset[[chr]] <- check_switch(lai_tmp_rfmix_subset)
  			#LAMPLD to determine the switch point
  			lai_tmp_rfmix_gai_subset <- read.table(paste("/net/topmed2/working/khlin/global_ancestry/temp/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
  			switch_point_chr_rfmix_gai_subset[[chr]] <- check_switch(lai_tmp_rfmix_gai_subset)
  		}
  		switch_point_data_family_rfmix[[family_member.idx]] <- switch_point_chr_rfmix
  		switch_point_data_family_rfmix_subset[[family_member.idx]] <- switch_point_chr_rfmix_subset
  		switch_point_data_family_rfmix_gai_subset[[family_member.idx]] <- switch_point_chr_rfmix_gai_subset
  	}
  	result_rfmix[[family_id]] <- count_switch(switch_point_data_family_rfmix, tolerence = tolerence, family_member_id= family_member_id)
  	result_rfmix_subset[[family_id]] <- count_switch(switch_point_data_family_rfmix_subset, tolerence = tolerence, family_member_id= family_member_id)
  	result_rfmix_gai_subset[[family_id]] <- count_switch(switch_point_data_family_rfmix_gai_subset, tolerence = tolerence, family_member_id= family_member_id)
  }
  (result_rfmix <- do.call(rbind, result_rfmix))
  (result_rfmix_subset <- do.call(rbind, result_rfmix_subset))
  (result_rfmix_gai_subset <- do.call(rbind, result_rfmix_gai_subset))
  (result <- rbind(data.frame(result_rfmix, method="rfmix", tolerence=tolerence),
               data.frame(result_rfmix_subset, method="rfmix_subset", tolerence=tolerence),
               data.frame(result_rfmix_gai_subset, method="rfmix_gai_subset", tolerence=tolerence)))
  return(result)
}
result_0 <- switch_analysis(tolerence = 0, n_trios=length(trios))
result_100 <- switch_analysis(tolerence = 100, n_trios=length(trios))

write.csv(rbind(result_0, result_100), paste0("result_0_100.csv"))

####################################################################
##check mendal error
##idea is to generate possible child genotype based on founder genotype
##write my own function to check

mendelian_check <- function(f1_geno, f2_geno, child_geno) {
	#there are 8 combinations considering the order matter
  geno_list <- mapply(function(x,y) paste0(x,y), f1_geno[,c(1,1,2,2)], f2_geno[,c(1,2,1,2)])
  child_geno_test <- mapply(function(x,y) paste0(x,y), child_geno[,c(1,2)], child_geno[,c(2,1)])
  result <- mapply(function(x,y) any(x %in% y), x=split(child_geno_test[,1:2],1:nrow(child_geno_test)), y=split(geno_list[,1:4],1:nrow(geno_list)))
  return(result)
}

#read in trio ID
trios <- readLines("/net/topmed2/working/khlin/infereed_trios.barnes.txt")

result_rfmix <- list()
result_rfmix_subset <- list()
result_rfmix_gai_subset <- list()
for(family_id in 1:length(trios)) {
  cat(paste0("family: ", family_id))
	family_member_id <- unlist(strsplit(trios[family_id], split="\t"))
	n_family_member <- length(family_member_id) - 1
	mendelian_error_trio_rfmix <- list()
	mendelian_error_trio_rfmix_subset <- list()
	mendelian_error_trio_rfmix_gai_subset <- list()
	for (family_member.idx in 3:(n_family_member)) {
	  cat(paste0("family_member: ", family_member.idx))
	  mendelian_error_rfmix <- list()
		mendelian_error_rfmix_subset <- list()
		mendelian_error_rfmix_gai_subset <- list()
		for(chr in 1:22) { 
		  cat(paste0("chr: ", chr))
			#founder 1 genotype
			id <- family_member_id[1+1] #since the first number is family id
			lai_f1_rfmix <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
			lai_f1_rfmix_subset <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id,"/",id,"_HGDP_subset_chr",chr,".0.Viterbi.txt", sep=""))
			lai_f1_rfmix_gai_subset <- read.table(paste("/net/topmed2/working/khlin/global_ancestry/temp/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
			#founder 2 genotype
			id <- family_member_id[2+1]
			lai_f2_rfmix <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
			lai_f2_rfmix_subset <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id,"/",id,"_HGDP_subset_chr",chr,".0.Viterbi.txt", sep=""))
			lai_f2_rfmix_gai_subset <- read.table(paste("/net/topmed2/working/khlin/global_ancestry/temp/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
			#child genotype
			id <- family_member_id[family_member.idx+1]
			lai_child_rfmix <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
			lai_child_rfmix_subset <- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id,"/",id,"_HGDP_subset_chr",chr,".0.Viterbi.txt", sep=""))
			lai_child_rfmix_gai_subset <- read.table(paste("/net/topmed2/working/khlin/global_ancestry/temp/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))

			#check mendelian error
			mendelian_error_rfmix[[chr]] <- mendelian_check(f1_geno=lai_f1_rfmix, f2_geno=lai_f2_rfmix, child_geno=lai_child_rfmix)
			mendelian_error_rfmix_subset[[chr]] <- mendelian_check(f1_geno=lai_f1_rfmix_subset, f2_geno=lai_f2_rfmix_subset, child_geno=lai_child_rfmix_subset)
			mendelian_error_rfmix_gai_subset[[chr]] <- mendelian_check(f1_geno=lai_f1_rfmix_gai_subset, f2_geno=lai_f2_rfmix_gai_subset, child_geno=lai_child_rfmix_gai_subset)
		}
		mendelian_error_trio_rfmix[[paste(family_member_id[1], family_member_id[family_member.idx+1], sep="_")]] <- mean(unlist(mendelian_error_rfmix))
		mendelian_error_trio_rfmix_subset[[paste(family_member_id[1], family_member_id[family_member.idx+1], sep="_")]] <- mean(unlist(mendelian_error_rfmix_subset))
		mendelian_error_trio_rfmix_gai_subset[[paste(family_member_id[1], family_member_id[family_member.idx+1], sep="_")]] <- mean(unlist(mendelian_error_rfmix_gai_subset))
	}
	result_rfmix[[family_id]] <- do.call(rbind, mendelian_error_trio_rfmix)
	result_rfmix_subset[[family_id]] <- do.call(rbind, mendelian_error_trio_rfmix_subset)
	result_rfmix_gai_subset[[family_id]] <- do.call(rbind, mendelian_error_trio_rfmix_gai_subset)
}
(result_rfmix <- do.call(rbind, result_rfmix))
(result_rfmix_subset <- do.call(rbind, result_rfmix_subset))
(result_rfmix_gai_subset <- do.call(rbind, result_rfmix_gai_subset))
(result <- rbind(data.frame(id=row.names(result_rfmix), mendelian_error=result_rfmix[,1], method="rfmix"),
                 data.frame(id=row.names(result_rfmix_subset), mendelian_error=result_rfmix_subset[,1], method="rfmix_subset"),
                 data.frame(id=row.names(result_rfmix_gai_subset), mendelian_error=result_rfmix_gai_subset[,1], method="rfmix_gai_subset")))

write.csv(result, paste0("result_mendelian_error.csv"), row.names = F)

#############################################################################
## global - summed local
############################################################################
trios <- readLines("/net/topmed2/working/khlin/infereed_trios.barnes.txt")
id <- unlist(sapply(strsplit(trios, split="\t"), function(x) x[-1]))

global.list <- list()
diff.list <- list()
for(i in 1:length(id)) {
  lai.list <- list() 
  for(chr in 1:22) {
    lai.list[[chr]]<- read.table(paste("/net/topmed2/working/khlin/global_ancestry/temp/",id[i],"/",id[i],"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
  }
  lai <- factor(unlist(lai.list), levels=c("1", "2","3","4","5","6","7"))
  (p <- round(prop.table(matrix(table(lai), 1), 1), 4))
  #save the result
  write.table(p, file = paste0("../output/GAI/",id[i],"/",id[i],"_global_lai_adaptive.txt"), quote = F, col.names = paste0("pop.", 1:7))
  
  
	id.global <- read.table(paste0("../output/GAI/",id[i],"/",id[i],"_global_lai.txt"), header=T)
	id.global_lai <- read.table(paste0("../output/GAI/",id[i],"/",id[i],"_global_lai_adaptive.txt"), header=T)
	diff.global <- id.global - id.global_lai
	
	diff.list[[i]] <- c(id[i], id.global, id.global_lai, diff.global)
}

diff <- do.call(rbind, diff.list)
write.table(diff, "global_diff.txt", quote = F, row.names = F)

