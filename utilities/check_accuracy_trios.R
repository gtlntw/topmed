#####################################################
##LAMPLD vs RFMix subset
#determine the switch point
check_switch <- function(lai_tmp) {
	hap1 <- diff(lai_tmp[,1])
	hap2 <- diff(lai_tmp[,2])
	switch_point <- which((hap1!=0 | hap2!=0)==T) + 1 #switching marker id
	return(switch_point)
}

#count the switch point in trios
count_switch <- function(switch_point_data_family, tolerence=0, id) {
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
trios <- readLines("../infereed_trios.barnes.txt")

tolerence <- 100 #search the nearby switch point within # markers
result_rfmix <- list()
result_rfmix_subset <- list()
result_lampld <- list()
for (family_id in 1:8) {
	family_member_id <- unlist(strsplit(trios[family_id], split="\t"))
	#read in trios file
	n_family_member <- length(family_member_id) - 1
	switch_point_data_family_rfmix <- list()
	switch_point_data_family_rfmix_subset <- list()
	switch_point_data_family_lampld <- list()
	for (family_member.idx in 1:(n_family_member)) {
		id <- family_member_id[family_member.idx + 1] #the id starts at the second element
		switch_point_chr_rfmix <- list()
		switch_point_chr_rfmix_subset <- list()
		switch_point_chr_lampld <- list()
		for (chr in 1:22) {
			#RFMix needs to determine the switchpoint
			lai_tmp_rfmix <- read.table(paste("../output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
			switch_point_chr_rfmix[[chr]] <- check_switch(lai_tmp_rfmix)
			#RFMix subset needs to determine the switchpoint
			lai_tmp_rfmix_subset <- read.table(paste("../output/LAI/",id,"/",id,"_HGDP_subset_chr",chr,".0.Viterbi.txt", sep=""))
			switch_point_chr_rfmix_subset[[chr]] <- check_switch(lai_tmp_rfmix_subset)
			#LAMPLD to determine the switch point
			lai_tmp_lampld <- read.table(paste("../output/LAI/",id,"/",id,"_chr",chr,"_lampped_plot.out", sep=""))
			switch_point_chr_lampld[[chr]] <- check_switch(lai_tmp_lampld)
		}
		switch_point_data_family_rfmix[[family_member.idx]] <- switch_point_chr_rfmix
		switch_point_data_family_rfmix_subset[[family_member.idx]] <- switch_point_chr_rfmix_subset
		switch_point_data_family_lampld[[family_member.idx]] <- switch_point_chr_lampld
	}
	result_rfmix[[family_id]] <- count_switch(switch_point_data_family_rfmix, tolerence = tolerence, id= family_member_id)
	result_rfmix_subset[[family_id]] <- count_switch(switch_point_data_family_rfmix_subset, tolerence = tolerence, id= family_member_id)
	result_lampld[[family_id]] <- count_switch(switch_point_data_family_lampld, tolerence = tolerence, id= family_member_id)
}
(result_rfmix <- do.call(rbind, result_rfmix))
(result_rfmix_subset <- do.call(rbind, result_rfmix_subset))
(result_lampld <- do.call(rbind, result_lampld))
write.csv(result_rfmix, paste0("result_rfmix_",tolerence,".csv"))
write.csv(result_rfmix_subset, paste0("result_rfmix_subset_",tolerence,".csv"))
write.csv(result_lampld, paste0("result_lampld_",tolerence,".csv"))


####################################################################
##check mendal error
##idea is to generate possible child genotype based on founder genotype
##write my own function to check

mendelian_check <- function(f1_geno, f2_geno, child_geno) {
	#there are 8 combinations considering the order matter
  geno_list <- paste(f1_geno[c(1,1,2,2)], f2_geno[c(1,2,1,2)], sep="")
  child_geno_test <- paste(child_geno[c(1,2)], child_geno[c(2,1)], sep="")
  return(any(child_geno_test %in% geno_list))
}

#read in trio ID
trios <- readLines("../infereed_trios.barnes.txt")

result_rfmix <- list()
result_rfmix_subset <- list()
result_lampld <- list()
for(family_id in 1:8) {
	family_member_id <- unlist(strsplit(trios[family_id], split="\t"))
	n_family_member <- length(family_member_id) - 1
	mendelian_error_trio_rfmix <- list()
	mendelian_error_trio_rfmix_subset <- list()
	mendelian_error_trio_lampld <- list()
	for (family_member.idx in 3:(n_family_member)) {
		mendelian_error_rfmix <- list()
		mendelian_error_rfmix_subset <- list()
		mendelian_error_lampld <- list()
		for(chr in 1:22) { 
			#founder 1 genotype
			id <- family_member_id[1+1] #since the first number is family id
			lai_f1_rfmix <- read.table(paste("../output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
			lai_f1_rfmix_subset <- read.table(paste("../output/LAI/",id,"/",id,"_HGDP_subset_chr",chr,".0.Viterbi.txt", sep=""))
			lai_f1_lampld <- read.table(paste("../output/LAI/",id,"/",id,"_chr",chr,"_lampped_plot.out", sep=""))
			#founder 2 genotype
			id <- family_member_id[2+1]
			lai_f2_rfmix <- read.table(paste("../output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
			lai_f2_rfmix_subset <- read.table(paste("../output/LAI/",id,"/",id,"_HGDP_subset_chr",chr,".0.Viterbi.txt", sep=""))
			lai_f2_lampld <- read.table(paste("../output/LAI/",id,"/",id,"_chr",chr,"_lampped_plot.out", sep=""))
			#child genotype
			id <- family_member_id[family_member.idx+1]
			lai_child_rfmix <- read.table(paste("../output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
			lai_child_rfmix_subset <- read.table(paste("../output/LAI/",id,"/",id,"_HGDP_subset_chr",chr,".0.Viterbi.txt", sep=""))
			lai_child_lampld <- read.table(paste("../output/LAI/",id,"/",id,"_chr",chr,"_lampped_plot.out", sep=""))

			#check mendelian error
			mendelian_error_rfmix[[chr]] <- sapply(1:nrow(lai_f1_rfmix), function(i) 
				mendelian_check(f1_geno=lai_f1_rfmix[i, ], f2_geno=lai_f2_rfmix[i, ], child_geno=lai_child_rfmix[i, ]))
			mendelian_error_rfmix_subset[[chr]] <- sapply(1:nrow(lai_f1_rfmix_subset), function(i) 
				mendelian_check(f1_geno=lai_f1_rfmix_subset[i, ], f2_geno=lai_f2_rfmix_subset[i, ], child_geno=lai_child_rfmix_subset[i, ]))
			mendelian_error_lampld[[chr]] <- sapply(1:nrow(lai_f1_lampld), function(i) 
				mendelian_check(f1_geno=lai_f1_lampld[i, ], f2_geno=lai_f2_lampld[i, ], child_geno=lai_child_lampld[i, ]))
		}
		mendelian_error_trio_rfmix[[paste(family_member_id[1], family_member_id[family_member.idx+1], sep="_")]] <- mean(unlist(mendelian_error_rfmix))
		mendelian_error_trio_rfmix_subset[[paste(family_member_id[1], family_member_id[family_member.idx+1], sep="_")]] <- mean(unlist(mendelian_error_rfmix_subset))
		mendelian_error_trio_lampld[[paste(family_member_id[1], family_member_id[family_member.idx+1], sep="_")]] <- mean(unlist(mendelian_error_lampld))
	}
	result_rfmix[[family_id]] <- do.call(rbind, mendelian_error_trio_rfmix)
	result_rfmix_subset[[family_id]] <- do.call(rbind, mendelian_error_trio_rfmix_subset)
	result_lampld[[family_id]] <- do.call(rbind, mendelian_error_trio_lampld)
}
(result_rfmix <- do.call(rbind, result_rfmix))
(result_rfmix_subset <- do.call(rbind, result_rfmix_subset))
(result_lampld <- do.call(rbind, result_lampld))
write.csv(result_rfmix, "result_rfmix_mendelian.csv")
write.csv(result_rfmix_subset, "result_rfmix_subset_mendelian.csv")
write.csv(result_lampld, "result_lampld_mendelian.csv")
