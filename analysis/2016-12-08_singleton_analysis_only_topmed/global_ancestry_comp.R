# id = c("NWD104907", "NWD536211")
###########################gai for 4318
id <- read.table("/net/topmed2/working/khlin/id_freeze3.txt", stringsAsFactors = F)$V1

global.list <- list()

for(i in 1: length(id)) {
	id.global <- tryCatch(read.table(paste0("/net/topmed2/working/khlin/output/LAI/",id[i],"/",id[i],"_global_lai.txt"), header=T), error=function(e) NULL)
	if(!is.null(id.global)) global.list[[i]] <- c(id[i], id.global)
}

global <- do.call(rbind, global.list)
write.table(global, "global.txt", quote = F, row.names=F, col.names=F)

############################################### take the difference between gai and summed lai in barbidos trios
# for(i in 1: length(id)) {
# 	id.global <- read.table(paste0("output/GAI/",id[i],"/",id[i],"_global.txt"), header=T)
# 	id.global_lai <- read.table(paste0("output/GAI/",id[i],"/",id[i],"_global_lai.txt"), header=T)
# 	diff.global <- id.global - id.global_lai
# 	
# 	diff.list[[i]] <- c(id[i], id.global, id.global_lai, diff.global)
# }
# 
# diff <- do.call(rbind, diff.list)
# write.table(diff, "global_diff.txt", quote = F)
