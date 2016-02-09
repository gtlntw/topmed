##parse parameters from command
parseCommandArgs <- function (evaluate = TRUE) 
{
  arglist <- list()
  args <- commandArgs()
  i <- which(args == "--args")
  if (length(i) == 0 || length(args) < 1) 
    return(invisible())
  args <- args[(i + 1):length(args)]
  for (i in seq(1, length(args), by = 2)) {
    value <- NA
    tryCatch(value <- as.double(args[i + 1]), warning = function(e) {
    })
    if (is.na(value)) {
      value <- args[i + 1]
      if (substr(value, 1, 2) == "c(") 
        value <- eval(parse(text = args[i + 1]))
    }
    if (evaluate) 
      assign(args[i], value, inherits = TRUE)
    arglist[[length(arglist) + 1]] <- value
    names(arglist)[length(arglist)] <- args[i]
  }
  return(arglist)
}
## passed in from the command prompt.
parseCommandArgs()
##print parameter
print(id)

global.data.list <- list()

for(chr in 1:22) {
	print(paste("Processing chr:", chr))
	id.data <- read.table(paste0("temp/",id,"/",id,"_filtered_setref_chr",chr,".frq.counts"), header=T)
	pop1.data <- read.table(paste0("HGDP/HGDP_chr",chr,"_pop_1.frq"), header=T, col.names=c("CHR","SNP","A1","A2","MAF.1","NCHROBS.1"))
	pop2.data <- read.table(paste0("HGDP/HGDP_chr",chr,"_pop_2.frq"), header=T, col.names=c("CHR","SNP","A1","A2","MAF.2","NCHROBS.2"))
	pop3.data <- read.table(paste0("HGDP/HGDP_chr",chr,"_pop_3.frq"), header=T, col.names=c("CHR","SNP","A1","A2","MAF.3","NCHROBS.3"))
	pop4.data <- read.table(paste0("HGDP/HGDP_chr",chr,"_pop_4.frq"), header=T, col.names=c("CHR","SNP","A1","A2","MAF.4","NCHROBS.4"))
	pop5.data <- read.table(paste0("HGDP/HGDP_chr",chr,"_pop_5.frq"), header=T, col.names=c("CHR","SNP","A1","A2","MAF.5","NCHROBS.5"))
	pop6.data <- read.table(paste0("HGDP/HGDP_chr",chr,"_pop_6.frq"), header=T, col.names=c("CHR","SNP","A1","A2","MAF.6","NCHROBS.6"))
	pop7.data <- read.table(paste0("HGDP/HGDP_chr",chr,"_pop_7.frq"), header=T, col.names=c("CHR","SNP","A1","A2","MAF.7","NCHROBS.7"))

	pop.data <- Reduce(function(...) merge(..., by=c("SNP","CHR","A1","A2")), list(pop1.data,pop2.data,pop3.data,pop4.data,pop5.data,pop6.data,pop7.data))
	global.data <- merge(id.data, pop.data, by=c("SNP","CHR","A1","A2"))
	global.data.list[[chr]] <- global.data
}
global.data <- do.call(rbind, global.data.list)

#when heterogeneous ancestral population
mdbinom <- function(x, f1, f2) {
	(x==0)*(1-f1)*(1-f2) + (x==1)*(f1*(1-f2)+(1-f1)*f2) + (x==2)*f1*f2
}

#precalculate likelihood at each marker given ancestral MAFs
for(i in 1:7) {
	for(j in i:7) {
		lk_exp <- paste0("lk",i,j," <- mdbinom(global.data$C1, global.data$MAF.",i,", global.data$MAF.",j,")")
		eval(parse(text=lk_exp))
	}
}

lkh <- function(c) { #likelihood function using precalculated data
	#if(!all(c>=0)) return(NA) #check if all positive probability
	sum(log(c[1]^2*lk11 + c[2]^2*lk22 +	c[3]^2*lk33 +	c[4]^2*lk44 +	c[5]^2*lk55 +	c[6]^2*lk66 +	c[7]^2*lk77 +
						2*c[1]*c[2]*lk12 + 2*c[1]*c[3]*lk13 + 2*c[1]*c[4]*lk14 + 2*c[1]*c[5]*lk15 + 2*c[1]*c[6]*lk16 + 2*c[1]*c[7]*lk17 +
						2*c[2]*c[3]*lk23 + 2*c[2]*c[4]*lk24 + 2*c[2]*c[5]*lk25 + 2*c[2]*c[6]*lk26 + 2*c[2]*c[7]*lk27 +
						2*c[3]*c[4]*lk34 + 2*c[3]*c[5]*lk35 + 2*c[3]*c[6]*lk36 + 2*c[3]*c[7]*lk37 +
						2*c[4]*c[5]*lk45 + 2*c[4]*c[6]*lk46 + 2*c[4]*c[7]*lk47 +
						2*c[5]*c[6]*lk56 + 2*c[5]*c[7]*lk57 + 
						2*c[6]*c[7]*lk67))
}

lkh_link <- function(x) { #use the concept from multinomial glm to limit the sum of probability equals to one
	sum_exp <- 1+ sum(exp(x))
	p <- c(exp(x)/sum_exp, 1/sum_exp)
	lkh(p)
}
lkh_link(rnorm(6))

#use 10 multiple random start and keep the best result
result.best.value = -Inf
for(i in 1:10) {
	print(paste("running iteration:", i))
	result <- optim(rnorm(6), lkh_link, method="BFGS", control = list(fnscale=-1, maxit=5000))
	if(result.best.value < result$value) {
		result.best.value <- result$value
		result.best <- result
	}
	#print(result)
	print(result$value)
	sum_exp <- 1+ sum(exp(result$par))
	(p <- round(c(exp(result$par)/sum_exp, 1/sum_exp), 4))
	print(p)
}
result.best$value
sum_exp <- 1+ sum(exp(result.best$par))
(p <- round(c(exp(result.best$par)/sum_exp, 1/sum_exp), 4))

#save the result
write.table(t(p), file = paste0("output/GAI/",id,"/",id,"_global.txt"), quote = F, col.names = paste0("pop.", 1:7))

##calculate Global based on result from LAI
# lai.list <- list() 
# for(chr in 1:22) {
# lai.list[[chr]]<- read.table(paste("../output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
# }
# lai <- factor(unlist(lai.list), levels=c("1", "2","3","4","5","6","7"))
# (p <- round(prop.table(matrix(table(lai), 1), 1), 4))
# #save the result
# write.table(p, file = paste0("output/GAI/",id,"/",id,"_global_lai.txt"), quote = F, col.names = paste0("pop.", 1:7))
# 

########################################################################
######11/12 calculation assuming only homogeneous ancestral population combination
########################################################################

# #precalculate likelihood at each marker given ancestral MAFs
# lk <- cbind(dbinom(global.data$C1,2, global.data$MAF.1),
# 						dbinom(global.data$C1,2, global.data$MAF.2),
# 						dbinom(global.data$C1,2, global.data$MAF.3),
# 						dbinom(global.data$C1,2, global.data$MAF.4),
# 						dbinom(global.data$C1,2, global.data$MAF.5),
# 						dbinom(global.data$C1,2, global.data$MAF.6),
# 						dbinom(global.data$C1,2, global.data$MAF.7))
# 
# lkh <- function(c) { #likelihood function using precalculated data
# 	#if(!all(c>=0)) return(NA) #check if all positive probability
# 	sum(log(c[1]*lk[,1] +
# 						c[2]*lk[,2] +
# 						c[3]*lk[,3] +
# 						c[4]*lk[,4] +
# 						c[5]*lk[,5] +
# 						c[6]*lk[,6] +
# 						c[7]*lk[,7]))
# }
# 
# lkh_link <- function(x) { #use the concept from multinomial glm to limit the sum of probability equals to one
# 	sum_exp <- 1+ sum(exp(x))
# 	p <- c(exp(x)/sum_exp, 1/sum_exp)
# 	lkh(p)
# }
# lkh_link(rep(0,6))
# 
# #use 10 multiple random start and keep the best result
# result.best.value = -Inf
# for(i in 1:10) {
# 	print(paste("running iteration:", i))
# 	result <- optim(rnorm(6), lkh_link, method="BFGS", control = list(fnscale=-1, maxit=5000))
# 	if(result.best.value < result$value) result.best <- result
# 	#print(result)
# }
# result.best$value
# sum_exp <- 1+ sum(exp(result$par))
# (p <- round(c(exp(result.best$par)/sum_exp, 1/sum_exp), 4))
# 
# #save the result
# write.table(t(p), file = paste0("output/GAI/",id,"/",id,"_global.txt"), quote = F, col.names = paste0("pop.", 1:7))
# 
# ##calculate Global based on result from LAI
# lai.list <- list() 
# for(chr in 1:22) {
# lai.list[[chr]]<- read.table(paste("../output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
# }
# lai <- factor(unlist(lai.list), levels=c("1", "2","3","4","5","6","7"))
# (p <- round(prop.table(matrix(table(lai), 1), 1), 4))
# #save the result
# write.table(p, file = paste0("output/GAI/",id,"/",id,"_global_lai.txt"), quote = F, col.names = paste0("pop.", 1:7))


######################################################
#####make HGDP_pop_id.txt
# iid_subpop = read.table(file = "/net/wonderland/home/gzajac/QC_genotyping/HGDP_938/HGDP_populations1.txt", sep="\t", header=T, stringsAsFactors=F)
# populations_table = read.table(file = "/net/wonderland/home/gzajac/QC_genotyping/HGDP_938/HGDP_populations.txt", sep="\t", header=T, stringsAsFactors=F)
# HGDP_populations = merge(iid_subpop, populations_table, by.x="Population", by.y="PopName1", all.x=T, all.y=F)
# 
# pop <- c("Africa","Central/South Asia","Eastern Asia","Europe","Native America","Oceania","Western Asia")
# Ancestry <- HGDP_populations$Ancestry
# Ancestry[which(Ancestry=="Africa2")] <- "Western Asia"
# pop_ancestry <- unlist(sapply(Ancestry, function(x) grep(x, pop)))
# 
# write.table(cbind(HGDP_populations[, 1:3],pop_ancestry), "HGDP_pop_id.txt",sep="\t", quote=F, row.names=F, col.names=F)
