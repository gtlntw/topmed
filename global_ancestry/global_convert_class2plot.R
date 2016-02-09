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
# print(id)

viterbi <- read.table(paste0("/net/topmed2/working/khlin/global_ancestry/temp/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt"), header=F)

#convert ex: 1,2,3 to 1,4,7
ref <- read.table(paste0("/net/topmed2/working/khlin/global_ancestry/temp/",id,"/",id,"_HGDP_ref.txt"), header = F)$V1
ref_n_class <- length(ref)
ref_class <- ref
names(ref_class) <- 1:ref_n_class
result <- cbind(ref_class[viterbi$V1], ref_class[viterbi$V2]) 

write.table(result, paste0("/net/topmed2/working/khlin/global_ancestry/temp/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt"), col.names=F, row.names=F, quote=F)
            