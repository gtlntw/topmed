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

##calculate Global based on result from LAI
lai.list <- list() 
for(chr in 1:22) {
lai.list[[chr]]<- read.table(paste("/net/topmed2/working/khlin/output/LAI/",id,"/",id,"_HGDP_chr",chr,".0.Viterbi.txt", sep=""))
}
lai <- factor(unlist(lai.list), levels=c("1", "2","3","4","5","6","7"))
(p <- round(prop.table(matrix(table(lai), 1), 1), 4))
#save the result
write.table(p, file = paste0("/net/topmed2/working/khlin/output/LAI/",id,"/",id,"_global_lai.txt"), quote = F, col.names = paste0("pop.", 1:7))

