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

# id = "NWD104907"
id.global <- tryCatch(read.table(paste0("/net/topmed2/working/khlin/global_ancestry/output/GAI/",id,"/",id,"_global.txt"), header=T), error=function(e) NULL)
HGDP_pop_id <- read.table("/net/topmed2/working/khlin/global_ancestry/HGDP_pop_id.txt", sep="\t", col.names=c("V1", "id", "V3", "pop"), stringsAsFactors = F)

#use 1% as the threshold to include the ancestral population in the reference panel
ref <- which(id.global>=0.01) 
ref_id <- with(HGDP_pop_id, id[which(pop %in% ref)])

write.table(ref, paste0("/net/topmed2/working/khlin/global_ancestry/temp/",id,"/",id,"_HGDP_ref.txt"), col.names=F, row.names=F, quote=F)
write.table(data.frame(ref_id), paste0("/net/topmed2/working/khlin/global_ancestry/temp/",id,"/",id,"_HGDP_ref_id.txt"), col.names=F, row.names=F, quote=F)
