args <- commandArgs(TRUE)
# print(args)
if(length(args) != 6){
  #print help
  cat("
Rscript utilities/eila.R i='\"1000g/NA19625/NA19625_chr20_lamp.012\"' eur='\"HGDP_938/LAMPLD/HGDP_europe_chr20.impute.hap\"' afr='\"HGDP_africa_chr20.impute.hap\"' nat='\"HGDP_native_america_chr20.impute.hap\"' out='\"test.out\"'      ")
  
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[i]))
  }
}

library(EILA, lib.loc = "/net/snowwhite/home/khlin/R/x86_64-pc-linux-gnu-library/2.13")
sampleData <- read.table(inp, header = F)
posData <- read.table(pos, header = F)$V1
eurData <- read.table(eur, header = F)
afrData <- read.table(afr, header = F)
natData <- read.table(nat, header = F)

sampleData <- sampleData[,c(1,1,1,1,1)]
res.eila <- eila(sampleData, posData, eurData, afrData, natData, rng.seed=12374607)

res.out <- res.eila$local.ancestry[, 1] 
write.table(res.out, out, quote = F, row.names = F, col.names = F)