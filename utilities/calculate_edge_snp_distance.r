edge_snp_pos <- list()
for(chr in 1:22) {
  print(paste("processing chr:",chr))
  temp <- read.table(paste("/net/topmed2/working/khlin/common_site/topmed_freeze2_chr",chr,"_HGDP_common.txt", sep=""), header = F, stringsAsFactors = F)[,2]
  first_hgdp_snp_pos <- temp[1]
  last_hgdp_snp_pos <- temp[length(temp)]
  
  temp <- read.table(paste("/net/topmed2/working/khlin/vcf/chr",chr,".freeze2.new.pos", sep=""), header = F, stringsAsFactors = F)[,2]
  temp <- temp[which(temp < first_hgdp_snp_pos | temp > last_hgdp_snp_pos)]
  dist <- ifelse(temp < first_hgdp_snp_pos, first_hgdp_snp_pos - temp, temp - last_hgdp_snp_pos) 
  edge_snp_pos[[chr]] <- data.frame(pos=temp, dist=dist, chr=chr)  
}

cal_dist <- function(x) {
  if(is.null(x)) return()
  dist <- x[,2]
  grp <- cut(dist, c(seq(0, 1.6*10^6, 2*10^5), Inf))
  tapply(dist, grp, length)
}

lapply(edge_snp_pos, cal_dist)

write.csv(do.call(rbind, edge_snp_pos), "calculate_edge_snp_distance.csv", row.names = F)

# result <- read.csv("C:/Users/khlin/Desktop/calculate_edge_snp_distance.csv", stringsAsFactors=FALSE)

# ggplot(result, aes(group=chr, x=chr, y=dist)) +
#   geom_boxplot() +
#   labs(title="distance to the first/last HGDP SNP by chromosome", x="chromosome", y="distance")
# 
# ggsave("calculate_edge_snp_distance.png", width=5, height=5)
