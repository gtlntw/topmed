#example: Rscript /net/wonderland/home/gzajac/ChromoPainter_test/make_HGDP_idfile_poplist.R /net/wonderland/home/gzajac/ChromoPainter_test/MGI_phase_1/HGDP_938_Illumina_20131217_chr1_sampleQC_genderQC_variantQC_imputed.impute.hap.indv MGI

argv <- commandArgs(TRUE)

impute.hap.indv.file = argv[1]
output = argv[2]
subset = argv[3]

impute.hap.indv.table = read.table(impute.hap.indv.file, sep="\t", header=F, stringsAsFactors=F, col.names = "IID")

iid_subpop = read.table(file = "/net/wonderland/home/gzajac/QC_genotyping/HGDP_938/HGDP_populations1.txt", sep="\t", header=T, stringsAsFactors=F)
populations_table = read.table(file = "/net/wonderland/home/gzajac/QC_genotyping/HGDP_938/HGDP_populations.txt", sep="\t", header=T, stringsAsFactors=F)

HGDP_populations = merge(iid_subpop, populations_table, by.x="Population", by.y="PopName1", all.x=T, all.y=F)
impute.hap.indv.HGDP_populations = merge(impute.hap.indv.table, merge(impute.hap.indv.table, HGDP_populations, by="IID", all.x=T, all.y=F, sort=F), by="IID", all=T, sort=F)

#impute.hap.indv.HGDP_populations$one = 1
impute.hap.indv.HGDP_populations$Ancestry = ifelse(is.na(impute.hap.indv.HGDP_populations$Ancestry), "0 0", 
												ifelse(impute.hap.indv.HGDP_populations$Ancestry == "Africa", "1 1",
												ifelse(impute.hap.indv.HGDP_populations$Ancestry == "Central/South Asia", "2 2",
												ifelse(impute.hap.indv.HGDP_populations$Ancestry == "Eastern Asia", "3 3",
												ifelse(impute.hap.indv.HGDP_populations$Ancestry == "Europe", "4 4",
												ifelse(impute.hap.indv.HGDP_populations$Ancestry == "Native America", "5 5",
												ifelse(impute.hap.indv.HGDP_populations$Ancestry == "Oceania", "6 6",
												ifelse(impute.hap.indv.HGDP_populations$Ancestry == "Western Asia" | impute.hap.indv.HGDP_populations$Ancestry == "Africa2", "7 7", "0 0"))))))))

#convert class from ex: 1,4,7 to 1,2,3
if(!is.na(subset)) {
  ref <- read.table(subset, header = F)$V1
  ref_n_class <- length(ref)
  ref_class <- c("0 0", sapply(1:ref_n_class, function(x) paste(x,x)))
  names(ref_class) <- c("0 0", sapply(ref, function(x) paste(x,x)))
  impute.hap.indv.HGDP_populations$Ancestry <- ref_class[impute.hap.indv.HGDP_populations$Ancestry]
}

cat(impute.hap.indv.HGDP_populations$Ancestry, file = output, sep = " ")

#impute.hap.indv.HGDP_populations$Ancestry = gsub(" ", "_", impute.hap.indv.HGDP_populations$Ancestry)



#idfile_table = (cbind(impute.hap.indv.HGDP_populations$IID, impute.hap.indv.HGDP_populations$Ancestry, impute.hap.indv.HGDP_populations$one))
#write.table(idfile_table, file = paste0(impute.hap.indv.file, ".idfile.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

#poplist = as.data.frame(table(HGDP_populations$Ancestry), stringsAsFactors=F)
#poplist$D = "D"
#poplist$neg_nine = -9
#write.table(cbind(poplist$Var1, poplist$D, 1/poplist$Freq, poplist$neg_nine), file = paste0(impute.hap.indv.file, ".poplist.txt"),quote = F, sep = "\t", row.names = F, col.names = F)
#cat(paste("Africa", "D", "1", "-9", sep="\t"), 1
#paste("Central/South_Asia", "D", "1", "-9", sep="\t"), 2
#paste("Eastern_Asia", "D", "1", "-9", sep="\t"), 3
#paste("Europe", "D", "1", "-9", sep="\t"), 4
#paste("Native_America", "D", "1", "-9", sep="\t"), 5
#paste("Oceania", "D", "1", "-9", sep="\t"), 6
#paste("Western_Asia/Northern_Africa", "D", "1", "-9", sep="\t"), 7
#paste(target_population_name,"R", sep="\t"),
#file = paste0(impute.hap.indv.file, ".poplist.txt"), sep="\n")

#HGDP_populations$one = 1
#HGDP_populations$Ancestry = ifelse(HGDP_populations$Ancestry == "Western Asia" | HGDP_populations$Ancestry == "Africa2", "Western Asia/Northern Africa", HGDP_populations$Ancestry)
#HGDP_populations$Ancestry = gsub(" ", "_", HGDP_populations$Ancestry)
#idfile_table = (cbind(HGDP_populations$IID, HGDP_populations$Ancestry, HGDP_populations$one))[order(HGDP_populations$IID),]
#write.table(idfile_table, file = "/net/wonderland/home/gzajac/ChromoPainter_test/HGDP_938.idfile.txt", quote = F, sep = "\t", row.names = F, col.names = F)

#poplist = as.data.frame(table(HGDP_populations$Ancestry), stringsAsFactors=F)
#poplist$D = "D"
#poplist$neg_nine = -9
#write.table(cbind(poplist$Var1, poplist$D, 1/poplist$Freq, poplist$neg_nine), file = "/net/wonderland/home/gzajac/ChromoPainter_test/HGDP_938.poplist.txt",quote = F, sep = "\t", row.names = F, col.names = F)
