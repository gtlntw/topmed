#!/usr/bin/env python
#!/usr/bin/python
#!python

import sys
import os
import re
import getopt
import numpy
from string import Template

def usage():
	print """
	create makefile for the pipeline of simulations

	 makefile [-h] [-i <input filename>] [-o <output filename>]

	 -h	  ask help

	 -l	  launch by local or slurm

	 -j	  job name

	 """

#option variables
opts = {}
opts["help"] = ""
opts["verbose"] = ""
opts["debug"] = ""
opts["outputDir"] = os.getcwd()
opts["makeFile"] = "makefile"
launchMethod = "local"
opts["id"] = ""
opts["jobName"] = "test"
idList = []

o, a = getopt.getopt(sys.argv[1:], "i:l:j:h")
for k,v in o:
	opts[k] = v
	if opts.has_key("-h"):
		usage(); sys.exit(0)
	elif k == "-l":
		launchMethod = v
		if launchMethod != "local" and launchMethod != "slurm":
			print "Launch method has to be local or slurm"; exit(1)
	elif k == "-j":
		opts["jobName"] = v
		opts["makeFile"] = "makefile_{jobName}".format(**opts)
	elif k == "-i":
		opts["id"] = v


##############
#print options
##############
print "Options"
print "sample ID	: {id}".format(**opts)
print "output directory : {outputDir}".format(**opts)
print "launch method	: {launchMethod}".format(launchMethod = launchMethod)
print "output jobName   : {jobName}".format(**opts)
print "makefile	 : {makeFile}".format(**opts)
print ""

#arrays for storing targets, dependencies and commands
tgts = []
deps = []
cmds = []
 
#temporary variables
tgt = ""
dep = ""
cmd = []
inputFiles = []
inputFilesOK = []
inputFile = ""
outputFile = ""
opts["param"] = "" #parameters for slurm when using shell script
opts["toolsDir"] = "/net/snowwhite/home/khlin/tools"

#create directory needed for slurm script 
if not os.path.exists(opts["outputDir"]): os.makedirs(opts["outputDir"])
opts["slurmScriptsDir"] = "{outputDir}/slurm_scripts".format(**opts)
if not os.path.exists(opts["slurmScriptsDir"]): os.makedirs(opts["slurmScriptsDir"])
opts["slurmScriptNo"] = 0
if not os.path.exists("output"): os.makedirs("output")
if not os.path.exists("log"): os.makedirs("log")

##########
#functions
##########
 
#run a job either locally or by slurm
def makeJob(method, tgt, dep, cmd):
	if method == "local":
		makeLocalStep(tgt, dep, cmd)
	elif method == "slurm":
		makeSlurm(tgt, dep, cmd)
 
#run slurm jobs
def makeSlurm(tgt, dep, cmd):
	tgts.append(tgt)
	deps.append(dep)
	cmd_tmp = []
	for c in cmd:
		opts.update({"command": c})
		#contains pipe or cd
		if re.search("\||^cd", c):
			opts["slurmScriptNo"] += 1
			opts["slurmScriptFile"] = "{slurmScriptsDir}/{slurmScriptNo}_{jobName}.sh".format(**opts)
			IN = open("{slurmScriptFile}".format(**opts), "w")
			IN.write("#!/bin/bash\n")
			IN.write("set pipefail; {command}".format(**opts))
			IN.close()
			os.chmod("{slurmScriptFile}".format(**opts), 0755)

			cmd_tmp.append("\tsrun -p nomosix,main -J {jobName} -D {outputDir} {param} {slurmScriptFile} \n".format(**opts))
		else:
			cmd_tmp.append("\tsrun -p nomosix,main -J {jobName} -D {outputDir} {param} {command} \n".format(**opts))
	cmd_tmp.append("\ttouch {tgt}\n".format(tgt=tgt))
	cmds.append(cmd_tmp)
 
#run a local job
def makeLocalStep(tgt, dep, cmd):
	tgts.append(tgt)
	deps.append(dep)
	cmd_tmp = []
	for c in cmd:
		cmd_tmp.append("\t{command}\n".format(command=c))
	cmd_tmp.append("\ttouch {tgt}\n".format(tgt=tgt))
	cmds.append(cmd_tmp)

######################
#read in sample id
######################
if opts["id"] == ".":
	infh = open("id.txt")
	for line in infh.xreadlines():
		item = line.strip()
		if not re.search("^#+", item):
			idList.append(item)
	infh.close()
else:
	idList = [opts["id"]] #use specified id
print "sample id:"
print idList


#############################################
# Start global ancestry inference
#############################################

# ######################
# #0.1. ld-pruning
# ######################
# for chr in xrange(1,23):
# 	opts["chr"] = chr
# 	tgt = "{outputDir}/temp/{id}/{id}_filtered_chr{chr}.plink.OK".format(**opts)
# 	dep = ""
# 	cmd = ["/net/snowwhite/home/khlin/tools/plink_1.9/plink --vcf /net/topmed2/working/khlin/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz \
# --indep-pairwise 50 5 0.1 -out {outputDir}/HGDP/HGDP_chr${chr}".format(**opts)]
# 	makeJob("local", tgt, dep, cmd)

# ######################
# #0.2. generate bim file for allele order and position
# ######################
# 	opts["chr"] = chr
# 	tgt = "{outputDir}/temp/{id}/{id}_filtered_chr{chr}.bim.OK".format(**opts)
# 	dep = "{outputDir}/temp/{id}/{id}_filtered_chr{chr}.plink.OK".format(**opts)
# 	cmd = ["/net/snowwhite/home/khlin/tools/plink_1.9/plink --vcf /net/topmed2/working/khlin/HGDP_938/HGDP_938_chr${chr}_filtered_phased.vcf.gz \
# --extract {outputDir}/HGDP/HGDP_chr${chr}.prune.in -make-bed -out {outputDir}/HGDP/HGDP_chr${chr}".format(**opts)]
# 	makeJob("local", tgt, dep, cmd)

# ######################
# #0.3. change rsid to chr:position
# ######################
# 	opts["chr"] = chr
# 	tgt = "{outputDir}/temp/{id}/{id}_filtered_setref_chr{chr}.vcf.gz.OK".format(**opts)
# 	dep = "{outputDir}/temp/{id}/{id}_filtered_chr{chr}.bim.OK".format(**opts)
# 	cmd = ["cat {outputDir}/HGDP/HGDP_chr${chr}.bim | awk '{{print $1\"\\t\"$1\":\"$4\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6}}' > {outputDir}/HGDP/HGDP_chr${chr}.bim.tmp && mv {outputDir}/HGDP/HGDP_chr${chr}.bim.tmp {outputDir}/HGDP/HGDP_chr${chr}.bim".format(**opts)]
# 	makeJob("local", tgt, dep, cmd)

# ######################
# #0.4 generate list of snps to keep for bcftools
# ######################
# 	opts["chr"] = chr
# 	tgt = "{outputDir}/temp/{id}/{id}_filtered_setref_chr{chr}.freq.counts.OK".format(**opts)
# 	inputFilesOK.append(tgt)
# 	dep = "{outputDir}/temp/{id}/{id}_filtered_setref_chr{chr}.vcf.gz.OK".format(**opts)
# 	cmd = ["cat {outputDir}/HGDP/HGDP_chr${chr}.bim | awk '{{print $1\"\\t\"$4}}' > {outputDir}/HGDP/HGDP_chr${chr}.pos".format(**opts)]
# 	makeJob("local", tgt, dep, cmd)

# ######################
# #0.5 global ancestry inference
# ######################
# 	for pop in xrange(1,8):
# 		tgt = "{outputDir}/temp/{id}/{id}_global.OK".format(**opts)
# 		dep = " ".join(inputFilesOK)
# 		cmd = ["cat HGDP_pop_id.txt | awk -F'\\t' '{{if($4=='${pop}') print $2\" \"$2}}' > HGDP_pop_${pop}_id.txt".format(**opts),
# 		"{toolsDir}/plink_1.9/plink --bfile {outputDir}/HGDP/HGDP_chr${chr} --keep {outputDir}/HGDP/HGDP_pop_${pop}_id.txt --freq --out {outputDir}/HGDP/HGDP_chr${chr}_pop_${pop} --keep-allele-order".format(**opts)]
# 		makeJob("local", tgt, dep, cmd)

######################
#1.0. log the start time
######################
tgt = "{outputDir}/log/start.runmake.gai.OK".format(**opts)
dep = ""
cmd = ["date | awk '{{print \"Global ancestry pipeline\\n\\nstart: \"$$0}}' > {outputDir}/log/runmake_gai_time.log".format(**opts)]
makeJob("local", tgt, dep, cmd)

opts["param"] = "--time=0-1" #indicate this is a quick job
for id in idList:
	opts["id"] = id
	if not os.path.exists("output/GAI/{id}".format(**opts)): os.makedirs("output/GAI/{id}".format(**opts)) #create GAI result folder
	if not os.path.exists("temp/{id}".format(**opts)): os.makedirs("temp/{id}".format(**opts)) #create temp folder
	######################
	#1.1. extract sample vcf and convert to bfile
	#freeze.1a /net/topmed2/working/gt-release/sftp-barnes/freeze.1a/topmed.freeze1.nhlbi.791.sftp-barnes.keep.chr{chr}.gtonly.vcf.gz
	#freeze.1b /net/topmed2/working/hmkang/snps/4318/filtered/chr{chr}.filtered.rehdr.gt2.vcf.gz
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_filtered_chr{chr}.plink.OK".format(**opts)
		dep = ""
		cmd = ["bcftools view /net/topmed2/working/hmkang/snps/4318/filtered/chr{chr}.filtered.rehdr.gt2.vcf.gz \
--types snps -M2 --exclude-uncalled -f PASS -s {id} -R {outputDir}/HGDP/HGDP_chr{chr}.pos --output-type z --output-file {outputDir}/temp/{id}/{id}_filtered_chr{chr}.vcf.gz && \
{toolsDir}/plink_1.9/plink --vcf {outputDir}/temp/{id}/{id}_filtered_chr{chr}.vcf.gz -make-bed --out {outputDir}/temp/{id}/{id}_filtered_chr{chr}".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#1.2. create variant id
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_filtered_chr{chr}.bim.OK".format(**opts)
		dep = "{outputDir}/temp/{id}/{id}_filtered_chr{chr}.plink.OK".format(**opts)
		cmd = ["cat {outputDir}/temp/{id}/{id}_filtered_chr{chr}.bim | awk '{{print $1\"\\t\"$1\":\"$4\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6}}' > {outputDir}/temp/{id}/{id}_filtered_chr{chr}.bim.tmp \
&& mv {outputDir}/temp/{id}/{id}_filtered_chr{chr}.bim.tmp {outputDir}/temp/{id}/{id}_filtered_chr{chr}.bim".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)
	
	######################
	#1.3. set ref allele
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_filtered_setref_chr{chr}.vcf.gz.OK".format(**opts)
		dep = "{outputDir}/temp/{id}/{id}_filtered_chr{chr}.bim.OK".format(**opts)
		cmd = ["{toolsDir}/plink_1.9/plink --bfile {outputDir}/temp/{id}/{id}_filtered_chr{chr} --a2-allele {outputDir}/HGDP/HGDP_chr{chr}.bim 6 2 --make-bed \
--out {outputDir}/temp/{id}/{id}_filtered_setref_chr{chr}".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#1.4 calculate the allele count for sample
	######################
	inputFilesOK = [] #clean up
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_filtered_setref_chr{chr}.freq.counts.OK".format(**opts)
		inputFilesOK.append(tgt)
		dep = "{outputDir}/temp/{id}/{id}_filtered_setref_chr{chr}.vcf.gz.OK".format(**opts)
		cmd = ["{toolsDir}/plink_1.9/plink --bfile {outputDir}/temp/{id}/{id}_filtered_setref_chr{chr} --freq counts --out {outputDir}/temp/{id}/{id}_filtered_setref_chr{chr} --keep-allele-order".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#1.5 global ancestry inference
	######################
	tgt = "{outputDir}/temp/{id}/{id}_global.OK".format(**opts)
	dep = " ".join(inputFilesOK)
	cmd = ["Rscript global_ancestry.R id {id}".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)


opts["outputDir"] = '/net/topmed2/working/khlin'
opts["outputDirGlobal"] = '/net/topmed2/working/khlin/global_ancestry'
opts["param"] = "--time=0-2  --exclude=sun08,sun05,sun02" #indicate this is a quick job
######################
#2.0. #RFMix local ancestry inference using all 7 ancestral populations filtered HGDP and keep only bi-allelic
######################
for chr in xrange(1,23):
	opts["chr"] = chr
	tgt = "{outputDir}/HGDP_938/HGDP_938_chr{chr}_filtered_phased.vcf.gz.OK".format(**opts)
	dep = ""
	cmd = ["bcftools view {outputDir}/HGDP_938/HGDP_938_chr{chr}_phased.vcf.gz \
--types snps -M2 --exclude-uncalled -f PASS \
--output-type z --output-file {outputDir}/HGDP_938/HGDP_938_chr{chr}_filtered_phased.vcf.gz \
bcftools index -t -f {outputDir}/HGDP_938/HGDP_938_chr{chr}_filtered_phased.vcf.gz".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)

######################
#2.1. ind. common snp between sample and HGDP
######################
for chr in xrange(1,23):
	opts["chr"] = chr
	tgt = "{outputDir}/common_site/topmed_chr{chr}_HGDP_common_global.txt.OK".format(**opts)
	dep = "{outputDir}/HGDP_938/HGDP_938_chr{chr}_filtered_phased.vcf.gz.OK".format(**opts)
	cmd = ["/net/snowwhite/home/khlin/bin/vcftools --remove-indels --remove-filtered-all --max-alleles 2 --gzvcf /net/topmed2/working/hmkang/snps/4318/filtered/chr{chr}.filtered.rehdr.gt2.vcf.gz --gzdiff {outputDir}/HGDP_938/HGDP_938_chr{chr}_filtered_phased.vcf.gz \
--diff-site --stdout | awk '{{if($4 == \"B\")  print $1 \"\\t\" $2}}' > {outputDir}/common_site/topmed_chr{chr}_HGDP_common_global.txt".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)

for id in idList:
	opts["id"] = id
	######################
	#3.0.  based on global ancestry, create adaptive ref and ref_id to subset HGDP
	######################
	tgt = "{outputDirGlobal}/temp/{id}/{id}_HGDP_ref_id.txt.OK".format(**opts)
	dep = ""
	cmd = ["Rscript {outputDirGlobal}/global_select_HGDP.R id {id}".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)


	######################
	#3.1. subset of HGDP based on global ancestry
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDirGlobal}/temp/{id}/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz.OK".format(**opts)
		dep = "{outputDirGlobal}/temp/{id}/{id}_HGDP_ref_id.txt.OK".format(**opts)
		cmd = ["Rscript {outputDirGlobal}/global_select_HGDP.R id {id} && bcftools view {outputDir}/HGDP_938/HGDP_938_chr{chr}_filtered_phased.vcf.gz \
--force-samples -S {outputDirGlobal}/temp/{id}/{id}_HGDP_ref_id.txt \
--output-type z --output-file {outputDirGlobal}/temp/{id}/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz && \
bcftools index -t -f {outputDirGlobal}/temp/{id}/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#3.2. RFMix: merge sample and HGDP (filtered HGDP and keep only bi-allelic)
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDirGlobal}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz.OK".format(**opts)
		dep = "{outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz.OK {outputDirGlobal}/temp/{id}/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz.OK {outputDir}/common_site/topmed_chr{chr}_HGDP_common_global.txt.OK".format(**opts)
		cmd = ["bcftools merge -O u -R {outputDir}/common_site/topmed_chr{chr}_HGDP_common_global.txt \
{outputDir}/output/LAI/{id}/{id}_filtered_phased_chr{chr}.vcf.gz {outputDirGlobal}/temp/{id}/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz | \
bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o {outputDirGlobal}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz && \
bcftools index -t -f {outputDirGlobal}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#3.3. make files needed for RFMix allele, location, classes and snp position files
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}_RFMix_files_prep.OK".format(**opts)
		dep = "{outputDirGlobal}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz.OK".format(**opts)
		cmd = ["bcftools view {outputDirGlobal}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz \
| bcftools query -f '[%GT]\\n' - | sed 's/|//g' > {outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}.alleles && \
bcftools view {outputDirGlobal}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz \
| bcftools query -f '%POS\\n' - | Rscript {outputDir}/utilities/generate_markerLocations_file.R stdin {outputDir}/genetic_map_GRCh37/genetic_map_chr{chr}_combined_b37.txt {outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}.locations && \
bcftools query -l {outputDirGlobal}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz | Rscript {outputDir}/utilities/make_classes_file.R \
stdin {outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}.classes {outputDirGlobal}/temp/{id}/{id}_HGDP_ref.txt && \
bcftools query -f '%CHROM\t%POS\\n' {outputDirGlobal}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz > {outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}.pos".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#3.4. run RFMix
	######################
	inputFiles = [] #clean up
	inputFilesOK = [] #clean up
	for chr in xrange(1,23):
		opts["chr"] = chr
		inputFiles.append("{outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}.0.Viterbi.txt".format(**opts))
		inputFiles.append("{outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}.pos".format(**opts))
		tgt = "{outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}_RFMix_run.OK".format(**opts)
		inputFilesOK.append(tgt)
		dep = "{outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}_RFMix_files_prep.OK".format(**opts)
		cmd = ["cd {toolsDir}/RFMix_v1.5.4 && \
python ./RunRFMix.py PopPhased {outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}.alleles \
{outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}.classes \
{outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}.locations \
-o {outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr} --forward-backward --num-threads 1 && \
Rscript {outputDirGlobal}/global_convert_class2plot.R id {id} chr {chr}".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#3.4. make local ancestry plot
	######################
	outputFile = "{outputDir}/output/plot/{id}_HGDP_subset_global_local_ancestry_RFMix.png".format(**opts)
	opts["outputFile"] = outputFile
	opts["inputFiles"] = " ".join(inputFiles)
	tgt = "{outputDir}/output/plot/{id}_HGDP_subset_global_local_ancestry_RFMix.png.OK".format(**opts)
	dep = " ".join(inputFilesOK)
	cmd = ["Rscript {outputDir}/utilities/plot_local_ancestry_pipeline.R {id}_HGDP_subset_global_RFMix {outputFile} {inputFiles} 40000".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)

	######################
	#3.5. copy results and clean up
	######################
	cmds_temp = []
	tgt = "{outputDir}/temp/{id}/{id}_cleanup.OK".format(**opts)
	dep = "{outputDir}/output/plot/{id}_HGDP_local_ancestry_RFMix.png.OK {outputDir}/output/plot/{id}_HGDP_subset_local_ancestry_RFMix.png.OK {outputDir}/output/plot/{id}_HGDP_subset_local_ancestry_LAMPLD.png.OK".format(**opts)
	for chr in xrange(1,23):
		opts["chr"] = chr
		cmds_temp.append("cp {outputDir}/temp/{id}/{id}_HGDP_chr{chr}.0.Viterbi.txt {outputDir}/output/LAI/{id}/ && \
cp {outputDir}/temp/{id}/{id}_HGDP_chr{chr}.pos {outputDir}/output/LAI/{id}/ && \
cp {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.0.Viterbi.txt {outputDir}/output/LAI/{id}/ && \
cp {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.pos {outputDir}/output/LAI/{id}/ && \
cp {outputDir}/temp/{id}/{id}_chr{chr}_lampped.out {outputDir}/output/LAI/{id}/ && \
cp {outputDir}/temp/{id}/{id}_chr{chr}_lampped_plot.out {outputDir}/output/LAI/{id}/ && \
cp {outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz {outputDir}/output/LAI/{id}/ && \
cp {outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz.tbi {outputDir}/output/LAI/{id}/".format(**opts))
	cmds_temp.append("find temp/{id} ! -name '*.OK' -type f -exec rm -f {{}} +".format(**opts))
	makeJob("local", tgt, dep, cmds_temp)
######################
#2.0. log end time
######################
tgt = "{outputDir}/log/end.runmake.gai.OK".format(**opts)
dep = "{outputDirGlobal}/temp/{id}/{id}_HGDP_chr{chr}_RFMix_run.OK".format(**opts)
cmd = ["date | awk '{{print \"\\nend: \"$$0}}' >> {outputDir}/log/runmake_gai_time.log".format(**opts)]
makeJob("local", tgt, dep, cmd)


#*******************
#Write out make file
#*******************
#create makefile named makefile so that we have easy access to make clean
MAK = open("makefile".format(**opts), "w")
MAK.write(".DELETE_ON_ERROR:\n\n")
MAK.write("all: {tgts}\n\n".format(tgts=" ".join(tgts)))
#target makefile
MAK_tgt = open("{makeFile}".format(**opts), "w")
MAK_tgt.write(".DELETE_ON_ERROR:\n\n")
MAK_tgt.write("all: {tgts}\n\n".format(tgts=" ".join(tgts))) 

#clean
tgts.append("clean")
deps.append("")
cmds.append("\trm -f temp/NWD*/*.* output/GAI/NWD*/*.*")

#clean_job
tgts.append("clean_job")
deps.append("")
cmds.append("\tps xu | grep python | grep runmake | awk '{{print $$2}}' | xargs --verbose kill; scancel -n {jobName}\n".format(**opts))
 
for tgt,dep,cmd in zip(tgts, deps, cmds):
	MAK.write("{tgt} : {dep}\n".format(tgt=tgt, dep=dep))
	MAK.writelines(cmd)
	MAK.write("\n")
	MAK_tgt.write("{tgt} : {dep}\n".format(tgt=tgt, dep=dep))
	MAK_tgt.writelines(cmd)
	MAK_tgt.write("\n")

MAK.close()
MAK_tgt.close()
