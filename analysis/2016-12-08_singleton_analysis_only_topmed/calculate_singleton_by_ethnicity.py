#!/usr/bin/env python
#!/usr/bin/python
#!python

import sys
import os
import re
import getopt
import numpy
import subprocess
from string import Template

######################
#Initialization
######################
idList=[]
unrelatedList=[]
ethnicityList=[] #unrelated samples only

#create temp and rseults directory
if not os.path.exists("temp"): os.makedirs("temp")
if not os.path.exists("result"): os.makedirs("result")

######################
#read in ethnicity ancestry and unrelated list
######################
##unrelated list
infh = open("unrelated_degree_2_freeze3a.txt")
for line in infh.xreadlines():
	item = line.strip().split()
	unrelatedList.append(item[0])
infh.close()

##global ancestry list
infh = open("ethnicity_max_topmed_only.txt")
for line in infh.xreadlines():
	sample_id, ethnicity = line.strip().split("\t")
	if sample_id in unrelatedList: #only unrelated
		ethnicityList.append([sample_id, ethnicity])
infh.close()
print ethnicityList
print len(ethnicityList)

######################
#output unrelated sample list by ancestry
######################
import random
#set seed, make sure results are tractable
random.seed(464351)
#shuffle the ethnicity list
random.shuffle(ethnicityList)
#output each ethnicity
for i in xrange(1,8):
	outfh = open("temp/unrelated_{ethnicity}.txt".format(ethnicity=i), 'w')
	for j in ethnicityList:
		if int(j[1]) == i:
			outfh.writelines("\t".join(j)+"\n")
	outfh.close()

##########
#functions for makefile pipeline
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
		if re.search("\||^cd|&&", c):
			opts["slurmScriptNo"] += 1
			opts["slurmScriptFile"] = "{slurmScriptsDir}/{slurmScriptNo}_{jobName}.sh".format(**opts)
			IN = open("{slurmScriptFile}".format(**opts), "w")
			IN.write("#!/bin/bash\n")
			IN.write("set pipefail; {command}".format(**opts))
			IN.close()
			os.chmod("{slurmScriptFile}".format(**opts), 0755)

			cmd_tmp.append("\tsrun -p topmed,nomosix -J {jobName} -D {outputDir} {param} {slurmScriptFile} \n".format(**opts))
		else:
			cmd_tmp.append("\tsrun -p topmed,nomosix -J {jobName} -D {outputDir} {param} {command} \n".format(**opts))
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

# ######################
# #count singleton
# ######################
# import subprocess
# import numpy
# resultList=[]
# idList=[]

# opts={}
# opts["ethnicityIdx"]=1
# opts["size"]=200
# opts["chr"]=22

# ##take the unrelated samples by ancestry
# infh = open("temp/unrelated_{ethnicityIdx}.txt".format(**opts))
# for line in infh:
# 	items = line.strip().split()
# 	idList.append(items[0])
# infh.close()
##create temp_id
# for i in numpy.arange(opts["size"], len(idList), step=opts["size"]):
# 	opts["cureent_size"] = i
# 	outfh = open("temp/id_{ethnicityIdx}_{cureent_size}_{chr}.txt".format(**opts), 'w')
# 	outfh.writelines('\n'.join(idList[0:i])) #first i samples
# 	outfh.close()

# 	cmd = "bcftools view --type snps -S temp/id_{ethnicityIdx}_{cureent_size}_{chr}.txt --force-samples -G /net/topmed2/working/khlin/topmed.freeze3.phased.hgdp.snp/topmed.freeze3.chr{chr}.phased.hgdp.snp.filtered.vcf.gz | bcftools query -i 'INFO/AC==1' -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' | wc -l | awk '{{print \"{ethnicityIdx}\t{cureent_size}\t{chr}\t\"$1}}' > result/count_{ethnicityIdx}_{cureent_size}_{chr}.txt".format(**opts)
	# subprocess.call(cmd, shell=True)
# cmd = "bcftools view --type snps -G -S temp/id_temp.txt /net/topmed3/working/hmkang/freeze2/10597.v2/final/topmed_freeze2_10597.chr{chr}.overlap_removed.svm_pass.genotypes.vcf.gz | bcftools query -i 'INFO/AC==1' -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' | wc -l".format(chr=chr)
 # result = subprocess.check_output(cmd, shell=True)

######################
##use makefile pipeline to be efficient
######################
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
opts={}

##additional parameters
launchMethod = "local" #need to use $$ if running local since make needs escape for $
opts["outputDir"] = os.getcwd()
opts["jobName"] = "singleton"
opts["makeFile"] = "makefile_{jobName}".format(**opts)
exclude = "--exclude=topmed,topmed5,topmed6,topmed7,topmed8" #,\"r[6301-6315]\""
opts["param"] = "{exclude} --time=0-12:0".format(exclude = exclude) #indicate this is a quick job

#create directory needed for slurm script 
opts["slurmScriptsDir"] = "{outputDir}/slurm_scripts".format(**opts)
if not os.path.exists(opts["slurmScriptsDir"]): os.makedirs(opts["slurmScriptsDir"])
opts["slurmScriptNo"] = 0

######################
#1.1. count singleton
######################
opts["size"]=200 ##increase the sample size by xxx samples
for a in [1,4,5]:
	opts["ethnicityIdx"] = a

	##take the unrelated samples by ancestry
	infh = open("temp/unrelated_{ethnicityIdx}.txt".format(**opts))
	idList=[]
	for line in infh:
		items = line.strip().split()
		idList.append(items[0])
	infh.close()
	#create sample id by ethnicity, current_size, chr
	for b in numpy.arange(opts["size"], len(idList), step=opts["size"]):
		opts["cureent_size"] = b
		for c in numpy.arange(1,23):
			opts["chr"] = c
			outfh = open("temp/id_{ethnicityIdx}_{cureent_size}_{chr}.txt".format(**opts), 'w')
			outfh.writelines('\n'.join(idList[0:opts["cureent_size"]])) #first i samples
			outfh.close()

			tgt = "result/count_{ethnicityIdx}_{cureent_size}_{chr}.txt.OK".format(**opts)
			inputFilesOK.append(tgt)
			dep = ""
			cmd = ["bcftools view --type snps -S temp/id_{ethnicityIdx}_{cureent_size}_{chr}.txt --force-samples -G /net/topmed4/working/hmkang/freeze3a/v4/release/passgt.minDP10/ALL.chr{chr}.freeze3a.pass.gtonly.genotypes.bcf | bcftools query -i 'INFO/AC==1' -f '%CHROM\\t%POS\\t%INFO/AC\\t%INFO/AN\\n' | wc -l | awk '{{print \"{ethnicityIdx}\\t{cureent_size}\\t{chr}\\t\"$$1}}' > {outputDir}/result/count_{ethnicityIdx}_{cureent_size}_{chr}.txt".format(**opts)]
			# cmd = ["bcftools view --type snps -S temp/id_{ethnicityIdx}_{cureent_size}_{chr}.txt --force-samples -G /net/topmed2/working/khlin/topmed.freeze3.phased.hgdp.snp/topmed.freeze3.chr{chr}.phased.hgdp.snp.filtered.vcf.gz | bcftools query -i 'INFO/AC==1' -f '%CHROM\\t%POS\\t%INFO/AC\\t%INFO/AN\\n' | wc -l | awk '{{print \"{ethnicityIdx}\\t{cureent_size}\\t{chr}\\t\"$$1}}' > {outputDir}/result/count_{ethnicityIdx}_{cureent_size}_{chr}.txt".format(**opts)]  ##use a smaller file for test run
			makeJob(launchMethod, tgt, dep, cmd)

#create all unrelated sample id by ecurrent_size, chr. 99=all unrelated individuals
#set seed, make sure results are tractable
random.seed(987362)
#shuffle the ethnicity list
random.shuffle(ethnicityList)
unrelated_ethnicityList = [x[0] for x in ethnicityList]
opts["size"]=200
for b in numpy.arange(opts["size"], len(unrelated_ethnicityList), step=opts["size"]):
	opts["cureent_size"] = b
	for c in numpy.arange(1,23):
		opts["chr"] = c
		outfh = open("temp/id_{cureent_size}_{chr}.txt".format(**opts), 'w')
		outfh.writelines('\n'.join(unrelated_ethnicityList[0:opts["cureent_size"]])) #first i samples
		outfh.close()

		tgt = "result/count_{cureent_size}_{chr}.txt.OK".format(**opts)
		inputFilesOK.append(tgt)
		dep = ""
		cmd = ["bcftools view --type snps -S temp/id_{cureent_size}_{chr}.txt --force-samples -G /net/topmed4/working/hmkang/freeze3a/v4/release/passgt.minDP10/ALL.chr{chr}.freeze3a.pass.gtonly.genotypes.bcf | bcftools query -i 'INFO/AC==1' -f '%CHROM\\t%POS\\t%INFO/AC\\t%INFO/AN\\n' | wc -l | awk '{{print \"99\\t{cureent_size}\\t{chr}\\t\"$$1}}' > {outputDir}/result/count_{cureent_size}_{chr}.txt".format(**opts)]
		# cmd = ["bcftools view --type snps -S temp/id_{ethnicityIdx}_{cureent_size}_{chr}.txt --force-samples -G /net/topmed2/working/khlin/topmed.freeze3.phased.hgdp.snp/topmed.freeze3.chr{chr}.phased.hgdp.snp.filtered.vcf.gz | bcftools query -i 'INFO/AC==1' -f '%CHROM\\t%POS\\t%INFO/AC\\t%INFO/AN\\n' | wc -l | awk '{{print \"{ethnicityIdx}\\t{cureent_size}\\t{chr}\\t\"$$1}}' > {outputDir}/result/count_{ethnicityIdx}_{cureent_size}_{chr}.txt".format(**opts)]  ##use a smaller file for test run
		makeJob(launchMethod, tgt, dep, cmd)


##combine into one final file for analysis
tgt = "{outputDir}/singleton_count.txt.OK".format(**opts)
dep = " ".join(inputFilesOK)
cmd = ["cat result/*.* > {outputDir}/singleton_count.txt".format(**opts)]
makeJob(launchMethod, tgt, dep, cmd)

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
#tgts.append("clean")
#deps.append("")
#cmds.append("\trm -f temp/NWD*/*.* output/LAI/NWD*/*.*")

#clean_job
tgts.append("clean_job")
deps.append("")
cmds.append("\tscancel -n {jobName}; ps xu | grep make | grep {jobName} | awk '{{print $$2}}' | xargs --verbose kill\n".format(**opts))
 
for tgt,dep,cmd in zip(tgts, deps, cmds):
	MAK.write("{tgt} : {dep}\n".format(tgt=tgt, dep=dep))
	MAK.writelines(cmd)
	MAK.write("\n")
	MAK_tgt.write("{tgt} : {dep}\n".format(tgt=tgt, dep=dep))
	MAK_tgt.writelines(cmd)
	MAK_tgt.write("\n")

MAK.close()
MAK_tgt.close()
