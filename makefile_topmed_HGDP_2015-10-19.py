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

			cmd_tmp.append("\tsrun -p nomosix,main --exclude=sun02 -J {jobName} -D {outputDir} {param} {slurmScriptFile} \n".format(**opts))
		else:
			cmd_tmp.append("\tsrun -p nomosix,main --exclude=sun02 -J {jobName} -D {outputDir} {param} {command} \n".format(**opts))
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

######################
#0.1. log the start time
######################
tgt = "{outputDir}/log/start.runmake.shapeit.OK".format(**opts)
dep = ""
cmd = ["date | awk '{{print \"Local ancestry pipeline\\n\\nstart: \"$$0}}' > {outputDir}/log/runmake_shapeit_time.log".format(**opts)]
makeJob("local", tgt, dep, cmd)

for id in idList:
	opts["id"] = id
	if not os.path.exists("output/LAI/{id}".format(**opts)): os.makedirs("output/LAI/{id}".format(**opts)) #create LAI result folder
	if not os.path.exists("temp/{id}".format(**opts)): os.makedirs("temp/{id}".format(**opts)) #create temp folder
	######################
	#0.2. extract chromosomes of sample with bi-allelic snps
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_filtered_chr{chr}.vcf.gz.OK".format(**opts)
		dep = ""
		cmd = ["bcftools view /net/topmed2/working/gt-release/sftp-barnes/freeze.1a/topmed.freeze1.nhlbi.791.sftp-barnes.keep.chr{chr}.gtonly.vcf.gz \
--types snps -M2 --exclude-uncalled -f PASS --regions {chr} --force-samples -s {id} \
--output-type z --output-file {outputDir}/temp/{id}/{id}_filtered_chr{chr}.vcf.gz && \
bcftools index -t -f {outputDir}/temp/{id}/{id}_filtered_chr{chr}.vcf.gz".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#0.3. run shapeit to phase and exclude missing snps if necessary
	######################
	opts["param"] = "--mem=20000" # make srun sh to specify memomry

	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_phased_chr{chr}.OK".format(**opts)
		dep = "{outputDir}/temp/{id}/{id}_filtered_chr{chr}.vcf.gz.OK".format(**opts)
		cmd = ["{toolsDir}/shapeit.v2.r837/bin/shapeit -phase \
-V {outputDir}/temp/{id}/{id}_filtered_chr{chr}.vcf.gz \
-M {outputDir}/genetic_map_GRCh37/genetic_map_chr{chr}_combined_b37.txt \
--input-ref {outputDir}/1000GP_Phase3/1000GP_Phase3_chr{chr}.hap.gz {outputDir}/1000GP_Phase3/1000GP_Phase3_chr{chr}.legend.gz {outputDir}/1000GP_Phase3/1000GP_Phase3.sample \
-O {outputDir}/temp/{id}/{id}_phased_chr{chr} \
--output-log {outputDir}/temp/{id}/{id}_phased_chr{chr}.Output || \
{toolsDir}/shapeit.v2.r837/bin/shapeit -phase \
-V {outputDir}/temp/{id}/{id}_filtered_chr{chr}.vcf.gz \
-M {outputDir}/genetic_map_GRCh37/genetic_map_chr{chr}_combined_b37.txt \
--input-ref {outputDir}/1000GP_Phase3/1000GP_Phase3_chr{chr}.hap.gz {outputDir}/1000GP_Phase3/1000GP_Phase3_chr{chr}.legend.gz {outputDir}/1000GP_Phase3/1000GP_Phase3.sample \
--exclude-snp {outputDir}/temp/{id}/{id}_phased_chr{chr}.Output.snp.strand.exclude \
-O {outputDir}/temp/{id}/{id}_phased_chr{chr} \
--output-log {outputDir}/temp/{id}/{id}_phased_chr{chr}.Output".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)
	
	opts["param"] = "" #clean param

	######################
	#0.4. convert back to vcf.gz
	######################
	inputFiles = [] #clean up
	inputFilesOK = [] #clean up

	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz.OK".format(**opts)
		inputFilesOK.append(tgt)
		dep = "{outputDir}/temp/{id}/{id}_phased_chr{chr}.OK".format(**opts)
		cmd = ["{toolsDir}/shapeit.v2.r837/bin/shapeit -convert \
--input-haps {outputDir}/temp/{id}/{id}_phased_chr{chr} \
--output-vcf {outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf \
--output-log {outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf && \
bgzip -f {outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf && \
tabix -f -p vcf {outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

######################
#0.5. log end time
######################
tgt = "{outputDir}/log/end.runmake.shapeit.OK".format(**opts)
dep = " ".join(inputFilesOK)
cmd = ["date | awk '{{print \"\\nend: \"$$0}}' >> {outputDir}/log/runmake_shapeit_time.log".format(**opts)]
makeJob("local", tgt, dep, cmd)


#############################################
#RFMix local ancestry inference using all 7 ancestral populations
######################
######################
#1.0. filtered HGDP and keep only bi-allelic
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
#1.1. ind. common snp between sample and HGDP
######################
opts["id"] = idList[0] # use the first sample to decide the common snps sites between topmed and HGDP
for chr in xrange(1,23):
	opts["chr"] = chr
	tgt = "{outputDir}/common_site/topmed_chr{chr}_HGDP_common.txt.OK".format(**opts)
	dep = "{outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz.OK {outputDir}/HGDP_938/HGDP_938_chr{chr}_filtered_phased.vcf.gz.OK".format(**opts)
	cmd = ["/net/snowwhite/home/khlin/bin/vcftools --gzvcf {outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz --gzdiff {outputDir}/HGDP_938/HGDP_938_chr{chr}_filtered_phased.vcf.gz \
--diff-site --stdout | awk '{{if($4 == \"B\")  print $1 \"\\t\" $2}}' > {outputDir}/common_site/topmed_chr{chr}_HGDP_common.txt".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)

######################
#1.2. log the start time
######################
tgt = "{outputDir}/log/start.runmake.HGDP_RFMix.OK".format(**opts)
dep = ""
cmd = ["date | awk '{{print \"Local ancestry pipeline\\n\\nstart: \"$$0}}' > {outputDir}/log/runmake_HGDP_RFMix_time.log".format(**opts)]
makeJob("local", tgt, dep, cmd)

for id in idList:
	opts["id"] = id
	######################
	#1.3. merge sample and HGDP (filtered HGDP and keep only bi-allelic)
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_HGDP_chr{chr}_filtered_phased.vcf.gz.OK".format(**opts)
		dep = "{outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz.OK {outputDir}/common_site/topmed_chr{chr}_HGDP_common.txt.OK".format(**opts)
		cmd = ["bcftools merge -O u -R {outputDir}/common_site/topmed_chr{chr}_HGDP_common.txt \
{outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz {outputDir}/HGDP_938/HGDP_938_chr{chr}_filtered_phased.vcf.gz | \
bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o {outputDir}/temp/{id}/{id}_HGDP_chr{chr}_filtered_phased.vcf.gz && \
bcftools index -t -f {outputDir}/temp/{id}/{id}_HGDP_chr{chr}_filtered_phased.vcf.gz".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#1.4. make files needed for RFMix allele, location, classes and snp position files
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_HGDP_chr{chr}_RFMix_files_prep.OK".format(**opts)
		dep = "{outputDir}/temp/{id}/{id}_HGDP_chr{chr}_filtered_phased.vcf.gz.OK".format(**opts)
		cmd = ["bcftools view {outputDir}/temp/{id}/{id}_HGDP_chr{chr}_filtered_phased.vcf.gz \
| bcftools query -f '[%GT]\\n' - | sed 's/|//g' > {outputDir}/temp/{id}/{id}_HGDP_chr{chr}.alleles && \
bcftools view {outputDir}/temp/{id}/{id}_HGDP_chr{chr}_filtered_phased.vcf.gz \
| bcftools query -f '%POS\\n' - | Rscript {outputDir}/utilities/generate_markerLocations_file.R stdin {outputDir}/genetic_map_GRCh37/genetic_map_chr{chr}_combined_b37.txt {outputDir}/temp/{id}/{id}_HGDP_chr{chr}.locations && \
bcftools query -l {outputDir}/temp/{id}/{id}_HGDP_chr{chr}_filtered_phased.vcf.gz | Rscript {outputDir}/utilities/make_classes_file.R \
stdin {outputDir}/temp/{id}/{id}_HGDP_chr{chr}.classes && \
bcftools query -f '%CHROM\t%POS\\n' {outputDir}/temp/{id}/{id}_HGDP_chr{chr}_filtered_phased.vcf.gz > {outputDir}/temp/{id}/{id}_HGDP_chr{chr}.pos".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#1.5. run RFMix
	######################
	inputFiles = [] #clean up
	inputFilesOK = [] #clean up
	for chr in xrange(1,23):
		opts["chr"] = chr
		inputFiles.append("{outputDir}/temp/{id}/{id}_HGDP_chr{chr}.0.Viterbi.txt".format(**opts))
		inputFiles.append("{outputDir}/temp/{id}/{id}_HGDP_chr{chr}.pos".format(**opts))
		tgt = "{outputDir}/temp/{id}/{id}_HGDP_chr{chr}_RFMix_run.OK".format(**opts)
		inputFilesOK.append(tgt)
		dep = "{outputDir}/temp/{id}/{id}_HGDP_chr{chr}_RFMix_files_prep.OK".format(**opts)
		cmd = ["cd {toolsDir}/RFMix_v1.5.4 && \
python ./RunRFMix.py PopPhased {outputDir}/temp/{id}/{id}_HGDP_chr{chr}.alleles \
{outputDir}/temp/{id}/{id}_HGDP_chr{chr}.classes \
{outputDir}/temp/{id}/{id}_HGDP_chr{chr}.locations \
-o {outputDir}/temp/{id}/{id}_HGDP_chr{chr} --forward-backward --num-threads 1".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#1.6. make local ancestry plot
	outputFile = "{outputDir}/output/plot/{id}_HGDP_local_ancestry_RFMix.png".format(**opts)
	opts["outputFile"] = outputFile
	opts["inputFiles"] = " ".join(inputFiles)
	tgt = "{outputDir}/output/plot/{id}_HGDP_local_ancestry_RFMix.png.OK".format(**opts)
	dep = " ".join(inputFilesOK)
	cmd = ["Rscript {outputDir}/utilities/plot_local_ancestry_pipeline.R {id}_HGDP_RFMix {outputFile} {inputFiles} 40000".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)

######################
#1.7. log end time
######################
tgt = "{outputDir}/log/end.runmake.HGDP_RFMix.OK".format(**opts)
dep = "{outputDir}/output/plot/{id}_HGDP_local_ancestry_RFMix.png.OK".format(**opts)
cmd = ["date | awk '{{print \"\\nend: \"$$0}}' >> {outputDir}/log/runmake_HGDP_RFMix_time.log".format(**opts)]
makeJob("local", tgt, dep, cmd)


#############################################
#RFMix local ancestry inference using only europe, africa, native american
#############################################


######################
#2.0. Europe, Africa, Native America subset of HGDP
######################
for chr in xrange(1,23):
	opts["chr"] = chr
	tgt = "{outputDir}/HGDP_938/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz.OK".format(**opts)
	dep = ""
	cmd = ["bcftools view {outputDir}/HGDP_938/HGDP_938_chr{chr}_filtered_phased.vcf.gz \
--force-samples -S {outputDir}/HGDP_938/HGDP_europe_africa_native_america.txt \
--output-type z --output-file {outputDir}/HGDP_938/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz && \
bcftools index -t -f {outputDir}/HGDP_938/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)

######################
#2.1. ind. common snp between sample and HGDP subset
######################
opts["id"] = idList[0] # use the first sample to decide the common snps sites between topmed and HGDP
for chr in xrange(1,23):
	opts["chr"] = chr
	tgt = "{outputDir}/common_site/topmed_chr{chr}_HGDP_subset_common.txt.OK".format(**opts)
	dep = "{outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz.OK {outputDir}/HGDP_938/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz.OK".format(**opts)
	cmd = ["/net/snowwhite/home/khlin/bin/vcftools --gzvcf {outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz --gzdiff {outputDir}/HGDP_938/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz \
--diff-site --stdout | awk '{{if($4 == \"B\")  print $1 \"\\t\" $2}}' > {outputDir}/common_site/topmed_chr{chr}_HGDP_subset_common.txt".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)

######################
#2.2. log the start time
######################
tgt = "{outputDir}/log/start.runmake.HGDP_subset_RFMix.OK".format(**opts)
dep = ""
cmd = ["date | awk '{{print \"Local ancestry pipeline\\n\\nstart: \"$$0}}' > {outputDir}/log/runmake_HGDP_subset_RFMix_time.log".format(**opts)]
makeJob("local", tgt, dep, cmd)

for id in idList:
	opts["id"] = id
	######################
	#2.3. merge sample and HGDP subset (filtered HGDP and keep only bi-allelic)
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz.OK".format(**opts)
		dep = "{outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz.OK {outputDir}/common_site/topmed_chr{chr}_HGDP_subset_common.txt.OK".format(**opts)
		cmd = ["bcftools merge -O u -R {outputDir}/common_site/topmed_chr{chr}_HGDP_subset_common.txt \
{outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz {outputDir}/HGDP_938/HGDP_938_chr{chr}_subset_filtered_phased.vcf.gz | \
bcftools view --types snps -M2 --exclude-uncalled -f PASS -O z -o {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz && \
bcftools index -t -f {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#2.4. make files needed for RFMix allele, location, classes and snp position files
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_RFMix_files_prep.OK".format(**opts)
		dep = "{outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz.OK".format(**opts)
		cmd = ["bcftools view {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz \
| bcftools query -f '[%GT]\\n' - | sed 's/|//g' > {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.alleles && \
bcftools view {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz \
| bcftools query -f '%POS\\n' - | Rscript {outputDir}/utilities/generate_markerLocations_file.R stdin {outputDir}/genetic_map_GRCh37/genetic_map_chr{chr}_combined_b37.txt {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.locations && \
bcftools query -l {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz | Rscript {outputDir}/utilities/make_classes_file.R \
stdin {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.classes && \
sed 's/4/2/g' {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.classes | sed 's/5/3/g' > {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.classes.tmp && \
mv {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.classes.tmp {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.classes && \
bcftools query -f '%CHROM\t%POS\\n' {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_filtered_phased.vcf.gz > {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.pos".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#2.5. run RFMix
	######################
	inputFiles = [] #clean up
	inputFilesOK = [] #clean up
	for chr in xrange(1,23):
		opts["chr"] = chr
		inputFiles.append("{outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.0.Viterbi.txt".format(**opts))
		inputFiles.append("{outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.pos".format(**opts))
		tgt = "{outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_RFMix_run.OK".format(**opts)
		inputFilesOK.append(tgt)
		dep = "{outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}_RFMix_files_prep.OK".format(**opts)
		cmd = ["cd {toolsDir}/RFMix_v1.5.4 && \
python ./RunRFMix.py PopPhased {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.alleles \
{outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.classes \
{outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.locations \
-o {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr} --forward-backward --num-threads 1 && \
sed 's/2/4/g' {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.0.Viterbi.txt | sed 's/3/5/g' > {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.0.Viterbi.txt.tmp && \
mv {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.0.Viterbi.txt.tmp {outputDir}/temp/{id}/{id}_HGDP_subset_chr{chr}.0.Viterbi.txt".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#2.6. make local ancest
	######################
	outputFile = "{outputDir}/output/plot/{id}_HGDP_subset_local_ancestry_RFMix.png".format(**opts)
	opts["outputFile"] = outputFile
	opts["inputFiles"] = " ".join(inputFiles)
	tgt = "{outputDir}/output/plot/{id}_HGDP_subset_local_ancestry_RFMix.png.OK".format(**opts)
	dep = " ".join(inputFilesOK)
	cmd = ["Rscript {outputDir}/utilities/plot_local_ancestry_pipeline.R {id}_HGDP_subset_RFMix {outputFile} {inputFiles} 40000".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)

######################
#2.7. log end time
######################
tgt = "{outputDir}/log/end.runmake.HGDP_subset_RFMix.OK".format(**opts)
dep = "{outputDir}/output/plot/{id}_HGDP_subset_local_ancestry_RFMix.png.OK".format(**opts)
cmd = ["date | awk '{{print \"\\nend: \"$$0}}' >> {outputDir}/log/runmake_HGDP_subset_RFMix_time.log".format(**opts)]
makeJob("local", tgt, dep, cmd)

# ##########################################################
# #local ancestry inference using LAMPLD
# #############################################


######################
#3.0. create subset reference panel by ancestral populations
######################
tgt = "{outputDir}/HGDP_938/LAMPLD/create_LAMP_ref_hap.OK".format(**opts)
dep = ""
cmd = ["./HGDP_938/create_LAMP_ref_hap.sh".format(**opts)]
makeJob("local", tgt, dep, cmd)

######################
#3.1. log the start time
######################
tgt = "{outputDir}/log/start.runmake.HGDP_subset_LAMPLD.OK".format(**opts)
dep = ""
cmd = ["date | awk '{{print \"Local ancestry pipeline using HGDP subset \\n\\nstart: \"$$0}}' > {outputDir}/log/runmake_HGDP_subset_LAMPLD_time.log".format(**opts)]
makeJob("local", tgt, dep, cmd)

######################
#3.2. convert sample vcf with common snp to HGDP to genotype dosage
######################
for id in idList:
	opts["id"] = id
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_chr{chr}_lamp.012.OK".format(**opts)
		dep = "{outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz.OK {outputDir}/common_site/topmed_chr{chr}_HGDP_subset_common.txt.OK".format(**opts)
		cmd = ["vcftools --gzvcf {outputDir}/temp/{id}/{id}_filtered_phased_chr{chr}.vcf.gz --012 --out {outputDir}/temp/{id}/{id}_chr{chr}_lamp --positions {outputDir}/common_site/topmed_chr{chr}_HGDP_subset_common.txt && \
rm -f {outputDir}/temp/{id}/{id}_chr{chr}_lamp.012.indv {outputDir}/temp/{id}/{id}_chr{chr}_lamp.012.pos {outputDir}/temp/{id}/{id}_chr{chr}_lamp.log && \
cut -f 2- {outputDir}/temp/{id}/{id}_chr{chr}_lamp.012 | sed \"s/\\t//g\" > {outputDir}/temp/{id}/{id}_chr{chr}_lamp.012.tmp && mv {outputDir}/temp/{id}/{id}_chr{chr}_lamp.012.tmp {outputDir}/temp/{id}/{id}_chr{chr}_lamp.012".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#3.3. run LAMPLD
	######################
	for chr in xrange(1,23):
		opts["chr"] = chr
		tgt = "{outputDir}/temp/{id}/{id}_chr{chr}_lampped.out.OK".format(**opts)
		dep = "{outputDir}/HGDP_938/LAMPLD/create_LAMP_ref_hap.OK {outputDir}/temp/{id}/{id}_chr{chr}_lamp.012.OK".format(**opts)
		cmd = ["perl {toolsDir}/LAMPLD-v1.1/run_LAMPLD.pl {outputDir}/HGDP_938/LAMPLD/chr{chr}.pos \
{outputDir}/HGDP_938/LAMPLD/HGDP_europe_chr{chr}.impute.hap \
{outputDir}/HGDP_938/LAMPLD/HGDP_native_america_chr{chr}.impute.hap \
{outputDir}/HGDP_938/LAMPLD/HGDP_africa_chr{chr}.impute.hap \
{outputDir}/temp/{id}/{id}_chr{chr}_lamp.012 \
{outputDir}/temp/{id}/{id}_chr{chr}_lampped.out".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd);

	######################
	#3.4. convert plot format 
	######################
	inputFiles = [] #clean up
	inputFilesOK = [] #clean up
	for chr in xrange(1,23):
		opts["chr"] = chr
		inputFiles.append("{outputDir}/temp/{id}/{id}_chr{chr}_lampped_plot.out".format(**opts))
		inputFiles.append("{outputDir}/common_site/topmed_chr{chr}_HGDP_subset_common.txt".format(**opts))
		tgt = "{outputDir}/temp/{id}/{id}_chr{chr}_lampped_plot.out.OK".format(**opts)
		inputFilesOK.append(tgt)
		dep = "{outputDir}/temp/{id}/{id}_chr{chr}_lampped.out.OK".format(**opts)
		cmd = ["python {outputDir}/utilities/lamped_out_2_plot.py -i {outputDir}/temp/{id}/{id}_chr{chr}_lampped.out \
-p {outputDir}/common_site/topmed_chr{chr}_HGDP_subset_common.txt -o {outputDir}/temp/{id}/{id}_chr{chr}_lampped_plot.out".format(**opts)]
		makeJob(launchMethod, tgt, dep, cmd)

	######################
	#3.5. make local ancestry plot
	######################
	outputFile = "{outputDir}/output/plot/{id}_HGDP_subset_local_ancestry_LAMPLD.png".format(**opts)
	opts["outputFile"] = outputFile
	opts["inputFiles"] = " ".join(inputFiles)
	tgt = "{outputDir}/output/plot/{id}_HGDP_subset_local_ancestry_LAMPLD.png.OK".format(**opts)
	dep = " ".join(inputFilesOK)
	cmd = ["Rscript {outputDir}/utilities/plot_local_ancestry_pipeline.R {id}_HGDP_subset_LAMPLD {outputFile} {inputFiles} 40000".format(**opts)]
	makeJob(launchMethod, tgt, dep, cmd)

######################
#3.6. log end time
######################
tgt = "{outputDir}/log/end.runmake.HGDP_subset_LAMPLD.OK".format(**opts)
dep = "{outputDir}/output/plot/{id}_HGDP_subset_local_ancestry_LAMPLD.png.OK".format(**opts)
cmd = ["date | awk '{{print \"\\nend: \"$$0}}' >> {outputDir}/log/runmake_HGDP_subset_LAMPLD_time.log".format(**opts)]
makeJob("local", tgt, dep, cmd)


######################
#4.0. copy results and clean up
######################
for id in idList:
	opts["id"] = id
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
cmds.append("\trm -f temp/NWD*/*.* output/LAI/NWD*/*.*")

#clean_job
tgts.append("clean_job")
deps.append("")
cmds.append("\tps xu | grep make | grep {jobName} | awk '{{print $$2}}' | xargs --verbose kill; scancel -n {jobName}\n".format(**opts))
 
for tgt,dep,cmd in zip(tgts, deps, cmds):
	MAK.write("{tgt} : {dep}\n".format(tgt=tgt, dep=dep))
	MAK.writelines(cmd)
	MAK.write("\n")
	MAK_tgt.write("{tgt} : {dep}\n".format(tgt=tgt, dep=dep))
	MAK_tgt.writelines(cmd)
	MAK_tgt.write("\n")

MAK.close()
MAK_tgt.close()
