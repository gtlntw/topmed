#!/bin/bash
#default: id=NA19625
id="."
launchMethod="slurm"
job_name="topmedlai"
jobNo=200 #number of jobs in parallel

while getopts i:l: opt; do
  case $opt in
  i)
      id=$OPTARG
      ;;
  l)
      launchMethod=$OPTARG
      ;;  esac
done

#clean start and end log
rm -rf log/start*.OK log/end*.OK slurm_scripts/*.*

#command to run
cmd='python makefile_topmed_HGDP.py -l ${launchMethod} -j ${job_name} -i ${id};
	 nohup make -f makefile_${job_name} -j ${jobNo} --keep-going > log/runmake_${job_name}.log 2>&1 &'

echo "running ${jobNo} parallel jobs"
eval $cmd