#!/bin/bash
job_name="topmedlai"
while getopts i:j:l: opt; do
  case $opt in
  i)
      ID=$OPTARG
      ;;
  j)
      job_name=$OPTARG
      ;;
  l)
      launchMethod=$OPTARG
      ;;  esac
done

#run the pipeline
[ -z "$ID" ] && echo need to give an id or a file of list \*.txt \"./runmake.sh -i xxxxxx\" && exit 1
nohup python runmake_topmed_HGDP.py -i ${ID} -j ${job_name}> log/runmake_${job_name}.log  2>&1 &
