#!/bin/bash
ID="."
job_name="topmedgai"
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
nohup python runmake_topmed_global.py -i ${ID} > log/runmake_${job_name}.log  2>&1 &
