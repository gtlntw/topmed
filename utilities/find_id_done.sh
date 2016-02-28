#!/bin/sh

#create id_done.txt
find ./temp -name '*_cleanup.OK' | cut -d / -f 3 | sort > id_done.txt
#generate id_to_run.txt
comm -23 id_10597.txt id_done.txt > id_to_run.txt