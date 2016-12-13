#!/bin/bash
cat log/runmake_topmedlai.log | grep runn
echo "# of touch files"
cat log/runmake_topmedlai.log | grep touch | wc -l
echo "# of Error files"
cat log/runmake_topmedlai.log | grep Error | wc -l
echo "# of slurm jobs in queue"
squeue -u khlin | grep topmed | grep PD | wc -l
echo "# of slurm running jobs"
squeue -u khlin | grep topmed | grep R | wc -l