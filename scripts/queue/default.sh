#!/bin/bash
#$ -S /bin/bash
#$ -wd WORKDIR
#$ -N JOBNAME
#$ -o SCRNAME.$JOB_ID.log
#$ -pe PENAME SLOTS -l vf=MEM
#$ -p PRIORITY

COMMAND
