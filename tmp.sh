#!/bin/bash
#$ -cwd
#$ -o Stdout.$JOB_ID
#$ -e Stderr.$JOB_ID
#$ -j n
#$ -pe shared 1
#$ -l h_data=4g,h_rt=8:00:00
#$ -m n
source /u/local/Modules/default/init/modules.sh
source /u/home/j/jshamsh1/.bashrc
conda activate r_env
which conda
which R
R --version
R -e "run = 1; source('sim_weak.R')"
