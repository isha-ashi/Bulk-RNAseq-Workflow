#!/bin/bash
#
#SBATCH -J qc    # job name
#SBATCH -o LOGPATHOLDER/fastqc-%A_%a.out.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1                # total number of mpi tasks requested
#SBATCH -N 1                # single node use
##SBATCH -p normal          # queue (partition) -- normal, development, etc.
#SBATCH -t hoursPlaceHolder:00:00          # run time (hh:mm:ss)
#SBATCH --mail-user=emailHolder
#SBATCH --mail-type=FAIL     # email me when the job fails
#SBATCH --array=0-arrayLengthHolder # job array index

echo "Started at:"
echo $(date)

module purge
module load fastqc/0.11.9

DATAPATH=dataPathPlaceHolder
FASTQNAMES=($(cat $DATAPATH/sampleList.txt))
 

echo "Doing fastqc analysis on:"
echo ${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}

QCOUTPATH=qcFolderNameHolder

## Make folder to save FASTQC results
mkdir -p $QCOUTPATH

## The job command:
jobCommandHolder1
jobCommandHolder2

echo "Finished at:"
echo $(date)
