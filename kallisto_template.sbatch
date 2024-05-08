#!/bin/bash
#
#SBATCH -J kallisto_align    # job name
#SBATCH -o LOGPATHOLDER/kallisto-%A_%a.out.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1                # total number of tasks requested
#SBATCH -N 1                # single node use
##SBATCH -p normal          # queue (partition) -- normal, development, etc.
#SBATCH -t hoursPlaceHolder:00:00          # run time (hh:mm:ss)
#SBATCH --mail-user=emailHolder
#SBATCH --mail-type=FAIL     # email me when the job fails
#SBATCH --array=0-arrayLengthHolder # job array index

echo "Started at:"
echo $(date)

## Loading kallisto module
module purge
module load kallisto/0.46.1

DATAPATH=dataPathPlaceHolder
OUTPATH=outPathPlaceHolder 
mkdir -p $OUTPATH

FASTQNAMES=($(cat $DATAPATH/sampleList.txt))
echo ${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}

kallisto quant -i TRANSCRIPTSINDEXPATH \
               -b 100 -o $OUTPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]} \
	             $DATAPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}_R1*.fastq.gz $DATAPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}_R2*.fastq.gz
	             
echo "Finished at:"
echo $(date)
