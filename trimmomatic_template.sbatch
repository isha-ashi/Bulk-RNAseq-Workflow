#!/bin/bash
#
#SBATCH -J trim    # job name
#SBATCH -o LOGPATHOLDER/trimmomatic-%A_%a.out.o%j       # output and error file name (%j expands to jobID)
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
module load trimmomatic/0.38

TRIMMOMATICDRI=`which trimmomatic-0.38.jar|xargs dirname`

DATAPATH=dataPathPlaceHolder
FASTQNAMES=($(cat $DATAPATH/sampleList.txt))
 
echo "Trimming fastq files:"
echo ${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}

INTERIMPATH=$DATAPATH/intermediate
TRIMOUTPATH=$DATAPATH/trimmed

## Make folder to save TRIMMOMATIC results
mkdir -p $INTERIMPATH
mkdir -p $TRIMOUTPATH

## The job command:
java -jar ${TRIMMOMATICDRI}/trimmomatic-0.38.jar PE -threads threadsHolder -phred33 \
          inputFileHolder \
          outFileHolder \
          ILLUMINACLIP:${TRIMMOMATICDRI}/adapters/fastaFileHolder:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50

echo "Finished at:"
echo $(date)
