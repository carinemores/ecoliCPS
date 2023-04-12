#$ -cwd
#$ -S /bin/bash
#$ -N flye
#$ -V
#$ -pe smp 32
#$ -l h_vmem=8G
#$ -e logs/flye.error.log
#$ -o logs/flye.out.log
#$ -t 1-24

IDS=(echo {01..24})
ID=barcode${IDS[$SGE_TASK_ID]}

ml Flye/2.9-foss-2018b-Python-2.7.15

mkdir -p scratch/flye
mkdir -p scratch/flye/$ID

flye --nano-raw scratch/basecall/$ID\_filtered.fastq -o scratch/flye/$ID/ -t 32
