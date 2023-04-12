#$ -cwd
#$ -S /bin/bash
#$ -N flye
#$ -V
#$ -pe smp 32
#$ -l h_vmem=8G
#$ -e logs/flye.error.log
#$ -o logs/flye.out.log

# PacBio HiFi reads assembly - Using Flye
# Version 1.0
set -e

########################
# Set working directory
########################
maindir="/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/data/raw/PacBio/fastq"

#####################
# Assembly using Flye
#####################
ml Flye/2.9-foss-2018b-Python-2.7.15


for i in $maindir/*.fastq; do
out=${i%%.hifi_reads.fastq}
flye --pacbio-hifi $i -o $out -t 32
done
