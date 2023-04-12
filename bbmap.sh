#$ -cwd
#$ -S /bin/bash
#$ -N bbmap
#$ -V
#$ -pe smp 32
#$ -l h_vmem=8G
#$ -e logs/bbmap.error.log
#$ -o logs/bbmap.out.log


maindir="/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN"
# Full path of the directory where the base-called Illumina read files are located
Illuminadata="${maindir}/data/raw/Illumina/NovaSeq"
# Output directory for QC'd iSEQ data
IlluminadataQC="${maindir}/data/processed/Illumina/qc_reads"


###############################
## QC Illumina data using BBMap
###############################
ml BBMap/38.87-foss-2018b

for r1 in $Illuminadata/*R1*fastq.gz; do
r2=${r1/R1/R2}

# qc1
ref1="/nfs/modules/modules/software/BBMap/38.26-foss-2018b/resources/adapters.fa" #Full path to the file with adapters (comes with BBMap)
out1_adapt=${r1%%.fastq.gz}_adapt.1.fq.gz
out2_adapt=${r2%%.fastq.gz}_adapt.2.fq.gz
output_adaptMatch=${r1%%.fastq.gz}_adaptMatch.fq.gz
output_adaptSingle=${r1%%.fastq.gz}_adaptSingle.fq.gz
output_adaptStats=${r1%%.fastq.gz}_adapt.stats
log_adaptLog=${r1%%.fastq.gz}.adapt.log

# qc2
ref2="/nfs/modules/modules/software/BBMap/38.26-foss-2018b/resources/phix174_ill.ref.fa.gz" #Full path to the file with Phi X (comes with BBMap)
out1_phix=${r1%%.fastq.gz}_adapt_phix.1.fq.gz
out2_phix=${r2%%.fastq.gz}_adapt_phix.2.fq.gz
output_phixMatch=${r1%%.fastq.gz}_phixMatch.fq.gz
output_phixSingle=${r1%%.fastq.gz}_phixSingle.fq.gz
output_phixStats=${r1%%.fastq.gz}_phix.stats
log_phixLog=${r1%%.fastq.gz}.phix.log

# qc3
output_qc1=${r1%%.fastq.gz}_qc.1.fq.gz
output_qc2=${r2%%.fastq.gz}_qc.2.fq.gz
output_single=${r1%%.fastq.gz}_single.fq.gz
output_fail=${r1%%.fastq.gz}_fail.fq.gz
output_stats=${r1%%.fastq.gz}_qc.stats
log_log=${r1%%.fastq.gz}.qc.log

## Remove adapters (qc1)
bbduk.sh -Xmx2G usejni=t threads=14 overwrite=t qin=33 in1=$r1 in2=$r2 ref=$ref1 ktrim=r k=23 mink=11 hdist=1 out1=$out1_adapt out2=$out2_adapt outm=$output_adaptMatch outs=$output_adaptSingle refstats=$output_adaptStats statscolumns=5 2> $log_adaptLog;
## Remove Phi X (qc2)
bbduk.sh -Xmx2G usejni=t threads=4 overwrite=t qin=33 in1=$out1_adapt in2=$out2_adapt out1=$out1_phix out2=$out2_phix outm=$output_phixMatch outs=$output_phixSingle ref=$ref2 k=31 hdist=1 refstats=$output_phixStats statscolumns=5 2> $log_phixLog;
## Quality trim (qc3)
bbduk.sh -Xmx2G usejni=t threads=14 overwrite=t qin=33 in1=$out1_phix in2=$out2_phix fastawrap=10000 out1=$output_qc1 out2=$output_qc2 outm=$output_fail outs=$output_single minlength=45 maq=20 maxns=0 overwrite=t stats=$output_stats statscolumns=5 trimq=14 qtrim=r 2>$log_log;
## Move QC'd reads to processed folder
for x in $output_qc1 $output_qc2 $output_single; do mv $x $IlluminadataQC/; done
done
