#$ -cwd
#$ -S /bin/bash
#$ -N pilon
#$ -V
#$ -pe smp 32
#$ -l h_vmem=8G
#$ -e logs/pilon2.error.log
#$ -o logs/pilon2.out.log
#$ -t 1-65


# Polish flye assemblies using pilon

ml BWA/0.7.17-foss-2018b
ml SAMtools/1.9-foss-2018b
ml Pilon/1.23-Java-1.8


########### remember to change the -t parameter (match the number of samples)

IIDS=("" "K3" "K4" "K5" "K6" "K7" "K8" "K9" "K10" "K11" "K12" "K13" "K14" "K16" "K17" "K18a" "K18ab" "K19" "K20" "K22" "K23" "K24" "K26" "K27" "K28" "K30" "K31" "K32" "K33" "K34" "K35" "K36" "K37" "K38" "K39" "K40" "K41" "K42" "K43" "K44" "K45" "K46" "K47" "K48" "K49" "K50" "K51" "K52" "K53" "K54" "K55" "K57" "K83" "K84" "K87" "K92" "K93" "K94" "K95" "K96" "K97" "K98" "K100" "K101" "K102" "K103")
IID=${IIDS[$SGE_TASK_ID]}

NIDS=("" "bc2004" "bc2005" "bc2006" "bc2007" "bc2008" "bc2009" "bc2010" "bc2011" "bc2012" "bc2013" "bc2014" "bc2015" "bc2016" "bc2017" "bc2018" "bc2019" "bc2020" "bc2021" "bc2022" "bc2023" "bc2024" "bc2025" "bc2026" "bc2027" "bc2028" "bc2029" "bc2030" "bc2031" "bc2032" "bc2033" "bc2034" "bc2035" "bc2036" "bc2037" "bc2038" "bc2039" "bc2040" "bc2041" "bc2042" "bc2043" "bc2044" "bc2045" "bc2046" "bc2047" "bc2048" "bc2049" "bc2050" "bc2051" "bc2052" "bc2053" "bc2054" "bc2055" "bc2056" "bc2057" "bc2058" "bc2059" "bc2060" "bc2061" "bc2062" "bc2063" "bc2064" "bc2065" "bc2066" "bc2067" "bc2068")
NID=${NIDS[$SGE_TASK_ID]}

# mkdir -p /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/scratch/processed/polish/pilon/$NID
# cp /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/data/processed/PacBio/assemblies/flye/all_assemblies/fasta/$NID\_assembly.fasta /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/scratch/processed/polish/pilon/$NID/contigs_0.fasta

for CYCLE in 0 1 2 3 4
do
ASSEMBLY=/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/scratch/processed/polish/pilon/$NID/contigs_$CYCLE.fasta
bwa index $ASSEMBLY
bwa mem -t 16 $ASSEMBLY /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/data/processed/Illumina/qc_reads/$IID\_qc.1.fq.gz /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/data/processed/Illumina/qc_reads/$IID\_qc.2.fq.gz | samtools view -bS - > /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/scratch/processed/polish/pilon/$NID/alignment_$CYCLE\.bam
samtools sort /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/scratch/processed/polish/pilon/$NID/alignment_$CYCLE\.bam > /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/scratch/processed/polish/pilon/$NID/alignment_$CYCLE\.sorted.bam
samtools index /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/scratch/processed/polish/pilon/$NID/alignment_$CYCLE\.sorted.bam

pilon -Xmx16G --genome $ASSEMBLY --frags /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/scratch/processed/polish/pilon/$NID/alignment_$CYCLE\.sorted.bam --output contigs_$(expr $CYCLE + 1) --outdir /nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/BRCCH_ECOLI_WGS_EAN/scratch/processed/polish/pilon/$NID/ --changes
done
