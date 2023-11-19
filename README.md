# ecoliCPS

## Introduction

`ecoliCPS` is a repository created to provide all scripts used to generate the genomes presented and described in the publication "Complete Genomes of E. coli with diverse K-antigens" by Mores et al. MRA (2023), including including setup, walkthroughs, example commands, and test data.


## 1. Installation Instructions and Dependencies
This project requires the installation of the following bioinformatics tools: Guppy v6.4.6, bbmap v38.87, Flye v2.9 and Pilon v1.23.

### Guppy v6.4.6
To install Guppy, follow the instructions provided on the [Oxford Nanopore Technologies website](https://nanoporetech.com). 

**Dependencies**: Compatible operating system, specific versions of Python may be needed. Requires hardware capable of supporting high-throughput sequencing data processing.

### bbmap v38.87
bbmap can also be installed via conda:
```
conda install -c bioconda bbmap
```
**Dependencies**: Java 7 or later. Sensitive to Java version and memory allocation.

### Flye v2.9
Flye can be installed using conda:
```bash
conda install -c bioconda flye
```
**Dependencies**: Python (version 2.7 or 3.5 and later) and a standard C++ compiler. Performance varies based on Python version and C++ compiler.

### Pilon v1.23
Pilon can be installed using conda as follows:
```
conda install -c bioconda pilon
```
**Dependencies**: Java 8 or later. Requires significant memory for processing large genomes.


* **Note: Ensure all necessary dependencies are installed before using these tools for optimal operation and accurate results.**



## 2. Processing the input data

### Illumina reads: 
In this project, Illumina reads were subjected to quality control and trimming using BBMap's BBDuk tool, streamlined into three key steps:

* **Adapter Trimming**: Removal of adapter sequences with specific parameters, including a kmer length of 23 and a minimum kmer length of 11.
* **PhiX Contamination Removal**: Elimination of PhiX control DNA using a kmer length of 31.
* **Quality Trimming and Filtering**: Trimming of low-quality bases and filtering out reads below a 45-base length threshold, with a minimum average quality of 20 and no ambiguous bases ('N's).

Define directories for raw and processed data: 
```
raw_data_dir="/path/to/raw/Illumina/data"
processed_data_dir="/path/to/processed/Illumina/data"
```
Define reference files for adapters and PhiX (comes with BBmap):
```
adapter_ref="/path/to/adapter/reference.fa"
phix_ref="/path/to/phix/reference.fa"
```
Loop through each pair of Illumina read files:
```
for r1 in ${raw_data_dir}/*R1*fastq.gz; do
    r2=${r1/R1/R2}

    #Adapter trimming
    bbduk.sh in1=$r1 in2=$r2 out1=${r1%.fastq.gz}_trimmed.1.fq.gz out2=${r2%.fastq.gz}_trimmed.2.fq.gz ref=$adapter_ref ktrim=r k=23 mink=11 hdist=1 stats=${r1%.fastq.gz}_adapter.stats

    #PhiX removal
    bbduk.sh in1=${r1%.fastq.gz}_trimmed.1.fq.gz in2=${r2%.fastq.gz}_trimmed.2.fq.gz out1=${r1%.fastq.gz}_phix_removed.1.fq.gz out2=${r2%.fastq.gz}_phix_removed.2.fq.gz ref=$phix_ref k=31 hdist=1 stats=${r1%.fastq.gz}_phix.stats

    #Quality trimming and filtering
    bbduk.sh in1=${r1%.fastq.gz}_phix_removed.1.fq.gz in2=${r2%.fastq.gz}_phix_removed.2.fq.gz out1=${r1%.fastq.gz}_final.1.fq.gz out2=${r2%.fastq.gz}_final.2.fq.gz trimq=14 qtrim=r minlength=45 maq=20 maxns=0 stats=${r1%.fastq.gz}_final.stats

done
```

### ONT reads: 
Guppy


### PacBio reads: 
In this project, the PacBio HiFi reads quality control was performed during the run using the Control SW Version 11.0.0.144466, rendering Q20+ reads. No further QC was performed.


## 3. Assembling the genomes

## 4. Polishing the assemblies
In this project, iterative polishing of genomic assemblies using Pilon with Illumina QC'd reads was perfomed for five cycles. In each cycle, it:
* Indexes the assembly with BWA.
* Aligns the quality-controlled Illumina reads to the assembly using BWA MEM, converting the output to BAM format with Samtools.
* Sorts and indexes the BAM file using Samtools.
* Runs Pilon for error correction, specifying input files and output directories. Pilon uses the sorted BAM file and the current assembly to produce a polished version.

```
for CYCLE in {0..4}; do

#Define the path for assembly file and output directory:
ASSEMBLY_PATH="<path_to_assembly_directory>/contigs_$CYCLE.fasta"
OUTPUT_DIR="<path_to_output_directory>"
ILLUMINA_READS_DIR="<path_to_illumina_qc_reads_directory>"

#Index the assembly with BWA:
bwa index $ASSEMBLY_PATH

#Align QC'd reads to the assembly and convert to BAM:
bwa mem -t <number_of_threads> $ASSEMBLY_PATH $ILLUMINA_READS_DIR/read1.fq.gz $ILLUMINA_READS_DIR/read2.fq.gz | samtools view -bS - > $OUTPUT_DIR/alignment_$CYCLE.bam

#Sort and index the BAM file using Samtools:
samtools sort $OUTPUT_DIR/alignment_$CYCLE.bam > $OUTPUT_DIR/alignment_$CYCLE.sorted.bam
samtools index $OUTPUT_DIR/alignment_$CYCLE.sorted.bam

#Polish the assembly using Pilon:
pilon --genome $ASSEMBLY_PATH --frags $OUTPUT_DIR/alignment_$CYCLE.sorted.bam --output contigs_$(($CYCLE + 1)) --outdir $OUTPUT_DIR --changes

done
```


## References

1.

