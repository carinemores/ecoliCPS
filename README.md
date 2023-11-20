# ecoliCPS

## Introduction

`ecoliCPS` is a repository created to provide all scripts used to generate the genomes presented and described in the publication "Complete Genomes of *E. coli* with diverse K-antigens" by Mores *et al.* MRA (2023), including including setup, walkthroughs and example commands.


## 1. Installation Instructions and Dependencies
The followig bioinformatic tools were required to generated the presented genomes: Guppy v6.4.6, bbmap v38.87, Flye v2.9, Pilon v1.23, BWA v0.7.17 and SAMtools v1.9.
Below you can find instructions on how to download each tool and its dependencies.

### bbmap v38.87
bbmap can be installed via conda:
```
conda install -c bioconda bbmap
```
**Dependencies**: Java 7 or later. Sensitive to Java version and memory allocation.

### Guppy v6.4.6
To install Guppy, follow the instructions provided on the [Oxford Nanopore Technologies website](https://nanoporetech.com). 

**Dependencies**: Compatible operating system, specific versions of Python may be needed. Requires hardware capable of supporting high-throughput sequencing data processing.

### Flye v2.9
Flye can be installed using conda:
```bash
conda install -c bioconda flye
```
**Dependencies**: Python (version 2.7 or 3.5 and later) and a standard C++ compiler. Performance varies based on Python version and C++ compiler.

### Pilon v1.23
Pilon can also be installed using conda:
```
conda install -c bioconda pilon
```
**Dependencies**: Java 8 or later. Requires significant memory for processing large genomes.

### Combined Installation Guide for BWA v0.7.17 and SAMtools v1.9 (required for the polishing step):

**Dependencies:**
* GCC compiler
* zlib (development libraries)
* ncurses (development libraries, only for SAMtools)
* bzip2 (development libraries, only for SAMtools)
* xz (development libraries, only for SAMtools)
  
**Installation Steps:**

Update Package List:
```
sudo apt-get update
```
Install Common Dependencies:
```
sudo apt-get install gcc zlib1g-dev
sudo apt-get install libncurses5-dev libbz2-dev liblzma-dev\
```
**Download and Install BWA:**
```
wget https://github.com/lh3/bwa/archive/refs/tags/v0.7.17.tar.gz
```
Extract the Tarball:
```
tar -xzf v0.7.17.tar.gz
```
Navigate to the BWA directory and compile:
```
cd bwa-0.7.17
make
```
**Download and Install SAMtools:**

Download SAMtools from its official website or GitHub repository:
```
wget https://github.com/samtools/samtools/archive/refs/tags/1.9.tar.gz
```
Extract the Tarball:
```
tar -xzf 1.9.tar.gz
```
Navigate to the SAMtools directory and compile:
```
cd samtools-1.9
./configure --prefix=/where/to/install
make
make install
```
* Replace /where/to/install with your desired installation path.


* **NOTE**: Ensure all necessary dependencies are installed before using these tools for optimal operation and accurate results.



## 2. Processing the input data

### Illumina reads: 
Illumina reads were subjected to quality control and trimming using BBMap's BBDuk tool, streamlined into three key steps:

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
* **Note**: Ensure that you replace the paths provided with your actual file paths before running the script. 

### ONT reads: 
This script is for running the Guppy basecaller, a tool used for basecalling Oxford Nanopore Technologies (ONT) sequencing data. The script specifies the input path for raw sequencing data, the output path for basecalled data, flowcell and kit types, and other relevant parameters for basecalling. It also includes options for trimming adapters and primers, barcode handling, and output compression.
```
#Set up Guppy basecaller with generic paths and parameters
GUPPY_PATH="path/to/guppy_basecaller"
INPUT_PATH="path/to/raw_fast5_files"
OUTPUT_PATH="path/to/basecalled_output"
FLOWCELL_TYPE="FLO-MIN114"  
KIT_TYPE="SQK-RBK114-24"       
BARCODE_KITS="SQK-RBK114-24"   

#Run Guppy basecaller
$GUPPY_PATH --input_path $INPUT_PATH \
            --save_path $OUTPUT_PATH \
            --flowcell $FLOWCELL_TYPE \
            --kit $KIT_TYPE \
            --barcode_kits $BARCODE_KITS \
            -x cuda:all \
            --trim_adapters \
            --trim_primers \
            --enable_trim_barcodes \
            --compress_fastq \
            --num_callers 8
```
* **Note**: Please replace `path/to/guppy_basecaller`, `path/to/raw_fast5_files`, and `path/to/basecalled_output` with your actual file paths. The `-x cuda:all` option is intended for systems equipped with NVIDIA CUDA-compatible GPUs and is used to enable GPU acceleration, which can significantly speed up the basecalling process. If the your system does not have NVIDIA GPUs, or if GPU acceleration is not desired, this option can be omitted. The basecaller will then use the CPU, which will be slower compared to GPU processing.


### PacBio reads: 
PacBio HiFi reads quality control was performed during the run using the Control SW Version 11.0.0.144466, rendering Q20+ reads. No further QC was performed.


## 3. Assembling the genomes
Flye is a *de novo* assembler for long and noisy reads, such as those produced by PacBio and Oxford Nanopore Technologies. 
Flye v2.9 offers specific modes for handling different types of long-read data:

* PacBio HiFi Reads: The `--pacbio-hifi` option is used for high-quality, low-error-rate PacBio reads. This mode is optimized for accuracy, taking advantage of the high fidelity of these reads.

* ONT Reads: The `--nano-raw` option is suitable for raw ONT reads. These reads are typically longer but have a higher error rate compared to PacBio HiFi reads. Flye employs algorithms to handle the unique error profiles of ONT data.

### Assembling PacBio reads:

Define variables for the input file and output directory:
```
PACBIO_INPUT="/path/to/pacbio_hifi_reads.fastq"
OUTPUT_DIR="/path/to/output/directory"
```
Run Flye for PacBio HiFi reads:
```
flye --pacbio-hifi $PACBIO_INPUT -o $OUTPUT_DIR -t <number_of_threads>
```

### Assembling ONT reads:

Define variables for the input file and output directory:
```
ONT_INPUT="/path/to/ont_reads.fastq"
OUTPUT_DIR="/path/to/output/directory"
```
Run Flye for ONT reads:
```
flye --nano-raw $ONT_INPUT -o $OUTPUT_DIR -t <number_of_threads>
```

* **Note**: In both scripts, ensure that you replace the paths provided with your actual file paths before running the script. Also, update `<number_of_threads>` with the appropriate value for your system. The `-o` flag specifies the output directory for the assembled genome.

## 4. Polishing the assemblies
To polish the presented genomes, iterative polishing was perfomed for five cycles. In each cycle, it:

* Indexes the assembly with BWA.
* Aligns the quality controlled Illumina reads to the assembly using BWA MEM, converting the output to BAM format with Samtools.
* Sorts and indexes the BAM file using Samtools.
* Runs Pilon for error correction, specifying input files and output directories. Pilon uses the sorted BAM file and the current assembly to produce a polished version.

Polishing Flye generated assemblies using Pilon:
```
for CYCLE in {0..4}; do
    
    #Define the path for assembly file and output directory
    ASSEMBLY_PATH="path/to/your/assembly_directory/contigs_$CYCLE.fasta"
    OUTPUT_DIR="path/to/output_directory"
    ILLUMINA_READS_DIR="path/to/your/illumina_qc_reads_directory"

    #Index the assembly with BWA
    bwa index $ASSEMBLY_PATH

    #Align QC'd reads to the assembly and convert to BAM
    bwa mem -t <number_of_threads> $ASSEMBLY_PATH $ILLUMINA_READS_DIR/read1.fq.gz $ILLUMINA_READS_DIR/read2.fq.gz | samtools view -bS - > $OUTPUT_DIR/alignment_$CYCLE.bam

    #Sort and index the BAM file
    samtools sort $OUTPUT_DIR/alignment_$CYCLE.bam > $OUTPUT_DIR/alignment_$CYCLE.sorted.bam
    samtools index $OUTPUT_DIR/alignment_$CYCLE.sorted.bam

    #Run Pilon for polishing
    pilon --genome $ASSEMBLY_PATH --frags $OUTPUT_DIR/alignment_$CYCLE.sorted.bam --output contigs_$(($CYCLE + 1)) --outdir $OUTPUT_DIR --changes

done

```
* **Note**: Replace the file paths according to your computing environment and dataset. The `<number_of_threads>` placeholder should be replaced with the desired number of threads for parallel processing.


## References

* Bushnell, B. (n.d.). BBMap: A Fast, Accurate, Splice-Aware Aligner. Lawrence Berkeley National Lab. (LBNL), Berkeley, CA (United States). Retrieved from https://sourceforge.net/projects/bbmap/.
* Oxford Nanopore Technologies. (n.d.). Guppy. Retrieved from https://nanoporetech.com/nanopore-sequencing-data-analysis.
* Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37, 540–546. https://doi:10.1038/s41587-019-0072-8.
* Walker, B. J., Abeel, T., Shea, T., Priest, M., Abouelliel, A., Sakthikumar, S., ... & Earl, A. M. (2014). Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. PLoS ONE, 9(11), e112963. https://doi:10.1371/journal.pone.0112963.
* Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics, 25(14), 1754-1760. https://doi:10.1093/bioinformatics/btp324.
* Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079. https://doi:10.1093/bioinformatics/btp352.


