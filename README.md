# ecoliCPS

## Introduction

`ecoliCPS` is a repository created to provide all scripts used to generate the genomes presented and described in the publication "Complete Genomes of *E. coli* with diverse K-antigens" by Mores *et al.* MRA (2023), including including setup, walkthroughs and example commands.


## 1. Installation Instructions and Dependencies

All tools can be installed via conda after enabling [bioconda channels](https://bioconda.github.io/)

```
conda create -n ecoliCPS bbmap=38.87 python=3.10 guppy=6.4.6 flye=2.9 pilon=1.23 bwa=0.7.17 samtools=1.9
```

## 2. Processing the input data

### Illumina reads: 
Illumina reads were were quality controlled using BBDuk:

* **Adapter Trimming**: Removal of adapter sequences with specific parameters, including a kmer length of 23 and a minimum kmer length of 11.
* **PhiX Contamination Removal**: Removal of PhiX control sequences using a kmer length of 31.
* **Quality Trimming and Filtering**: Trimming of low-quality bases and filtering out reads below a 45-base length threshold, with a minimum average quality of 20. Remaining sequences are not allowed to contain ambiguous bases ('N's).

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
    base_name=$(basename $r1 .fastq.gz)

    # Adapter trimming
    bbduk.sh in1=$r1 in2=$r2 out1=${processed_data_dir}/${base_name}_trimmed.1.fq.gz out2=${processed_data_dir}/${base_name}_trimmed.2.fq.gz ref=$adapter_ref ktrim=r k=23 mink=11 hdist=1 stats=${processed_data_dir}/${base_name}_adapter.stats

    # PhiX removal
    bbduk.sh in1=${processed_data_dir}/${base_name}_trimmed.1.fq.gz in2=${processed_data_dir}/${base_name}_trimmed.2.fq.gz out1=${processed_data_dir}/${base_name}_phix_removed.1.fq.gz out2=${processed_data_dir}/${base_name}_phix_removed.2.fq.gz ref=$phix_ref k=31 hdist=1 stats=${processed_data_dir}/${base_name}_phix.stats

    # Quality trimming and filtering
    bbduk.sh in1=${processed_data_dir}/${base_name}_phix_removed.1.fq.gz in2=${processed_data_dir}/${base_name}_phix_removed.2.fq.gz out1=${processed_data_dir}/${base_name}_final.1.fq.gz out2=${processed_data_dir}/${base_name}_final.2.fq.gz trimq=14 qtrim=r minlength=45 maq=20 maxns=0 stats=${processed_data_dir}/${base_name}_final.stats

done
```

* **NOTES**:
  * Ensure that you replace the paths provided (e.g. `/path/to/raw/Illumina/data`) with your actual file paths before running the script.
  * The `basename` command is used to extract the base name of the file (i.e., file name without its path and extension) and it is utilized to generate consistent and recognizable names for the output files. This ensures that the output files have predictable names, facilitating their use in subsequent steps of the workflow.

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
            --num_callers <number_of_callers>
```
* **NOTES**:
  * Please replace `path/to/guppy_basecaller`, `path/to/raw_fast5_files`, and `path/to/basecalled_output` with your actual file paths.
  * Replace `<number_of_callers>` with the number of caller threads suitable for your system.
  * The `-x cuda:all` option is intended for systems equipped with NVIDIA CUDA-compatible GPUs and is used to enable GPU acceleration, which can significantly speed up the basecalling process. If the your system does not have NVIDIA GPUs, or if GPU acceleration is not desired, this option can be omitted. The basecaller will then use the CPU, which will be slower compared to GPU processing. 


### PacBio reads: 
PacBio HiFi reads quality control was performed during the run using the Control SW Version `11.0.0.144466`, rendering Q20+ reads. No further QC was performed.


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

* **NOTES**:
  * In both scripts, ensure that you replace the paths provided with your actual file paths.
  * Please replace `<number_of_threads>` with the appropriate value for your system.
  * The `-o` flag specifies the output directory for the assembled genome.

## 4. Polishing the assemblies
Five cycles of iterative polishing was performed on the assembled genomes. In each cycle, it:

* Indexes the assembly with `BWA`.
* Aligns the quality controlled Illumina reads to the assembly using `BWA MEM`, converting the output to `BAM` format with Samtools.
* Sorts and indexes the `BAM` file using `SAMtools`.
* Runs `Pilon` for error correction, specifying input files and output directories. Pilon uses the sorted `BAM` file and the current assembly to produce a polished version.

Polishing `Flye` generated assemblies using `Pilon`:

```
for CYCLE in {0..4}; do

    #Define the path for assembly file and output directory
    ASSEMBLY_PATH="path/to/your/assembly_directory/contigs_$CYCLE.fasta"
    OUTPUT_DIR="path/to/output_directory"
    ILLUMINA_READS_DIR="path/to/your/illumina_qc_reads_directory"

    # Index the assembly with BWA
    bwa index $ASSEMBLY_PATH

    # Align QC'd reads to the assembly and convert to BAM
    for qc_read in ${ILLUMINA_READS_DIR}/*_final.1.fq.gz; do
        r1=$qc_read
        r2=${qc_read%.1.fq.gz}.2.fq.gz
        base_name=$(basename $qc_read _final.1.fq.gz)

        bwa mem -t <number_of_threads> $ASSEMBLY_PATH $r1 $r2 | samtools view -bS - > $OUTPUT_DIR/${base_name}_alignment_$CYCLE.bam

        # Sort and index the BAM file
        samtools sort $OUTPUT_DIR/${base_name}_alignment_$CYCLE.bam > $OUTPUT_DIR/${base_name}_alignment_$CYCLE.sorted.bam
        samtools index $OUTPUT_DIR/${base_name}_alignment_$CYCLE.sorted.bam

        # Run Pilon for polishing
        pilon --genome $ASSEMBLY_PATH --frags $OUTPUT_DIR/${base_name}_alignment_$CYCLE.sorted.bam --output contigs_$(($CYCLE + 1))_${base_name} --outdir $OUTPUT_DIR --changes

    done
done
```

* **NOTES**:
  * Replace the file paths according to your computing environment and dataset.
  * The `<number_of_threads>` placeholder should be replaced with the number of threads for parallel processing.
  * The `basename` command is used to extract the base name of the file (i.e., file name without its path and extension) and it is utilized to generate consistent and recognizable names for the output files. 


## References

* Bushnell, B. (n.d.). BBMap: A Fast, Accurate, Splice-Aware Aligner. Lawrence Berkeley National Lab. (LBNL), Berkeley, CA (United States). Retrieved from https://sourceforge.net/projects/bbmap/.
* Oxford Nanopore Technologies. (n.d.). Guppy. Retrieved from https://nanoporetech.com/nanopore-sequencing-data-analysis.
* Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37, 540–546. https://doi:10.1038/s41587-019-0072-8.
* Walker, B. J., Abeel, T., Shea, T., Priest, M., Abouelliel, A., Sakthikumar, S., ... & Earl, A. M. (2014). Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. PLoS ONE, 9(11), e112963. https://doi:10.1371/journal.pone.0112963.
* Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics, 25(14), 1754-1760. https://doi:10.1093/bioinformatics/btp324.
* Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079. https://doi.org/10.1093/bioinformatics/btp352.


