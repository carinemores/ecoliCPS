#Define variables for the input file and output directory:
ONT_INPUT="/path/to/ont_reads.fastq"
OUTPUT_DIR="/path/to/output/directory"

#Run Flye for ONT reads:
flye --nano-raw $ONT_INPUT -o $OUTPUT_DIR -t <number_of_threads>
