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
