#Define directories for raw and processed data:
raw_data_dir="/path/to/raw/Illumina/data"
processed_data_dir="/path/to/processed/Illumina/data"

#Define reference files for adapters and PhiX (comes with BBmap):
adapter_ref="/path/to/adapter/reference.fa"
phix_ref="/path/to/phix/reference.fa"

#Loop through each pair of Illumina read files:
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
