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
