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
