#! /bin/bash
#SBATCH --job-name=guppy
#SBATCH --output=guppy_23CR10_out.log
#SBATCH --error=guppy_23CR10_error.log
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gpus-per-node=8
#SBATCH --tmp=1024
#SBATCH --time=24:00:00

/cluster/project/sunagawa/fieldc/ont-guppy/bin/guppy_basecaller --input_path /cluster/work/sunagawa/fieldc/23CR10/fast5 --save_path /cluster/work/sunagawa/fieldc/23CR10/basecall --flowcell FLO-MIN114 --kit SQK-RBK114-24 --barcode_kits SQK-RBK114-24 -x cuda:all --trim_adapters --trim_primers --enable_trim_barcodes --compress_fastq --num_callers 8
