#!/bin/bash
#SBATCH --job-name=Salmon_quant
#SBATCH --cpus-per-task=8
#SBATCH --mem=24000
#SBATCH --output=outfile_salmon.%j
#SBATCH --error=errfile_salmon.%j



cd /mnt/beegfs/workshop/A_demo_files/files_for_differentialExpression/fastq_files

for file in *_1.fastq.gz; do
 

  filename=$(basename "$file" _1.fastq.gz)  ## create a prefix called filename
  
 salmon quant -i /mnt/beegfs/workshop/A_reference_files/salmon_index \
 -l A \
 -1 ${filename}_1.fastq.gz \
 -2 ${filename}_2.fastq.gz \
 -o /mnt/beegfs/workshop/<SSO>/results/salmon_quants/${filename}_quant \
 --seqBias \
 --gcBias \
 -p 8

done 
 
