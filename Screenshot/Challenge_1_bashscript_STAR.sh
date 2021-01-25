#!/usr/bin/bash
#if u use hpc or cluster server which softwares are saved as module.
module purge
module load subread-1.6.2-gcc-5.4.0-7zywp5u
module load samtools-1.9-gcc-5.4.0-vf6vvem
module load star/2.5.0a

datadir=/home/training/Practical_day2/SS2_Data/Fastq
refdir=/home/training/Practical_day2/refdata-gex-GRCh38-2020-A/star
gtfdir=/home/training/Practical_day2/refdata-gex-GRCh38-2020-A/genes

cd $datadir


for fastqs in $(ls *.fastq.gz | cut -f 1 -d "."| uniq);
do
  STAR --runThreadN 8 --genomeDir ${refdir} \
       --readFilesCommand zcat \
       --readFilesIn ${fastqs}.1.fastq.gz  ${fastqs}.2.fastq.gz \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix ${fastqs}_ ;

  samtools index ${fastqs}_Aligned.sortedByCoord.out.bam

done

## more options you could add for filtering or quanlity control
# --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
# --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
# --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \

featureCounts -a  ${gtfdir}/genes.gtf \
              -t exon -g gene_id \
              -o STAR_QuantSeq_featureCounts.txt \
                 24087_7_73_Aligned.sortedByCoord.out.bam \
                 24087_5_146_Aligned.sortedByCoord.out.bam \
                 23728_8_119_Aligned.sortedByCoord.out.bam

cut -f1,7-9 STAR_QuantSeq_featureCounts.txt > SS2_STAR_EM.txt
