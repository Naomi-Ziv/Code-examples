#!/bin/sh

#comments in shell scripts are defined by a hashtag

conda activate SeqA
cd /data/naomi/matts_analysis/files/

samtools faidx C_auris_B11221_chromosomes.fasta
samtools faidx C_auris_B8441_version_s01-m01-r17_chromosomes.fasta

for name in 1_S1 3_S3 4_S4 7_S7 8_S8 9_S9 10_S10 11_S11 12_S12 13_S13 14_S14 32_S25 35_S26 52_S29 53_S30
do
    echo ${name}
    samtools view -@ 4 -bt C_auris_B11221_chromosomes.fasta.fai strain${name}_q.sam -o strain${name}.bam
    samtools sort -@ 4 strain${name}.bam -o strain${name}_sorted.bam
    samtools index strain${name}_sorted.bam

    bcftools mpileup -Ou -d 1000 -f C_auris_B11221_chromosomes.fasta strain${name}_sorted.bam | bcftools call -mv -Oz --ploidy 1 -o strain${name}_sorted.vcf
done

for name in 2_S2 5_S5 6_S6 16_S16 17_S17 18_S18 23_S23 24_S23 25_S24 38_S27 41_S28 54_S31 55_S32
do
    echo ${name}
    samtools view -@ 4 -bt C_auris_B8441_version_s01-m01-r17_chromosomes.fasta.fai strain${name}_q.sam -o strain${name}.bam
    samtools sort -@ 4 strain${name}.bam -o strain${name}_sorted.bam
    samtools index strain${name}_sorted.bam

    bcftools mpileup -Ou -d 1000 -f C_auris_B8441_version_s01-m01-r17_chromosomes.fasta strain${name}_sorted.bam | bcftools call -mv -Oz --ploidy 1 -o strain${name}_sorted.vcf
done

