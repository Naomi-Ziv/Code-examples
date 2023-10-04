#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=05:00:00
#PBS -l mem=4GB
#PBS -N NZ_snps
#PBS -M nz375@nyu.edu

module purge
module load bwa/gnu/0.7.8
module load samtools/intel/0.1.19

REF=/scratch/nz375/Ref/Ref.SGD020311.fasta

cd /scratch/nz375/

#The following example is for sample C1S2, originally the 'name' variable iterated on all samples 
for name in C1S2
do
    READ=`ls /scratch/nz375/reads/$name*`
    bwa samse $REF /scratch/nz375/${name}.sai $READ > ${name}.sam
    samtools view -bS ${name}.sam | samtools sort - ${name}.sort
    samtools index ${name}.sort.bam
    samtools mpileup -I -uf $REF ${name}.sort.bam > ${name}.sort.bcf
    bcftools index ${name}.sort.bcf
    bcftools view -bvcg ${name}.sort.bcf | bcftools view - > /scratch/nz375/${name}.sort.snp
done

exit 0;
