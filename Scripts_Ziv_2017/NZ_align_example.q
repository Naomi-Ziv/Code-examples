#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=03:00:00
#PBS -N NZ_align
#PBS -M nz375@nyu.edu

#source /etc/profile.d/env-modules.sh
module purge
module load bwa/gnu/0.7.8
module load samtools/intel/0.1.19

REF=/scratch/nz375/Ref/Ref.SGD020311.fasta

bwa index -a bwtsw $REF

cd /scratch/nz375/

#The following example is for sample C1S2, originally the script included all samples
READ=/scratch/nz375/reads/C1S2_1_CGATGT_L004_R1_001.fastq
bwa aln -t 12 $REF $READ > /scratch/nz375/C1S2.sai

exit 0;
