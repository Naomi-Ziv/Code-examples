Analysis for:
Broad susceptibility of Candida auris strains to 8-hydroxyquinolines and mechanisms of resistance
mBio, 2023, DOI: 10.1128/mbio.01376-23

Attaching some scripts I used for this project - I'll just say that they are completely uncommented (I usually make scripts I publish a bit more readable but I'm just sending it as is at this point - if you want me to go through and spruce up or double check functionality - let me know).

MattAuris.R - A R script that organizes the sequencing data. Input: sampleinfo.txt (also attached), the gff files (also attached the versions I used) and the sorted vcf files, Output: VarDataset.tsv   

mpileup.sh - A shell script that goes through calling variants. Input: sam files, Output: sorted vcf files  