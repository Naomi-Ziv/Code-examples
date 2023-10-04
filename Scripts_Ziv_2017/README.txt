Resolving the complex genetic basis of phenotypic variation and variability of cellular growth 
Naomi Ziv, Bentley M. Shuster, Mark L. Siegal, David Gresham
2017 
#############################################################
For any questions, explanations, raw image data and information on production of figures, please contact Naomi Ziv, nz375@nyu.edu   
#############################################################

Explanation of files:
1) QTL mapping 
2) Other analysis files
3) Datasets in .txt format
#######################################################################################################################################################################################

1) QTL mapping 
#############################################################

	A) HPC.QTL.R - R script used for QTL mapping (one and two dimensional QTL scans), run on high performance computing cluster. Two dimensional scan permutations are computationally intensive. Requires R packages ‘snow’ and ‘QTL’. Uses as input: L.geno.csv and H.geno.csv. Will create output: QTLmapping.Rdata    

	B) Ziv2017QTLAnalysis.R - Main R script used for analysis of QTL (calling loci and fitting multi-QTL models). Requires R packages ‘QTL’ and ‘plyr’. Can be edited slightly to assess different thresholds. Uses as input: QTLmapping.Rdata and RenameQTLMarkers.txt. Will create output hkQTLworkspace0.05.Rdata 

	C) L.geno.csv - Data for QTL mapping. Comma-delimited text file containing phenotype and genotype data for QTL mapping. Compatible with RQTL. Growth limiting glucose concentration (0.22mM glucose). Phenotypes are: GR - mean growth rate, GRsd - growth rate standard deviation, GRcount - number of colonies, lag - median lag duration, lagmed - lag MAD, lag count - number of colonies, GR.PC - plate corrected mean growth rate, GRsd.PC - plate corrected growth rate standard deviation (NOT USED), lag.PC - plate corrected median lag duration, lagmad.PC - plate corrected lag duration MAD (NOT USED), GR.GRsd.RES - loess residuals growth rate standard deviation on mean, lag.lagmad.RES - loess residuals lag duration MAD on median.   

	D) H.geno.csv - Data for QTL mapping. Comma-delimited text file containing phenotype and genotype data for QTL mapping. Compatible with RQTL. Non-growth limiting glucose concentration (4.44mM glucose). Phenotypes are: GR - mean growth rate, GRsd - growth rate standard deviation, GRcount - number of colonies, lag - median lag duration, lagmed - lag MAD, lag count - number of colonies, GR.PC - plate corrected mean growth rate, GRsd.PC - plate corrected growth rate standard deviation (NOT USED), lag.PC - plate corrected median lag duration, lagmad.PC - plate corrected lag duration MAD (NOT USED), GR.GRsd.RES - loess residuals growth rate standard deviation on mean, lag.lagmad.RES - loess residuals lag duration MAD on median. 

	E) RenameQTLMarkers.txt - Tab delimited text file describing the locations of markers used for genotyping, from (Gerke et al. 2009). 

	F) QTLmapping.Rdata - R data file containing R variables created with the script HPC.QTL.R, loaded by Ziv2017QTLAnalysis.R for further QTL analysis. Contains: cross.H - RQTL cross object for high glucose, cross.L - RQTL cross object for low glucose, operm.hk.H - permutations for one dimensional scan for high glucose, operm.hk.L - permutations for one dimensional scan for low glucose, operm2.hk.H - permutations for two dimensional scan for high glucose, operm2.hk.L - permutations for two dimensional scan for low glucose, out.hk.H - one dimensional scan for high glucose, out.hk.L - one dimensional scan for low glucose, out2.hk.H - two dimensional scan for high glucose, out2.hk.L - two dimensional scan for low glucose.  

	G) hkQTLworkspace0.05.Rdata - R data file containing R variables created with the script Ziv2017QTLAnalysis.R. Contains all of QTLmapping.Rdata variables plus: Afits.H - additive loci models for high glucose, Afits.L - additive loci models for low glucose, findunique - function for collapsing close loci, fitadd - function for fitting additive loci models, fitint - function for fitting additive and interaction loci models, getpys - function for getting approximate physical position, getQTL - function for identifying QTL, Ifits.H - additive and interaction loci models for high glucose, Ifits.L - additive and interaction loci models for low glucose, pys - approximate physical positions, qtls.H - QTLs identified for high glucose, qtls.L - QTLs identified for low glucose.   

	H) Segdata.Rdata - R data file containing R variables. Contains: Segdata - microcolony growth data, see also SegregantGrowth.txt in section 3, agg - growth summary statistics.       

2) Other analysis files
#############################################################

	A) Ziv2017SEQAnalysis.R - R script used for analysis. Requires R packages ‘QTL’ and ‘plyr’. Uses as input: SEQdata.Rdata, Multipool.Rdata. Also uses QTLmapping.Rdata and RenameQTLMarkers.txt which are found in section 1. 

	B) doMultipool.sh - shell script for running MULTIPOOL on many samples with different parameter settings. 

	C) NZ_align_example.q - example of a file used for aligning fastq read files, script represents an example for a single sample, originally scripts iterated on all samples.

	D) NZ_snps_example.q - example of a file used for calling SNPs, script represents an example for a single sample, originally scripts iterated on all samples.

	E) snp.py - python script used for extracting and reorganizing information about SNPs. Creates .snp files. 

	F) Oak.Vineyard.multipool.datagen.R - R script used for organizing sequencing data, in general and for multipool analysis. Uses as input: .snp files, creates SEQdata.Rdata, txt files per chromosome and sample with snp positions and read counts - used for multipool analysis and Multipool.Rdata. 

	G) SEQdata.Rdata - R data file containing R variables. Contains: SEQdata - results from all sequencing analysis, see also SNPsummary.txt in section 3.

	H) Multipool.Rdata - R data file containing R variables. Contains: Multipool - results from all MULTIPOOL analysis, see also Multipool.txt in section 3.

	I) HXTdata.Rdata - R data file containing R variables. Contains: AllRepdata - Microcolony growth data for allele replacement strains, RecHemdata - Microcolony growth data for Reciprocol hemizygote strains, see also AlleleReplaceGrowth.txt and ReciprocalHemizygoteGrowth.txt in section 3.


3) Datasets in .txt format
#############################################################

	A) SegregantGrowth.txt - Sergeant growth data, Tab delimited text file containing data for 1,304,734 micro-colonies. Columns are: GR - growth rate, lag - lag duration, colony - unique colony identifier (within plate), Isize - initial size, Genotype - strain, Environment - media, well, plate, field. Strains are: BC241 - Vineyard, BC248 - Oak, BC251 - F1, BC252 - F1 and 477 F2s with unique identifiers (i.e: F2-12-A1).

	B) AlleleReplaceGrowth.txt - Allele replacement strains growth data, Tab delimited text file containing data for 298,066 micro-colonies. Columns are: GR - growth rate, lag - lag duration, colony - unique colony identifier (within plate), Isize - initial size, Genotype - strain, Environment - media, well, plate, field. Strains are: DGY279 - Oak, DGY1059 - Oak with Vineyard HXT7 whole locus, DGY1150 - Oak with Vineyard HXT7 single aa, DGY278 - Vineyard, DGY1055 - Vineyard with Oak HXT7 whole locus, DGY1077 - Vineyard with Oak HXT7 single aa, DGY419 - F1, DGY1197 - F1 with Oak HXT7 whole locus, DGY1195 - F1 with Vineyard HXT7 whole locus.

	C) ReciprocalHemizygoteGrowth.txt - Reciprocol hemizygote strains growth data, Tab delimited text file containing data for 92,922 micro-colonies. Columns are: GR - growth rate, lag - lag duration, colony - unique colony identifier (within plate), Isize - initial size, Genotype - strain, Environment - media, well, plate, field. Strains are: F1-O6-KO - F1 with oak HXT6 knocked out, F1-W6-KO - F1 with vineyard HXT6 knocked out, F1-O7-KO - F1 with oak HXT7 knocked out, F1-W7-KO - F1 with vineyard HXT7 knocked out.

	D) SNPsummary.txt - Summary of all sequence data, Tab delimited text file summarizing SNPs identified by genomic sequencing. Contains 1,895,619 rows and columns: car - chromosome, pos - position, ref - reference allele, alt - alternate allele, qual - quality, RefF - number of reference forward reads, RefR - number of reference reverse reads, AltF - number of alternate forward reads, AltR - number of alternate reverse reads, Rdepth - read depth, Sample - sample identifier, set - identifier that groups samples, dilution - chemostat dilution rate (if applicable_, Id - SNP identifier (0 - not found in parents, 1 - oak, 2 - vineyard, 3 - found in both parents, 1.1/2.2 - positions where each parent has different allele from reference, 1.2/2.2 - alleles with frequency < 0.75 in parents).   

	E) Multipool.txt - Summary of multipool results, Tab delimited file summarizing MULTIPOOL reseals. Contains 5,742,416 rows and columns: bin - bin start position, AF - allele frequency, LOD - LOD score, tag - identifier, car - chromosome, RF - MULTIPOOL recombination rate parameter, N - MULTIPOOL number of individuals parameter.       