#Used for paper: Resolving the complex genetic basis of phenotypic variation and variability of cellular growth, Genetics, 2017
#Naomi Ziv, nz375@nyu.edu

#QTL mapping using HPC

#set directory
setwd('/scratch/nz375/')
#setwd('/scratch/nz375/QTLpermScan2')

#libraries
require(snow)
#install.packages('/home/nz375/qtl_1.35-3.tar','/home/nz375/Rpackages/',NULL)
require(qtl,'/home/nz375/Rpackages/')

#read in data as backcross and change to doubled haploids
cross.L<-read.cross(format='csv',file='L.geno.csv',genotypes=c('A','B'))
class(cross.L)[1]<-'dh'
cross.H<-read.cross(format='csv',file='H.geno.csv',genotypes=c('A','B'))
class(cross.H)[1]<-'dh'

#calculate QTL genotype probabilities
cross.L<-calc.genoprob(cross.L,step=1,error.prob=0.01)
cross.H<-calc.genoprob(cross.H,step=1,error.prob=0.01)
#single QTL mapping
out.hk.L<-scanone(cross.L,pheno.col=1:12,method='hk')
out.hk.H<-scanone(cross.H,pheno.col=1:12,method='hk')
#threshold by permutation
operm.hk.L<-scanone(cross.L,n.perm=10000,pheno.col=c('GR.PC','GR.GRsd.RES','lag.PC','lag.lagmad.RES'),verbose=T,n.cluster=12,method='hk')
operm.hk.H<-scanone(cross.H,n.perm=10000,pheno.col=c('GR.PC','GR.GRsd.RES','lag.PC','lag.lagmad.RES'),verbose=T,n.cluster=12,method='hk')
#save workspace
#save.image(file='workspace1.Rdata')

#corser grid
#calculate QTL genotype probabilities
cross.L<-calc.genoprob(cross.L,step=5,error.prob=0.01)
cross.H<-calc.genoprob(cross.H,step=5,error.prob=0.01)

#double QTL mapping
out2.hk.L<-scantwo(cross.L,pheno.col=c('GR.PC','GR.GRsd.RES','lag.PC','lag.lagmad.RES'),method='hk')

#save.image(file='workspace2.Rdata')

out2.hk.H<-scantwo(cross.H,pheno.col=c('GR.PC','GR.GRsd.RES'),method='hk')

#save workspace
#save.image(file='workspace2.Rdata')

#threshold by permutation
operm2.hk.L<-scantwo(cross.L,n.perm=1000,pheno.col=c('GR.PC','GR.GRsd.RES','lag.PC','lag.lagmad.RES'),verbose=T,n.cluster=20,method='hk')

#save.image(file='workspace3.Rdata')

operm2.hk.H<-scantwo(cross.H,n.perm=1000,pheno.col=c('GR.PC','GR.GRsd.RES'),verbose=T,n.cluster=20,method='hk')

#save workspace
save.image(file='QTLmapping.Rdata')
#Later renamed to QTLmapping.Rdata


#####scantwo permutations
#This was done once to assess the genome wide distrbution of LOD scores on a permutated dataset. 
#first read a permuted dataset and perform scantwo

#read in data as backcross and change to doubled haploids
#cross.L.P<-read.cross(format='csv',file='L.geno.Perm1.csv',genotypes=c('A','B'))
#class(cross.L.P)[1]<-'dh'

#calculate QTL genotype probabilities
#cross.L.P<-calc.genoprob(cross.L.P,step=1,error.prob=0.01)
#cross.H<-calc.genoprob(cross.H,step=1,error.prob=0.01)

#double QTL mapping
#out2.hk.L.P<-scantwo(cross.L.P,pheno.col=c('GR.PC','GR.GRsd.RES','lag.PC','lag.lagmad.RES'),method='hk')

#save.image(file='scantwo.P1.Rdata')





