#Used for paper: Resolving the complex genetic basis of phenotypic variation and variability of cellular growth, Genetics, 2017
#Naomi Ziv, nz375@nyu.edu

#Set working directory
setwd('~/Desktop/QTLpaper')

#Additional code:
#####Organizing growth rate data - Oak.Vineyard.QTL.datagen.R
#input is microcolony growth assay output folders, output is L.phenotypes, H.phenotypes and Segdata.Rdata
#####Formating data for RQTL (was run twice, once for each envirionment) - QTLsetup.R
#input is L.phenotypes and H.phenotypes, output is L.geno and H.geno
#####Running RQTL (HPC) - HPC.QTL.R:
#input is L.geno and H.geno, output is QTLmapping.Rdata

#load mapping data
load('./Data/QTLmapping.Rdata')
#cross.L, cross.H, out.hk.L, out.hk.H, out2.hk.L, out2.hk.H, operm.hk.L, operm.hk.H, operm2.hk.L, operm2.hk.H

#libraries
library(qtl)
library(plyr)

##############################################################
#Custom functions            #################################
##############################################################
##############################################################
#Identifing unique loci
findunique<-function(chrlist,collapse=30){
	#chrlist is a dataframe with loci position ($pos) and placeholder for unique identifier ($unique, all should be 0), for one chromosome
	while(sum(chrlist$unique==0)!=0){
		temp<-which(abs((chrlist$pos-chrlist$pos[which(chrlist$unique==0)[1]]))<=collapse)
		if(sum(chrlist$unique[temp]!=0)>0){print('warning')}
		#warning means a loci is within collapse distence of two differnt loci (which are not within collapse distence) - look at positions and determine best collapse
		chrlist$unique[temp]<-max(chrlist$unique)+1
		}
		return(chrlist)}
##############################################################
##############################################################
#Calling QTL
getQTL<-function(lodcolumn,out,operm,out2,operm2,alpha,collapse=30){
	#out,operm,out2 and operm2 are result of scanone and scantwo, lodcolumn is number, alpha is a vector of 2 (for additive and interactions)
	#example:lodcolumn=1;out=out.hk.L;operm=operm.hk.L;out2=out2.hk.L;operm2=operm2.hk.L;alpha=c(0.1,0.1);cross=cross.L;collapse=30
	require(qtl)
	require(plyr)
	#assumes findunique function 
	
	#single QTL
	temp1<-summary(out,perms=operm,alpha=alpha[1],lodcolumn=lodcolumn)
	#addtional QTL on same chromosome from 2d scan
	temp2<-summary(out2,perms=operm2,alpha=c(0,0,0,alpha[1],alpha[1]),lodcolumn=lodcolumn,what='add',allpairs=F)
	#organize - get chr and pos
	data<-data.frame(chr=c(temp1$chr,temp2$chr1,temp2$chr2),pos=c(temp1$pos,temp2$pos1,temp2$pos2))
	
	#organize - collapse loci (loci found from 2d with loci found from 1d, take mean pos)
	if(nrow(data)>0){
		data<-ldply(split(cbind(data,unique=0),data$chr),findunique,collapse=collapse)
		data<-aggregate(pos~chr+unique,data,mean)
		data<-cbind(data[,c('chr','pos')],name=paste('Q',1:nrow(data),sep=''))
		data$name<-as.character(data$name)}
	
	intbound<-nrow(data)
	#interacting QTL
	temp3<-summary(out2,perms=operm2,alpha=c(1,0,alpha[2],0,0),lodcolumn=lodcolumn,what='int')
	int<-matrix(,nrow(temp3),3)
	if(nrow(temp3)>0){
			for(i in 1:nrow(temp3)){
				#is 1st interacting loci already a QTL
				q1<-which(data$chr==temp3$chr1[i] & data$pos<(temp3$pos1[i]+collapse) & data$pos>(temp3$pos1[i]-collapse))
				if(length(q1)==0){
					#1st interacting loci is new
					data<-rbind(data,data.frame(chr=temp3$chr1[i],pos=temp3$pos1[i],name=paste('QI',nrow(data)+1,sep='')))
					q1<-nrow(data)
				}
				#is 2st interacting loci already a QTL
				q2<-which(data$chr==temp3$chr2[i] & data$pos<(temp3$pos2[i]+collapse) & data$pos>(temp3$pos2[i]-collapse))		
				if(length(q2)==0){
					#2st interacting loci is new
					data<-rbind(data,data.frame(chr=temp3$chr2[i],pos=temp3$pos2[i],name=paste('QI',nrow(data)+1,sep='')))
					q2<-nrow(data)
				}
				#if two interactions share loci, take mean pos
				if(q1>intbound){data$pos[q1]<-mean(c(data$pos[q1],temp3$pos1[i]))}
				if(q2>intbound){data$pos[q2]<-mean(c(data$pos[q2],temp3$pos2[i]))}
				#info on which loci interact
				int[i,]<-c(q1,q2,temp3$lod.int[i])
			}
		}
	
	attr(data,'trait')<-colnames(out)[lodcolumn+2]
	return(list(data,int))
	
	}
##############################################################
##############################################################
#Fitting multi-QTL models - additive
fitadd<-function(qtl,cross){
	#qtl is the result of running getQTL, cross is the cross object
	require(qtl)
	#only qtl list
	qtl<-qtl[[1]]
	trait<-attr(qtl,'trait')
	#additive loci
	sub<-grep(paste('Q','[1-9]',sep=''),qtl$name)
	#fit
	if(length(sub)>0){
		mak<-makeqtl(cross,qtl$chr[sub],qtl$pos[sub],qtl$name[sub],what='prob')
		return(fitqtl(cross,trait,mak,method='hk',get.ests=T))
	}
}
##############################################################
##############################################################
#Fitting multi-QTL models - additive+interactions
fitint<-function(qtl,cross){
	#qtl is the result of running getQTL, cross is the cross object
	require(qtl)
	#only qtl list
	qtllist<-qtl[[1]]
	trait<-attr(qtllist,'trait')	
	#interactions
	int<-matrix(qtl[[2]],ncol=3)
	#fit
	if(nrow(qtllist)>0){
		mak<-makeqtl(cross,qtllist$chr,qtllist$pos,qtllist$name,what='prob')
		#formula
		sub<-grep(paste('Q','[1-9]',sep=''),qtllist$name)
		formula<-paste('y ~ ',Reduce(paste,paste(qtllist$name[sub],'+')))
		if(length(sub)==0){formula<-substr(formula,1,nchar(formula)-1)}
		if(nrow(int)>0){
			formula<-paste(formula,Reduce(paste,paste(mak$altname[mak$name%in%qtllist$name[int[,1]]],'*',mak$altname[mak$name%in%qtllist$name[int[,2]]],'+')))
		}
		formula<-substr(formula,1,nchar(formula)-1)
		#fit
		return(fitqtl(cross,trait,mak,formula=formula,method='hk',get.ests=T))
	}
}
##############################################################
##############################################################
##Calculate approximate physical position of pseudomarkers
getpys<-function(QTLscan, markerfile){
	#Marker file should have no header (history), V1=name, V2=chromosome, V3=position
	#Marker names in first column should match row names in QTLscan object ('L11')
	#read in marker information
	markers<-read.table(markerfile)
	#placeholder for physical position
	pyspos<-numeric(nrow(QTLscan))
	#find markers in QTLscan object
	rows<-match(markers$V1,rownames(QTLscan))
	#create physical position of each psudo marker based on position of markers
	for(i in 1:(nrow(markers)-1)){
		if(markers$V2[i]==markers$V2[i+1]){
			dis<-(markers$V3[i+1]-markers$V3[i])/(QTLscan$pos[rows[i+1]]-QTLscan$pos[rows[i]])
			pyspos[rows[i]:rows[i+1]]<-markers$V3[i]+(QTLscan$pos[rows[i]:rows[i+1]]-QTLscan$pos[rows[i]])*dis
		}
	}	
	return(pyspos)
}


##############################################################
##Analysis		             #################################
##############################################################
##############################################################

#Find QTL
qtls.L<-mlply(.fun=getQTL,1:4,out.hk.L,operm.hk.L,out2.hk.L,operm2.hk.L,alpha=c(0.05,0.05))
qtls.H<-mlply(.fun=getQTL,1:2,out.hk.H,operm.hk.H,out2.hk.H,operm2.hk.H,alpha=c(0.05,0.05))

#different threshold
#qtls.L<-mlply(.fun=getQTL,1:4,out.hk.L,operm.hk.L,out2.hk.L,operm2.hk.L,alpha=c(0.1,0.1))
#qtls.H<-mlply(.fun=getQTL,1:2,out.hk.H,operm.hk.H,out2.hk.H,operm2.hk.H,alpha=c(0.1,0.1))

#Fit multi-QTL models
Afits.L<-llply(qtls.L,fitadd,cross=cross.L)
Afits.H<-llply(qtls.H,fitadd,cross=cross.H)
Ifits.L<-llply(qtls.L,fitint,cross=cross.L)
Ifits.H<-llply(qtls.H,fitint,cross=cross.H)
#Calculate pysical position
pys<-getpys(out.hk.L,'./Files/RenameQTLMarkers.txt')

#save workspace
save.image(file='./Data/hkQTLworkspace0.05.Rdata')
#save.image(file='./Data/hkQTLworkspace.Rdata')


		
