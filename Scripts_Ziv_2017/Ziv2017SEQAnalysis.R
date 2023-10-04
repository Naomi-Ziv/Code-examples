#Used for paper: Resolving the complex genetic basis of phenotypic variation and variability of cellular growth, Genetics, 2017
#Naomi Ziv, nz375@nyu.edu

#Set working directory
setwd('~/Desktop/QTLpaper')

#Additional code:
#####Formating sequencing data - NZ_align_example.q, NZ_snps_example.q and snp.py
#input is .fastq read files, output is .snp file 
#.q (HPC) scripts represent examples for a single sample, originally scripts iterated on all samples
#####Organizing sequencing data, in general and for multipool analysis - Oak.Vineyard.multipool.datagen.R
#input is .snp files, 
#output is 1) SEQdata.Rdata, 2) txt files per chromosome and sample with snp positions and read counts - used for multipool analysis and 3) Multipool.Rdata
#####Multipool analysis on multiple samples - doMultipool.sh
#input is txt files created by Oak.Vineyard.multipool.datagen.R and output are multipool output files
#uses mp_inference.py  
#####HXT growth rate - organized by - Oak.Vineyard.HXT.datagen.R
#input is microcolony growth output files
#output is HXTdata.Rdata (no addtional analysis performed, used for figures) 

#load SNP data 
load('./Data/SEQdata.Rdata')
#load mapping data (needed to count recombination events in F2s)
load('./Data/QTLmapping.Rdata')
#cross.L, cross.H, out.hk.L, out.hk.H, out2.hk.L, out2.hk.H, operm.hk.L, operm.hk.H, operm2.hk.L, operm2.hk.H
#load Multipool data
load('./Data/Multipool.Rdata')

#libraries
library(qtl)
library(plyr)

##############################################################
#Custom functions            #################################
##############################################################
##############################################################

#function to count recombination events
recombin<-function(data,out=T){
	#data is filtered data frame per chr (and sample) with column Id
	#or vector of Id values (numeric)
	#function will always give number of transitions
	#if out=T, output will also have number of transitions excluding transitions based on one data point
	
	#get only vector of Ids and get rid of NAs
	if(!is.null(dim(data))){
		data<-data$Id
	}
	data<-data[!is.na(data)]
	
	#get transitions (did the Id change)
	trans<-data[2:length(data)]-data[1:(length(data)-1)]
	#count transitions - always output
	out1<-sum(trans!=0)
	
	#disregard transitions based on one data point
	#two subsequent transitions are disregarded (because you switch and switch back)
	#three subsequent transitions count as a single transition
	out2<-out1
	#get positions of transitions
	temp<-which(trans!=0)
	while(length(temp)>1){
		if(temp[2]-temp[1]==1){
			if(length(temp==2)){
				temp<-NA
			}else{
				temp<-temp[3:length(temp)]
			}
			out2<-out2-2
		}else{
			temp<-temp[2:length(temp)]
		}	
	}
	
	#return
	if(out){return(c(out1,out2))}else{return(out1)}
	#return(which(trans!=0))
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

#########################################################
#####     Advanced intercross population           ######
#########################################################

#subset data - only defined oak and vineyard alleles
SEQdata<-subset(SEQdata,Id==1 | Id==2)
#allele relative to oak
SEQdata<-cbind(SEQdata,allele=NA)
SEQdata$allele[SEQdata$Id==1]<-(SEQdata$AltF+SEQdata$AltR)[SEQdata$Id==1]
SEQdata$allele[SEQdata$Id==2]<-(SEQdata$RefF+SEQdata$RefR)[SEQdata$Id==2]

#create continous position across genome
SEQdata<-cbind(SEQdata,cpos=0)
SEQdata$cpos[which(SEQdata$chr=='chr01')]<-subset(SEQdata,chr=='chr01')$pos
for (i in 1:17){
	SEQdata$cpos[which(SEQdata$chr==levels(SEQdata$chr)[i+1])]<-subset(SEQdata,chr==levels(SEQdata$chr)[i+1])$pos+max(SEQdata$cpos[which(SEQdata$chr==levels(SEQdata$chr)[i])])
}

#####number of snps with minor allele frequency>10% in F12 compared to F1
nrow(subset(SEQdata,Sample=='F12'&(AltF+AltR)/Rdepth<0.9))/nrow(subset(SEQdata,Sample=='F1'&(AltF+AltR)/Rdepth<0.9))

#####Calling the number of recombination events (transtions between oak and vineyard snps) in the 3 F12 clones
#filter for read depth and quality 
#also filter ~0.3% of snps (91689/92028) with allele frequency between 0.25-0.75 that may be due to duplications or hetrozygosity 
clonedata<-droplevels(subset(SEQdata,set=='clones'&Rdepth<300&Rdepth>20&(allele/Rdepth<0.25|allele/Rdepth>0.75)&chr!='chrM'&qual>100))
summary(clonedata$Sample)
#F12c1 F12c2 F12c3 
#30861 29204 31624
#also filter snps that are not segregating (~3%)
clonedata<-clonedata[!clonedata$cpos%in%subset(SEQdata,Sample=='F1'&(allele/Rdepth>0.75|allele/Rdepth<0.25))$cpos,]
summary(clonedata$Sample)
#F12c1 F12c2 F12c3 
#29950 28446 30543 

#split data according to sample and chromosome  
clonesplit<-split(clonedata,list(clonedata$Sample,clonedata$chr))

#count recombination events for each chromosome
Revents<-sapply(clonesplit,recombin)

#events per clone 
rowSums(Revents[,seq(from=1,by=3,length.out=16)])
rowSums(Revents[,seq(from=2,by=3,length.out=16)])
rowSums(Revents[,seq(from=3,by=3,length.out=16)])

#####Calling the number of recombination events (transtions between oak and vineyard snps) in the 374 F2s
F2s<-matrix(,16,374)
for(i in 1:16){
	F2s[i,]<-apply(cross.L[[1]][[i]][[1]],1,recombin,out=F)
}

#mean and range
mean(colSums(F2s))
#31.9385
range(colSums(F2s))
#13 55


#####Compare interval mapping and multipool by looking at maximum LOD score in intervals defined by markers 
#get marker locations
markers<-read.table('./Files/RenameQTLMarkers.txt')

#get intervals between markers and positions
maxintlod<-cbind(	ldply(split(which(rownames(out.hk.L)%in%markers$V1),markers$V2),function(x){data.frame(x1=x[1:(length(x)-1)],x2=x[2:length(x)])}),
					ldply(split(markers,markers$V2),function(x){data.frame(name1=x$V1[1:(nrow(x)-1)],name2=x$V1[2:nrow(x)],pys1=x$V3[1:(nrow(x)-1)],pys2=x$V3[2:nrow(x)])}))
#get LOD score max in each interval based on QTL mapping, for low and high glucose
maxintlod<-cbind(maxintlod,mdply(maxintlod[,2:3],function(x1,x2){data.frame(GR.PC.L=max(out.hk.L[x1:x2,'GR.PC']),GR.PC.H=max(out.hk.H[x1:x2,'GR.PC']))}))
#get LOD score max in each interval based on Multipool
temp<-droplevels(subset(Multipool,tag=='C5'|tag=='C6'|tag=='ComLow'|tag=='ComHigh'))
maxintpool<-ldply(split(temp,list(temp$tag,temp$RF,temp$N)),
	function(data){
		maxs<-numeric(nrow(maxintlod))
		for(i in 1:nrow(maxintlod)){
			maxs[i]<-max(subset(data,data$chr==levels(Multipool$chr)[as.numeric(maxintlod$.id[i])]&data$bin>maxintlod$pys1[i]&data$bin<maxintlod$pys2[i])$LOD)
		}
		return(data.frame(maxs=maxs))
		})
#organize
maxintlod<-cbind(maxintlod,Reduce(cbind,split(maxintpool$maxs,maxintpool$.id)))
colnames(maxintlod)[13:28]<-levels(factor(maxintpool$.id))

#number of intervals
nrow(maxintlod)
#size of intervals
range(with(maxintlod,pys2-pys1))
#14103 139320
mean(with(maxintlod,pys2-pys1))
#50853.37
median(with(maxintlod,pys2-pys1))
#49653

#what is shared between specific conditions (all intervals with LOD>3 reported in table S3)
#advanced intercross and interval mapping
maxintlod[(which(maxintlod$GR.PC.L>3&maxintlod$ComLow.1000.1000>3)),1:8]
maxintlod[(which(maxintlod$GR.PC.H>3&maxintlod$ComHigh.1000.1000>3)),1:8]
#bulk F2 and interval mapping
maxintlod[(which(maxintlod$GR.PC.L>3&maxintlod$C5.1000.1000>3)),1:8]
maxintlod[(which(maxintlod$GR.PC.H>3&maxintlod$C6.1000.1000>3)),1:8]
#only in F2s
maxintlod[(which(maxintlod$GR.PC.L<3&maxintlod$ComLow.1000.1000<3&maxintlod$C5.1000.1000>3)),1:8]
maxintlod[(which(maxintlod$GR.PC.H<3&maxintlod$ComHigh.1000.1000<3&maxintlod$C6.1000.1000>3)),1:8]


#####Compare 2-lod drop intervals for QTL on chromosome 4 between methods

#Calculate pysical position for interval mapping
pys<-getpys(out.hk.L,'./Files/RenameQTLMarkers.txt')

#interval mapping, low glucose - (67.49597 kb)
pys[with(out.hk.L,min(which(chr==4&pos>max(out.hk.L,chr=4)$pos&GR.PC<(max(out.hk.L,chr=4)$GR.PC-2))))]-pys[with(out.hk.L,max(which(chr==4&pos<max(out.hk.L,chr=4)$pos&GR.PC<(max(out.hk.L,chr=4)$GR.PC-2))))]

#interval mapping, high glucose - (245.8889 kb)
pys[with(out.hk.H,min(which(chr==4&pos>max(out.hk.H,chr=4)$pos&GR.PC<(max(out.hk.H,chr=4)$GR.PC-2))))]-pys[with(out.hk.H,max(which(chr==4&pos<max(out.hk.H,chr=4)$pos&GR.PC<(max(out.hk.H,chr=4)$GR.PC-2))))]

#bulk segregant mapping
#function to calculate interval
getint<-function(chromo,sample,pop,rf){
	#get data on peak
	temp<-subset(Multipool,tag==sample&N==pop&RF==rf&chr==chromo)[with(subset(Multipool,tag==sample&N==pop&RF==rf&chr==chromo),which.max(LOD)),]
	#calculate interval
	min(subset(Multipool,chr==chromo&tag==sample&LOD<(temp[,'LOD']-2)&bin>temp[,'bin']&N==pop&RF==rf)$bin)-max(subset(Multipool,chr==chromo&tag==sample&LOD<(temp[,'LOD']-2)&bin<temp[,'bin']&N==pop&RF==rf)$bin)
}

#advanced intercross
getint('chr04','ComLow',1000,1000)
getint('chr04','ComHigh',1000,1000)
#pooled F2s
getint('chr04','C5',1000,1000)
getint('chr04','C6',1000,1000)
getint('chr04','C5',200,2500)
getint('chr04','C6',200,2500)

