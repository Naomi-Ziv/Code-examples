#Oak-Vineyard QTL paper - data generation
#Organizes sequence data on Oak/Vineyard segregants

#directory with .snp files 
setwd('/Volumes/Aristotle/NYU/Sequencing/SnpFilesOakVinyard/')

#needed for organizing multipool data at the end
library(plyr)

#########################################################
#####     Organizing sequencing data - general    #######
#########################################################
#3 levels of grouping: sample, set and group
SEQdata<-rbind(
	cbind(read.table('F1.snp'),'Rdepth'=NA,sample='F1',set='AIP',dilution=NA,group='setup',time=NA),
	cbind(read.table('F6.snp'),'Rdepth'=NA,sample='F6',set='AIP',dilution=NA,group='setup',time=NA),
	cbind(read.table('F12.snp'),'Rdepth'=NA,sample='F12',set='AIP',dilution=NA,group='setup',time=NA),
	
	cbind(read.table('F12_clone1.snp'),'Rdepth'=NA,sample='F12c1',set='clones',dilution=NA,group='setup',time=NA),
	cbind(read.table('F12_clone2.snp'),'Rdepth'=NA,sample='F12c2',set='clones',dilution=NA,group='setup',time=NA),
	cbind(read.table('F12_clone3.snp'),'Rdepth'=NA,sample='F12c3',set='clones',dilution=NA,group='setup',time=NA),
	
	cbind(read.table('C1S2.snp'),'Rdepth'=NA,sample='C1S2',set='C1',dilution='low',group='F12',time=1650),
	cbind(read.table('C1S4.snp'),'Rdepth'=NA,sample='C1S4',set='C1',dilution='low',group='F12',time=2820),
	cbind(read.table('C1S6.snp'),'Rdepth'=NA,sample='C1S6',set='C1',dilution='low',group='F12',time=4080),
	
	cbind(read.table('C2S2.snp'),'Rdepth'=NA,sample='C2S2',set='C2',dilution='low',group='F12',time=1650),
	cbind(read.table('C2S4.snp'),'Rdepth'=NA,sample='C2S4',set='C2',dilution='low',group='F12',time=2820),
	cbind(read.table('C2S6.snp'),'Rdepth'=NA,sample='C2S6',set='C2',dilution='low',group='F12',time=4080),
	
	cbind(read.table('C3S2.snp'),'Rdepth'=NA,sample='C3S2',set='C3',dilution='High',group='F12',time=1650),
	cbind(read.table('C3S4.snp'),'Rdepth'=NA,sample='C3S4',set='C3',dilution='High',group='F12',time=2820),
	cbind(read.table('C3S6.snp'),'Rdepth'=NA,sample='C3S6',set='C3',dilution='High',group='F12',time=4080),
	
	cbind(read.table('C4S1.snp'),'Rdepth'=NA,sample='C4S1',set='C4',dilution='High',group='F12',time=1410),
	cbind(read.table('C4S4.snp'),'Rdepth'=NA,sample='C4S4',set='C4',dilution='High',group='F12',time=2820),
	cbind(read.table('C4S6.snp'),'Rdepth'=NA,sample='C4S6',set='C4',dilution='High',group='F12',time=4080),
	
	cbind(read.table('C5S0.snp'),'Rdepth'=NA,sample='C5S0',set='C5',dilution='low',group='F2',time=1245),
	cbind(read.table('C5S4.snp'),'Rdepth'=NA,sample='C5S4',set='C5',dilution='low',group='F2',time=3045),
	cbind(read.table('C5S7.snp'),'Rdepth'=NA,sample='C5S7',set='C5',dilution='low',group='F2',time=4485),
	
	cbind(read.table('C6S0.snp'),'Rdepth'=NA,sample='C6S0',set='C6',dilution='High',group='F2',time=1245),
	cbind(read.table('C6S4.snp'),'Rdepth'=NA,sample='C6S4',set='C6',dilution='High',group='F2',time=3045),
	cbind(read.table('C6S6.snp'),'Rdepth'=NA,sample='C6S6',set='C6',dilution='High',group='F2',time=4245),
	
	cbind(read.table('F2SegPool.snp'),'Rdepth'=NA,sample='F2SegPool',set='F2SegPool',dilution=NA,group='F2',time=NA)
	
	)
	
colnames(SEQdata)<-c('chr','pos','ref','alt','qual','RefF','RefR','AltF','AltR','Rdepth','Sample','set','dilution','group','time')
#read depth = total number of quaility reads mapped to base
SEQdata$Rdepth<-rowSums(SEQdata[,c('RefF','RefR','AltF','AltR')])
#remove snps called with more than one alternate allele (12087/1907706)
SEQdata<-subset(SEQdata,sapply(SEQdata$alt,function(x){length(unlist(strsplit(as.character(x),'')))})==1)

#get parental data
vine<-read.table('Vineyard.snp')
oak<-read.table('Oak.snp')
colnames(oak)<- colnames(vine)<-c('chr','pos','ref','alt','qual','RefF','RefR','AltF','AltR')
vine<-cbind(vine,Rdepth=rowSums(vine[,c('RefF','RefR','AltF','AltR')]))
oak<-cbind(oak,Rdepth=rowSums(oak[,c('RefF','RefR','AltF','AltR')]))
#remove snps called with more than one alternate allele
#832/67055 for oak, 624/50202 for vine
oak<-subset(oak,sapply(oak$alt,function(x){length(unlist(strsplit(as.character(x),'')))})==1)
vine<-subset(vine,sapply(vine$alt,function(x){length(unlist(strsplit(as.character(x),'')))})==1)

#define snps and index snps in data
vinepos<-paste(vine$chr,vine$pos,vine$alt)
oakpos<-paste(oak$chr,oak$pos,oak$alt)
snpspos<-paste(SEQdata$chr,SEQdata$pos,SEQdata$alt)
#placeholder for index
SEQdata=cbind(SEQdata,Id=0)
#oak=1,vineyard=2,both=3
SEQdata$Id[snpspos%in%setdiff(oakpos,vinepos)]<-1
SEQdata$Id[snpspos%in%setdiff(vinepos,oakpos)]<-2
SEQdata$Id[snpspos%in%intersect(vinepos,oakpos)]<-3
#correction for positions where oak and vineyard have different alleles from reference (100 positions)
SEQdata$Id[snpspos%in%oakpos[!(oakpos%in%vinepos)&paste(oak$chr,oak$pos)%in%(paste(vine$chr,vine$pos))]]<-1.1
SEQdata$Id[snpspos%in%vinepos[!(vinepos%in%oakpos)&paste(vine$chr,vine$pos)%in%(paste(oak$chr,oak$pos))]]<-2.1
#correction for positions where oak and vineyard have allele frequency<0.75 
#usually indicates repetitive or duplicated region and will confound results in data 
#4714 in snps in oak, 2005 in vineyard
SEQdata$Id[snpspos%in%oakpos[!(oakpos%in%vinepos)&with(oak,((AltF+AltR)/Rdepth)<0.75)]]<-1.2
SEQdata$Id[snpspos%in%vinepos[!(vinepos%in%oakpos)&with(vine,((AltF+AltR)/Rdepth)<0.75)]]<-2.2

#save data
save(SEQdata,file='~/Desktop/QTLpaper/Data/SEQdata.Rdata')


#########################################################
#####     Formating sequencing data - multipool    ######
#########################################################

#subset data - only defined oak and vineyard alleles
SEQdata<-subset(SEQdata,Id==1 | Id==2)
#allele relative to oak
SEQdata<-cbind(SEQdata,allele=NA)
SEQdata$allele[SEQdata$Id==1]<-(SEQdata$AltF+SEQdata$AltR)[SEQdata$Id==1]
SEQdata$allele[SEQdata$Id==2]<-(SEQdata$RefF+SEQdata$RefR)[SEQdata$Id==2]

#For multipool analysis, contrasted samples should have same snps
#chrtxt function that will output files used by multipool - per chromosome and sample 
#three columns: snp position and read counts for ref and alt alleles
chrtxt<-function(data, rename=NA){
	#data id filtered data for coupled samples
	snps<-lapply(split(data,data$Sample),function(data){paste(data$chr,data$pos)})
	if(length(snps)>1){
		inter<-Reduce(intersect,snps)
		data<-data[paste(data$chr,data$pos)%in%inter,]
	}
	if(!is.na(rename)){data$Sample<-factor(data$Sample,labels=rename)}
	for(i in levels(data$Sample)){
		dir.create(paste('~/multipool-master/',i,sep=''))
		temp<-subset(data,Sample==i)
		for(j in levels(data$chr)[1:16]){
			temp3<-subset(temp,chr==j)
			temp3<-cbind(temp3$pos,temp3$allele,temp3$Rdepth-temp3$allele)
			write.table(temp3,file=paste('~/multipool-master/',i,'/',j,'.txt',sep=''),row.names=F,col.names=F)
			rm(temp3)
			}
	rm(temp)
	}
	}
	
#create files for multipool, each time define temp as the subseted data frame and use the chrtxt function
temp<-droplevels(subset(SEQdata,(Sample=='C1S2'|Sample=='C1S6') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp)
temp<-droplevels(subset(SEQdata,(Sample=='C2S2'|Sample=='C2S6') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp)
temp<-droplevels(subset(SEQdata,(Sample=='C3S2'|Sample=='C3S6') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp)
temp<-droplevels(subset(SEQdata,(Sample=='C4S1'|Sample=='C4S6') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp)
temp<-droplevels(subset(SEQdata,(Sample=='C5S0'|Sample=='C5S7') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp)
temp<-droplevels(subset(SEQdata,(Sample=='C6S0'|Sample=='C6S6') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp)
temp<-droplevels(subset(SEQdata,(Sample=='F1') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp)
temp<-droplevels(subset(SEQdata,(Sample=='F2SegPool') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp)
temp<-droplevels(subset(SEQdata,(Sample=='F12') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp)
temp<-droplevels(subset(SEQdata,(Sample=='C1S2'|Sample=='C2S2') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp,c('NullLow1','NullLow2'))
temp<-droplevels(subset(SEQdata,(Sample=='C3S2'|Sample=='C4S1') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp,c('NullHigh1','NullHigh2'))
temp<-droplevels(subset(SEQdata,(Sample=='C5S0'|Sample=='C6S0') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
chrtxt(temp,c('NullF21','NullF22'))

#Combining replicate data - low glucose condition
#make sure contrasted samples have same snps
low<-droplevels(subset(SEQdata,(Sample=='C1S2'|Sample=='C1S6'|Sample=='C2S2'|Sample=='C2S6') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
snps<-lapply(split(low,low$Sample),function(data){paste(data$chr,data$pos)})
inter<-Reduce(intersect,snps)
low<-low[paste(low$chr,low$pos)%in%inter,]
for(i in c('ComLow1','ComLow2')){
	dir.create(paste('~/multipool-master/',i,sep=''))
	for(j in levels(low$chr)[1:16]){
		if(i=='ComLow1'){
			temp1<-subset(low,Sample=='C1S2'&chr==j)
			temp2<-subset(low,Sample=='C2S2'&chr==j)
			
		}
		if(i=='ComLow2'){
			temp1<-subset(low,Sample=='C1S6'&chr==j)
			temp2<-subset(low,Sample=='C2S6'&chr==j)
		}
		temp3<-cbind(temp1$pos,temp1$allele+temp2$allele,(temp1$Rdepth-temp1$allele)+(temp2$Rdepth-temp2$allele))
		write.table(temp3,file=paste('~/multipool-master/',i,'/',j,'.txt',sep=''),row.names=F,col.names=F)
		rm(temp1,temp2,temp3)
	}
}

#Combining replicate data - high glucose condition
#make sure contrasted samples have same snps
high<-droplevels(subset(SEQdata,(Sample=='C3S2'|Sample=='C3S6'|Sample=='C4S1'|Sample=='C4S6') & qual>100 & Rdepth>20 & allele/Rdepth>0.1 & allele/Rdepth<0.9))
snps<-lapply(split(high,high$Sample),function(data){paste(data$chr,data$pos)})
inter<-Reduce(intersect,snps)
high<-high[paste(high$chr,high$pos)%in%inter,]
for(i in c('ComHigh1','ComHigh2')){
	dir.create(paste('~/multipool-master/',i,sep=''))
	for(j in levels(high$chr)[1:16]){
		if(i=='ComHigh1'){
			temp1<-subset(high,Sample=='C3S2'&chr==j)
			temp2<-subset(high,Sample=='C4S1'&chr==j)
			
		}
		if(i=='ComHigh2'){
			temp1<-subset(high,Sample=='C3S6'&chr==j)
			temp2<-subset(high,Sample=='C4S6'&chr==j)
		}
		temp3<-cbind(temp1$pos,temp1$allele+temp2$allele,(temp1$Rdepth-temp1$allele)+(temp2$Rdepth-temp2$allele))
		write.table(temp3,file=paste('~/multipool-master/',i,'/',j,'.txt',sep=''),row.names=F,col.names=F)
		rm(temp1,temp2,temp3)
	}
}


###############################################################
#####     Formating multipool output data 				 ######
#####	  After multipool was run using doMultipool.sh   ######
###############################################################

#Set up combinations
pargrid<-expand.grid(tag=c('C1','C2','C3','C4','C5','C6','F12','ComLow','ComHigh','NullLow','NullHigh','NullF2'),chr=levels(SEQdata$chr)[1:16],RF=c(1000,2500),N=c(200,1000))
#Get data
Multipool<-mdply(pargrid,
	function(tag,chr,RF,N){
		temp<-read.table(paste('~/multipool-master/out/',tag,'.',chr,'.txt.',RF,'.',N,sep=''),header=T,sep='\t')
		temp<-cbind(temp,tag=tag,chr=chr,RF=RF,N=N)
		return(temp)})
#Clean up column names		
colnames(Multipool)<-c('bin','AF','LOD','tag','chr','RF','N')
#multipool$chr<-as.numeric(multipool$chr)
#save data
save(Multipool,file='~/Desktop/QTLpaper/Data/Multipool.Rdata')


