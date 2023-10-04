#Analysis for paper: C.albicans white-opaque Sap manuscript
#Naomi Ziv 6/2019

#Organization  
setwd('~/Desktop/MattLohseSap/')
library(flowCore)
library(ggcyto)
#Careful with loading order as want tidyverse (dplyr) filter function by default 
library(tidyverse)
#Also need Biobase installed
	
##############################################################################	
##############################################################################	
#Function for combining plate ID files and formatting for flowCore
#This has a bunch of very specifc formatting matching...
#Not very general
plate<-function(name,prefix='FlowWOcomp',home='~/Desktop/MattLohseSap/'){
	require(tidyverse)
	
	#add backspace for path
	path<-paste(name,'/',sep='')
	
	#loop on plate ID files 
	plate<-data.frame()
	for(i in grep(prefix,list.files(paste(home,path,sep='')),value=T)){
		#get plate information
		pID<-read.csv(paste(home,path,i,sep=''),header=T)   
		pID<-filter(pID,!is.na(Genotype))  
		pID<-separate(pID,Comments, into=c('Time','Dilution','Volume'),sep=';',remove=T,convert=T)
		#get plate number
		num<-str_extract(i,'(?<=\\.)[:digit:]')
		#combine
		pID<-cbind(pID,plate=num)
		plate<-rbind(plate,pID)
	}
	#some more formatting
	plate<-unite(plate,'row',plate,Well,remove=F)
	rownames(plate)<-paste(plate$row,'.fcs',sep='')
	
	return(plate)
	
}	

##############################################################################	
#Function for combining sample ID files and formatting for flowCore
#This has a bunch of very specifc formatting matching...
#Not very general
samples<-function(name,prefix='FlowWOrep',home='~/Desktop/MattLohseSap/'){
	require(tidyverse)
	
	#add backspace for path
	path<-paste(name,'/',sep='')
	
	#loop on plate ID files 
	plate<-data.frame()
	for(i in grep(prefix,list.files(paste(home,path,sep='')),value=T)){
		#get plate information
		pID<-read.csv(paste(home,path,i,sep=''),header=T)   
		pID<-filter(pID,!is.na(Genotype))  
		#get plate number
		num<-str_extract(i,'(?<=\\.)[:digit:]')
		#combine
		pID<-cbind(pID,plate=num)
		plate<-rbind(plate,pID)
	}
	#some more formatting
	plate<-unite(plate,'row',plate,Sample,remove=F,sep='')
	rownames(plate)<-paste(plate$row,'.fcs',sep='')
	
	return(plate)
	
}	


##############################################################################	
#Function for calculating growth rates / lag duration, useful for flow data
#apply per group in data frame - expects data with time/total columns
grcalc<-function(data, size, rsquare=0.9,fold=4){
	#data is data frame with Time/Total
	#size is the size of the sliding window (total size is size+1)
	
	#intialize matrix and output
	info<-matrix(0,nrow(data)-size,5)
	output<-data.frame(intercept=NA,slope=NA,rsquare=NA,time=NA,lag=NA)
	#sliding window
	for (i in 1:(nrow(data)-size)){
		if(sum(!is.na(data$Total[i:(i+size)]))>(size-2)){
			reg<-lm(log(data$Total)[i:(i+size)]~data$Time[i:(i+size)])
			info[i,1]<-coef(reg)[1]
			info[i,2]<-coef(reg)[2]
			info[i,3]<-summary(reg)$adj.r.squared
			info[i,4]<-data$Time[i]
			info[i,5]<-data$Total[i+size]/data$Total[i]
			}else{
				info[i,]<-NA
			}
		}
		
	#find highest slope with decent r^2	and fold
	pass<-info[which(info[,3]>rsquare&info[,5]>fold),]
	pass<-matrix(pass,ncol=5)
	if(nrow(pass)>0){
		temp<-which.max(pass[,2])
		lag<-(log(data$Total[!is.na(data$Total)][1])-pass[temp,1])/pass[temp,2]
		output[1,]<-c(pass[temp,1:4],lag)
	}
		
		return(output)

} 

##############################################################################	
#Little function for use in plotting reporter boxplots with specific whiskers
quants<-function(x){
	result<-quantile(x,probs=c(0.05,0.25,0.5,0.75,0.95))
	names(result)<-c('ymin','lower','middle','upper','ymax')
	return(result)
}	

##############################################################################
##############################################################################
#Filters

#cell filter - for debris - important for Myoglobin
rectGate<-rectangleGate(filterId='cells',"FSC.A"=c(200000,Inf))	

#fluorescence filter - for accuri errors - 0 value 
rectGate2<-rectangleGate(filterId='cells2',"FL1.A"=c(1,Inf))

#for plotting
GRlimit1<-c(2.5*10^5,1.1*10^9)
min_b<-c((1:10)*10^4,(1:10)*10^5,(1:10)*10^6,(1:10)*10^7,(1:10)*10^8)
min_b2<-c((1:10)*10^(-3),(1:10)*10^(-2),(1:10)*10^(-1),(1:10)*10^0,(1:10)*10^1,(1:10)*10^2)

##############################################################################
##############################################################################
#######################							##############################
#######################		Growth Curves		##############################
#######################							##############################
##############################################################################
##############################################################################

#Experiment - 20180612 -
######################### 
######################### 
#White and Opaque cell, different isolates, different protein as nitrogen source 
#Doing Myoglobin separately to use different FSC.A threshold

name<-'20180612'

#get info
info612<-plate(name)
path<-paste(name,'/',sep='')

#get data (not myoglobin samples)
data612<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info612[info612$Environment!='Myoglobin',],data.frame(labelDescription=colnames(info612))),alter.names=T)
#get data (myoglobin samples)
dataMyo<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info612[info612$Environment=='Myoglobin',],data.frame(labelDescription=colnames(info612))),alter.names=T)

#using a threshold of 200000 on FSC.A for filtering out debris (a problem for myoglobin samples)
#ggcyto(data612,aes(x=FSC.A))+geom_histogram()+facet_wrap(~Environment)+scale_x_log10()
#ggcyto(dataMyo,aes(x=FSC.A))+geom_histogram()+facet_wrap(~Genotype)+geom_vline(xintercept=200000)+scale_x_log10()

#filter
dataMyo<-Subset(dataMyo,rectGate)	
	
#count cells and add to info	
counts<-rbind(fsApply(data612,each_col,length),fsApply(dataMyo,each_col,length))
counts<-data.frame(rowname=rownames(counts),Count=counts[,1])
info612<-rownames_to_column(info612) 
info612<-inner_join(info612,counts)	

#new columns for plotting
info612<-mutate(info612,Total=((Count*Dilution)/Volume)*1000,WO=grepl('White',info612$Genotype),Isolate=fct_collapse(Genotype,SC5314=c('Opaque.1','Opaque.2','White.1','White.2'),L26=c('Opaque.L26','White.L26'),P37005=c('Opauqe.P37005','White.P37005')),Strain=str_sub(Genotype,,-3),exp=name)


#Experiment - 20190304 -
######################### 
######################### 
#Opaque cells, Sap deletion strains 

name<-'20190304'

#get info
info304<-plate(name)
path<-paste(name,'/',sep='')

#get data
data304<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info304,data.frame(labelDescription=colnames(info304))),alter.names=T)

#no obvious problem
#ggcyto(data304,aes(x=FSC.A))+geom_histogram()+facet_wrap(~Genotype)+scale_x_log10()

#count cells	
info304<-cbind(info304,Count=fsApply(data304,each_col,length)[,1])
info304<-info304%>%
		rownames_to_column()%>%
		mutate(Total=((Count*Dilution)/Volume)*1000,WO=grepl('White',info304$Genotype),Isolate='SC5314',Strain=str_sub(Genotype,,-3),exp=name)	
	
	
#Experiment - 20190311 -
######################### 
######################### 
#Opaque cells, TF deletions 	

name<-'20190311'

#get info
info311<-plate(name)
path<-paste(name,'/',sep='')

#get data
data311<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info311,data.frame(labelDescription=colnames(info311))),alter.names=T)

#no obvious problem - EFG1 KO strains have high variability
#ggcyto(data311,aes(x=FSC.A))+geom_histogram()+facet_wrap(~Genotype)+scale_x_log10()

#count cells	
info311<-cbind(info311,Count=fsApply(data311,each_col,length)[,1])
info311<-info311%>%
		rownames_to_column()%>%
		mutate(Total=((Count*Dilution)/Volume)*1000,WO=grepl('White',info311$Genotype),Isolate='SC5314',Strain=str_sub(Genotype,,-3),exp=name)	
		
		
#Experiment - 20190327 -
######################### 
######################### 
#White cells, TF deletions 	

name<-'20190327'

#get info
info327<-plate(name)
path<-paste(name,'/',sep='')

#get data
data327<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info327,data.frame(labelDescription=colnames(info327))),alter.names=T)

#no obvious problem - EFG1 KO strains have high variability (less in one of the replicates)
#ggcyto(data327,aes(x=FSC.A))+geom_histogram()+facet_wrap(~Genotype)+scale_x_log10()

#count cells	
info327<-cbind(info327,Count=fsApply(data327,each_col,length)[,1])
info327<-info327%>%
		rownames_to_column()%>%
		mutate(Total=((Count*Dilution)/Volume)*1000,WO=grepl('White',info327$Genotype),Isolate='SC5314',Strain=str_sub(Genotype,,-3),exp=name)	
		
		
#Experiment - 20190321 -
######################### 
######################### 
#opaque cells, STP1 deletion 	

name<-'20190321'

#get info
info321<-plate(name)
path<-paste(name,'/',sep='')

#get data
data321<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info321,data.frame(labelDescription=colnames(info321))),alter.names=T)

#no obvious problem 
#ggcyto(data321,aes(x=FSC.A))+geom_histogram()+facet_wrap(~Genotype+Environment)+scale_x_log10()

#count cells	
info321<-cbind(info321,Count=fsApply(data321,each_col,length)[,1])
info321<-info321%>%
		rownames_to_column()%>%
		mutate(Total=((Count*Dilution)/Volume)*1000,WO=grepl('White',info321$Genotype),Isolate='SC5314',Strain=str_sub(Genotype,,-3),exp=name)	
		
#Experiment - 20190219 -
######################### 
######################### 
#White and Opaque cells, 25 or 37 degrees, also STP1 deletion 
#There also seems to be a white/opaque mix - not used

name<-'20190219'

#get info
info219<-plate(name)
path<-paste(name,'/',sep='')

#get data 
data219<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info219,data.frame(labelDescription=colnames(info219))),alter.names=T)

#no obvious problem 
#ggcyto(data219,aes(x=FSC.A))+geom_histogram()+facet_wrap(~Genotype+Environment)+scale_x_log10()
	
#count cells	
info219<-cbind(info219,Count=fsApply(data219,each_col,length)[,1])
info219<-info219%>%
		rownames_to_column()%>%
		mutate(Total=((Count*Dilution)/Volume)*1000,WO=grepl('White',info219$Genotype),Isolate='SC5314',Strain=str_sub(Genotype,,-3),exp=name)


#Experiment - 20190429 -
######################### 
######################### 
#White and Opaque cells, no glucose, either BSA or BSA+AS, also SAP 5KO deletion 

name<-'20190429'

#get info
info429<-plate(name)
path<-paste(name,'/',sep='')

#get data 
data429<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info429,data.frame(labelDescription=colnames(info429))),alter.names=T)

#no obvious problem 
#ggcyto(data429,aes(x=FSC.A))+geom_histogram()+facet_wrap(~Genotype+Environment)+scale_x_log10()
	
#count cells	
info429<-cbind(info429,Count=fsApply(data429,each_col,length)[,1])
info429<-info429%>%
		rownames_to_column()%>%
		mutate(Total=((Count*Dilution)/Volume)*1000,WO=grepl('White',info429$Genotype),Isolate='SC5314',Strain=str_sub(Genotype,,-3),exp=name)


##############################
#######  Growth rates  #######
##############################

#only looking at BSA/AS 25 degrees
#myoglobin (both White/opaque) seems to have two phases (similar rates) with a pause
#White Hemoglobin seems to have two phases (different rates)  
infoGR<-bind_rows(filter(info612,Time<70),filter(info304,Time<100),info311,info327,info321) %>%
		filter(!Environment%in%c('Stock','ConA','NoN'))%>%
		filter(!Strain%in%c('WhiteEFG1'))%>%
		group_by(Environment,Genotype,Strain,exp)

infoGR2<-group_modify(infoGR,~grcalc(.x,size=4))%>%
			ungroup()

########################
#######  Table   #######
########################

#############Log growth
##########################
GRtable1<- infoGR2%>%
			select(-Genotype,-exp,-intercept,-rsquare,-time)%>%
			group_by(Environment,Strain)%>%
			summarise(GR=mean(slope,na.rm=T),GRsd=sd(slope,na.rm=T),lagsd=sd(lag,na.rm=T),lag=mean(lag,na.rm=T),count=sum(!is.na(slope)),total.count=n())%>%
			mutate(GR.95CI_lower=GR-1.96*((GRsd)/sqrt(count)),GR.95CI_upper=GR+1.96*((GRsd)/sqrt(count)),lag.95CI_lower=lag-1.96*((lagsd)/sqrt(count)),lag.95CI_upper=lag+1.96*((lagsd)/sqrt(count)))%>%
			select(Environment,Strain,GR,GR.95CI_lower,GR.95CI_upper,lag,lag.95CI_lower,lag.95CI_upper,count,total.count)%>%
			mutate_at(3:8,round,2)

write_tsv(GRtable1,'~/dropbox/MattNaomiSAP/GRtable1.txt')



#######################
#######  Plots  #######
#######################
#plot growth curves


pdf('~/Dropbox/MattNaomiSAP/GrowthCurves.1.pdf',width=3.5,height=4.2,colormodel='cmyk')
#############AS verse BSA
##########################
ggplot(filter(info612,Environment%in%c('BSA','BSA.ND','AS')),aes(x=Time,y=Total,color=WO,fill=WO,shape=grepl('BSA',Environment),linetype=Isolate,group=paste(Genotype,Environment)))+
	geom_line(size=0.8)+
	geom_point(size=2,color=1)+
	scale_y_continuous(expression(paste('Cell density (',ml^-1,')')),trans='log',breaks=c(10^5,10^6,10^7,10^8,10^9),limits=GRlimit1,minor_breaks=min_b)+
	scale_x_continuous('Time (hr)',limits=c(0,70))+
	scale_color_manual('Cell type',values=c('#3A54A3','#EB2627'),labels=c('Opaque','White'))+
	scale_fill_manual(values=c('#3A54A3','#EB2627'),labels=c('Opaque','White'),guide='none')+
	scale_shape_manual('Nitrogen',values=c(22,21),labels=c('Ammonium\nSulfate','BSA'))+
 	theme_bw(8)+theme(legend.spacing=unit(0.3,'lines'),legend.position='bottom',legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'))
 	
 
 
#############other proteins
########################## 
ggplot(filter(info612,Environment%in%c("Hemoglobin","HSA","Myoglobin")),aes(x=Time,y=Total,color=WO,fill=WO,shape=factor(Environment,c("HSA","Myoglobin","Hemoglobin")),group=paste(Genotype,Environment)))+
	geom_line(size=0.8)+
	geom_point(size=2,color=1)+
	scale_y_continuous(expression(paste('Cell density (',ml^-1,')')),trans='log',breaks=c(10^5,10^6,10^7,10^8,10^9),limits=GRlimit1,minor_breaks=min_b)+
	scale_x_continuous('Time (hr)',limits=c(0,70))+
	scale_color_manual('Cell type',values=c('#3A54A3','#EB2627'),labels=c('Opaque','White'))+
	scale_fill_manual(values=c('#3A54A3','#EB2627'),labels=c('Opaque','White'),guide='none')+
	scale_shape_manual('Nitrogen',values=c(23,24,25),labels=c("HSA","Myoglobin","Hemoglobin"))+
 	theme_bw(8)+theme(legend.position='bottom',legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'))

dev.off()

pdf('~/Dropbox/MattNaomiSAP/GrowthCurves.1-2.pdf',width=4.9,height=3.3,colormodel='cmyk')
#############Sap deletions
##########################
ggplot(info304,aes(x=Time,y=Total,color=Strain,fill=Strain,group=Genotype))+
	geom_line(size=0.8)+
	geom_point(size=2,color=1,pch=21)+
	geom_line(,filter(info304,Genotype=='Opaque5KO.3'),size=0.8)+
	geom_point(,filter(info304,Genotype=='Opaque5KO.3'),size=2,color=1,pch=21)+
	scale_y_continuous(expression(paste('Cell density (',ml^-1,')')),trans='log',breaks=c(10^5,10^6,10^7,10^8,10^9),limits=GRlimit1,minor_breaks=min_b)+
	scale_x_continuous('Time (hr)',limits=c(0,100),breaks=c(0,20,40,60,80,100))+
	scale_color_manual('Genotype',labels=c('Wild Type',expression(paste(Delta,italic('sap1/2/3/8/99'))), expression(paste(Delta,italic('sap2/3/8/99'))), expression(paste(Delta,italic('sap1/3/8/99'))), expression(paste(Delta,italic('sap1/2/8/99'))), expression(paste(Delta,italic('sap1/2/3/99'))), expression(paste(Delta,italic('sap1/2/3/8')))), values=c('#3A54A3','grey20','#9470b3','#d95f02','#66a61e','#e6ab02','grey90'))+
	scale_fill_manual(values=c('#3A54A3','grey20','#9470b3','#d95f02','#66a61e','#e6ab02','grey80'),guide='none')+
 	theme_bw(8)+theme(legend.text.align=0,legend.direction='vertical',legend.justification='left',legend.key.width=unit(1.6,'lines'),legend.title.align=0.5,legend.key.height=unit(0.7,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'))
 	
dev.off()


pdf('~/dropbox/MattNaomiSAP/GrowthCurves.2.pdf',width=7.2,height=5,colormodel='cmyk')
#############TF deletions
##########################
temp<-filter(rbind(info311,info327,info321),Environment=='BSA'&Strain!="OpaqueWOR3CSR1"&Strain!="WhiteEFG1")
temp$Strain<-factor(temp$Strain)
temp$Strain<-factor(temp$Strain,levels(temp$Strain)[c(1:4,6,5,7:9,11,10)])

ggplot(temp,aes(x=Time,y=Total,color=Strain,fill=Strain,group=paste(Genotype,Environment,exp)))+
	geom_line()+
	geom_point(size=2,color=1,pch=21)+
	scale_y_continuous(expression(paste('Cell density (',ml^-1,')')),trans='log',breaks=c(10^5,10^6,10^7,10^8,10^9),limits=GRlimit1,minor_breaks=min_b)+
	scale_x_continuous('Time (hr)',breaks=seq(0,350,50),limits=range(info321$Time))+
	scale_color_manual('Genotype',labels=c('Opaque - Wild Type',expression(paste(Delta,italic('csr1'))), expression(paste(Delta,italic('efg1'))), expression(paste(Delta,italic('ssy1'))), expression(paste(Delta,italic('wor3'))), expression(paste(Delta,italic('stp1'))), 'White - Wild Type', expression(paste(Delta,italic('csr1'))), expression(paste(Delta,italic('ssy1'))), expression(paste(Delta,italic('wor3'))), expression(paste(Delta,italic('stp1')))), values=c('#3A54A3',"#b3e2cd", '#73dbf0' ,'#edb28a',"#f4cae4",'grey','#EB2627',"#b3e2cd" ,'#edb28a',"#f4cae4",'grey'))+
	scale_fill_manual(values=c('#3A54A3',"#b3e2cd", '#73dbf0' ,'#edb28a',"#f4cae4",'grey','#EB2627',"#b3e2cd" ,'#edb28a',"#f4cae4",'grey'),guide='none')+
	facet_wrap(~factor(WO,,c('Opaque cells','White cells')),ncol=1,scales='free_x')+
 	theme_bw(8)+theme(legend.title.align=0.5,legend.text.align=0,legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8))


dev.off()


pdf('~/Dropbox/MattNaomiSAP/GrowthCurves.3.pdf',width=3.5,height=4.2,colormodel='cmyk')

#############37 degrees
##########################

ggplot(filter(info219,Environment%in%c('BSA.37','BSA')&Strain%in%c('White','Opaque')),aes(x=Time,y=Total,color=Strain,fill=Strain,group=paste(Genotype,Environment),shape=Environment))+
	geom_line(size=0.8)+
	geom_point(size=2,color=1)+
	scale_y_continuous(expression(paste('Cell density (',ml^-1,')')),trans='log',breaks=c(10^5,10^6,10^7,10^8,10^9),limits=GRlimit1,minor_breaks=min_b)+
	scale_x_continuous('Time (hr)',breaks=seq(0,140,20))+
	scale_color_manual('Cell type',values=c('#3A54A3','#EB2627'),labels=c('Opaque','White'))+
	scale_fill_manual(values=c('#3A54A3','#EB2627'),labels=c('Opaque','White'),guide='none')+
	scale_shape_manual('Temperature',labels=c(expression(paste(25*degree,' C')),expression(paste(37*degree,' C'))),values=c(21,22))+	theme_bw(8)+theme(legend.title.align=0.5,legend.position='bottom',legend.direction='vertical',legend.justification='center',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'))
 	
#############no carbon source
##########################

ggplot(filter(info429,Strain%in%c('White','Opaque')),aes(x=Time,y=Total,color=Strain,fill=Strain,group=paste(Genotype,Environment),shape=Environment))+
	geom_line(size=0.8)+
	geom_point(size=2,color=1)+
	scale_y_continuous(expression(paste('Cell density (',ml^-1,')')),trans='log',breaks=c(10^5,10^6,10^7,10^8,10^9),limits=GRlimit1,minor_breaks=min_b)+
	scale_x_continuous('Time (hr)',breaks=seq(0,180,20))+
	scale_color_manual('Cell type',values=c('#3A54A3','#EB2627'),labels=c('Opaque','White'))+
	scale_fill_manual(values=c('#3A54A3','#EB2627'),labels=c('Opaque','White'),guide='none')+
	scale_shape_manual('Media',labels=c('No carbon source\n(has AmS and BSA)','No carbon source\n(has BSA)'),values=c(24,25))+	
	theme_bw(8)+theme(legend.title.align=0.5,legend.position='bottom',legend.direction='vertical',legend.justification='center',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(1.25,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'))
 

dev.off()




##############################################################################
##############################################################################
#######################							##############################
#######################		Mixed cultures		##############################
#######################							##############################
##############################################################################
##############################################################################

#Experiment - 20170918 -
######################### 
######################### 
#Opaque & white cells, fluor swaps,  

name<-'20170918'

#get info
info918<-plate(name)
path<-paste(name,'/',sep='')

#seperate mixes
info918<-separate(info918,Genotype, into=c('p1','Genotype1','p2','Genotype2'),sep='\\.',remove=
F,convert=T)

#get data
#look at subset for testing thresholds
#data918<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info918[(info918$p1==100|info918$p2==100|info918$p1==50)&info918$Environment=='BSA',],data.frame(labelDescription=colnames(info918))),alter.names=T)
#all the data
data918<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info918,data.frame(labelDescription=colnames(info918))),alter.names=T)

#filter and calculate ratio
data918<-Subset(data918,rectGate2)	
data918<-transform(data918,SSCratio=`FL1.A`/`SSC.A`)

#look at flourecence
#ggcyto(data918,aes(x=SSCratio))+geom_histogram(binwidth=0.002)+facet_wrap(~Genotype,scales='free')+ggcyto_par_set(limits=list(x=c(0,0.2),y=c(0,7000)))+geom_vline(xintercept=0.040)

#Using 0.04 as threshold on FL1.A/SSC.A, this does not work for the no nitrogen (NoN) samples

#count cells	
info918<-cbind(info918,Count=fsApply(data918,each_col,length)[,1],GFPos=fsApply(data918,each_col,function(x){sum(x>0.04)})[,'SSCratio'])

#reorganizing
temp<-info918
info918<-info918%>%
		rownames_to_column()%>%
		gather('Which','Geno',Genotype1,Genotype2)%>%
		mutate(	Total=((Count*Dilution)/Volume)*1000,
				gTotal=NA,
				p=NA,
				exp=name)%>%
		select(-p1,-p2)

info918[info918$Which=='Genotype1','gTotal']<-(((temp$Count-temp$GFPos)*temp$Dilution)/temp$Volume)*1000
info918[info918$Which=='Genotype2','gTotal']<-((temp$GFPos*temp$Dilution)/temp$Volume)*1000
info918[info918$Which=='Genotype1','p']<-temp$p1
info918[info918$Which=='Genotype2','p']<-temp$p2
info918<-mutate(info918,WO=grepl('White',info918$Geno)|grepl('Aa',info918$Geno))	
rm(temp)	
				
#Experiment - 20170822 -
######################### 
######################### 
#Opaque & white cells, intial experiment, resembles 918 set but not close enough to put on same graph
#White 80% population is acting strange (not growing well)  

name<-'20170822'

#get info
info822<-plate(name)
path<-paste(name,'/',sep='')

#seperate mixes
info822<-separate(info822,Genotype, into=c('p1','Genotype1','p2','Genotype2'),sep='\\.',remove=
F,convert=T)

#get data
data822<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info822,data.frame(labelDescription=colnames(info822))),alter.names=T)

#filter and calculate ratio
data822<-Subset(data822,rectGate2)	
data822<-transform(data822,SSCratio=`FL1.A`/`SSC.A`)

#look at flourecence
#ggcyto(data822,aes(x=SSCratio))+geom_histogram(binwidth=0.002)+facet_wrap(~Genotype,scales='free')+ggcyto_par_set(limits=list(x=c(0,0.2),y=c(0,7000)))+geom_vline(xintercept=0.040)

#count cells	
info822<-cbind(info822,Count=fsApply(data822,each_col,length)[,1],GFPos=fsApply(data822,each_col,function(x){sum(x>0.04)})[,'SSCratio'])

#reorganizing
temp<-info822
info822<-info822%>%
		rownames_to_column()%>%
		gather('Which','Geno',Genotype1,Genotype2)%>%
		mutate(	Total=((Count*Dilution)/Volume)*1000,
				gTotal=NA,
				p=NA,
				exp=name)%>%
		select(-p1,-p2)

info822[info822$Which=='Genotype1','gTotal']<-(((temp$Count-temp$GFPos)*temp$Dilution)/temp$Volume)*1000
info822[info822$Which=='Genotype2','gTotal']<-((temp$GFPos*temp$Dilution)/temp$Volume)*1000
info822[info822$Which=='Genotype1','p']<-temp$p1
info822[info822$Which=='Genotype2','p']<-temp$p2
info822<-mutate(info822,WO=grepl('White',info822$Geno)|grepl('Aa',info822$Geno))	
rm(temp)	
		
		
#Experiment - 20180605 -
######################### 
######################### 
#Opaque & white cells, also with SAP/STP1 deletions  

name<-'20180605'

#get info
info605<-plate(name)
path<-paste(name,'/',sep='')

#seperate mixes
info605<-separate(info605,Genotype, into=c('p1','Genotype1','p2','Genotype2'),sep='\\.',remove=
F,convert=T)

#get data
data605<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info605,data.frame(labelDescription=colnames(info605))),alter.names=T)

#filter and calculate ratio
data605<-Subset(data605,rectGate2)	
data605<-transform(data605,SSCratio=`FL1.A`/`SSC.A`)

#look at flourecence
#ggcyto(data605,aes(x=SSCratio))+geom_histogram(binwidth=0.002)+facet_wrap(~Genotype,scales='free')+ggcyto_par_set(limits=list(x=c(0,0.2),y=c(0,7000)))+geom_vline(xintercept=0.040)

#count cells	
info605<-cbind(info605,Count=fsApply(data605,each_col,length)[,1],GFPos=fsApply(data605,each_col,function(x){sum(x>0.04)})[,'SSCratio'])

#reorganizing
temp<-info605
info605<-info605%>%
		rownames_to_column()%>%
		gather('Which','Geno',Genotype1,Genotype2)%>%
		mutate(	Total=((Count*Dilution)/Volume)*1000,
				gTotal=NA,
				p=NA,
				exp=name)%>%
		select(-p1,-p2)

info605[info605$Which=='Genotype1','gTotal']<-(((temp$Count-temp$GFPos)*temp$Dilution)/temp$Volume)*1000
info605[info605$Which=='Genotype2','gTotal']<-((temp$GFPos*temp$Dilution)/temp$Volume)*1000
info605[info605$Which=='Genotype1','p']<-temp$p1
info605[info605$Which=='Genotype2','p']<-temp$p2
info605<-mutate(info605,WO=grepl('White',info605$Geno)|grepl('Aa',info605$Geno))	
rm(temp)	

# #reorganizing

# info605<-info605%>%
		# rownames_to_column()%>%
		# mutate(	Total=((Count*Dilution)/Volume)*1000,
				# Total1=(((Count-GFPos)*Dilution)/Volume)*1000,
				# Total2=((GFPos*Dilution)/Volume)*1000,
				# WO1=grepl('White',info605$Genotype1),
				# WO2=grepl('White',info605$Genotype2),
				# exp=name)
	
			
		
#Experiment - 20180924 -
######################### 
######################### 
#Opaque & other candida species 

name<-'20180924'

#get info
info924<-plate(name)
path<-paste(name,'/',sep='')

#seperate mixes
info924<-separate(info924,Genotype, into=c('p1','Genotype1','p2','Genotype2'),sep='\\.',remove=
F,convert=T)

#get data
data924<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info924,data.frame(labelDescription=colnames(info924))),alter.names=T)

#filter and calculate ratio
data924<-Subset(data924,rectGate2)	
data924<-transform(data924,SSCratio=`FL1.A`/`SSC.A`)

#look at flourecence
#ggcyto(data924,aes(x=SSCratio))+geom_histogram(binwidth=0.002)+facet_wrap(~Genotype,scales='free')+ggcyto_par_set(limits=list(x=c(0,0.2),y=c(0,7000)))+geom_vline(xintercept=0.040)

#count cells	
info924<-cbind(info924,Count=fsApply(data924,each_col,length)[,1],GFPos=fsApply(data924,each_col,function(x){sum(x>0.04)})[,'SSCratio'])

#reorganizing
temp<-info924
info924<-info924%>%
		rownames_to_column()%>%
		gather('Which','Geno',Genotype1,Genotype2)%>%
		mutate(	Total=((Count*Dilution)/Volume)*1000,
				gTotal=NA,
				p=NA,
				exp=name)%>%
		select(-p1,-p2)

info924[info924$Which=='Genotype1','gTotal']<-(((temp$Count-temp$GFPos)*temp$Dilution)/temp$Volume)*1000
info924[info924$Which=='Genotype2','gTotal']<-((temp$GFPos*temp$Dilution)/temp$Volume)*1000
info924[info924$Which=='Genotype1','p']<-temp$p1
info924[info924$Which=='Genotype2','p']<-temp$p2
info924<-mutate(info924,WO=grepl('White',info924$Geno)|grepl('Aa',info924$Geno))	
rm(temp)	
				
				
##############################
#######  Growth rates  #######
##############################

M.infoGR<-bind_rows(	filter(info605,Time<75&p!=0&!grepl('White612',info605$Genotype)),
						filter(info918,Environment=='BSA'&Time<140&p!=0&Geno!='Opaque3KO'&Geno!='OpaqueSAP2'&Geno!='White3KO'&Geno!='WhiteSAP2')) %>%
						rename(o.Total=Total,Total=gTotal) %>%
						group_by(Environment,Genotype,Geno,p,exp)

M.infoGR2<-group_modify(M.infoGR,~grcalc(.x,size=4))%>%
			ungroup()

########################
#######  Table   #######
########################

#############Log growth
##########################
GRtable2<- M.infoGR2%>%
			select(-Genotype,-exp,-intercept,-rsquare,-time)%>%
			group_by(Environment,Geno,p)%>%
			summarise(GR=mean(slope,na.rm=T),GRsd=sd(slope,na.rm=T),lagsd=sd(lag,na.rm=T),lag=mean(lag,na.rm=T),count=sum(!is.na(slope)),total.count=n())%>%
			mutate(GR.95CI_lower=GR-1.96*((GRsd)/sqrt(count)),GR.95CI_upper=GR+1.96*((GRsd)/sqrt(count)),lag.95CI_lower=lag-1.96*((lagsd)/sqrt(count)),lag.95CI_upper=lag+1.96*((lagsd)/sqrt(count)))%>%
			select(Environment,Geno,p,GR,GR.95CI_lower,GR.95CI_upper,lag,lag.95CI_lower,lag.95CI_upper,count,total.count)%>%
			mutate_at(4:9,round,2)

write_tsv(GRtable2,'~/dropbox/MattNaomiSAP/GRtable2.txt')
					
		
		
#######################
#######  Plots  #######
#######################		

##Main White verses opaque

pdf('~/dropbox/MattNaomiSAP/Mix.1.pdf',width=7,height=5,colormodel='cmyk')

ggplot(filter(info918,Environment=='BSA'&p!=0&Geno!='Opaque3KO'&Geno!='OpaqueSAP2'&Geno!='White3KO'&Geno!='WhiteSAP2'),aes(x=Time,y=gTotal,color=paste(WO,p),fill=paste(WO,p),group=paste(Genotype,Environment,exp)))+
	geom_line(size=0.8)+
	geom_point(size=2,shape=21,col=1)+
	geom_point(aes(x=Inf,y=Inf,color='NA'))+
	geom_point(aes(x=Inf,y=Inf,color='NA2'))+
	geom_label(aes(x=1,y=4*10^8,label=c('White cells','Opaque cells')),data.frame(WO=c(TRUE,FALSE),Genotype=NA,Environment=NA,p=NA,exp=NA),size=2.5,color=1,fill='white',hjust=0)+
	scale_y_continuous(expression(paste('Cell density (',ml^-1,')')),trans='log',breaks=c(10^4,10^5,10^6,10^7,10^8),minor_breaks=min_b)+
	scale_x_continuous('Time (hr)',breaks=seq(0,140,20),limits=c(0,140))+
	scale_color_manual(expression(bold('Opaque\t\t\t\t\t\t\t\t\t\t\tWhite')),breaks=c('FALSE 100','FALSE 80','FALSE 50','FALSE 20','FALSE 5','NA','NA2','TRUE 20','TRUE 50','TRUE 80','TRUE 95','TRUE 100'),values=c('#3A54A3','#b2df8a','#6a3d9a','#6baed6','#33a02c','white','white','#EB2627','#b15928','#fc9272','#ffff99','#f768a1'),labels=c('100%','80%','50%','20%','5%','0%','0%','20%','50%','80%','95%','100%'))+
	scale_fill_manual(,breaks=c('FALSE 100','FALSE 80','FALSE 50','FALSE 20','FALSE 5','TRUE 20','TRUE 50','TRUE 80','TRUE 95','TRUE 100'),values=c('#3A54A3','#b2df8a','#6a3d9a','#6baed6','#33a02c','#EB2627','#b15928','#fc9272','#ffff99','#f768a1'),guide='none')+
	#scale_color_manual(expression(bold('Opaque\t\t\t\t\t\t\t\t\t\t\tWhite')),breaks=c('FALSE 100','FALSE 80','FALSE 50','FALSE 20','FALSE 5','NA','NA2','TRUE 20','TRUE 50','TRUE 80','TRUE 95','TRUE 100'),values=c('#3A54A3','#9ecae1','#c6dbef','#6baed6','#3182bd','white','white','#EB2627','#fcbba1','#fc9272','#fb6a4a','#d53e4f'),labels=c('100%','80%','50%','20%','5%','0%','0%','20%','50%','80%','95%','100%'))+
	#scale_fill_manual(,breaks=c('FALSE 100','FALSE 80','FALSE 50','FALSE 20','FALSE 5','TRUE 20','TRUE 50','TRUE 80','TRUE 95','TRUE 100'),values=c('#3A54A3','#9ecae1','#c6dbef','#6baed6','#3182bd','#EB2627','#fcbba1','#fc9272','#fb6a4a','#d53e4f'),guide='none')+
	facet_grid(rows=vars(WO),scales='free',space='free_y')+
theme_bw(8)+theme(legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8,color='white'))	+guides(color=guide_legend(ncol=2))

dev.off()		
		
#colors not changed
#822 - not sure we want to do anything with this
ggplot(filter(info822,Environment=='BSA'&p!=0),aes(x=Time,y=gTotal,color=paste(WO,p),fill=paste(WO,p),group=paste(Genotype,Environment,exp)))+
	geom_line(size=0.8)+
	geom_point(size=2,shape=21,col=1)+
	scale_y_continuous(expression(paste('Cell density (',ml^-1,')')),trans='log',breaks=c(10^4,10^5,10^6,10^7,10^8),minor_breaks=min_b)+
	scale_x_continuous('Time (hr)',breaks=seq(0,140,20),limits=c(0,140))+
	scale_color_manual('Starting proportion',breaks=c('FALSE 100','FALSE 80','FALSE 50','FALSE 20','FALSE 5','TRUE 100','TRUE 95','TRUE 80','TRUE 50','TRUE 20'),values=c('#3A54A3','#9ecae1','#c6dbef','#6baed6','#3182bd','#EB2627','#fcbba1','#fc9272','#fb6a4a','#d53e4f'),labels=c('100% Opaque','80%','50%','20%','5%','100% White','95%','80%','50%','20%'))+
	scale_fill_manual(,breaks=c('FALSE 100','FALSE 80','FALSE 50','FALSE 20','FALSE 5','TRUE 100','TRUE 95','TRUE 80','TRUE 50','TRUE 20'),values=c('#3A54A3','#9ecae1','#c6dbef','#6baed6','#3182bd','#EB2627','#fcbba1','#fc9272','#fb6a4a','#d53e4f'),guide='none')+
	facet_grid(rows=vars(WO),scales='free',space='free_y')+
theme_bw(8)+theme(legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8,color='white'))	


##Opaque mutants
info605<-info605%>%
 			separate(Genotype,into=c(NA,'col1',NA,'col2'),remove=F)%>%
 			unite(Mix,col1,col2)%>%
 			mutate(Mix2=factor(Mix,c("White613_Opaque629","OpaqueSTP1_Opaque629","Opaque5KO_Opaque629",  "Opaque5KO_White612"),c('White / Opaque',expression(paste(Delta,italic('stp1'),' Opaque / Opaque')),expression(paste(Delta,italic('sap1/2/3/8/99'),' Opaque / Opaque')),expression(paste(Delta,italic('sap1/2/3/8/99'),' Opaque / White')))))	
 
 
pdf('~/dropbox/MattNaomiSAP/Mix.2.pdf',width=7.8,height=6,colormodel='cmyk')

 			
 			ggplot(filter(info605,Environment=='BSA'&p!=0),aes(x=Time,y=gTotal,color=paste(Geno,p),fill=paste(Geno,p),group=paste(Genotype,Geno,Environment,exp)))+
	geom_line(size=0.8)+
	geom_point(size=2,shape=21,col=1)+
		geom_point(aes(x=Inf,y=Inf,color='NA'))+
		geom_point(aes(x=Inf,y=Inf,color='NA2'))+
		geom_point(aes(x=Inf,y=Inf,color='NA3'))+
		geom_point(aes(x=Inf,y=Inf,color='NA4'))+
	scale_y_continuous(expression(paste('Cell density (',ml^-1,')')),trans='log',breaks=c(10^4,10^5,10^6,10^7,10^8),minor_breaks=min_b)+
	scale_x_continuous('Time (hr)',breaks=seq(0,80,20),limits=c(0,75))+
	scale_color_manual('Starting proportion',breaks=c('White613 100','White613 95','White613 80','White613 50','NA','NA2','OpaqueSTP1 100','OpaqueSTP1 95','OpaqueSTP1 80','OpaqueSTP1 50','NA3','NA4','Opaque5KO 100','Opaque5KO 95','Opaque5KO 80','Opaque5KO 50','Opaque5KO 20','Opaque5KO 5', 'Opaque629 100','Opaque629 50','Opaque629 20','Opaque629 5'),values=c(rep('white',4),'#252525','#bdbdbd','#d9d9d9','#969696','#737373','#525252','#3A54A3','#b2df8a','#6a3d9a','#6baed6','#252525','#969696','#737373','#525252','#EB2627','#fc9272','#ffff99','#f768a1','#EB2627','#fc9272','#ffff99','#f768a1'),labels=c('100% White','95%','80%','50%','','',expression(paste('100% ',Delta,italic('stp1'))),'95%','80%','50%','','',c(expression(paste('100% ',Delta,italic('sap1/2/3/8/99'))),'95%','80%','50%','20%','5%','100% Opaque','50%','20%','5%','100% White\n(fluorescent)','95%','80%','50%')))+
	scale_fill_manual(,breaks=c('White613 100','White613 95','White613 80','White613 50','OpaqueSTP1 100','OpaqueSTP1 95','OpaqueSTP1 80','OpaqueSTP1 50','Opaque5KO 100','Opaque5KO 95','Opaque5KO 80','Opaque5KO 50','Opaque5KO 20','Opaque5KO 5', 'Opaque629 100','Opaque629 50','Opaque629 20','Opaque629 5','White612 100','White612 95','White612 80','White612 50'),values=c('#252525','#bdbdbd','#d9d9d9','#969696','#737373','#525252','#3A54A3','#b2df8a','#6a3d9a','#6baed6','#252525','#969696','#737373','#525252','#EB2627','#fc9272','#ffff99','#f768a1','#EB2627','#fc9272','#ffff99','#f768a1'),guide='none')+
	facet_wrap(~Mix2+factor(Geno,unique(info605$Geno),c('White',expression(paste(Delta,italic('stp1'),' Opaque')),expression(paste(Delta,italic('sap1/2/3/8/99'),' Opaque')),'Opaque  (fluorescent)','White  (fluorescent)')),ncol=4,dir='v',labeller=label_parsed)+
	theme_bw(8)+guides(color=guide_legend(ncol=1))+theme(legend.title.align=0.5,legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),legend.text.align=0,legend.position='bottom')	+guides(color=guide_legend(ncol=4,nrow=6))

dev.off()
	
##Mixes with other species
pdf('~/dropbox/MattNaomiSAP/Mix.3.pdf',width=6,height=9,colormodel='cmyk')

ggplot(filter(mutate(info924,sh=factor(grepl('50',info924$Genotype)+grepl('White612',info924$Genotype))),Environment=='BSA'&p!=0&Geno%in%c('CdSB2','CEM180','CtOpaque1592','CtWhite1582','CtSB113')),aes(x=Time,y=gTotal,shape=sh),group=paste(Genotype,Environment,exp))+
	geom_line(size=0.8)+
	geom_point(size=2,col=1,fill='white')+
	scale_y_continuous(expression(paste('Cell density (',ml^-1,')')),trans='log',breaks=c(10^4,10^5,10^6,10^7,10^8),minor_breaks=min_b)+
	scale_x_continuous('Time (hr)',breaks=seq(0,80,20))+
	scale_shape_manual('Condition',values=c(21,22,23),labels=c('Alone (100%)','Mixed with Opaque cells (50%)','Mixed with White cells (50%)'))+
	facet_wrap(~factor(Geno,c("CdSB2","CEM180","CtSB113","CtWhite1582","CtOpaque1592"),c(expression(italic('C.dubliniensis')),expression(italic('C.parapsilosis')),expression(italic('C.tropicalis')),expression(paste(italic('C.tropicalis'),' - white')),expression(paste(italic('C.tropicalis'),' - opaque')))),ncol=1,labeller=label_parsed)+
theme_bw(8)+theme(legend.title.align=0.5,legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),legend.text.align=0)

dev.off()


	
	
##############################################################################
##############################################################################
#######################							##############################
#######################		Reporters			##############################
#######################							##############################
##############################################################################
##############################################################################


#Experiment - 20190501 -
######################### 
######################### 
#Response to carbon starvation + pilot of mixes 

name<-'20190501'

#get info
info501<-samples(name)
path<-paste(name,'/',sep='')

#seperate mixes
info501<-separate(info501,Genotype, into=c('p1','Genotype1','p2','Genotype2'),sep='\\.',remove=
F,convert=T)

#get data
data501<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info501,data.frame(labelDescription=colnames(info501))),alter.names=T)

#calculate ratio
data501<-transform(data501,SSCratio=`BB515.A`/`SSC.A`)

#look at flourecence
#ggcyto(data501,aes(x=SSCratio))+geom_histogram(binwidth=0.002)+facet_wrap(~Genotype,scales='free')+ggcyto_par_set(limits=list(x=c(0,3),y=c(0,100)))+geom_vline(xintercept=0.5)

#get actual data and bind to other information
DF501<-do.call('rbind',fsApply(data501,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF501<-select(DF501,-Time)
DF501<-info501%>%
		rownames_to_column()%>%
		left_join(DF501)%>%
		mutate( GFPpos=SSCratio>0.5,
				rep=1,
				exp=name)

#summarise data
sum501<-DF501%>%
			group_by(Environment,Genotype,Time,rep,GFPpos,exp,Genotype1)%>%
			summarise(avg=median(PE.CF594.A/SSC.A),count=n())


#Experiment - 20190507 -
######################### 
######################### 
#Mixes

name<-'20190507'

#get info
info507<-samples(name)
path<-paste(name,'/',sep='')

#seperate mixes
info507<-separate(info507,Genotype, into=c('p1','Genotype1','p2','Genotype2'),sep='\\.',remove=
F,convert=T)
info507<-separate(info507,Genotype, into=c('Genotype','rep'),sep='_',remove=
T,convert=T)

#get data
data507<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info507,data.frame(labelDescription=colnames(info507))),alter.names=T)

#calculate ratio
data507<-transform(data507,SSCratio=`BB515.A`/`SSC.A`)

#look at flourecence
#ggcyto(data507,aes(x=SSCratio))+geom_histogram(binwidth=0.002)+facet_wrap(~Genotype,scales='free')+ggcyto_par_set(limits=list(x=c(0,3),y=c(0,100)))+geom_vline(xintercept=0.5)

#get actual data and bind to other information
DF507<-do.call('rbind',fsApply(data507,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF507<-select(DF507,-Time)
DF507<-info507%>%
		rownames_to_column()%>%
		left_join(DF507)%>%
		mutate( GFPpos=SSCratio>0.5,
				exp=name)
				
#summarise data
sum507<-DF507%>%
			group_by(Environment,Genotype,Time,rep,GFPpos,exp,Genotype1)%>%
			summarise(avg=median(PE.CF594.A/SSC.A),count=n())


#Experiment - 20190513 -
######################### 
######################### 
#Opaque cells

name<-'20190513'

#get info
info513<-samples(name)
path<-paste(name,'/',sep='')

#seperate mixes
info513<-separate(info513,Genotype, into=c('p1','Genotype1','p2','Genotype2'),sep='\\.',remove=
F,convert=T)
info513<-separate(info513,Genotype, into=c('Genotype','rep'),sep='_',remove=
T,convert=T)

#get data
data513<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info513,data.frame(labelDescription=colnames(info513))),alter.names=T)

#calculate ratio
data513<-transform(data513,SSCratio=`BB515.A`/`SSC.A`)

#look at flourecence
#ggcyto(data513,aes(x=SSCratio))+geom_histogram(binwidth=0.002)+facet_wrap(~Genotype,scales='free')+ggcyto_par_set(limits=list(x=c(0,3),y=c(0,100)))+geom_vline(xintercept=0.5)

#get actual data and bind to other information
DF513<-do.call('rbind',fsApply(data513,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF513<-select(DF513,-Time)
DF513<-info513%>%
		rownames_to_column()%>%
		left_join(DF513)%>%
		mutate( GFPpos=SSCratio>0.5,
				exp=name)
				
#summarise data
sum513<-DF513%>%
			group_by(Environment,Genotype,Time,rep,GFPpos,exp,Genotype1)%>%
			summarise(avg=median(PE.CF594.A/SSC.A),count=n())



#Experiment - 20190522 -
######################### 
######################### 
#Mixes

name<-'20190522'

#get info
info522<-samples(name)
path<-paste(name,'/',sep='')

#seperate mixes
info522<-separate(info522,Genotype, into=c('p1','Genotype1','p2','Genotype2'),sep='\\.',remove=
F,convert=T)
info522<-separate(info522,Genotype, into=c('Genotype','rep'),sep='_',remove=
T,convert=T)

#get data
data522<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info522,data.frame(labelDescription=colnames(info522))),alter.names=T)

#calculate ratio
data522<-transform(data522,SSCratio=`BB515.A`/`SSC.A`)

#look at flourecence
#ggcyto(data522,aes(x=SSCratio))+geom_histogram(binwidth=0.002)+facet_wrap(~Genotype+Environment,scales='free')+ggcyto_par_set(limits=list(x=c(0,3),y=c(0,100)))+geom_vline(xintercept=0.5)

#get actual data and bind to other information
DF522<-do.call('rbind',fsApply(data522,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF522<-select(DF522,-Time)
DF522<-info522%>%
		rownames_to_column()%>%
		left_join(DF522)%>%
		mutate( GFPpos=SSCratio>0.5,
				exp=name)
				
#summarise data
sum522<-DF522%>%
			group_by(Environment,Genotype,Time,rep,GFPpos,exp,Genotype1)%>%
			summarise(avg=median(PE.CF594.A/SSC.A),count=n())



#######################
#######  Plots  #######
#######################		

#SAP2/uga4
pdf('~/dropbox/MattNaomiSAP/Reporters.1.pdf',width=7,height=3.5,colormodel='cmyk')

ggplot(filter(mutate(ungroup(sum507),se=grepl('50',sum507$Genotype)),Environment=='BSA'&!GFPpos&Genotype1%in%c('WhiteSAP2','WhiteUGA4')),aes(x=Time,y=avg,color=se,fill=se,group=paste(Genotype,rep)))+
	stat_summary(aes(y=PE.CF594.A/SSC.A,group=paste(Genotype,rep,Time),fill=se),data=filter(mutate(DF507,se=grepl('50',DF507$Genotype)),Environment=='BSA'&!GFPpos&Genotype1%in%c('WhiteSAP2','WhiteUGA4')),color=1,position=position_dodge(width=4),width=5,geom='boxplot',fun.data=quants)+
	geom_line(size=0.4,position=position_dodge(width=4))+
	geom_point(size=0.8,shape=21,col=1,position=position_dodge(width=4))+
	scale_y_continuous('Reporter expression (AU)',trans='log',breaks=c(10^(-3),10^(-2),10^(-1),10^(0),10^1,10^2),minor_breaks=min_b2)+
	scale_x_continuous('Time (hr)',breaks=seq(0,140,20),limits=c(-3,142))+
	scale_color_manual('Starting proportion',labels=c('100% White','50% White'),values=c('#EB2627','#fc9272'))+
	scale_fill_manual(values=c('#EB2627','#fc9272'),guide='none')+
	facet_wrap(~factor(Genotype1,,c(expression(paste(italic('SAP2'),' reporter')),expression(paste(italic('UGA4'),' reporter')))),labeller=label_parsed)+
theme_bw(8)+theme(legend.title.align=0.5,legend.background=(element_rect(color=1)),legend.position=c(0.9,0.87),legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8))	
	
dev.off()

#other reporters
pdf('~/dropbox/MattNaomiSAP/Reporters.2.pdf',width=7,height=3.5,colormodel='cmyk')

ggplot(filter(mutate(ungroup(sum507),se=grepl('50',sum507$Genotype)),Environment=='BSA'&!GFPpos&Genotype1%in%c('WhiteOPT1','WhiteOPT2')),aes(x=Time,y=avg,color=se,fill=se,group=paste(Genotype,rep)))+
	stat_summary(aes(y=PE.CF594.A/SSC.A,group=paste(Genotype,rep,Time),fill=se),data=filter(mutate(DF507,se=grepl('50',DF507$Genotype)),Environment=='BSA'&!GFPpos&Genotype1%in%c('WhiteOPT1','WhiteOPT2')),color=1,position=position_dodge(width=4),width=5,geom='boxplot',fun.data=quants)+
	geom_line(size=0.4,position=position_dodge(width=4))+
	geom_point(size=0.8,shape=21,col=1,position=position_dodge(width=4))+
	scale_y_continuous('Reporter expression (AU)',trans='log',breaks=c(10^(-3),10^(-2),10^(-1),10^(0),10^1,10^2),minor_breaks=min_b2)+
	scale_x_continuous('Time (hr)',breaks=seq(0,140,20),limits=c(-3,142))+
	scale_color_manual('Starting proportion',labels=c('100% White','50% White'),values=c('#EB2627','#fc9272'))+
	scale_fill_manual(values=c('#EB2627','#fc9272'),guide='none')+
	facet_wrap(~factor(Genotype1,,c(expression(paste(italic('OPT1'),' reporter')),expression(paste(italic('OPT2'),' reporter')))),labeller=label_parsed)+
theme_bw(8)+theme(legend.title.align=0.5,legend.background=(element_rect(color=1)),legend.position=c(0.9,0.87),legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8))	
	
dev.off()

#carbon limitation
pdf('~/dropbox/MattNaomiSAP/Reporters.3.pdf',width=6,height=5,colormodel='cmyk')

sumRep<-bind_rows(
			filter(ungroup(sum501),Environment=='AS.NoGlu'&!GFPpos&Genotype%in%c('100.WhiteSAP2.0.NA','100.WhiteUGA4.0.NA','100.OpaqueSAP2.0.NA','100.OpaqueUGA4.0.NA')),
			filter(ungroup(sum507),!GFPpos&Genotype%in%c('100.WhiteSAP2.0.NA','100.WhiteUGA4.0.NA')),
			filter(ungroup(sum513),!GFPpos&Genotype%in%c('100.OpaqueSAP2.0.NA','100.OpaqueUGA4.0.NA')))

ggplot(sumRep,aes(x=Time,y=avg,group=paste(Genotype,Environment,rep),shape=Environment,fill=Genotype1))+
	#geom_boxplot(aes(y=PE.CF594.A/SSC.A,group=paste(Genotype,Time),fill=Genotype1),data=filter(DF501,Environment=='AS.NoGlu'&!GFPpos&Genotype1%in%c('WhiteSAP2','WhiteUGA4','OpaqueSAP2','OpaqueUGA4')),outlier.shape=NA,color=1,width=5)+
	geom_line(size=0.5)+
	geom_point(size=2,col=1)+
	scale_y_continuous('Reporter expression (AU)',trans='log',breaks=c(10^(-3),10^(-2),10^(-1),10^(0),10^1,10^2),minor_breaks=min_b2,limits=c(9.775171e-04, 1.130679e+02))+
	scale_x_continuous('Time (hr)',breaks=seq(0,140,20),limits=c(-5,45))+
	scale_shape_manual('Media',labels=c('No carbon source\n(has AmS)','BSA\n(has glucose)','No nitrogen source\n(has glucose)'),values=c(22,21,23))+
	scale_fill_manual(values=rep(c('#3A54A3','#EB2627'),each=2),guide='none')+
	facet_wrap(~factor(Genotype1,,c(expression(paste(italic('SAP2'),' reporter - Opaque cells')),expression(paste(italic('UGA4'),' reporter - Opaque cells')),expression(paste(italic('SAP2'),' reporter - White cells')),expression(paste(italic('UGA4'),' reporter - White cells')))),ncol=2,labeller=label_parsed)+	
theme_bw(8)+theme(legend.title.align=0.5,legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(1.6,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),plot.title=element_text(hjust=0.5))	

dev.off()

#proteinese k
pdf('~/dropbox/MattNaomiSAP/Reporters.4.pdf',width=6,height=5,colormodel='cmyk')

ggplot(filter(ungroup(sum522),!GFPpos),aes(x=Time,y=avg,linetype=Environment,group=paste(Genotype,Environment,rep)))+
	stat_summary(aes(y=PE.CF594.A/SSC.A,group=paste(Genotype,Environment,Time,rep),fill=Genotype1),data=filter(DF522,!GFPpos),color=1,width=5,position=position_dodge(width=4),geom='boxplot',fun.data=quants)+
	geom_line(size=0.5,position=position_dodge(width=4))+
	geom_point(size=1,shape=16,col=1,position=position_dodge(width=4))+
	scale_y_continuous('Reporter expression (AU)',trans='log',breaks=c(10^(-3),10^(-2),10^(-1),10^(0),10^1,10^2),minor_breaks=min_b2,limits=c(9.775171e-04, 1.130679e+02))+
	scale_x_continuous('Time (hr)',breaks=seq(0,140,20))+
	scale_fill_manual(values=rep('#EB2627',4),guide='none')+
	scale_linetype_discrete('Media',labels=c('BSA','BSA + Proteinase K'))+
	facet_wrap(~factor(Genotype1,,c(expression(paste(italic('OPT1'),' reporter - White cells')),expression(paste(italic('OPT2'),' reporter - White cells')),expression(paste(italic('SAP2'),' reporter - White cells')),expression(paste(italic('UGA4'),' reporter - White cells')))),ncol=2,labeller=label_parsed)+	
theme_bw(8)+theme(legend.title.align=0.5,legend.direction='vertical',legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),axis.text=element_text(size=8,color=1),legend.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),plot.title=element_text(hjust=0.5))	

dev.off()