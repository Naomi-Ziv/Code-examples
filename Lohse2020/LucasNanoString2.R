#Analysis for paper: C.albicans white-opaque Sap manuscript
#Naomi Ziv 4/2019, updated 8/2019, 2/2020

#Organization  
setwd('~/Dropbox/Brenes.et.al.sap/LucasNanostring/')
#setwd('/Volumes/Dodgson/Papers/Brenes.et.al.sap/LucasNanostring/')
#Libraries
library(tidyverse)

#also uses 'squish' from 'scales'

#Read samples
samples<-read_csv('LucasSamples.csv',col_types=cols(.default=col_factor(NULL)))

#Read data
x=read_tsv('NanostringRaw.txt')
xx<-gather(x,Sample,reading,-MessageCat,-Name,-Gene,factor_key=T)
data<-cbind(xx,samples[match(xx$Sample,samples$Sample),2:4])

#Filter data to reflect only relevent samples:
data<-data%>%
			#all coculture samples are 16hr - not well controlled experiment, later used reporters
			filter(Time!='16hr')%>%
			#not using various media conditions, some have only single repeats (removing PBS is important for normalization as HTA1 is changing drastically)
			filter(!Media%in%c('PBS','AS+BSA','SD.AA.uri','AS.PepA.5','AS.DMSO','BSA.PepA.5','BSA.DMSO','BSA.notfilter','BSA.filter.notused','BSA.dialyzed.notused'))%>%
			#not using early 'BSA' conditions, afterwards everything BSA.filter (except for Stp1.T at log and TF052(stp1) at 2hr)
			filter(Media!='BSA'|(Strain%in%c('Stp1.Tw','Stp1.To','TF052w','TF052o')&Media=='BSA'))%>%
			#not looking at EFG1 mutations, SAP combinations (other than quintuple) and WOR3/CSR1 double
			filter(!Strain%in%c('Sap.3.10.99w','Sap.3.10.99o','Sap.2.3.10.99w','Sap.1.2.3.99w','Sap.1.2.3.99o','Sap.2.3.10.99o','EFG1.T208Ew','EFG1.T208Eo','EFG1.T208Aw','EFG1.T208Ao','WOR3del.TF106.w','WOR3del.TF106.o'))%>%
			#addtional removal of samples that ended up not being used in the paper (3/2020)
			filter(!Strain%in%c('Sap.1.2.3.8.99w','TF132o','TF140o'))%>%
			filter(!(Strain=='TF052o'&Media%in%c('BSA.filter.noUri','NoNitro.noUri')))%>%
			filter(!(Strain=='TF052w'&Media=='BSA'))%>%
			droplevels()
			

#calculate geometric mean for each gene based on all samples
data<-data%>%
		group_by(Name)%>%
		summarise(gmean=mean(log(reading)))%>%
		full_join(data)%>%
		ungroup()
		
#look at positive controls
ggplot(filter(data,MessageCat=='Positive'&Name!='POS_F'),aes(x=gmean,y=reading,group=Sample))+geom_point(alpha=0.1)+geom_line(alpha=.2)+scale_y_continuous(trans='log')+geom_abline(intercept=0,slope=1,color='red')

#look at reference controls
ggplot(filter(data,Name=='PGA59'|Name=='MTLa1'|Name=='TBP1'|Name=='HTA1'),aes(x=gmean,y=reading,group=Sample))+geom_point(alpha=0.1)+geom_line(alpha=.2)+scale_y_continuous(trans='log',limits=range(data$reading),breaks=c(10,100,1000,10000,100000))+geom_abline(intercept=0,slope=1,color='red')

ggplot(filter(data,Name=='PGA59'|Name=='MTLa1'|Name=='TBP1'|Name=='HTA1'),aes(x=gmean,y=reading,group=Sample))+geom_point(alpha=0.1)+geom_line(alpha=.2)+scale_y_continuous(trans='log',limits=range(data$reading))+geom_abline(intercept=0,slope=1,color='red')+facet_wrap(~Media+Strain)

#normalize - subtract diffrence in reference mean per sample to grand reference mean 
data<-data%>%
		group_by(Sample)%>%
		filter(Name=='PGA59'|Name=='MTLa1'|Name=='TBP1'|Name=='HTA1')%>%
		summarise(ref=mean(filter(data,Name=='PGA59'|Name=='MTLa1'|Name=='TBP1'|Name=='HTA1')$gmean)-mean(log(reading)))%>%
		full_join(data)%>%
		ungroup()%>%
		mutate(readingN=exp(log(reading)+ref))

#look at normalized reference controls
ggplot(filter(data,Name=='PGA59'|Name=='MTLa1'|Name=='TBP1'|Name=='HTA1'),aes(x=gmean,y=readingN,group=Sample))+geom_point(alpha=0.1)+geom_line(alpha=.2)+scale_y_continuous(trans='log',limits=range(data$reading))+geom_abline(intercept=0,slope=1,color='red')

#remove controls
data2<-filter(data,MessageCat=='Endogenous')

#Summarise replicates - 
summ<-data2%>%
		group_by(Name,Strain,Media,Time)%>%
		summarise(Mean=mean(readingN),SD=sd(readingN),N=n(),Min=min(readingN),Max=max(readingN))%>%
		ungroup()	

summ2<-data2%>%
		mutate(time=Time,strain=Strain,media=Media)%>%
		group_by(Name,Strain,Media,Time)%>%
		nest()
		
#Calculate stats for WT AS vrs BSA 
#fold change BSA/AS
FC<- summ%>%
		select(-SD,-N,-Min,-Max)%>%
		filter(Strain%in%c('WTo','WTw')&Media%in%c('AS','BSA.filter'))%>%
		spread(Media,Mean)%>%
		mutate(B.A=BSA.filter/AS)
			
NstatsWT<- data2%>% 
			filter(Strain%in%c('WTw','WTo')&Media%in%c('AS','BSA.filter'))%>%
			group_by(Name,Strain,Time)%>%
			summarise(pval=t.test(log(readingN)~Media,var.equal=T)$p.val)%>%
			ungroup()%>%
			unite(Strain,Time,col=panel,remove=F)%>%
			arrange(pval)%>%
			mutate(FDR.5=0.05*((1:nrow(FC))/nrow(FC)),FDR.1=0.01*((1:nrow(FC))/nrow(FC)))%>%
			mutate(FDR2.5=pval<FDR.5,FDR2.1=pval<FDR.1)%>%
			full_join(FC)			
					
###########################################################
#######  Exporting data - all genes and conditions  #######
###########################################################
#Get Matt names
Mnames<-read_tsv('~/dropbox/MattNaomiSAP/SampleNames.txt')
Strain1<-c("WTw","WTo","Stp1.Tw","Stp1.To","TF052w","TF052o","EFG1del.Stp1act.w","EFG1del.Stp1act.o","WOR3del.Stp1act.w","WOR3del.Stp1act.o","TF156w","TF156o","TF170w","TF170o","Sap.1.2.3.8.99o","TF106o","Ssy1o","TF106.Stp1act.o","TF106w","TF106.Stp1act.w")
Strain2<-c("WT.w","WT.o","Stp1act.w","Stp1act.o","dStp1.w","dStp1.o","dEfg1.Stp1act.w","dEfg1.Stp1act.o","dWor3.Stp1act.w","dWor3.Stp1act.o","dEfg1.w","dEfg1.o","dWor3.w","dWor3.o","dSap.1.2.3.8.99.o","dCsr1.o","dSsy1.o","dCsr1.Stp1act.o","dCsr1.w","dCsr1.Stp1act.w")

#Genes as rows, all conditions as columns - raw data
Mattraw<-data2%>%
		select(Name,Gene,reading,Strain,Media,Time,Sample)%>%
		unite(Strain, Media, Time,col='Condition',sep='_')%>%
		mutate(Condition2=factor(Condition,Mnames$list,Mnames$name))%>%
		unite(Condition2,Sample,col='Condition3',sep=' ')%>%
		select(-Condition)%>%
		spread(Condition3,reading)
		
write_tsv(Mattraw,'~/dropbox/MattNaomiSAP/WiderRaw.txt')


Mattnorm<-data2%>%
		select(Name,Gene,readingN,Strain,Media,Time,Sample)%>%
		mutate(readingN=round(readingN,digits=1))%>%
		unite(Strain, Media, Time,col='Condition',sep='_')%>%
		mutate(Condition2=factor(Condition,Mnames$list,Mnames$name))%>%
		unite(Condition2,Sample,col='Condition3',sep=' ')%>%
		select(-Condition)%>%
		spread(Condition3,readingN)
		
write_tsv(Mattnorm,'~/dropbox/MattNaomiSAP/WiderNorm.txt')


Mattrawrep<-summ%>%
		select(Name,Mean,Strain,Media,Time)%>%
		mutate(Mean=round(Mean,digits=1))%>%
		unite(Strain, Media, Time,col='Condition',sep='_')%>%
		mutate(Condition2=factor(Condition,Mnames$list,Mnames$name))%>%
		select(-Condition)%>%
		spread(Condition2,Mean)
		
write_tsv(Mattrawrep,'~/dropbox/MattNaomiSAP/WideNorm.txt')

Naomi<-data2%>%
		select(Name,Gene,Strain,Media,Time,Sample,reading,readingN)%>%
		mutate(Strain=factor(Strain,Strain1,Strain2))%>%
		mutate(readingN=round(readingN,digits=1))

write_tsv(Naomi,'~/dropbox/MattNaomiSAP/LongNorm.txt')

##################################################
#######  Comparing AS/BSA WT White/Opaque  #######
##################################################


#############
################

#############################################
#######  Genes and colors and ranges  #######
#############################################

#Range for all-gene plots (Wildtypes in AS/BSA)
ran1<-range(c(filter(summ,Strain%in%c('WTo','WTw'))$Min,filter(summ,Strain%in%c('WTo','WTw'))$Max))

#Genes for main color bar plots
genes1<-c('SAP2','SAP3','OPT2','PTR2','OPT5','SAP8','OPT7','OPT1','SAP1','SAP99','OPT4')
#Genes for expanded colorbar plot
genes2<-c('SAP2','SAP3','OPT2','PTR2','OPT5','SAP8','OPT7','OPT1','SAP1','SAP99','OPT4','UGA4')
#genes2<-c('SAP2','SAP3','OPT2','PTR2','OPT5','SAP8','OPT7','OPT1','SAP1','SAP99','OPT4','UGA4','DAL5','MEP2','SSY1','STP1','SAP98')

#Colors for protein and nitrogen depletion plots 
#colors<-c(SAP1='#33a02c', SAP2='#6a3d9a', SAP3='#ff7f00', SAP8='#1f78b4', OPT1='#e31a1c', OPT2='#b2df8a', OPT5='#ffff99', OPT7='#a6cee3', PTR2='#fb9a99', UGA4='#b15928', DAL5='#cab2d6', MEP2='#fdbf6f')
#colors<-c(SAP1='#33a02c', SAP2='#6a3d9a', SAP3='#ff7f00', SAP8='#1f78b4', OPT1='#e31a1c', OPT2='#b2df8a', OPT5='#ffff99', OPT7='#a6cee3', PTR2='#fb9a99', UGA4='#b15928', OPT3='#cab2d6')
#colors=brewer.pal(12,'Set3')[c(3:8,10:12,1:2)]
colors=c(rainbow(9)[c(1:8)],'brown','grey')
names(colors)=c('OPT1','OPT2','OPT5','OPT7','PTR2','SAP1','SAP2','SAP3','SAP8','UGA4')

####################################
#######  Tables - all genes  #######
####################################

#############Log growth
##########################
logtable<-summ%>%
		select(-SD,-N,-Min,-Max)%>%
		filter(Strain%in%c('WTo','WTw')&Media%in%c('AS','BSA.filter')&Time=='log')%>%
		spread(Strain,Mean)%>%
		nest(WTw,WTo)%>%
		spread(Media,data)%>%
		unnest()%>%
		mutate(OW.AS=log2(WTo/WTw),OW.BSA=log2(WTo1/WTw1),BA.W=log2(WTw1/WTw),BA.O=log2(WTo1/WTo))%>%
		mutate(OW.AS=sign(OW.AS)*2^(abs(OW.AS)),OW.BSA=sign(OW.BSA)*2^(abs(OW.BSA)),BA.W=sign(BA.W)*2^(abs(BA.W)),BA.O=sign(BA.O)*2^(abs(BA.O)))	

colnames(logtable)<-c('Name','Time','WhiteAmS','OpaqueAmS','WhiteBSA','OpaqueBSA','OWAms','OWBSA','BSAAmsW','BSAAmsO')

write_tsv(logtable,'~/dropbox/MattNaomiSAP/logtable.txt')


#############2 hr
##########################
hr2table<-summ%>%
		select(-SD,-N,-Min,-Max)%>%
		filter(Strain%in%c('WTo','WTw')&Media%in%c('AS','BSA.filter','NoNitro.withUri')&Time=='2hr')%>%
		spread(Strain,Mean)%>%
		nest(WTw,WTo)%>%
		spread(Media,data)%>%
		unnest()%>%
		mutate(BA.W=log2(WTw2/WTw),BA.O=log2(WTo2/WTo),NA.W=log2(WTw1/WTw),NA.O=log2(WTo1/WTo))%>%
		mutate(BA.W=sign(BA.W)*2^(abs(BA.W)),BA.O=sign(BA.O)*2^(abs(BA.O)),NA.W=sign(NA.W)*2^(abs(NA.W)),NA.O=sign(NA.O)*2^(abs(NA.O)))	

colnames(hr2table)<-c('Name','Time','WhiteAmS','OpaqueAmS','WhiteNoN','OpaqueNoN','WhiteBSA','OpaqueBSA','BSAAmsW','BSAAmsO','NoNAmsW','NoNAmsO')

write_tsv(hr2table,'~/dropbox/MattNaomiSAP/hr2table.txt')

####################################
#######  Stats - volcano  ##########
####################################


#############
pdf('~/dropbox/MattNaomiSAP/Stats.1.pdf',width=7,height=7,colormodel='cmyk')	
#############
lbs<-setNames(c('Opaque cells','Opaque cells','White cells','White cells'),unique(NstatsWT$panel))
bre<-c((1:10)*10^(-5),(1:10)*10^(-4),(1:10)*10^(-3),(1:10)*10^(-2),(1:10)*10^(-1))

temp<-filter(NstatsWT,abs(log2(B.A))>log2(5))
nudgesy<-rep(1,nrow(temp))
nudgesx<-rep(0,nrow(temp))
nudgesy[temp$Name=='SAP8'&temp$Time=='2hr']<-(-1)
nudgesy[temp$Name=='UGA4'&temp$Time=='2hr']<-(0)
nudgesx[temp$Name=='UGA4'&temp$Time=='2hr']<-(1)
nudgesy[temp$Name%in%c('OPT3','SAP1','GAP1')&temp$Time=='log']<-(0)
nudgesx[temp$Name%in%c('OPT3','SAP1','GAP1')&temp$Time=='log']<-(1)

ggplot(NstatsWT,aes(x=B.A,y=pval))+
		geom_vline(xintercept=c(5,1/5),color=2)+
		geom_hline(yintercept=c(max(filter(NstatsWT,FDR2.5)$pval)),color=2)+
		geom_point()+
		geom_text(aes(label=Name),temp,size=2,nudge_y=.15*nudgesy,nudge_x=.9*nudgesx)+
		geom_label(aes(label=anno,x=x,y=y),data.frame(x=rep(1/32,4),y=rep(10^-5,4),panel=c('WTo_log','WTw_log','WTo_2hr','WTw_2hr'),anno=c('Log growth','Log growth','2hr','2hr')),size=2.5)+
		scale_x_continuous('Fold change (BSA/AmS)',trans='log2',breaks=c(1/2,1/8,1/32,2^(seq(1,9,2))),labels=c('1/2','1/8','1/32','2','8','32','128','512'))+
		scale_y_continuous('P-value',trans=scales::trans_new('-log10',function(x){-log10(x)},function(x){1/10^(x)}),minor_breaks=bre,breaks=c(10^-1,10^-2,10^-3,10^-4,10^-5))+
		facet_wrap(~panel,labeller=labeller(panel=lbs))+
		
				theme_bw(8)+theme(axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color='grey80'),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8))

#############
dev.off()
#############


###################################
#######  Plots - all genes  #######
###################################
		
#############						
pdf('~/dropbox/MattNaomiSAP/Nanostring.1.pdf',width=7.8,height=4,colormodel='cmyk')						
#############

#############White verses opaque
##########################
temp1<-summ%>%
		select(-SD,-N)%>%
		filter(Strain%in%c('WTo','WTw')&Media%in%c('AS','BSA.filter'))%>%
		nest(Mean,Min,Max)%>%
		spread(Strain,data)%>%
		unnest()%>%
		mutate(Ratio=log2(Mean1/Mean),Media=factor(Media,,c('Ammonium Sulfate','BSA')))		
#Mean is white, Mean1 is opaque

#############White verses opaque - log
##########################
temp<-filter(temp1,Time=='log',abs(Ratio)>log2(5))
nudgesx<-sign(temp$Ratio)
nudgesy<-sign(temp$Ratio)
nudgesx[temp$Name=='SAP8']<-(-1)
nudgesy[temp$Name=='SAP99'&temp$Media=='BSA']<-(-1)


ggplot(filter(temp1,Time=='log'),aes(y=Mean1,x=Mean,ymax=Max1,ymin=Min1,xmax=Max,xmin=Min))+
		geom_errorbar(color='grey60',size=.25)+
		geom_errorbarh(color='grey60',size=.25)+
		geom_point(size=1.5)+
		geom_text(aes(label=Name),temp,size=2,nudge_y=0.15*nudgesy,nudge_x=-0.2*nudgesx)+
		geom_abline(intercept=0,slope=1)+
		scale_x_log10('Gene expression in White cells (normalized count)',limits=ran1)+
		scale_y_log10('Gene expression in Opaque cells (normalized count)',limits=ran1)+
		facet_wrap(~Media)+
		annotate('label',x=2,y=250000,label='Log growth',size=2.5)+
		theme_bw(8)+theme(axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8))

#############White verses opaque - 2hr
##########################

temp<-filter(temp1,Time=='2hr',abs(Ratio)>log2(5))
nudgesx<-sign(temp$Ratio)
nudgesy<-sign(temp$Ratio)
nudgesx[temp$Name=='SAP3'&temp$Media=='BSA']<-(-1)
nudgesy[temp$Name=='SAP1'&temp$Media=='BSA']<-(-1)
nudgesx[temp$Name=='OPT5'&temp$Media=='BSA']<-(-1)
nudgesy[temp$Name=='OPT5'&temp$Media=='BSA']<-(-1)
nudgesy[temp$Name=='OPT4'&temp$Media=='BSA']<-(-1)
nudgesx[temp$Name=='OP4'&temp$Media=='BSA']<-(-1)
nudgesy[temp$Name=='PTR2'&temp$Media=='BSA']<-(-1)
nudgesy[temp$Name=='PTR22'&temp$Media=='Ammonium Sulfate']<-(1)
nudgesy[temp$Name=='OPT7'&temp$Media=='Ammonium Sulfate']<-(1)

ggplot(filter(temp1,Time=='2hr'),aes(y=Mean1,x=Mean,ymax=Max1,ymin=Min1,xmax=Max,xmin=Min))+
		geom_errorbar(color='grey60',size=.25)+
		geom_errorbarh(color='grey60',size=.25)+
		geom_point(size=1.5)+
		geom_text(aes(label=Name),temp,size=2,nudge_y=0.15*nudgesy,nudge_x=-0.2*nudgesx)+
		geom_abline(intercept=0,slope=1)+
		scale_x_log10('Gene expression in White cells (normalized count)',limits=ran1)+
		scale_y_log10('Gene expression in Opaque cells (normalized count)',limits=ran1)+
		facet_wrap(~Media)+
		annotate('label',x=2,y=250000,label='2 hours',size=2.5)+
		theme_bw(8)+theme(axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8))

#############
dev.off()
#############

#############
pdf('~/dropbox/MattNaomiSAP/Nanostring.2.pdf',width=4,height=4,colormodel='cmyk')						
#############

#############AS verses BSA
##########################
temp2<-summ%>%
		select(-SD,-N)%>%
		filter(Strain%in%c('WTo','WTw')&Media%in%c('AS','BSA.filter'))%>%
		nest(Mean,Min,Max)%>%
		spread(Media,data)%>%
		unnest()%>%
		mutate(Ratio=log2(Mean1/Mean))				
#Mean is AS, Mean1 is BSA		

#############AS verses BSA -log
##########################
temp<-filter(temp2,Time=='log',abs(Ratio)>log2(5))
nudgesx<-sign(temp$Ratio)
nudgesy<-sign(temp$Ratio)
nudgesy[temp$Name=='OPT2'&temp$Strain=='WTw']<-(-1)
nudgesx[temp$Name=='OPT2'&temp$Strain=='WTo']<-(-1)
nudgesy[temp$Name=='PTR2'&temp$Strain=='WTo']<-(-1)
nudgesy[temp$Name=='SAP3'&temp$Strain=='WTo']<-(-1)
nudgesx[temp$Name=='SAP3'&temp$Strain=='WTo']<-(-1)
nudgesx[temp$Name=='SAP2'&temp$Strain=='WTo']<-(-1)

nudgesy[temp$Name=='OPT3'&temp$Strain=='WTo']<-(-1)
nudgesy[temp$Name=='WH11'&temp$Strain=='WTw']<-(-1)

ggplot(filter(temp2,Time=='log'),aes(y=Mean1,x=Mean,ymax=Max1,ymin=Min1,xmax=Max,xmin=Min))+
		geom_errorbar(color='grey60',size=.25)+
		geom_errorbarh(color='grey60',size=.25)+
		geom_point(size=1.5)+
		geom_text(aes(label=Name,color=Strain),temp,size=2,nudge_y=0.15*nudgesy,nudge_x=-0.2*nudgesx)+
		geom_abline(intercept=0,slope=1)+
		scale_x_log10('Gene expression in Ammonium Sulfate (normalized count)',limits=ran1)+
		scale_y_log10('Gene expression in BSA (normalized count)',limits=ran1)+
		scale_color_manual('Cell type',labels=c('White','Opaque'),values=c('#EB2627','#3A54A3'))+
		
		annotate('label',x=2,y=250000,label='Log growth',size=2.5)+
		theme_bw(8)+theme(legend.background=element_rect(color=1),legend.position=c(0.82,0.1),legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,legend.text=element_text(size=8),axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8))
		
#############AS verses BSA -2hr
##########################
		
temp<-filter(temp2,Time=='2hr',abs(Ratio)>log2(5))
nudgesx<-sign(temp$Ratio)
nudgesy<-sign(temp$Ratio)
nudgesy[temp$Name=='GAP1'&temp$Strain=='WTw']<-(1)
nudgesx[temp$Name=='OPT2'&temp$Strain=='WTw']<-(-1)
nudgesy[temp$Name=='UGA4'&temp$Strain=='WTw']<-(-1)
nudgesy[temp$Name=='AOX2'&temp$Strain=='WTw']<-(-1)

ggplot(filter(temp2,Time=='2hr'),aes(y=Mean1,x=Mean,ymax=Max1,ymin=Min1,xmax=Max,xmin=Min))+
		geom_errorbar(color='grey60',size=.25)+
		geom_errorbarh(color='grey60',size=.25)+
		geom_point(size=1.5)+
		geom_text(aes(label=Name,color=Strain),temp,size=2,nudge_y=0.15*nudgesy,nudge_x=-0.2*nudgesx)+
		geom_abline(intercept=0,slope=1)+
		scale_x_log10('Gene expression in Ammonium Sulfate (normalized count)',limits=ran1)+
		scale_y_log10('Gene expression in BSA (normalized count)',limits=ran1)+
		scale_color_manual('Cell type',labels=c('White','Opaque'),values=c('#EB2627','#3A54A3'))+
		
		annotate('label',x=2,y=250000,label='2 hours',size=2.5)+
		theme_bw(8)+theme(legend.background=element_rect(color=1),legend.position=c(0.82,0.1),legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,legend.text=element_text(size=8),axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8))		

#############
dev.off()
#############

#############
# #pdf('~/dropbox/MattNaomiSAP/Nanostring.3.pdf',width=6,height=4,colormodel='cmyk')
# #############

# #############AS verses BSA/hemoglobin/myoglobin -log - opaque
# ##########################
# temp3<-summ%>%
		# select(-SD,-N)%>%
		# filter(Time=='log'&Strain%in%c('WTo')&Media%in%c('AS','BSA.filter','BSA.dialyzed','Hemoglobin.dialyzed','Myoglobin.dialyzed'))%>%
		# nest(Mean,Min,Max)%>%
		# spread(Media,data)%>%
		# gather('Protein',,c('BSA.filter','BSA.dialyzed','Hemoglobin.dialyzed','Myoglobin.dialyzed'))%>%
		# unnest()	
# #Mean is AS, Mean1 is other proteins

# #temp<-filter(temp3,Name%in%c('SAP1','SAP2','SAP3','SAP8','OPT1','OPT2','OPT5','OPT7','PTR2'))
# temp<-filter(temp3,Name%in%c('SAP1','SAP2','SAP3','SAP8','OPT1','OPT2','OPT5','OPT7','PTR2','UGA4'))
		
# ggplot(temp3,aes(x=Mean,y=Mean1,xmin=Min,xmax=Max,ymin=Min1,ymax=Max1,shape=Protein))+
		# geom_errorbar(color='grey60',size=.25)+
		# geom_errorbarh(color='grey60',size=.25)+
		# geom_point(size=1.5,fill=1)+
		# geom_point(aes(fill=Name),temp,size=2.5,color=1)+
		# geom_abline(intercept=0,slope=1)+
		# scale_fill_manual('Genes',values=colors[levels(factor(temp$Name))])+
		# scale_x_log10('Gene expression in Ammonium Sulfate (normalized count)',limits=ran1)+
		# scale_y_log10('Gene expression in protein (normalized count)',limits=ran1)+
		# scale_shape_manual('Nitrogen source',values=c(22,21,25,24),breaks=c('BSA.filter','BSA.dialyzed','Myoglobin.dialyzed','Hemoglobin.dialyzed'),labels=c('BSA','BSA (dialyzed)','Myoglobin (dialyzed)','Hemoglobin (dialyzed)'))+
		
		# annotate('label',x=2,y=250000,label='Log growth',size=2.5)+
		# ggtitle('Opaque cells')+
		# theme_bw(8)+theme(legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,legend.text=element_text(size=8),axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),plot.title=element_text(hjust=0.5))+
		# guides(fill=guide_legend(override.aes=list(shape=21),order=1))


# #############AS verses BSA/hemoglobin/myoglobin -log - White
# ##########################
# temp4<-summ%>%
		# select(-SD,-N)%>%
		# filter(Time=='log'&Strain%in%c('WTw')&Media%in%c('AS','BSA.filter','Hemoglobin.dialyzed','Myoglobin.dialyzed'))%>%
		# nest(Mean,Min,Max)%>%
		# spread(Media,data)%>%
		# gather('Protein',,c('BSA.filter','Hemoglobin.dialyzed','Myoglobin.dialyzed'))%>%
		# unnest()

# #temp<-filter(temp4,Name%in%c('SAP1','SAP2','SAP3','SAP8','OPT1','OPT2','OPT5','OPT7','PTR2'))
# temp<-filter(temp4,Name%in%c('SAP1','SAP2','SAP3','SAP8','OPT1','OPT2','OPT5','OPT7','PTR2','UGA4'))
		
# ggplot(temp4,aes(x=Mean,y=Mean1,xmin=Min,xmax=Max,ymin=Min1,ymax=Max1,shape=Protein))+
		# geom_errorbar(color='grey60',size=.25)+
		# geom_errorbarh(color='grey60',size=.25)+
		# geom_point(size=1.5,fill=1)+
		# geom_point(aes(fill=Name),temp,size=2.5,color=1)+
		# geom_abline(intercept=0,slope=1)+
		# scale_fill_manual('Genes',values=colors[levels(factor(temp$Name))])+
		# scale_x_log10('Gene expression in Ammonium Sulfate (normalized count)',limits=ran1)+
		# scale_y_log10('Gene expression in protein (normalized count)',limits=ran1)+
		# scale_shape_manual('Nitrogen source',values=c(21,25,24),breaks=c('BSA.filter','Myoglobin.dialyzed','Hemoglobin.dialyzed'),labels=c('BSA','Myoglobin (dialyzed)','Hemoglobin (dialyzed)'))+
		
		# annotate('label',x=2,y=250000,label='Log growth',size=2.5)+
		# ggtitle('White cells')+
		# theme_bw(8)+theme(legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,legend.text=element_text(size=8),axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),plot.title=element_text(hjust=0.5))+
		# guides(fill=guide_legend(override.aes=list(shape=21),order=1))



# #############AS verses BSA/no nitro -2hr  - opaque
# ##########################
# temp5<-summ%>%
		# select(-SD,-N)%>%
		# filter(Time=='2hr'&Strain%in%c('WTo')&Media%in%c('AS','BSA.filter',"NoNitro.withUri","BSA.filter.noUri","NoNitro.noUri"))%>%
		# nest(Mean,Min,Max)%>%
		# spread(Media,data)%>%
		# gather('Protein',,c('BSA.filter',"NoNitro.withUri","BSA.filter.noUri","NoNitro.noUri"))%>%
		# unnest()	
# #Mean is AS, Mean1 is other proteins

# #temp<-temp5%>%
# #		filter(Name%in%c('SAP3','UGA4','OPT5','OPT2','SAP2','DAL5','MEP2','SAP8','PTR2','OPT7'))
# temp<-filter(temp5,Name%in%c('SAP1','SAP2','SAP3','SAP8','OPT1','OPT2','OPT5','OPT7','PTR2','UGA4'))
		
		
# ggplot(temp5,aes(x=Mean,y=Mean1,xmin=Min,xmax=Max,ymin=Min1,ymax=Max1,shape=Protein))+
		# geom_errorbar(color='grey60',size=.25)+
		# geom_errorbarh(color='grey60',size=.25)+
		# geom_point(size=1.5,fill=1)+
		# geom_point(aes(fill=Name),temp,size=2.5,color=1)+
		# geom_abline(intercept=0,slope=1)+
		# scale_fill_manual('Genes',values=colors[levels(factor(temp$Name))])+
		# scale_x_log10('Gene expression in Ammonium Sulfate (normalized count)',limits=ran1)+
		# scale_y_log10('Gene expression in BSA or Nitrogen depletion (normalized count)',limits=ran1)+
		# scale_shape_manual('Nitrogen source',values=c(21,22,24,25),labels=c('BSA','BSA (No uridine)','Nitrogen depletion','Nitrogen depletion\n(No uridine)'))+
		
		# annotate('label',x=2,y=250000,label='2 hours',size=2.5)+
		# ggtitle('Opaque cells')+
		# theme_bw(8)+theme(legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,legend.text=element_text(size=8),axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),plot.title=element_text(hjust=0.5))+
		# guides(fill=guide_legend(override.aes=list(shape=21),order=1))


# #############AS verses BSA/no nitro -2hr - White
# ##########################
# temp6<-summ%>%
		# select(-SD,-N)%>%
		# filter(Time=='2hr'&Strain%in%c('WTw')&Media%in%c('AS','BSA.filter',"NoNitro.withUri"))%>%
		# nest(Mean,Min,Max)%>%
		# spread(Media,data)%>%
		# gather('Protein',,c('BSA.filter',"NoNitro.withUri"))%>%
		# unnest()	
# #Mean is AS, Mean1 is other proteins

# #temp<-temp6%>%
# #		filter(Name%in%c('UGA4','OPT2','SAP2','OPT7'))
# temp<-filter(temp6,Name%in%c('SAP1','SAP2','SAP3','SAP8','OPT1','OPT2','OPT5','OPT7','PTR2','UGA4'))
		
# ggplot(temp6,aes(x=Mean,y=Mean1,xmin=Min,xmax=Max,ymin=Min1,ymax=Max1,shape=Protein))+
		# geom_errorbar(color='grey60',size=.25)+
		# geom_errorbarh(color='grey60',size=.25)+
		# geom_point(size=1.5,fill=1)+
		# geom_point(aes(fill=Name),temp,size=2.5,color=1)+
		# geom_abline(intercept=0,slope=1)+
		# scale_fill_manual('Genes',values=colors[levels(factor(temp$Name))])+
		# scale_x_log10('Gene expression in Ammonium Sulfate (normalized count)',limits=ran1)+
		# scale_y_log10('Gene expression in BSA or Nitrogen depletion (normalized count)',limits=ran1)+
		# scale_shape_manual('Nitrogen source',values=c(21,24),labels=c('BSA','Nitrogen depletion'))+
		
		# annotate('label',x=2,y=250000,label='2 hours',size=2.5)+
		# ggtitle('White cells')+
		# theme_bw(8)+theme(legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,legend.text=element_text(size=8),axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),plot.title=element_text(hjust=0.5))+
		# guides(fill=guide_legend(override.aes=list(shape=21),order=1))

##############################################################################
##############################################################################

pdf('~/dropbox/MattNaomiSAP/Nanostring.3.pdf',width=7.8,height=4,colormodel='cmyk')
#############

#############AS verses BSA/hemoglobin/myoglobin -log - opaque
##########################
temp34<-rbind(temp3,temp4)
temp34$Strain<-factor(temp34$Strain,c('WTo','WTw'),c('Opaque Cells','White Cells'))

temp<-filter(temp34,Name%in%c('SAP1','SAP2','SAP3','SAP8','OPT1','OPT2','OPT5','OPT7','PTR2','UGA4'))
		
ggplot(temp34,aes(x=Mean,y=Mean1,xmin=Min,xmax=Max,ymin=Min1,ymax=Max1,shape=Protein))+
		geom_errorbar(color='grey60',size=.25)+
		geom_errorbarh(color='grey60',size=.25)+
		geom_point(size=1.5,fill=1)+
		geom_point(aes(fill=Name),temp,size=2.5,color=1)+
		geom_abline(intercept=0,slope=1)+
		scale_fill_manual('Genes',values=colors[levels(factor(temp$Name))])+
		scale_x_log10('Gene expression in Ammonium Sulfate (normalized count)',limits=ran1)+
		scale_y_log10('Gene expression in protein (normalized count)',limits=ran1)+
		scale_shape_manual('Nitrogen source',values=c(22,21,25,24),breaks=c('BSA.filter','BSA.dialyzed','Myoglobin.dialyzed','Hemoglobin.dialyzed'),labels=c('BSA','BSA (dialyzed)','Myoglobin (dialyzed)','Hemoglobin (dialyzed)'))+
		facet_wrap(~Strain)+
		annotate('label',x=2.1,y=250000,label='Log growth',size=2.5)+
		theme_bw(8)+theme(legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,legend.text=element_text(size=8),axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),plot.title=element_text(hjust=0.5),legend.margin=margin(r=0,l=-8))+
		guides(fill=guide_legend(override.aes=list(shape=21),order=1))

#########################################################

temp56<-rbind(temp5,temp6)
temp56$Strain<-factor(temp56$Strain,c('WTo','WTw'),c('Opaque Cells','White Cells'))

temp<-filter(temp56,Name%in%c('SAP1','SAP2','SAP3','SAP8','OPT1','OPT2','OPT5','OPT7','PTR2','UGA4'))
		
		
ggplot(temp56,aes(x=Mean,y=Mean1,xmin=Min,xmax=Max,ymin=Min1,ymax=Max1,shape=Protein))+
		geom_errorbar(color='grey60',size=.25)+
		geom_errorbarh(color='grey60',size=.25)+
		geom_point(size=1.5,fill=1)+
		geom_point(aes(fill=Name),temp,size=2.5,color=1)+
		geom_abline(intercept=0,slope=1)+
		scale_fill_manual('Genes',values=colors[levels(factor(temp$Name))])+
		scale_x_log10('Gene expression in Ammonium Sulfate (normalized count)',limits=ran1)+
		scale_y_log10('Gene expression in BSA or Nitrogen depletion (normalized count)',limits=ran1)+
		scale_shape_manual('Nitrogen source',values=c(21,22,24,25),labels=c('BSA','BSA (No uridine)','Nitrogen depletion','Nitrogen depletion       \n(No uridine)'))+
		facet_wrap(~Strain)+
		annotate('label',x=2,y=250000,label='2 hours',size=2.5)+
		theme_bw(8)+theme(legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,legend.text=element_text(size=8),axis.text=element_text(size=8,color=1),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),plot.title=element_text(hjust=0.5) ,legend.margin=margin(r=0,l=-8))+
		guides(fill=guide_legend(override.aes=list(shape=21),order=1))


#############
dev.off()
#############


########################################
#######  Plots - Subset - Color  #######
########################################

#############
pdf('~/dropbox/MattNaomiSAP/Nanostring.4.pdf',width=7.8,height=3,colormodel='cmyk')	
#############

#############Regulators - opaque
##########################
#156-efg1,170-wor3,106-csr1

#New figure after reviewers

temp7<-summ%>%
		select(-SD,-N)%>%
		filter(Time=='log'&	((Strain%in%c('WTo','TF156o','TF170o','TF106o','Ssy1o','TF052o','EFG1del.Stp1act.o','TF106.Stp1act.o','WOR3del.Stp1act.o')&Media%in%c('AS','BSA.filter'))|
							(Strain=='Stp1.To'&Media%in%c('AS','BSA'))))%>%
		filter(Name%in%genes1)%>%
		mutate(cond=factor(paste(Strain,Media),
								c('WTo AS','TF156o AS','TF170o AS','TF106o AS', 'TF052o AS', 'WTo BSA.filter',  'TF156o BSA.filter','TF170o BSA.filter','TF106o BSA.filter', 'TF052o BSA.filter','Ssy1o BSA.filter','Stp1.To AS', 'Stp1.To BSA','EFG1del.Stp1act.o AS','WOR3del.Stp1act.o AS','TF106.Stp1act.o AS', 'EFG1del.Stp1act.o BSA.filter','WOR3del.Stp1act.o BSA.filter','TF106.Stp1act.o BSA.filter')))								
								
ggplot(temp7,aes(x=cond,y=factor(Name,genes1[length(genes1):1]),fill=Mean))+
	geom_raster()+
	geom_text(aes(x=1:19,y=12,label=1:19,fill=NA),data.frame(),size=2.5)+
	geom_hline(yintercept=c(3.5,4.5,6.5),color='white',size=1)+
	geom_vline(xintercept=c(1.5,5.5,6.5,11.5,13.5,16.5),color='white',size=c(1,2,1,2,1,1))+
	scale_y_discrete('\n\n\n\n\n\nShared\nprotein response\n\n\n\n\n\nOpaque specific\nprotein response\n\nWhite specific\nprotein response\n\n\nBasal Opaque\nprogram\n\n\n\n',expand=expand_scale(add=c(0,1.3)))+
	scale_x_discrete('',labels=c(
							expression('AmS; Wildtype'),  
							expression(paste('AmS; ',Delta,italic('efg1'))),
							expression(paste('AmS; ',Delta,italic('wor3'))), 
							expression(paste('AmS; ',Delta,italic('csr1'))),
							expression(paste('AmS; ',Delta,italic('stp1'))), 
							expression('BSA; Wildtype'),
							expression(paste('BSA; ',Delta,italic('efg1'))),
							expression(paste('BSA; ',Delta,italic('wor3'))),
							expression(paste('BSA; ',Delta,italic('csr1'))),
							expression(paste('BSA; ',Delta,italic('stp1'))), 
							expression(paste('BSA; ',Delta,italic('ssy1'))),
							expression(paste('AmS; ',italic('STP1'),Delta,'2-61')),
							expression(paste('BSA; ',italic('STP1'),Delta,'2-61')),
							expression(paste('AmS; ',Delta,italic('efg1'),',',italic(' STP1'),Delta,'2-61')),
							expression(paste('AmS; ',Delta,italic('wor3'),',',italic(' STP1'),Delta,'2-61')),
							expression(paste('AmS; ',Delta,italic('csr1'),',',italic(' STP1'),Delta,'2-61')),
							expression(paste('BSA; ',Delta,italic('efg1'),',',italic(' STP1'),Delta,'2-61')),
							expression(paste('BSA; ',Delta,italic('wor3'),',',italic(' STP1'),Delta,'2-61')),
							expression(paste('BSA; ',Delta,italic('csr1'),',',italic(' STP1'),Delta,'2-61'))))+
	scale_fill_gradientn('Gene expression\n(Normalized count)',colors=c('#ffffcc','#41b6c4','#0c2c84'),trans='log10',breaks=c(100,1000,10000,100000),limits=c(10,200000),oob=scales::squish)+
	ggtitle('Opaque cells - Log growth')+
	theme(text=element_text(size=8,color=1),rect=element_blank(),axis.text=element_text(color=1),axis.text.x=element_text(angle=30,vjust=1,hjust=1), plot.title=element_text(hjust=0.5,size=8),legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,axis.title.y=element_text(angle=0,size=6,vjust=0.5),legend.margin=margin(r=0,l=-8))	


dev.off()


#old figures
# temp7<-summ%>%
		# select(-SD,-N)%>%
		# filter(Time=='log'&	((Strain%in%c('WTo','TF156o','TF170o','TF106o')&Media%in%c('AS','BSA.filter'))|
							# (Strain=='TF052o'&Media=='AS')|
							# (Strain=='Stp1.To'&Media%in%c('AS'))))%>%
		# filter(Name%in%genes1)%>%
		# mutate(cond=factor(paste(Strain,Media),
								# c('WTo AS', 'TF156o AS','TF170o AS','TF106o AS', 'TF052o AS', 'WTo BSA.filter', 'Stp1.To AS', 'TF156o BSA.filter','TF170o BSA.filter','TF106o BSA.filter')))
								
								
# ggplot(temp7,aes(x=cond,y=factor(Name,genes1[length(genes1):1]),fill=Mean))+
	# geom_raster()+
	# geom_hline(yintercept=c(3.5,4.5,6.5),color='white',size=1)+
	# geom_vline(xintercept=c(1.5,5.5,7.5),color='white',size=c(1,2,1))+
	# scale_y_discrete('\n\n\n\n\n\nShared\nprotein response\n\n\n\n\n\nOpaque specific\nprotein response\n\nWhite specific\nprotein response\n\n\nBasal Opaque\nprogram\n\n\n\n')+
	# scale_x_discrete('',labels=c(
							# expression('AmS; Wildtype'),
							# expression(paste('AmS; ',Delta,italic('efg1'))),
							# expression(paste('AmS; ',Delta,italic('wor3'))), 
							# expression(paste('AmS; ',Delta,italic('csr1'))), 
							# expression(paste('AmS; ',Delta,italic('stp1'))), 
							# expression('BSA; Wildtype'),
							# expression(paste('AmS; ',italic('STP1'),'active')),
							# expression(paste('BSA; ',Delta,italic('efg1'))),
							# expression(paste('BSA; ',Delta,italic('wor3'))),
							# expression(paste('BSA; ',Delta,italic('csr1')))))+
	# scale_fill_gradientn('Gene expression\n(Normalized count)',colors=c('#ffffcc','#41b6c4','#0c2c84'),trans='log10',breaks=c(100,1000,10000,100000),limits=c(10,200000),oob=scales::squish)+
	# ggtitle('Opaque cells - Log growth')+
	# theme(text=element_text(size=8,color=1),rect=element_blank(),axis.text=element_text(color=1),axis.text.x=element_text(angle=30,vjust=1,hjust=1), plot.title=element_text(hjust=0.5,size=8),legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,axis.title.y=element_text(angle=0,size=6,vjust=0.5))	


# #############Regulators - opaque - STP1 interactions
# ##########################

# temp8<-summ%>%
		# select(-SD,-N)%>%
		# filter(Time=='log'&	((Strain%in%c('WTo','EFG1del.Stp1act.o','TF106.Stp1act.o','WOR3del.Stp1act.o')&Media%in%c('AS','BSA.filter'))|
							# (Strain=='Stp1.To'&Media%in%c('AS','BSA'))))%>%
		# filter(Name%in%genes1)%>%
		# mutate(cond=factor(paste(Strain,Media),
								# c('WTo AS', 'WTo BSA.filter', 'Stp1.To AS', 'Stp1.To BSA','EFG1del.Stp1act.o AS','WOR3del.Stp1act.o AS','TF106.Stp1act.o AS', 'EFG1del.Stp1act.o BSA.filter','WOR3del.Stp1act.o BSA.filter','TF106.Stp1act.o BSA.filter')))
								
# ggplot(temp8,aes(x=cond,y=factor(Name,genes1[length(genes1):1]),fill=Mean))+
	# geom_raster()+
	# geom_hline(yintercept=c(3.5,4.5,6.5),color='white',size=1)+
	# geom_vline(xintercept=c(2.5,4.5,7.5),color='white',size=1)+
	# scale_y_discrete('\n\n\n\n\n\nShared\nprotein response\n\n\n\n\n\nOpaque specific\nprotein response\n\nWhite specific\nprotein response\n\n\nBasal Opaque\nprogram\n\n\n\n')+
	# scale_x_discrete('',labels=c(
							# expression('AmS; Wildtype'), 
							# expression('BSA; Wildtype'),
							# expression(paste('AmS; ',italic('STP1'),'active')),
							# expression(paste('BSA; ',italic('STP1'),'active')),
							# expression(paste('AmS; ',Delta,italic('efg1'),',',italic(' STP1'),'active')),
							# expression(paste('AmS; ',Delta,italic('wor3'),',',italic(' STP1'),'active')),
							# expression(paste('AmS; ',Delta,italic('csr1'),',',italic(' STP1'),'active')),
							# expression(paste('BSA; ',Delta,italic('efg1'),',',italic(' STP1'),'active')),
							# expression(paste('BSA; ',Delta,italic('wor3'),',',italic(' STP1'),'active')),
							# expression(paste('BSA; ',Delta,italic('csr1'),',',italic(' STP1'),'active'))))+
	# scale_fill_gradientn('Gene expression\n(Normalized count)',colors=c('#ffffcc','#41b6c4','#0c2c84'),trans='log10',breaks=c(100,1000,10000,100000),limits=c(10,200000),oob=scales::squish)+
	# ggtitle('Opaque cells - Log growth')+
	# theme_bw()+
	# theme(text=element_text(size=8,color=1),rect=element_blank(),axis.text=element_text(color=1),axis.text.x=element_text(angle=30,vjust=1,hjust=1), plot.title=element_text(hjust=0.5,size=8),legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,axis.title.y=element_text(angle=0,size=6,vjust=0.5))	

# #############
# dev.off()
# #############


#############
pdf('~/dropbox/MattNaomiSAP/Nanostring.5.pdf',width=6,height=3,colormodel='cmyk')	
#############

#############Regulators - white
##########################

temp9<-summ%>%
		select(-SD,-N)%>%
		filter(Time=='log'&	((Strain%in%c('WTw','EFG1del.Stp1act.w','WOR3del.Stp1act.w')&Media%in%c('AS','BSA.filter'))|
							(Strain%in%c('TF106.Stp1act.w','TF052w','TF156w','TF170w','TF106w')&Media=='AS')|
							(Strain=='Stp1.Tw'&Media%in%c('AS','BSA'))))%>%
		filter(Name%in%genes1)%>%
		mutate(cond=factor(paste(Strain,Media),
								c('WTw AS', 'TF156w AS','TF170w AS','TF106w AS', 'TF052w AS', 'WTw BSA.filter', 'Stp1.Tw AS', 'Stp1.Tw BSA','EFG1del.Stp1act.w AS','WOR3del.Stp1act.w AS','TF106.Stp1act.w AS', 'EFG1del.Stp1act.w BSA.filter','WOR3del.Stp1act.w BSA.filter')))

	
ggplot(temp9,aes(x=cond,y=factor(Name,genes1[length(genes1):1]),fill=Mean))+
	geom_raster()+
	geom_text(aes(x=1:13,y=12,label=1:13,fill=NA),data.frame(),size=2.5)+
	geom_hline(yintercept=c(3.5,4.5,6.5),color='white',size=1)+
	geom_vline(xintercept=c(1.5,5.5,6.5,8.5,11.5),color='white',size=c(1,2,1,1,1))+
	scale_y_discrete('\n\n\n\n\n\nShared\nprotein response\n\n\n\n\n\nOpaque specific\nprotein response\n\nWhite specific\nprotein response\n\n\nBasal Opaque\nprogram\n\n\n\n',expand=expand_scale(add=c(0,1.3)))+
	scale_x_discrete('',labels=c(
					expression('AmS; Wildtype'), 
					expression(paste('AmS; ',Delta,italic('efg1'))),
					expression(paste('AmS; ',Delta,italic('wor3'))), 
					expression(paste('AmS; ',Delta,italic('csr1'))), 
					expression(paste('AmS; ',Delta,italic('stp1'))), 
					expression('BSA; Wildtype'),
					expression(paste('AmS; ',italic('STP1'),Delta,'2-61')),
					expression(paste('BSA; ',italic('STP1'),Delta,'2-61')),
					expression(paste('AmS; ',Delta,italic('efg1'),',',italic(' STP1'),Delta,'2-61')),
					expression(paste('AmS; ',Delta,italic('wor3'),',',italic(' STP1'),Delta,'2-61')),
					expression(paste('AmS; ',Delta,italic('csr1'),',',italic(' STP1'),Delta,'2-61')),
					expression(paste('BSA; ',Delta,italic('efg1'),',',italic(' STP1'),Delta,'2-61')),
					expression(paste('BSA; ',Delta,italic('wor3'),',',italic(' STP1'),Delta,'2-61'))
					))+
	scale_fill_gradientn('Gene expression\n(Normalized count)',colors=c('#ffffcc','#41b6c4','#0c2c84'),trans='log10',breaks=c(100,1000,10000,100000),limits=c(10,200000),oob=scales::squish)+
	ggtitle('White cells - Log growth')+
	theme(text=element_text(size=8,color=1),rect=element_blank(),axis.text=element_text(color=1),axis.text.x=element_text(angle=30,vjust=1,hjust=1), plot.title=element_text(hjust=0.5,size=8),legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,axis.title.y=element_text(angle=0,size=6,vjust=0.5),legend.margin=margin(r=0,l=-8))	
	
#############
dev.off()
#############

#############
pdf('~/dropbox/MattNaomiSAP/Nanostring.6.pdf',width=7.8,height=3,colormodel='cmyk')	
#############

#############STP1/SSY/SAP response
##########################

# temp10<-summ%>%
		# select(-SD,-N)%>%
		# filter(				((Strain%in%c('WTo')&Media%in%c('AS','BSA.filter','NoNitro.withUri'))|
							# (Strain%in%c('Ssy1o')&Media=='BSA.filter')|
							# (Strain%in%c('Sap.1.2.3.8.99o')&Media=='BSA.filter')|
							# (Strain%in%c('TF052o')&Media%in%c('AS','BSA.filter','NoNitro.withUri','BSA'))))%>%
		# filter(Name%in%genes1)%>%
		# mutate(cond=factor(paste(Strain,Media,Time),c('WTo AS log','TF052o AS log', 'WTo BSA.filter log', 'TF052o BSA.filter log' ,'Ssy1o BSA.filter log', 'WTo AS 2hr','TF052o AS 2hr', 'WTo BSA.filter 2hr','WTo NoNitro.withUri 2hr', 'TF052o BSA 2hr', 'TF052o NoNitro.withUri 2hr','Sap.1.2.3.8.99o BSA.filter 2hr')))

temp10<-summ%>%
		select(-SD,-N)%>%
		filter(Time=='2hr'&		((Strain%in%c('WTo')&Media%in%c('AS','BSA.filter','NoNitro.withUri'))|
							(Strain%in%c('Sap.1.2.3.8.99o')&Media=='BSA.filter')|
							(Strain%in%c('TF052o')&Media%in%c('AS','BSA.filter','NoNitro.withUri','BSA'))))%>%
		filter(Name%in%genes2)%>%
		mutate(cond=factor(paste(Strain,Media,Time),c( 'WTo AS 2hr','TF052o AS 2hr', 'WTo BSA.filter 2hr','WTo NoNitro.withUri 2hr', 'TF052o BSA 2hr', 'TF052o NoNitro.withUri 2hr','Sap.1.2.3.8.99o BSA.filter 2hr')))		
		
ggplot(temp10,aes(x=cond,y=factor(Name,genes2[length(genes2):1]),fill=Mean))+
	geom_raster()+
	geom_text(aes(x=1:7,y=13,label=1:7,fill=NA),data.frame(),size=2.5)+
	geom_hline(yintercept=c(1.5,4.5,5.5,7.5),color='white',size=1)+
	#geom_hline(yintercept=c(3.5,6.5,9.5,10.5,12.5),color='white',size=c(1,2,1,1,1))+
	geom_vline(xintercept=c(1.5,2.5,4.5,6.5),color='white',size=c(1,2,1,1))+
	scale_y_discrete('\n\n\n\n\n\n\n\n\n\n\n\nShared\nprotein response\n\n\n\n\nOpaque specific\nprotein response\n\nWhite specific\nprotein response\n\nBasal Opaque\nprogram\n\nNitrogen\nstarvation response\n\n\n\n\n\n\n\n\n',expand=expand_scale(add=c(0,1.3)))+
	#Orginal scale for this figure, if you want to change back need to change to genes2 in ggplot call also need a height of 4
	#scale_y_discrete('\n\n\n\n\n\nShared\nprotein response\n\n\n\n\n\nOpaque specific\nprotein response\n\nWhite specific\nprotein response\n\n\nBasal Opaque\nprogram\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')+
	scale_x_discrete('',labels=c(
					#expression('AmS; Wildtype; log'),  
					#expression(paste('AmS; ',Delta,italic('stp1'),'; log')), 
					#expression('BSA; Wildtype; log'),
					#expression(paste('BSA; ',Delta,italic('stp1'),'; log')),
					#expression(paste('BSA; ',Delta,italic('ssy1'),'; log')),
					expression('AmS; Wildtype; 2hr'),  
					expression(paste('AmS; ',Delta,italic('stp1'),'; 2hr')), 
					expression('BSA; Wildtype; 2hr'),
					expression('Nitrogen depletion; Wildtype; 2hr'),
					expression(paste('BSA; ',Delta,italic('stp1'),'; 2hr')),
					expression(paste('Nitrogen depletion; ',Delta,italic('stp1'),'; 2hr')),
					expression(paste('BSA; ',Delta,italic('sap1/2/3/8/99'),'; 2hr'))
					))+
	scale_fill_gradientn('Gene expression\n(Normalized count)',colors=c('#ffffcc','#41b6c4','#0c2c84'),trans='log10',breaks=c(100,1000,10000,100000),limits=c(10,200000),oob=scales::squish)+
annotate('text',label='*',color='red',x=c(6.6,6.6,6.6,6.6,6.6),y=c(3,4,7,11,12))+	#annotate('text',label='*',color='red',x=c(1.6,3.6,6.6,9.6,10.6,4.6,11.6,11.6,11.6,11.6,11.6),y=c(2,2,2,2,2,3,8,9,12,16,17))+
	ggtitle('Opaque cells - 2 hours')+
	theme(text=element_text(size=8,color=1),rect=element_blank(),axis.text=element_text(color=1),axis.text.x=element_text(angle=30,vjust=1,hjust=1), plot.title=element_text(hjust=0.5,size=8),legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,axis.title.y=element_text(angle=0,size=6,vjust=0.5))

#############
dev.off()
#############

######################################
#######  Some addtional plots  #######
######################################

#############Sap deletions on AmS 
##########################
#############
pdf('~/dropbox/MattNaomiSAP/Nanostring.7.pdf',width=4,height=4,colormodel='cmyk')						
#############
	

temp11<-summ%>%
		select(-SD,-N)%>%
		filter(Time=='log'&Strain%in%c('Sap.1.2.3.8.99o','WTo')&Media%in%c('AS'))%>%
		nest(Mean,Min,Max)%>%
		spread(Strain,data)%>%
		unnest()%>%
		mutate(Ratio=log2(Mean1/Mean))

temp<-filter(temp11,grepl('SAP',Name))
nudgesx<-sign(temp$Ratio)
nudgesy<-sign(temp$Ratio)
nudgesy[temp$Name=='SAP2']<-(1)
nudgesy[temp$Name=='SAP7']<-(-1.3)
nudgesx[temp$Name=='SAP7']<-(0.8)
nudgesx[temp$Name=='SAP4']<-(0.8)
nudgesy[temp$Name=='SAP4']<-(-0.9)

		
ggplot(temp11,aes(x=Mean,y=Mean1,xmin=Min,xmax=Max,ymin=Min1,ymax=Max1))+
		geom_errorbar(color='grey60',size=.25)+
		geom_errorbarh(color='grey60',size=.25)+
		geom_point(size=1,fill=1)+
		geom_text(aes(label=Name),temp,size=1.75,color=1,nudge_y=0.1*nudgesy,nudge_x=-0.2*nudgesx)+
		geom_abline(intercept=0,slope=1)+
		#scale_fill_manual('Genes',values=colors[levels(factor(temp$Name))])+
		scale_x_log10('Gene expression in Wildtype (normalized count)',limits=ran1)+
		scale_y_log10(expression(paste('Gene expression in ',Delta,italic('sap1/2/3/8/99'),' (normalized count)')),limits=ran1)+
		scale_shape_manual('Nitrogen source',values=c(21,24),labels=c('BSA','Nitrogen depletion'))+
		
		annotate('label',x=2,y=250000,label='Log growth',size=2.5)+
		ggtitle('Opaque cells - AmS')+
		theme_bw(8)+theme(legend.key.width=unit(1.6,'lines'),legend.key.height=unit(0.7,'lines'),legend.title.align=0.5,legend.text=element_text(size=8),axis.text=element_text(size=8),panel.grid.major=element_line(color='grey70'),panel.grid.minor=element_line(color=NA),strip.background=element_rect(color=NA,fill=NA),strip.text=element_text(size=8),plot.title=element_text(hjust=0.5))+
		guides(fill=guide_legend(override.aes=list(shape=21),order=1))

#############

#############
dev.off()
#############







