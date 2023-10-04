#Analysis for paper: Multiple molecular events underlie stochastic switching between two heritable cell states in a eukaryotic system
#Plos Biology, 2022
#Naomi Ziv 11/2020, update 7/2021, 2/2022

#Organization  
library(flowCore)
library(ggcyto)
#Careful with loading order as want tidyverse (dplyr) filter function by default 
library(tidyverse)
library(scales)
#Also need Biobase installed
library(flowViz)
#For fitting
library(drc)


##############################################################################	
##############################################################################	
#Function for combining plate ID files and formatting for flowCore
#This has a bunch of very specifc formatting matching...
#Not very general

#file names were changed using example code:
# setwd('~/Desktop/')
# name='S.cer2'
# path<-paste(name,'/',sep='')
# #rename file, only once
  # temp<-list.files(path=paste(path,'24Hours_CA',sep=''))
  # temp<-separate(data.frame(temp),temp,c(NA,NA,'well',NA,NA,'fcs'))
  # temp<-unite(temp,name,well,fcs,sep='.')
  # temp$name<-paste('5',temp$name,sep='_')
  # setwd('~/Desktop/S.cer2/24Hours_CA/')
  # file.rename(list.files(),temp$name)

plate<-function(name,prefix='Scer',home='~/Desktop/Paper2021/'){
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
	plate<-unite(plate,'row',plate,Well,remove=F)
	rownames(plate)<-paste(plate$row,'.fcs',sep='')

	return(plate)
	
}	

##############################################################################	
##############################################################################	
##############################################################################	
##############################################################################	
#Function for getting the cell asic data from outputed files 

get.cellasic<-function(j,exp='190404',channels=c('gfp'),home='/Volumes/Newton/'){
	
	name<-paste(exp,j,sep='.')

	areas<-read.table(paste(home,exp,'.2/files/Areas.',name,'.txt',sep=''))
	tops<-read.table(paste(home,exp,'.2/files/Tops.',name,'.txt',sep=''))
	tops<-split(tops,rep(channels,each=nrow(areas)))
	means<-read.table(paste(home,exp,'.2/files/Means.',name,'.txt',sep=''))
	means<-split(means,rep(channels,each=nrow(areas)))
	cells<-read_lines(paste(home,exp,'.2/files/Cells.',name,'.txt',sep=''))
	
	
	data=cbind(time=rep(1:nrow(areas),ncol(areas)),gather(areas,key='cell',value='area'))
	
	for(i in 1:length(channels)){
		
		if(channels[i]=='gfp'){
			temp<-gather(tops[[i]])
			data<-cbind(data,temp$value)
			colnames(data)[ncol(data)]<-channels[i]		
		}else{
			temp<-gather(means[[i]])
			data<-cbind(data,temp$value)
			colnames(data)[ncol(data)]<-channels[i]	
		}
						
	}
	
			
	data$cell<-factor(data$cell,unique(data$cell),labels=cells)
	
	
	#needed to update regexp in 2/2022 "\\.(?=[:digit:]$)" "\\.(?=[[:digit:]]$)"
	data<-data%>%
			filter(area>0)%>%
			separate(cell,c('M','D'),"\\.(?=[[:digit:]]$)",remove=F,fill='left')%>%
			mutate(field=j,M=replace_na(M,0),exp=exp)
	
	return(data)
		
}

##############################################################################	
##############################################################################
##############################################################################	
##############################################################################	
#Function for getting order for pedigrees - cell asic

get.vec<-function(agg){
	
	agg$M[agg$cell=='0']<-NA
	agg<-agg%>%
		arrange(M,desc(D))
		
	vec<-'0'
	place=1
	ds<-filter(agg,M==vec[place])%>%pull(cell)
	vec<-c(vec,as.character(ds))
	
	while(length(vec)<nrow(agg)){
		place=place+1
		ds<-filter(agg,M==vec[place])%>%pull(cell)
		oldvec<-vec
		vec<-c(vec[1:place],as.character(ds))
		if(place<length(oldvec)){vec<-c(vec,oldvec[(place+1):length(oldvec)])}
		}
		
	vec=vec[length(vec):1]
	return(vec)

}

##############################################################################	
##############################################################################



#####################################	Scer - Data 	##################################
##########################################################################################

setwd('~/Desktop/Paper2021/')

#repeat 1
name<-'S.cer1'

#get info
info<-plate(name)
path<-paste(name,'/',sep='')

#get data
data<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info,data.frame(labelDescription=colnames(info))),alter.names=T)
#create ratios
data<-transform(data,GFP=`BB515.A`/`FSC.A`,mCherry=`PE.CF594.A`/`FSC.A`)

#repeat 2
name<-'S.cer2'

#get info
info2<-plate(name)
path<-paste(name,'/',sep='')

#get data
data2<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(info2,data.frame(labelDescription=colnames(info2))),alter.names=T)
#create ratios
data2<-transform(data2,GFP=`BB515.A`/`FSC.A`,mCherry=`PE.CF594.A`/`FSC.A`)

#####################################	Gates	##################################
##################################################################################

#gates were defined on repeat 1,
#gates were saved as S.cer.Gates.Rdata

# #get data
# DF<-do.call('rbind',fsApply(data,function(x){data.frame(rowname=identifier(x),exprs(x))}))
# DF<-info%>%rownames_to_column()%>%left_join(DF)

# # #Cell population
# temp<-sample_n(DF,100000)
# plot(temp$FSC.A,temp$SSC.A,xlim=c(0,200000))
# p1<-locator(8)
# p1m<-matrix(,8,2)
# colnames(p1m)<-c('FSC.A','SSC.A')
# p1m[,1]<-p1$x
# p1m[,2]<-p1$y
# p1m
          # # FSC.A       SSC.A
# # [1,]   4027.4871  -3831.834
# # [2,]    294.4983  13233.620
# # [3,]   7829.3304  47837.185
# # [4,]   7061.1312  95419.829
# # [5,]  58353.6132 163930.939
# # [6,] 130514.3804 152521.395
# # [7,] 136028.1332  53768.313
# # [8,]  31354.5856  -5004.498

# Pgate<-polygonGate(.gate=p1m)
# #xyplot(`SSC.A`~`FSC.A`,data=data,filter=Pgate,smooth=F)


# #GFP/BV650.A

# #For high GFP
# temp<-sample_n(filter(DF,group=='G'),20000)
# plot(log10(temp$GFP),log10(temp$BV650.A))
# p1<-locator(8)
# p1m<-matrix(,8,2)
# colnames(p1m)<-c('GFP','BV650.A')
# p1m[,1]<-10^p1$x
# p1m[,2]<-10^p1$y
# p1m
             # # GFP    BV650.A
# # [1,] 0.001658140  322.064273
# # [2,] 0.010692819  369.252499
# # [3,] 3.170792579 2320.626877
# # [4,] 1.695873273   39.283129
# # [5,] 0.518261598    9.979876
# # [6,] 0.030005243   14.521956
# # [7,] 0.003134042   23.976754
# # [8,] 0.002072906   44.509897
# Ggate<-polygonGate(.gate=p1m)

# #For medium GFP
# temp<-sample_n(filter(DF,group=='WG'),20000)
# plot(log10(temp$GFP),log10(temp$BV650.A))
# p1<-locator(8)
# p1m<-matrix(,8,2)
# colnames(p1m)<-c('GFP','BV650.A')
# p1m[,1]<-10^p1$x
# p1m[,2]<-10^p1$y
# p1m
             # # GFP    BV650.A
# # [1,] 0.001742340  412.12703
# # [2,] 0.012613808  546.37085
# # [3,] 0.106727269 3369.60284
# # [4,] 0.250855966 4519.67853
# # [5,] 0.323146685 1066.51633
# # [6,] 0.096233774   14.05245
# # [7,] 0.013876336   11.81769
# # [8,] 0.001643924   37.51769
# WGgate<-polygonGate(.gate=p1m)

# #For no GFP
# temp<-sample_n(filter(DF,group=='NG'),20000)
# plot(log10(temp$GFP),log10(temp$BV650.A))
# p1<-locator(8)
# p1m<-matrix(,8,2)
# colnames(p1m)<-c('GFP','BV650.A')
# p1m[,1]<-10^p1$x
# p1m[,2]<-10^p1$y
# p1m
             # # GFP    BV650.A
# # [1,] 0.0005588709  210.836585
# # [2,] 0.0081230196 1263.377319
# # [3,] 0.0122662685 1365.801432
# # [4,] 0.0208991005  588.255482
# # [5,] 0.0047373595    6.296307
# # [6,] 0.0020188278    6.006317
# # [7,] 0.0010087515   16.331419
# # [8,] 0.0006739103   52.133858
# Ngate<-polygonGate(.gate=p1m)


#Filtering
load('S.cer.Gates.Rdata')

#Flow errors - should have FSC threshold of 5000
#Also taking olny positive values for GFP/mCherry - only issue is cuts off background florecence for mCherry
rectGate<-rectangleGate(filterId='cells',"FSC.A"=c(5000,Inf),"BB515.A"=c(1,Inf),"PE.CF594.A"=c(1,Inf))	
#filter
data<-Subset(data,rectGate)	
data2<-Subset(data2,rectGate)	
#Filter cell population
data<-Subset(data,Pgate)
data2<-Subset(data2,Pgate)
#Split
data<-split(data,factor(pData(data)$group))
data2<-split(data2,factor(pData(data2)$group))
#filter
G<-Subset(data[[1]],Ggate)
NG<-Subset(data[[2]],Ngate)
WG<-Subset(data[[3]],WGgate)
G2<-Subset(data2[[1]],Ggate)
NG2<-Subset(data2[[2]],Ngate)
WG2<-Subset(data2[[3]],WGgate)

#get filtered data
DF.g<-do.call('rbind',fsApply(G,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF.ng<-do.call('rbind',fsApply(NG,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF.wg<-do.call('rbind',fsApply(WG,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF<-info%>%rownames_to_column()%>%left_join(rbind(DF.g,DF.ng,DF.wg))

DF.g<-do.call('rbind',fsApply(G2,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF.ng<-do.call('rbind',fsApply(NG2,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF.wg<-do.call('rbind',fsApply(WG2,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF.ca<-do.call('rbind',fsApply(data2[[4]],function(x){data.frame(rowname=identifier(x),exprs(x))}))
DF2<-info2%>%rownames_to_column()%>%left_join(rbind(DF.g,DF.ng,DF.wg,DF.ca))

DFall<-full_join(DF,DF2)

#save
save(data,data2,Ggate,Ngate,WGgate,Pgate,DFall,file='Scer.Rdata')


#####################################	Fits	##################################
##################################################################################

#summarize	
agg<-DFall%>%group_by(Well,Environment,Genotype,Comments,time,type,date,plate,group)%>%
			summarize(GFP=median(GFP),mCherry=median(mCherry),count=n())%>%
			filter(group!='CA')
			
agg2<-agg%>%group_by(Comments,Environment,time,group)%>%
			summarize(GFPm=mean(GFP),GFPsd=sd(GFP),mCm=mean(mCherry),mCsd=sd(mCherry),count=n())%>%
			mutate(GFPmin=GFPm-1.96*(GFPsd/sqrt(count)), GFPmax=GFPm+1.96*(GFPsd/sqrt(count)), mCmin=mCm-1.96*(mCsd/sqrt(count)),mCmax=mCm+1.96*(mCsd/sqrt(count)))%>%
			filter(group!='CA')

aggtile<-DFall%>%
		group_by(Genotype)%>%
		mutate(GGFP=ntile(log10(GFP),100))%>%
		group_by(Genotype,Comments,time,type,date,plate,GGFP)%>%
		summarise(MedGFP=median(GFP),MedCHE=median(mCherry))%>%ungroup()

#model hills
		
llmodels<-filter(DFall, group!='CA')%>%
				filter(	(Comments%in%c('WOR1-GFP','WOR1-GFP-MY')&Environment<2)|
						Comments%in%c('WOR1-GFP-206','WOR1-GFP-206-MY'))%>%
				group_by(Genotype,Comments,time,type,date,group,plate)%>%
				nest()%>%
				mutate(fit = map(data, ~ drm(log10(mCherry)~GFP,data=.,fct=LL2.4(),robust='median')))%>%
				mutate(pre = map(fit, predict), coef = map(fit, coef))
				
forplot<-dplyr::select(llmodels,-fit,-coef)%>%unnest()%>%ungroup()
forplot2<-dplyr::select(llmodels,-fit,-data, -pre)%>%
				unnest()%>%
				mutate(names = rep(letters[2:5],32))%>%
				ungroup()%>%
				spread(names,coef)%>%
				mutate(rowname=1:32,Ee=exp(e))

hill<-function(b,c,d,e,x=seq(min(forplot$GFP),max(forplot$GFP),length.out=1000)){return(data.frame(x=x,y=c+(d-c)/(1+exp(b*(log(x)-e)))))}

forplot3<-pmap_dfr(dplyr::select(forplot2,b,c,d,e),hill,.id='rowname')%>%
				mutate(rowname=as.numeric(rowname))%>%
				left_join(forplot2)

#some t.tests for comparing parameters
temp<-filter(forplot2,time=='24')
t.test(filter(temp,!grepl('206',Comments))$b*-1,filter(temp,grepl('206',Comments))$b*-1)
t.test(exp(filter(temp,!grepl('206',Comments))$e),exp(filter(temp,grepl('206',Comments))$e))

#save
save(data,data2,Ggate,Ngate,WGgate,Pgate,DFall,llmodels,agg,agg2,forplot,forplot2,forplot3,file='Scer.Rdata')
#load('Scer.Rdata')

		
#####################################	Data	##################################
##################################################################################
#FULL verses No-UTR

# #file names were changed using example code:
# setwd('~/Desktop/')
# name='20160906'
# path<-paste(name,'/',sep='')
# #rename file, only once
  # temp<-list.files(path=paste(path,'plate2',sep=''))
  # temp<-separate(data.frame(temp),temp,c(NA,NA,'well',NA,'fcs'))
  # temp<-unite(temp,name,well,fcs,sep='.')
  # temp$name<-paste('2',temp$name,sep='_')
  # setwd('~/Desktop/20160906/plate2/')
  # file.rename(list.files(),temp$name)

setwd('~/Desktop/Paper2021/')

#FULL verses No-UTR
name<-'20160906'

#get info
repinfo<-plate(name,'Flow')
path<-paste(name,'/',sep='')

#get data
rep<-read.flowSet(path=paste(path,'files',sep=''),phenoData=Biobase::AnnotatedDataFrame(repinfo,data.frame(labelDescription=colnames(repinfo))),alter.names=T)
#create ratios
rep<-transform(rep,GFP=`FITC.A`/`FSC.A`)

#get data
DFrep<-do.call('rbind',fsApply(rep,function(x){data.frame(rowname=identifier(x),exprs(x))}))
DFrep<-repinfo%>%rownames_to_column()%>%left_join(DFrep)

DFrep<-separate(DFrep,Environment,sep='_',into=c(NA,'Environment'))%>%
		filter(Comments%in%c('NOGFP','FULL','FULL.MY','NOUTR.MY','NOUTR'))
		
#save
save(DFrep,file='ScerSup.Rdata')
		


########################################################################################	
########################################################################################	
#####################################	Cell asic data	####################################
########################################################################################
setwd('~/Desktop/Paper2021/')

#####################################	Pedigrees	##################################
######################################################################################

#WOR1-GFP-206 201019 - 63
pedi1<-get.cellasic(63,'201019',channels=c('gfp'))
pedi1A<-aggregate(time~cell+M+D,pedi1,min)
pedi1$cell<-factor(pedi1$cell,get.vec(pedi1A))

sub1<-filter(pedi1, cell%in%filter(pedi1A,time<114)$cell)
daughters1<-mutate(read_csv('201019.63.celldivsion.csv'),time=125)

#WOR1-GFP 190404 155
pedi2<-get.cellasic(155,'190404',channels=c('gfp'))
pedi2<-filter(pedi2,time<100)
pedi2A<-aggregate(time~cell+M+D,pedi2,min)
pedi2$cell<-factor(pedi2$cell,get.vec(pedi2A))

sub2<-filter(pedi2, cell%in%filter(pedi2A,time<70)$cell)
daughters2<-mutate(read_csv('190404.155.celldivsion.csv'),time=105)


#####################################	Switching rates	##################################
######################################################################################

setwd('~/Desktop/Paper2021/')

switch<-read_delim('Cellasic.txt',delim='\t')%>%
		filter(Genotype%in%c('WOR1-GFP','WOR1-GFP-206'))%>%
		filter(Full)%>%
		filter(Experiment!='190404')%>%
		dplyr::select(Experiment,Genotype,Field,Switch12.24,RemoveMD, RemoveR)

agg<-switch%>%
		gather(,,Switch12.24,RemoveMD,RemoveR)%>%
		group_by(Experiment, Genotype,key,value)%>%
		summarize(count=n())
		
	
poisexp<-function(data){
	data=data[[1]]
	tot<-sum(data$count)
	zero<-filter(data,value==0)$count
	lamb<--log(zero/tot)
	return(data.frame(exp=dpois(data$value,lamb)*tot))
}
agg2<-agg%>%group_by(Experiment, Genotype,key)%>%
				nest(value,count)%>%
				
				mutate(exp=lmap(data,.f=poisexp))%>%unnest()

agg3<-agg2%>%spread(key,count)%>%gather(,'count',Switch12.24,RemoveMD, RemoveR,exp)
	
			
# ggplot(filter(agg2,grepl('Remove',key)),aes(x=value,y=count,fill=factor(key,c('RemoveMD','RemoveR'))))+
			# geom_col(position='dodge')+
			# geom_point(aes(y=exp))+
			# facet_wrap(~Genotype+Experiment,scales='free_y')

temp=filter(agg3,Genotype=='WOR1-GFP-206')
temp$count[is.na(temp$count)]<-0
#make small correction so pois expectation for 3 is actually 3 or more - sould add up to 168,
#change exp for 3 to 0.274 from 0.258 
temp2<-aggregate(count~key,temp,sum)
temp$count[temp$key=='exp'&temp$value==3]<-temp$count[temp$key=='exp'&temp$value==3]+((temp2$count[temp2$key=='RemoveMD'])-(temp2$count[temp2$key=='exp']))

#chisq tests, simutalted p.values because small expected values
chisq.test(filter(temp,key=='Switch12.24')$count,p=filter(temp,key=='exp')$count,rescale.p=T,simulate.p.value=T)
chisq.test(filter(temp,key=='RemoveMD')$count,p=filter(temp,key=='exp')$count,rescale.p=T,simulate.p.value=T)
chisq.test(filter(temp,key=='RemoveR')$count,p=filter(temp,key=='exp')$count,rescale.p=T,simulate.p.value=T)


#save
save(pedi1,pedi1A,pedi2,pedi2A,daughters1,daughters2,sub1,sub2,temp,file='pedi.Rdata')



#####################################	Traces	##################################
##################################################################################

#############trace plots
smallbud<-function(time){temp=min(time);time[time<(temp+3)]<-NA;return(time)}

setwd('~/Desktop/Paper2021/')

#For activation of mothers
sync<-rbind(	get.cellasic(96,'190404',channels=c('gfp')),
				get.cellasic(103,'190404',channels=c('gfp')),
				get.cellasic(106,'190404',channels=c('gfp')),
				get.cellasic(108,'190404',channels=c('gfp')),
				get.cellasic(114,'190404',channels=c('gfp')),
				get.cellasic(123,'190404',channels=c('gfp')),
				get.cellasic(126,'190404',channels=c('gfp')),
				get.cellasic(144,'190404',channels=c('gfp')),
				get.cellasic(155,'190404',channels=c('gfp')),
				get.cellasic(33,'190320',channels=c('gfp')),
				get.cellasic(37,'190320',channels=c('gfp')),
				get.cellasic(40,'190320',channels=c('gfp')),
				get.cellasic(67,'190320',channels=c('gfp')),
				get.cellasic(69,'190320',channels=c('gfp')),
				get.cellasic(63,'201019',channels=c('gfp')),
				get.cellasic(41,'201019',channels=c('gfp')),
				get.cellasic(156,'201019',channels=c('gfp')),
				get.cellasic(184,'201019',channels=c('gfp')),
				get.cellasic(123,'201019',channels=c('gfp')),
				get.cellasic(142,'201019',channels=c('gfp')),
				get.cellasic(127,'201019',channels=c('gfp')))

						
syncsub<-filter(sync, 	(cell%in%c('1')&field==96)|
						(cell%in%c('0')&field==106&time<87)|
						(cell%in%c('0')&field==103&time<51)|
						(cell%in%c('1.1')&field==108&time<57)|
						(cell%in%c('1')&field==114&time<71&time>4)|
						(cell%in%c('1.1')&field==123&exp==190404&time<69)|
						(cell%in%c('1.1.2')&field==126&time<83)|
						(cell%in%c('2')&field==144&time<64)|
						(cell%in%c('1.1')&field==155&time<69)|
						(cell%in%c('1.1.2')&field==33&time<57)|
						(cell%in%c('1.3')&field==37&time<58)|
						(cell%in%c('1.1.1')&field==40&time<71)|
						(cell%in%c('1.1')&field==67&time<75)|
						(cell%in%c('0')&field==69&time<53)|
						(cell%in%c('2')&field==123&exp==201019&time<112)|
						(cell%in%c('1.1')&field==142&time<107)|
						(cell%in%c('0')&field==127&time<117))
						
syncsub2<-filter(sync, 	(cell%in%c('2')&field==63&time<110)|
						(cell%in%c('1.1.1')&field==41&time<95)|
						(cell%in%c('2')&field==156&time<109)|
						(cell%in%c('3.2')&field==184))

						

syncsub<-group_by(syncsub,cell,field)%>%
			mutate(time=smallbud(time),Geno='GFP')%>%
			filter(!is.na(time))
			
syncsub2<-group_by(syncsub2,cell,field)%>%
			mutate(time=smallbud(time),Geno='mGFP')%>%
			filter(!is.na(time))
			

osc<-rbind(get.cellasic(4,'190320',channels=c('gfp')),get.cellasic(106,'190404',channels=c('gfp')),get.cellasic(80,'190320',channels=c('gfp')))

oscsub<-filter(osc, 	(cell=='0'&field==4&time<51)|
						(cell=='0'&field==106&time>45)|
						(cell=='f'&field==106&time<50)|
						(cell=='1.1.1'&field==80)&time>21)
						
sync2<-filter(sync,field==127&time<111)
sync2<-group_by(sync2,cell,field)%>%
			mutate(time=smallbud(time))%>%
			filter(!is.na(time))

#ggplot(syncsub,aes(x=time,y=gfp,group=paste(cell, field),color=cell))+geom_point()+geom_line()

#save
#also added re (defined below)
save(osc,oscsub,sync,sync2,syncsub,syncsub2,re,file='trace.Rdata')



########################################################################################	
########################################################################################	
#####################################	Other data	####################################
########################################################################################
				
#####################################	C.albicans - mother-daughter Length distributions 	########################
####################################################################################################################
	
re<-read_csv('190320LengthResults.csv')
re2<-read_csv('190404LengthResults.csv')
re3<-read_csv('201019LengthResults.csv')					

re<-rbind(cbind(re,exp='190320'),cbind(re2,exp='190404'),cbind(re3,exp='201019'))
			
			
#####################################	C.albicans - switching rates on plates 	####################################
####################################################################################################################		
plateSw<-rbind(read_delim('NZ-201125.txt','\t'),read_delim('NZ-201225.txt','\t'))
plateSw<-mutate(plateSw,Environment='GLU')
plateSw$Environment[grepl('glc',plateSw$Plate)]<-'GLC'

plateSw<-filter(plateSw,Environment=='GLU'&Genotype%in%c('WT','3xmotif','WOR1-GFP','WOR1-GFP-206'))

PSagg<-group_by(plateSw,Genotype,Cell,Batch,Environment)%>%
		summarize(Colonies=sum(Colonies),Sectors=sum(Sectors),count=n())%>%
		mutate(Percent=100*Sectors/Colonies)

PSagg3<-summarize(group_by(plateSw,Genotype,Cell,Environment),Colonies=sum(Colonies),Sectors=sum(Sectors))%>%
		mutate(Percent=100*Sectors/Colonies,SE=sqrt((Percent*(100-Percent))/(Colonies)),ymin=Percent-1.96*SE,ymax=Percent+1.96*SE)


