# Naomi Ziv 
# nz375@nyu.edu
# Analysis for paper: 
# Genetic and non-genetic determinants of cell-growth variation assessed by high-throughput microscopy

# This script fits mixed-effect models to diffrent phenotypes (growth rate, lag duration, growth rate deviation, protein fluorescence) in four different datasets
# It also performs various forms of data analysis (fitting the Monod equation to growth rates, calculating correlation between growth rates and CIT1 expression)
# It creates many variables (mainly in the form of data-frames) which organize all data produced from the anaylsis
# This script does not make any plots
# A list of variables used for plotting can be found at the end of the script 
# Set working directory to directory containing data files

####################################
# Packages 
####################################

library(lme4)
library(lmodel2)

####################################
# Custom functions
####################################

#####################################################################
#####################################################################
#####################################################################

##################
#combpara.R
##################

#This function combines the estimates and standard errors of fixed effect parameters
#It assumes a model/data with Strain and Nutrient.concentration (+interaction)
#It returns a data frame 

#response should match the response variable as it appears in the model/data

combpara<-function(model,data,response='GR'){
	
grid <- expand.grid(
	Strain = levels(data$Strain), 
	Nutrient.Concentration = levels(data$Nutrient.Concentration),
	response=0)
colnames(grid)[3]<-response

mm<-model.matrix(terms(model),grid)
#creating a matrix of 1 and 0 specifing which parameters apply to which conditions 
 
grid[,3]=mm %*% fixef(model)
#matrix multiplication, just sums the individual parameters estimates for each condition

pvar1<-diag(mm %*% tcrossprod(vcov(model),mm))
#vcov is the variance-covariance matrix of main parameters 
#tcrossprod is the crossproduct of a transpose - same as vcov(model) %*% t(mm),
#The result is the combined variance for each condition (For example: variance parameter 1 + variance parameter 2 + 2*co-variance between parameters 1 and 2)

grid<-transform(grid,low=grid[,3]-1.96*sqrt(pvar1),high=grid[,3]+1.96*sqrt(pvar1))

return(grid)	
}

##################
#subtrand.R
##################

#This function subtracts random effects (conditional means) from data
#It assumes a model/data with random effects
#It returns a named 1 column matrix 

#response should match the response variable as it appears in the model/data

subtrand<-function(model,data,response='GR',trans=function(x){return(x)}){
	
	#find response column
	col<-which(colnames(data)==response)
	
	#intialize 
	mat<-matrix(data[,col],nrow(data),1)
	colnames(mat)<-paste('T.',response,sep='')
	
	#random effects
	rand<-ranef(model,drop=T,postVar=T)
	split<-strsplit(names(rand),split=':')
	
	#loop on the random effect terms
	for (i in 1:length(rand)){
		
		#if random term is not an interaction
		#####################################
		if (length(split[[i]])==1){
			
			#find random term column
			indata<-which(colnames(data)==split[[i]])
			
			#loop on random effect levels
			for (j in names(rand[[i]])){
				
				#find rows in data
				rows<-which(data[,indata]==j)
				# subtract effect
				mat[rows]<-mat[rows]-trans(rand[[i]][j])
				
				rm(rows)
			}
			rm(indata,j)	
		}
		
		#if random term is an interaction
		#####################################
		if (length(split[[i]])==2){
			
			#find random term column
			indata1<-which(colnames(data)==split[[i]][1])
			indata2<-which(colnames(data)==split[[i]][2])
			
			#loop on random effect levels
			for (j in names(rand[[i]])){
				
				#split interaction into individual levels
				level<-strsplit(j,':')
				#find rows in data
				rows<-which(data[,indata1]==level[[1]][1] & data[,indata2]==level[[1]][2])
				# subtract effect
				mat[rows]<-mat[rows]-trans(rand[[i]][j])
				
				rm(level,rows)
			}
			rm(indata1,indata2,j)
		}
		
	}
	
	return(mat)
	
}

##################
#bestmodel.R
##################

#This function fits the monod model (or alternative models) of substrate-limited growth 
#It assumes a data with Strain, Nutrient.concentration 
#It assumes Nutrient.Concentration is a factor and transforms it to a numeric value

#It was used on transformed data (after subtraction of plate and well effects) 

#It returns a list with a data frame (predict - contains predicted values) and a list (list - contains coefficients, comparison to basic monod model and AIC/BIC values).

#conc is the vector of concentrations to predict

#response should match the response variable as it appears in the data

bestmodel<-function(data,model='response~((Nutrient*max)/(k+Nutrient))',start=c(max=1,k=1),response='T.GR',conc=seq(0,6,0.01)){
	
	col<-which(colnames(data)==response)
	
	#data
	monod<-data.frame(
			response=data[,col],
			Strain=data$Strain,
			Nutrient=as.numeric(levels(data$Nutrient.Concentration)[data$Nutrient.Concentration]))
	
	Predict<-NA
	list<-vector('list',length(levels(data$Strain)))
	names(list)<-levels(data$Strain)
	
	co=0
	for (i in levels(data$Strain)){
		co=co+1
		
		list[[co]]<-vector('list',3)
		
		#monod
		nls0<-nls(response~((Nutrient*max)/(k+Nutrient)),
			data=subset(monod,Strain==i),start=c(max=1,k=1))
			
		nls1<-nls(model,
			data=subset(monod,Strain==i),start=start)
			
		Predict<-c(Predict,predict(nls1,list(Nutrient=conc)))
		list[[co]][[1]]<-summary(nls1)$coef
		list[[co]][[2]]<-anova(nls0,nls1)
		list[[co]][[3]]<-matrix(c(AIC(nls0),AIC(nls1),BIC(nls0),BIC(nls1)),2,2,dimnames=list(c('Monod','Alternative'),c('AIC','BIC')))
			
	}

predict<-data.frame(
			Strain=rep(levels(data$Strain),each=length(conc)),
			Nutrient=rep(conc,length(levels(data$Strain))),
			Predict=Predict[2:length(Predict)])
			 
return(list(predict=predict,list=list))
}
 
#####################################################################
#####################################################################
#####################################################################

####################################
# Data
# Extreme values were removed during data generation
# Lag duration was corrected for well position (LWO, WILD and SOIL)
####################################

data<-read.table('Ziv.et.al.2013.data.txt',header=T)
data$plate<-factor(data$plate)
data$Nutrient.Concentration<-factor(data$Nutrient.Concentration)
data<-split(data,data$dataset)
LWOdata<-droplevels(data$LWO)
WILDdata<-droplevels(data$WILD)
SOILdata<-droplevels(data$SOIL)
rm(data)
CIT1data<-read.table('Ziv.et.al.2013.CIT1data.txt',header=T)
CIT1data$plate<-factor(CIT1data$plate)
CIT1data$Nutrient.Concentration<-factor(CIT1data$Nutrient.Concentration)

####################################
# Proportion of lagging cells in dataset LWO
# Proportion taken before correction for well position (lag.old)
####################################

nolag<-function(x){
	data<-x[!is.na(x)] 
	sum(data > (-1) & data < 1)/length(data)
	}
	
laggersLWO<-aggregate(lag.old~Strain+Nutrient.Concentration,data=LWOdata,FUN=nolag)

#data for testing models only for "laggers"
#lagdataLWO<-LWOdata[which(LWOdata$lag.old>=1),]

####################################
# LWO and WILD - main datasets
# Growth rate, growth rate deviation and Lag duration
####################################

##############################
##### Model construction #####  
##############################

# ###### With REML=F for testing models

# GR1<-lmer(GR~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate),data=LWOdata,REML=F)

# L1log<-lmer(log(lag+1)~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate),data=LWOdata,REML=F)

# GR2<-lmer(GR~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate),data=WILDdata,REML=F)

# L2log<-lmer(log(lag+1)~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate),data=WILDdata,REML=F)

# ###### Test terms 

# anova(update(GR1,.~.-(1|plate)),GR1) #3.969e-14 ***
# anova(update(GR1,.~.-(1|well)),GR1) #2.2e-16 ***
# anova(update(GR1,.~.-(Strain:Nutrient.Concentration)),GR1) #1.586e-07 ***
# anova(GR1,update(GR1,.~.+(1|Strain:plate))) #0.9969
# anova(GR1,update(GR1,.~.+(1|Strain:well))) #0.9979
# anova(GR1,update(GR1,.~.+(1|Nutrient.Concentration:plate))) #2.2e-16 ***
# anova(GR1,update(GR1,.~.+(1|Nutrient.Concentration:well))) #0.9979
# #need addtional Nutrient.Concentration:plate

# anova(update(L1log,.~.-(1|plate)),L1log) #2.2e-16 ***
# anova(update(L1log,.~.-(1|well)),L1log) #2.2e-16 ***
# anova(update(L1log,.~.-(Strain:Nutrient.Concentration)),L1log) #2.2e-16 ***
# anova(L1log,update(L1log,.~.+(1|Strain:plate))) #0.03298 *
# anova(L1log,update(L1log,.~.+(1|Strain:well))) #1
# anova(L1log,update(L1log,.~.+(1|Nutrient.Concentration:plate))) #7.386e-13 ***
# anova(L1log,update(L1log,.~.+(1|Nutrient.Concentration:well))) #1
# #need addtional Nutrient.Concentration:plate
# #need addtional Strain:plate

# anova(update(GR2,.~.-(1|plate)),GR2) #2.2e-16 ***
# anova(update(GR2,.~.-(1|well)),GR2) #2.2e-16 ***
# anova(update(GR2,.~.-(Strain:Nutrient.Concentration)),GR2) #2.2e-16 ***
# anova(GR2,update(GR2,.~.+(1|Strain:plate))) #0.02117 *
# anova(GR2,update(GR2,.~.+(1|Strain:well))) #1
# anova(GR2,update(GR2,.~.+(1|Nutrient.Concentration:plate))) #1.285e-13 ***
# anova(GR2,update(GR2,.~.+(1|Nutrient.Concentration:well))) #1
# #need addtional Nutrient.Concentration:plate
# #need addtional Strain:plate

# anova(update(L2log,.~.-(1|plate)),L2log) #2.2e-16 ***
# anova(update(L2log,.~.-(1|well)),L2log) #2.2e-16 ***
# anova(update(L2log,.~.-(Strain:Nutrient.Concentration)),L2log) #2.2e-16 ***
# anova(L2log,update(L2log,.~.+(1|Strain:plate))) #1.275e-06 ***
# anova(L2log,update(L2log,.~.+(1|Strain:well))) #1
# anova(L2log,update(L2log,.~.+(1|Nutrient.Concentration:plate))) #1.493e-06 ***
# anova(L2log,update(L2log,.~.+(1|Nutrient.Concentration:well))) #1
# #need addtional Nutrient.Concentration:plate
# #need addtional Strain:plate

######Final models - fit with REML=T

GR1<-lmer(GR~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate)+(1|Nutrient.Concentration:plate),data=LWOdata)

L1log<-lmer(log(lag+1)~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate)+(1|Nutrient.Concentration:plate)+(1|Strain:plate),data=LWOdata)

#For showing why we do the log transformation - plot residuals verses fitted values
#L1<-lmer(lag~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate)+(1|Nutrient.Concentration:plate)+(1|Strain:plate),data=LWOdata)

GR2<-lmer(GR~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate)+(1|Nutrient.Concentration:plate)+(1|Strain:plate),data=WILDdata)

L2log<-lmer(log(lag+1)~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate)+(1|Nutrient.Concentration:plate)+(1|Strain:plate),data=WILDdata)

###############################
##### Parameter estimates #####  
###############################

###### Combine parameter estimates and standard errors
###### Only considers fixed effects (Strain+Nutrient.Concentration)
###### transform back estimates for lag

combGR1<-combpara(GR1,LWOdata)
combL1<-combpara(L1log,LWOdata,'lag')
combL1<-cbind(combL1[,1:2],exp(combL1[,3:5])-1)

combGR2<-combpara(GR2,WILDdata)
combL2<-combpara(L2log,WILDdata,'lag')
combL2<-cbind(combL2[,1:2],exp(combL2[,3:5])-1)

###### Random effects

randGR1<-ranef(GR1,drop=T,postVar=T)
randL1<-ranef(L1log,drop=T,postVar=T)

randGR2<-ranef(GR2,drop=T,postVar=T)
randL2<-ranef(L2log,drop=T,postVar=T)

###### Add growth rate residuals for later estimation of deviations

LWOdata<-cbind(LWOdata,Res.GR=NA)
LWOdata$Res.GR[!is.na(LWOdata$GR)]<-resid(GR1)

WILDdata<-cbind(WILDdata,Res.GR=NA)
WILDdata$Res.GR[!is.na(WILDdata$GR)]<-resid(GR2)

#######################################################
##### Modeling the variance by modeling residuals #####  
#######################################################

# ###### With REML=F for testing models

# GRD1<-lmer(log(abs(Res.GR))~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate), data=LWOdata, REML=F)

# GRD2<-lmer(log(abs(Res.GR))~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate), data=WILDdata, REML=F)

# ###### Test terms 

# anova(update(GRD1,.~.-(1|plate)),GRD1) # 2.231e-06 ***
# anova(update(GRD1,.~.-(1|well)),GRD1) #2.2e-16 ***
# anova(update(GRD1,.~.-(Strain:Nutrient.Concentration)),GRD1) #4.277e-05 ***
# anova(GRD1,update(GRD1,.~.+(1|Strain:plate))) # 8.365e-05 ***
# anova(GRD1,update(GRD1,.~.+(1|Strain:well))) #0.9926
# anova(GRD1,update(GRD1,.~.+(1|Nutrient.Concentration:plate))) #1
# anova(GRD1,update(GRD1,.~.+(1|Nutrient.Concentration:well))) # 0.9926
# # #need addtional Strain:plate

# anova(update(GRD2,.~.-(1|plate)),GRD2) #2.091e-10 ***
# anova(update(GRD2,.~.-(1|well)),GRD2) #2.2e-16 ***
# anova(update(GRD2,.~.-(Strain:Nutrient.Concentration)),GRD2) #1.973e-07 ***
# anova(GRD2,update(GRD2,.~.+(1|Strain:plate))) #3.014e-05 ***
# anova(GRD2,update(GRD2,.~.+(1|Strain:well))) #0.9914
# anova(GRD2,update(GRD2,.~.+(1|Nutrient.Concentration:plate))) #6.325e-11 ***
# anova(GRD2,update(GRD2,.~.+(1|Nutrient.Concentration:well))) #0.9914
# # #need addtional Nutrient.Concentration:plate
# # #need addtional Strain:plate

######Final models - fit with REML=T

GRD1<-lmer(log(abs(Res.GR))~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate)+(1|Strain:plate), data=LWOdata)

GRD2<-lmer(log(abs(Res.GR))~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate)+(1|Nutrient.Concentration:plate)+(1|Strain:plate), data=WILDdata)

###########################################
##### Parameter estimates - residuals #####  
###########################################

combGRD1<-combpara(GRD1,LWOdata,'Res.GR')
combGRD1<-cbind(combGRD1[,1:2],exp(combGRD1[,3:5]))

combGRD2<-combpara(GRD2,WILDdata,'Res.GR')
combGRD2<-cbind(combGRD2[,1:2],exp(combGRD2[,3:5]))

randGRD1<-ranef(GRD1,drop=T,postVar=T)
randGRD2<-ranef(GRD2,drop=T,postVar=T)

###############################
##### Data transformation #####  
###############################

###### Transformation of original data by subtracting plate and well (and their interactions) conditional means

LWOdata<-cbind(LWOdata,laglog=log(LWOdata$lag+1),Reslog=log(abs(LWOdata$Res.GR)))
LWOdata<-cbind(LWOdata, subtrand(GR1,LWOdata), subtrand(L1log,LWOdata,'laglog'),subtrand(GRD1,LWOdata,'Reslog'))
LWOdata<-cbind(LWOdata,T.lag=exp(LWOdata$T.laglog)-1,T.Res=exp(LWOdata$T.Reslog))

WILDdata<-cbind(WILDdata,laglog=log(WILDdata$lag+1),Reslog=log(abs(WILDdata$Res.GR)))
WILDdata<-cbind(WILDdata, subtrand(GR2,WILDdata), subtrand(L2log,WILDdata,'laglog'),subtrand(GRD2,WILDdata,'Reslog'))
WILDdata<-cbind(WILDdata,T.lag=exp(WILDdata$T.laglog)-1,T.Res=exp(WILDdata$T.Reslog))

###################################
##### Monod for Lab/Wine/Oak ###### 
###################################

monod0<-bestmodel(LWOdata)

monod2<-bestmodel(LWOdata,model='response~(max*(smin-Nutrient))/(smin-Nutrient-k)',start=c(max=1,k=1,smin=.01))

monod5<-bestmodel(LWOdata,model='response~(k*a-Nutrient*max+Nutrient*a)/(-Nutrient-k)',start=c(max=1,k=1,a=.01))

monod3<-bestmodel(LWOdata,model='response~((-Nutrient*max) + (k*a))/(-Nutrient-k)',start=c(max=1,k=1,a=.01))

wester<-bestmodel(LWOdata,model='response~log(Nutrient)*y+x',start=c(y=1,x=1))

#The first element has predicted values for a range of concentrations
#The second element has coefficients, errors and model comparisons
#Standard errors are based on linearization (from summary.nls)
#However similar standerd errors were also estimated by bootstrapping

#########################
##### Data summary ###### 
#########################

aggcount<-function(x){sum(!is.na(x))}

#####Summary by well

aggwellLWO<-cbind(
	aggregate(GR~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=mean,na.rm=T),
	lag=aggregate(lag~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=median,na.rm=T)$lag,
	GR.SD=aggregate(GR~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=sd,na.rm=T)$GR,
	lag.mad=aggregate(lag~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=mad,na.rm=T,constant=1)$lag,
	Res.GR=aggregate(abs(Res.GR)~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=mean,na.rm=T)$'abs(Res.GR)',
	T.GR=aggregate(T.GR~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=mean,na.rm=T)$T.GR,
	T.lag=aggregate(T.lag~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=median,na.rm=T)$T.lag,
	T.Res=aggregate(T.Res~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=mean,na.rm=T)$T.Res,
	count=aggregate(GR~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=aggcount)$GR,
	count.lag=aggregate(lag~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=aggcount)$lag
)

aggwellWILD<-cbind(
	aggregate(GR~Strain+Nutrient.Concentration+plate+well,data=WILDdata,FUN=mean,na.rm=T),
	lag=aggregate(lag~Strain+Nutrient.Concentration+plate+well,data=WILDdata,FUN=median,na.rm=T)$lag,
	GR.SD=aggregate(GR~Strain+Nutrient.Concentration+plate+well,data=WILDdata,FUN=sd,na.rm=T)$GR,
	lag.mad=aggregate(lag~Strain+Nutrient.Concentration+plate+well,data=WILDdata,FUN=mad,na.rm=T,constant=1)$lag,
	Res.GR=aggregate(abs(Res.GR)~Strain+Nutrient.Concentration+plate+well,data=WILDdata,FUN=mean,na.rm=T)$'abs(Res.GR)',
	T.GR=aggregate(T.GR~Strain+Nutrient.Concentration+plate+well,data=WILDdata,FUN=mean,na.rm=T)$T.GR,
	T.lag=aggregate(T.lag~Strain+Nutrient.Concentration+plate+well,data=WILDdata,FUN=median,na.rm=T)$T.lag,
	T.Res=aggregate(T.Res~Strain+Nutrient.Concentration+plate+well,data=WILDdata,FUN=mean,na.rm=T)$T.Res,
	count=aggregate(GR~Strain+Nutrient.Concentration+plate+well,data=WILDdata,FUN=aggcount)$GR,
	count.lag=aggregate(lag~Strain+Nutrient.Concentration+plate+well,data=WILDdata,FUN=aggcount)$lag
)

#########################
##### Lag variance ###### 
#########################

#####median/mad loess regression
#####excluding wells with very low counts (3 out of 336)
#####doesn't change conclusion (p-value: 0.185 verses 0.193)

lagvar<-cbind(
	aggregate(lag~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=median,na.rm=T),
	lag.mad=aggregate(lag~Strain+Nutrient.Concentration+plate+well,data=LWOdata,FUN=mad,na.rm=T,constant=1)$lag,
	residuals=NA)[aggwellLWO$count.lag>50,]

lagvar$residuals=loess(lag.mad~lag,data=lagvar)$residuals

#summary(aov(residuals~Strain,data=lagvar))

####################################
# CIT1
# Growth rate and CIT1 fluorescence
# Anaylsis of lag duration and growth rate deviation is similar to LWO
# Lag duration has not been corrected for well position  
####################################

##############################
##### Model construction #####  
##############################

# ###### With REML=F for testing models

# GR3<-lmer(GR~Nutrient.Concentration+(1|well)+(1|plate),data=CIT1data,REML=F)

# pix1<-lmer(log(perpix)~Nutrient.Concentration+(1|well)+(1|plate),data=CIT1data,REML=F)
# Using log transformation because of heteroscedasticity

# ###### Test terms 

# anova(update(GR3,.~.-(1|plate)),GR3) #4.157e-06 ***
# anova(update(GR3,.~.-(1|well)),GR3) #2.2e-16 ***
# anova(GR3,update(GR3,.~.+(1|Nutrient.Concentration:plate))) #4.994e-07 ***
# anova(GR3,update(GR3,.~.+(1|Nutrient.Concentration:well))) #0.9972
# #need addtional Nutrient.Concentration:plate

# anova(update(pix1,.~.-(1|plate)),pix1) #2.2e-16 ***
# anova(update(pix1,.~.-(1|well)),pix1) #2.2e-16 ***
# anova(pix1,update(pix1,.~.+(1|Nutrient.Concentration:plate))) #3.632e-08 ***
# anova(pix1,update(pix1,.~.+(1|Nutrient.Concentration:well))) #0.9983
# #need addtional Nutrient.Concentration:plate

######Final models - fit with REML=T

GR3<-lmer(GR~Nutrient.Concentration+(1|well)+(1|plate)+(1|Nutrient.Concentration:plate),data=CIT1data)

pix1<-lmer(log(perpix)~Nutrient.Concentration+(1|well)+(1|plate)+(1|Nutrient.Concentration:plate),data=CIT1data)

###############################
##### Parameter estimates #####  
##### Data transformation #####
##### Data summary        #####  
###############################

combGR3<-combpara(GR3,CIT1data)
combpix1<-combpara(pix1,CIT1data,'perpix')

randGR3<-ranef(GR3,drop=T,postVar=T)
randpix1<-ranef(pix1,drop=T,postVar=T)

CIT1data<-cbind(CIT1data,pixlog=log(CIT1data$perpix))
CIT1data<-cbind(CIT1data, subtrand(GR3,CIT1data), subtrand(pix1,CIT1data,'pixlog'))

#####Summary by well

aggwellCIT1<-cbind(
	aggregate(GR~Nutrient.Concentration+plate+well,data=CIT1data,FUN=mean,na.rm=T),
	pix=aggregate(perpix~Nutrient.Concentration+plate+well,data=CIT1data,FUN=mean,na.rm=T)$perpix,
	pix.SD=aggregate(perpix~Nutrient.Concentration+plate+well,data=CIT1data,FUN=sd,na.rm=T)$perpix,
	pixlog=aggregate(pixlog~Nutrient.Concentration+plate+well,data=CIT1data,FUN=mean,na.rm=T)$pixlog,
	GR.SD=aggregate(GR~Nutrient.Concentration+plate+well,data=CIT1data,FUN=sd,na.rm=T)$GR,
	pixlog.SD=aggregate(pixlog~Nutrient.Concentration+plate+well,data=CIT1data,FUN=sd,na.rm=T)$pixlog,
	T.GR=aggregate(T.GR~Nutrient.Concentration+plate+well,data=CIT1data,FUN=mean,na.rm=T)$T.GR,
	T.pixlog=aggregate(T.pixlog~Nutrient.Concentration+plate+well,data=CIT1data,FUN=mean,na.rm=T)$T.pixlog,
	count=aggregate(GR~Nutrient.Concentration+plate+well,data=CIT1data,FUN=aggcount)$GR
	)

##########################################
##### Growth rate / CIT1 correlation #####  
##########################################

#correlation of all data
#cor.test(CIT1data$T.pixlog,CIT1data$T.GR)
#pearson=(-0.8285306), p-value < 2.2e-16,  -0.8312943 -0.8257260 
#cor.test(CIT1data$T.pixlog,CIT1data$T.GR,method='s')
#spearman=(-0.8352062), p-value < 2.2e-16

#spearman correlation within conditions
correlation<-data.frame(Nutrient.Concentration=levels(CIT1data$Nutrient.Concentration),estimate=NA,p.value=NA,bootlow=NA,boothigh=NA,bootSE=NA)
co=0
for (i in levels(CIT1data$Nutrient.Concentration)){
	co=co+1
	temp<-cor.test(
		subset(CIT1data,Nutrient.Concentration==i)$T.pixlog,
		subset(CIT1data,Nutrient.Concentration==i)$T.GR,method='s')
	correlation$estimate[co]<-temp$estimate
	correlation$p.value[co]<-temp$p.value
	rm(temp)
}

#function for calculating confidence intervals by bootstrapping
bootspear<-function(x,y,samp=1000){
	index<-which(!is.na(x))
	mat<-matrix(sample(index,length(index)*samp,replace=T),length(index),samp)
	vec<-numeric(samp)
	for(i in 1:samp){
		vec[i]<-cor(x[mat[,i]],y[mat[,i]],method='s')
	}
	return(c(quantile(vec,c(0.025,0.975)),sd(vec)))
}

#spearman correlation confidence intervals within conditions
#1000 gives same results as 10000
co=0
for (i in levels(CIT1data$Nutrient.Concentration)){
	co=co+1
	temp<-bootspear(
		subset(CIT1data,Nutrient.Concentration==i)$T.pixlog,
		subset(CIT1data,Nutrient.Concentration==i)$T.GR)
	correlation[co,c('bootlow','boothigh','bootSE')]<-temp
	rm(temp)
}

#Type 2 regression - use RMA for plotting
#regressIImod will have all the model objects to access if needed
#regressII will only have parameters for plotting
regressIImod<-as.list(levels(CIT1data$Nutrient.Concentration)) 
names(regressIImod)<-levels(CIT1data$Nutrient.Concentration)
regressII<-data.frame(Nutrient.Concentration=levels(CIT1data$Nutrient.Concentration),slope=NA,intercept=NA,x1=NA,x2=NA,y1=NA,y2=NA)
co=0
for (i in levels(CIT1data$Nutrient.Concentration)){
	co=co+1
	regressIImod[[co]]<-lmodel2(T.pixlog~T.GR,data=subset(CIT1data,Nutrient.Concentration==i),range.y='relative',range.x='relative')
	regressII$slope[co]<-regressIImod[[co]]$regression.results[4,3] #4th row is RMA
	regressII$intercept[co]<-regressIImod[[co]]$regression.results[4,2] #4th row is RMA
}
#calculation of line segments for each condition
regressII$x1<-aggregate(T.GR~Nutrient.Concentration,data=CIT1data,FUN=quantile,probs=0.01)$T.GR
regressII$x2<-aggregate(T.GR~Nutrient.Concentration,data=CIT1data,FUN=quantile,probs=0.99)$T.GR
regressII$y1<-regressII$slope*regressII$x1+regressII$intercept
regressII$y2<-regressII$slope*regressII$x2+regressII$intercept 

####################################
# SOIL
# Growth rate and growth rate deviation
# For the purpose of assessing reproducability in growth rate deviation
####################################

##############################
##### Model construction #####  
##############################

# ###### With REML=F for testing models

# GR4<-lmer(GR~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate),data=SOILdata,REML=F)

# anova(update(GR4,.~.-(1|well)),GR4) #2.2e-16 ***
# anova(update(GR4,.~.-(1|plate)),GR4) #1
# anova(update(GR4,.~.-(Strain:Nutrient.Concentration)),GR4) #0.001315 **
# anova(GR4,update(GR4,.~.+(1|Strain:plate))) # 2.654e-14 ***
# anova(GR4,update(GR4,.~.+(1|Strain:well))) #1
# anova(GR4,update(GR4,.~.+(1|Nutrient.Concentration:plate))) #0.8506
# anova(GR4,update(GR4,.~.+(1|Nutrient.Concentration:well))) # 1
## Despite plate not being significant, strainXplate is

######Final models - fit with REML=T

GR4<-lmer(GR~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|plate)+(1|Strain:plate)+(1|well),data=SOILdata)

#######################################################
##### Modeling the variance by modeling residuals #####  
#######################################################

SOILdata<-cbind(SOILdata,Res.GR=NA)
SOILdata$Res.GR[!is.na(SOILdata$GR)]<-resid(GR4)

# ###### With REML=F for testing models

# GRD4<-lmer(log(abs(Res.GR))~Strain+Nutrient.Concentration+Strain:Nutrient.Concentration+(1|well)+(1|plate),data=SOILdata,REML=F)

# anova(update(GRD4,.~.-(1|well)),GRD4) #2.2e-16 ***
# anova(update(GRD4,.~.-(1|plate)),GRD4) #4.802e-06 ***
# anova(update(GRD4,.~.-(Strain:Nutrient.Concentration)),GRD4) #0.9797
# anova(GRD4,update(GRD4,.~.+(1|Strain:plate))) # 0.06673 
# anova(GRD4,update(GRD4,.~.+(1|Strain:well))) #0.9938
# anova(GRD4,update(GRD4,.~.+(1|Nutrient.Concentration:plate))) #0.3945
# anova(GRD4,update(GRD4,.~.+(1|Nutrient.Concentration:well))) # 0.9938

# More minimal model
# GRD4<-lmer(log(abs(Res.GR))~Strain+Nutrient.Concentration+(1|well)+(1|plate),data=SOILdata,REML=F)

# anova(update(GRD4,.~.-(1|well)),GRD4) #2.2e-16 ***
# anova(update(GRD4,.~.-(1|plate)),GRD4) #4.803e-06 ***
# anova(update(GRD4,.~.-(Strain)),GRD4) #2.2e-16 ***
# anova(update(GRD4,.~.-(Nutrient.Concentration)),GRD4) #0.8469
# No need for Nutrient.Concentration

######Final models - fit with REML=T

GRD4<-lmer(log(abs(Res.GR))~Strain+(1|well)+(1|plate),data=SOILdata)

###############################
##### Parameter estimates #####  
##### Data transformation #####
##### Data summary        #####  
###############################

combGR4<-combpara(GR4,SOILdata)

combGRD4<-combpara(GRD4,SOILdata,'Res.GR')
combGRD4<-cbind(combGRD4[,1:2],exp(combGRD4[,3:5]))

randGR4<-ranef(GR4,drop=T,postVar=T)

randGRD4<-ranef(GRD4,drop=T,postVar=T)

SOILdata<-cbind(SOILdata,Reslog=log(abs(SOILdata$Res.GR)))
SOILdata<-cbind(SOILdata, subtrand(GR4,SOILdata),subtrand(GRD4,SOILdata,'Reslog'))
SOILdata<-cbind(SOILdata,T.Res=exp(SOILdata$T.Reslog))

#####Summary by well

aggwellSOIL<-cbind(
	aggregate(GR~Strain+Nutrient.Concentration+plate+well,data=SOILdata,FUN=mean,na.rm=T),
	GR.SD=aggregate(GR~Strain+Nutrient.Concentration+plate+well,data=SOILdata,FUN=sd,na.rm=T)$GR,
	Res.GR=aggregate(abs(Res.GR)~Strain+Nutrient.Concentration+plate+well,data=SOILdata,FUN=mean,na.rm=T)$'abs(Res.GR)',
	T.GR=aggregate(T.GR~Strain+Nutrient.Concentration+plate+well,data=SOILdata,FUN=mean,na.rm=T)$T.GR,
	T.Res=aggregate(T.Res~Strain+Nutrient.Concentration+plate+well,data=SOILdata,FUN=mean,na.rm=T)$T.Res,
	count=aggregate(GR~Strain+Nutrient.Concentration+plate+well,data=SOILdata,FUN=aggcount)$GR
)

####################################
# Save subset of workspace for plotting
####################################

save(LWOdata,WILDdata,laggersLWO,combGR1,combL1,combGR2,combL2,randGR1,randL1,randGR2,randL2,combGRD1,combGRD2,randGRD1,randGRD2,monod0,monod2,monod3,monod5,wester,aggwellLWO,aggwellWILD,CIT1data,combGR3,combpix1,randGR3,randpix1,correlation,regressII,aggwellCIT1,SOILdata,combGR4,randGR4,combGRD4,randGRD4,aggwellSOIL,lagvar,file='Ziv.et.el.2013.Analysis.Rdata')