library(tidyverse)

setwd('~/Dropbox/MattAuris/')

sampleinfo<-read_tsv('sampleinfo.txt')

gff_1<-read_tsv('~/Dropbox/MattAuris/C_auris_B11221_features.gff',comment='#',col_names=c('CHROM','SOURCE','FEATURE','START','END','SCORE','STRAND','FRAME','INFO'))
gff_2<-read_tsv('~/Dropbox/MattAuris/C_auris_B8441_version_s01-m01-r17_features.gff',comment='#',col_names=c('CHROM','SOURCE','FEATURE','START','END','SCORE','STRAND','FRAME','INFO'))
gff<-rbind(gff_1,gff_2)

names<-read_tsv('Names.txt')


findgene<-function(chrom, pos, gff){

    genedist=function(startdist,enddist,strand){
        if(strand=='+'){return(startdist)}
        if(strand=='-'){return(enddist)}
    }

	temp<-filter(gff,CHROM==chrom)
	temp2<-temp3<-pos-temp$START
	temp2[temp2<0]=NA
	temp3[temp3>0]=NA
	gene1<-temp[which.min(temp2),]
	gene2<-temp[which.max(temp3),]
	if(nrow(gene1)==1 && pos<=gene1$END){
            return(data.frame(GENE1=str_extract(gene1$INFO,'(?<=ID=)([:alnum:]|_)+(?=;)'),
                              GENE2=NA,
                              NOTE='Gene',
                              STRAND1=gene1$STRAND,
                              STRAND2=NA,
                              DISTANCE1=genedist(pos-gene1$START,gene1$END-pos,gene1$STRAND),
                              DISTANCE2=NA))
	}else if(nrow(gene1)==1 && pos>=gene1$END & nrow(gene2)==1){
            return(data.frame(GENE1=str_extract(gene1$INFO,'(?<=ID=)([:alnum:]|_)+(?=;)'),
                              GENE2=str_extract(gene2$INFO,'(?<=ID=)([:alnum:]|_)+(?=;)'),
                              NOTE='Intergenic',
                              STRAND1=gene1$STRAND,
                              STRAND2=gene2$STRAND,
                              DISTANCE1=pos-gene1$END,
                              DISTANCE2=gene2$START-pos))
	}else if(nrow(gene1)==0 && nrow(gene2)==1){
            return(data.frame(GENE1=NA,
                              GENE2=str_extract(gene2$INFO,'(?<=ID=)([:alnum:]|_)+(?=;)'),
                              NOTE='Intergenic_CHRstart',
                              STRAND1=NA,
                              STRAND2=gene2$STRAND,
                              DISTANCE1=NA,
                              DISTANCE2=gene2$START-pos))
	}else if(nrow(gene1)==1 && nrow(gene2)==0){
            return(data.frame(GENE1=str_extract(gene1$INFO,'(?<=ID=)([:alnum:]|_)+(?=;)'),
                              GENE2=NA,
                              NOTE='Intergenic_CHRend',
                              STRAND1=gene1$STRAND,
                              STRAND2=NA,
                              DISTANCE1=pos-gene1$END,
                              DISTANCE2=NA))
	}else{return(data.frame(GENE1=NA,GENE2=NA,NOTE=NA,STRAND1=NA,STRAND2=NA,DISTANCE1=NA,DISTANCE2=NA))}
	rm(temp,temp2,temp3,gene1,gene2)
}

findname<-function(gene1,gene2,names){
	temp<-filter(names,B8441==gene1|B11221==gene1)
	if(nrow(temp)==1){
		name1=paste(temp$Cauris,temp$Calb,temp$Scer,sep='_')}else{
			name1='NotInFile'}
	temp2<-filter(names,B8441==gene2|B11221==gene2)
	if(nrow(temp2)==1){
		name2=paste(temp2$Cauris,temp2$Calb,temp2$Scer,sep='_')}else{
			name2='NotInFile'}
	return(data.frame(NAME1=name1,NAME2=name2))
}

dataset=data.frame()

for(i in list.files(,'.vcf')){
	
	sample<-str_extract(i,'(?<=strain)[:graph:]+(?=_sorted)')
	
	data<-read_tsv(i,comment='##')%>%
		mutate(CHROM=`#CHROM`)%>%
		select(CHROM,POS:INFO,-ID,-FILTER)%>%
		mutate(INDEL=str_detect(INFO,'INDEL'))%>%
		mutate(SIZE=str_length(ALT)-str_length(REF))%>%
		mutate(READS=as.numeric(str_extract(INFO,'(?<=DP=)[:digit:]+')))%>%
		mutate(READS_2=str_extract(INFO,'(?<=DP4=)[:graph:]+(?=;)'))%>%
		separate(READS_2,c('REF.F','REF.R','ALT.F','ALT.R'),sep=',',convert=T)%>%
		mutate(READS_2=REF.F+REF.R+ALT.F+ALT.R,REF_2=REF.F+REF.R,ALT_2=ALT.F+ALT.R)%>%
		select(-INFO)%>%
		mutate(genes=map2(CHROM,POS,findgene,gff=filter(gff,FEATURE=='gene')))%>%
		unnest()%>%
		mutate(names=map2(GENE1,GENE2,findname,names=names))%>%
		unnest()%>%
		mutate(Samples=sample)
	
	dataset<-rbind(dataset,data)
	
}

dataset<-left_join(dataset,sampleinfo)%>%unite(varID,CHROM:ALT,remove=F)

par<-filter(dataset,Note=='Parent')
dmso<-filter(dataset,Note=='DMSO')

dataset<-dataset%>%mutate(PAR_VAR=varID%in%par$varID,DMSO_VAR=varID%in%dmso$varID)%>%select(-varID)


dataset


write_tsv(dataset,'VarDataset.tsv')




par_sum<-summarize(group_by(par,Genome), count=n(),indels=sum(INDEL),snps=sum(!INDEL),gene=sum(NOTE=='Gene'),intergenic=sum(NOTE!='Gene'))%>%pivot_longer(-Genome)
par_sum_q<-summarize(group_by(filter(par,QUAL>100&ALT_2/READS_2>0.75),Genome), count=n(),indels=sum(INDEL),snps=sum(!INDEL),gene=sum(NOTE=='Gene'),intergenic=sum(NOTE!='Gene'))%>%pivot_longer(-Genome)


ggplot(par_sum,aes(x=factor(name,c('count','indels','snps','gene','intergenic')),y=value,fill=Genome))+
			geom_col(position='dodge')+
			scale_x_discrete('')+
			theme_bw()
			
ggplot(par_sum_q,aes(x=factor(name,c('count','indels','snps','gene','intergenic')),y=value,fill=Genome))+
			geom_col(position='dodge')+
			scale_x_discrete('')+
			ggtitle('Filtered: QUAL>100 & ALT_2/READS_2>0.75')+
			theme_bw()
