rm(list=ls())
library(reshape2)
# Should be where text files were downloaded from 
# http://snp.cshl.org/downloads/datafreeze.html.en
# I used data from http://snp.cshl.org/downloads/genotypes/2005-06_16c.1_phaseI/full/non-redundant/
setwd("U:/CF/Ethnic")


# this is from table S1 of Serre et al (2008)
# http://www.plosone.org/article/info:doi/10.1371/journal.pone.0001382#pone.0001382.s008
locations=read.csv("locations.csv",header=TRUE,stringsAsFactors=FALSE)

files=dir()
files=files[grep(".txt",files)]

f=1
for (f in 1:length(files)) {
  print(f)
  chdata=read.delim(files[f],sep=" ",header=TRUE,stringsAsFactors=FALSE)
  colnames(chdata)[1]="rsid"
  which.rs=chdata$rsid %in% locations$RSID
  if (sum(which.rs)>0) {
    ethnicity=substring(files[f],nchar(files[f])-6,nchar(files[f])-4)
    
    chdata=chdata[which.rs,]
    chdata.m=melt(chdata[,c(1,12:ncol(chdata))],id=1)
    chdata.m$variable <- as.character(chdata.m$variable)
    if(f==1) {
      snpdata=data.frame(Ethnicity=ethnicity,chdata.m,stringsAsFactors=FALSE)
    } else {
      snpdata=rbind(snpdata,data.frame(Ethnicity=ethnicity,chdata.m,stringsAsFactors=FALSE))
    }    
  }
}
colnames(snpdata)[3:4]=c("Subject","SNP")
write.csv(snpdata,"RawCombinedData.csv",row.names=FALSE)
table(snpdata$Ethnicity)
table(snpdata$Subject)
table(snpdata$SNP)


#######
snp.mat=dcast(snpdata, Ethnicity + Subject ~ rsid)
write.csv(snp.mat,"RawSNPMatrix.csv",row.names=FALSE)
dim(snp.mat)

table(snp.mat$Ethnicity)


#######
is.na(snp.mat)<-(snp.mat=="NN")
ggmatplot(as.matrix(is.na(snp.mat[3:100])),bw=F)
ggmatplot(snp.mat[,sample(3:ncol(snp.mat),100)],bw=F)

snp.bin.mat=matrix(0,nrow(snp.mat),ncol(snp.mat)-2)
colnames(snp.bin.mat)=colnames(snp.mat)[3:ncol(snp.mat)]
rownames(snp.bin.mat)=paste0(ifelse(snp.mat$Ethnicity=="CHB" | snp.mat$Ethnicity=="JPT","ASN",snp.mat$Ethnicity),
                             snp.mat$Subject)
is.na(snp.bin.mat)=is.na(snp.mat[,3:ncol(snp.mat)])

cols.rm=c()
for (c in 3:ncol(snp.mat)) {
  print(c)
  missing.ethn=table(substring(rownames(snp.bin.mat),1,3),is.na(snp.mat[,c]))
  if (any(missing.ethn[,1]==0)) {
    cols.rm=c(cols.rm,c-2)
  } else {
    snptbl=table(snp.mat[,c])
    homogeneous=substring(names(snptbl),1,1)==substring(names(snptbl),2,2)
    which.max.cat=names(which.max(snptbl[homogeneous]))
    snp.bin.mat[snp.mat[,c]!=which.max.cat,c-2]=1
  }
}
snp.bin.mat=snp.bin.mat[,-cols.rm]
# snp.bin.mat=snp.bin.mat[,colSums(snp.bin.mat,na.rm=TRUE)>0]
ggmatplot(snp.bin.mat[,sample(1:ncol(snp.bin.mat),400)],rownames=T)+coord_equal()
mean(is.na(snp.bin.mat)) # 0.527%
dim(snp.bin.mat) # 269x1322
write.csv(snp.bin.mat,"SNPBinaryMatrix.csv",row.names=TRUE)
# snp.bin.mat=read.csv("SNPBinaryMatrix.csv",row.names=1)

snp.bin.mat.na=snp.bin.mat
snp.bin.mat.na[is.na(snp.bin.mat.na)]=0.5
pca=prcomp(snp.bin.mat.na)
qplot(PC1,PC2,data=data.frame(pca$x),size=I(3),colour=as.factor(substring(rownames(snp.bin.mat),1,3)))+
  labs(colour="Ethnicity")+coord_equal()
plot(cumsum(pca$sdev^2)/sum(pca$sdev^2),ylim=c(0,1),type='b')
plot(pca$sdev^2)
