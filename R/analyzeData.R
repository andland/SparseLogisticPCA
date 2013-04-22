rm(list=ls())
library(ggplot2)
setwd("C:/Users/Andrew/SparseLogisticPCA")
source("R/sparse_logistic_pca.R")
source("C:/Users/Andrew/Dropbox/Stat/Research/MedData/genPCAfun.R")
# dat=read.csv("C:/Users/Andrew/Dropbox/Stat/Research/MedData/MedData_for_Analysis_PU.csv")
# meddata.nums=meddata[,5:ncol(meddata)]
# dat=meddata.nums[sort(sample(nrow(meddata.nums),100)),]
# dat=dat[rowSums(dat)>0,colSums(dat)>0]
# dat=read.csv("C:/Users/Andrew/Dropbox/Blog/HOF Factor Analysis/HOF votes.csv",row.names=1)
# dat=dat[rowSums(dat)>0,colSums(dat)>0]

# Load Data
locations=read.csv("Data/locations.csv")
# snpdata=read.csv("Raw Data/RawCombinedData.csv")

snp.mat=read.csv("Data/RawSNPMatrix.csv")
is.na(snp.mat)<-(snp.mat=="NN")

snp.bin.mat=read.csv("Data/SNPBinaryMatrix.csv",row.names=1)
snp.bin.mat=snp.bin.mat[,colSums(snp.bin.mat,na.rm=TRUE)>0]
snp.bin.mat.na=snp.bin.mat
snp.bin.mat.na[is.na(snp.bin.mat.na)]=0.5

ethn=as.factor(substring(rownames(snp.bin.mat),1,3))

# Visualze Data
ggmatplot(as.matrix(is.na(snp.mat[3:102])),bw=F)
ggmatplot(snp.mat[,3:502],bw=F)+coord_equal()
ggsave("Plots/SampleSNPData.pdf",ggmatplot(snp.mat[,3:402],bw=F)+coord_equal(),
       width=1.5*10,height=10)
ggmatplot(snp.bin.mat[,sample(1:ncol(snp.bin.mat),400)],rownames=F)+coord_equal()
udv=svd(scale(snp.bin.mat.na,T,F))
ggsave("Plots/AllBinaryData.pdf",ggmatplot(snp.bin.mat)+coord_equal()+theme(legend.position="none"),
       width=4.33*10,height=10)
ggsave("Plots/AllBinaryDataRand.pdf",ggmatplot(snp.bin.mat[,sample(ncol(snp.bin.mat))])+coord_equal()+theme(legend.position="none"),
       width=4.33*10,height=10)
ggsave("Plots/AllBinaryDataPC1.pdf",ggmatplot(snp.bin.mat[,order(udv$v[,1])])+coord_equal()+theme(legend.position="none"),
       width=4.33*10,height=10)
ggsave("Plots/AllBinaryDataPC2.pdf",ggmatplot(snp.bin.mat[,order(udv$v[,2])])+coord_equal()+theme(legend.position="none"),
       width=4.33*10,height=10)


ggmatplot(snp.bin.mat[,order(udv$v[,1])])+coord_equal()+theme(legend.position="none")

# Standard PCA
pca=prcomp(snp.bin.mat.na)
qplot(PC1,PC2,data=data.frame(pca$x),size=I(2),colour=ethn)+
  labs(colour="Ethnicity")+coord_equal()+labs(title="Standard PCA")
ggsave("Plots/StandardPCA.pdf")
pdf("Plots/StandardPCA Diagnostics.pdf")
plot(cumsum(pca$sdev^2)/sum(pca$sdev^2),ylim=c(0,1))
barplot(pca$sdev[1:20]^2)
dev.off()

# LDA
library(MASS)
snp.lda=lda(snp.bin.mat.na,ethn)
plot(snp.lda,col=(1:3)[ethn])

# Sparse Logistic PCA
lambdas=c(0,1.5^(-18:-10))
lambdas=seq(0,0.0006,length=10)
BICs=numeric(length(lambdas))
zeros=numeric(length(lambdas))

for (lambda in lambdas) {
  print(which(lambda==lambdas))
  
  slpca=sparse.logistic.pca(snp.bin.mat,lambda=lambda,k=2,quiet=TRUE,max.iters=100)
  
  zeros[lambda==lambdas]=slpca$zeros
  BICs[lambda==lambdas]=slpca$BIC
}
plot(lambdas[1:10],BICs[1:10],type='b')
plot(lambdas,zeros,type='b')
plot(zeros,BICs,type='b')
lambda=lambdas[which.min(BICs)]

slpca=sparse.logistic.pca(snp.bin.mat,lambda=lambda,k=10,quiet=TRUE,max.iters=100)

slpca=sparse.logistic.pca(snp.bin.mat,lambda=0.0015,k=10,quiet=FALSE,max.iters=100,randstart=T)
lpca=sparse.logistic.pca(snp.bin.mat,lambda=0,k=10,quiet=FALSE,max.iters=100,randstart=F)
pca=svd(scale(snp.bin.mat.na,T,F))

qplot(X1,X2,data=data.frame(slpca$A),colour=ethn)+labs(title=paste("Sparse Logistic PCA, k =",ncol(slpca$A)))
ggsave(paste0("Plots/Sparse Logistic PCA k=",ncol(slpca$A)," - PCA Start.pdf"))
qplot(X1,X2,data=data.frame(lpca$A),colour=ethn)+labs(title=paste("Logistic PCA, k =",ncol(lpca$A)))
ggsave(paste0("Plots/Logistic PCA k=",ncol(lpca$A)," - PCA Start.pdf"))

pairs(slpca$A,col=(1:3)[ethn])
pairs(lpca$A,col=(1:3)[ethn])

plot(slpca$B[,1],lpca$B[,1]); abline(v=0); abline(0,1,lty=2)
plot(slpca$B[,2],lpca$B[,2]); abline(v=0); abline(0,1,lty=2)