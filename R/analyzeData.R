library(ggplot2)
setwd("C:/Users/Andrew/SparseLogisticPCA")
source("R/sparse_logistic_pca.R")
ggmatplot <- function(mat,bw=TRUE,rotate.y=90,rownames=NULL) {
  require(ggplot2)
  if (is.null(rownames))
    rownames=(nrow(mat)<100 & any(!is.null(rownames(mat))))
  matm=melt(as.matrix(mat))
  matm$Var1=as.factor(matm$Var1)
  matm$Var2=as.factor(matm$Var2)
  if (is.null(rownames(mat))) 
    rownames(mat)=1:nrow(mat)
  if (is.null(colnames(mat))) 
    colnames(mat)=1:ncol(mat)
  
  
  plt<-ggplot(matm,aes(Var2,Var1))+geom_tile(aes(fill=value))+
    scale_y_discrete(limits = rev(rownames(mat)))+
    scale_x_discrete(limits = colnames(mat))+labs(x=NULL,y=NULL)
  if (bw) {
    plt<-plt+scale_fill_continuous(low="white",high="black")
  }
  if (rotate.y!=0) {
    plt<-plt+theme(axis.text.x = element_text(angle=rotate.y,vjust=0)) 
  }
  if (!rownames) {
    plt<-plt+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  return(plt)
}

# Load Data
locations=read.csv("Data/locations.csv")
# snpdata=read.csv("Data/RawCombinedData.csv")

snp.mat=read.csv("Data/RawSNPMatrix.csv")
is.na(snp.mat)<-(snp.mat=="NN")

snp.bin.mat=read.csv("Data/SNPBinaryMatrix.csv",row.names=1)
snp.bin.mat=snp.bin.mat[,colSums(snp.bin.mat,na.rm=TRUE)>0]
snp.bin.mat.na=snp.bin.mat
snp.bin.mat.na[is.na(snp.bin.mat.na)]=0.5

ethn=as.factor(substring(rownames(snp.bin.mat),1,3))

# Visualze Data
ggmatplot(as.matrix(is.na(snp.mat[3:102])),bw=F)
ggmatplot(snp.mat[,3:102],bw=F)
ggmatplot(snp.bin.mat[,sample(1:ncol(snp.bin.mat),400)],rownames=F)+coord_equal()
ggsave("Plots/AllBinaryData.pdf",ggmatplot(snp.bin.mat)+coord_equal()+theme(legend.position="none"),
       width=4.33*10,height=10)


# Standard PCA
pca=prcomp(snp.bin.mat.na)
qplot(PC1,PC2,data=data.frame(pca$x),size=I(3),colour=ethn)+
  labs(colour="Ethnicity")+coord_equal()
plot(cumsum(pca$sdev^2)/sum(pca$sdev^2),ylim=c(0,1))
barplot(pca$sdev[1:20]^2)

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

slpca=sparse.logistic.pca(snp.bin.mat,lambda=0.0015,k=10,quiet=FALSE,max.iters=1000,randstart=TRUE)
lpca=sparse.logistic.pca(snp.bin.mat,lambda=0,k=10,quiet=FALSE,max.iters=1000,randstart=TRUE)
pca=svd(scale(snp.bin.mat.na,T,F))

qplot(X1,X2,data=data.frame(slpca$A),colour=ethn)
qplot(X1,X2,data=data.frame(lpca$A),colour=ethn)

pairs(slpca$A,col=(1:3)[ethn])
pairs(lpca$A,col=(1:3)[ethn])

plot(slpca$B[,1],lpca$B[,1]); abline(v=0); abline(0,1,lty=2)
plot(slpca$B[,2],lpca$B[,2]); abline(v=0); abline(0,1,lty=2)