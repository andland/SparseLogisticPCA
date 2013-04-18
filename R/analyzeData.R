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

snpdata=read.csv("Data/RawCombinedData.csv")
snp.mat=read.csv("Data/RawSNPMatrix.csv")
is.na(snp.mat)<-(snp.mat=="NN")
snp.bin.mat=read.csv("Data/SNPBinaryMatrix.csv",row.names=1)
# snp.bin.mat=snp.bin.mat[,colSums(snp.bin.mat,na.rm=TRUE)>0]

# Visualze Data
ggmatplot(as.matrix(is.na(snp.mat[3:102])),bw=F)
ggmatplot(snp.mat[,3:103],bw=F)
ggmatplot(snp.bin.mat[,sample(1:ncol(snp.bin.mat),400)],rownames=F)+coord_equal()
ggmatplot(snp.bin.mat)+coord_equal()


# Standard PCA
snp.bin.mat.na=snp.bin.mat
snp.bin.mat.na[is.na(snp.bin.mat.na)]=0.5
pca=prcomp(snp.bin.mat.na)
qplot(PC1,PC2,data=data.frame(pca$x),size=I(3),colour=as.factor(substring(rownames(snp.bin.mat),1,3)))+
  labs(colour="Ethnicity")+coord_equal()
plot(cumsum(pca$sdev^2)/sum(pca$sdev^2),ylim=c(0,1),type='b')
plot(pca$sdev^2)

