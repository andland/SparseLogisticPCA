sparse.var.pca <- function(dat,lambda=0,k=2,quiet=TRUE,max.iters=100) {
  require(grplasso)
  n=nrow(dat)
  p=ncol(dat)
  
  # Initialize #
  ##################
  mu=colMeans(dat)
  X=as.matrix(scale(dat,center=TRUE,scale=FALSE))
  X.vec=as.numeric(X)
  udv=svd(X)
  
  A=matrix(udv$u[,1:k],n,k)
  B=matrix(udv$v[,1:k],p,k) %*% matrix(diag(udv$d[1:k]),k,k)
  # row.names(A)=row.names(dat); row.names(B)=colnames(dat)
  loss.trace=numeric(max.iters)
  
  for (m in 1:max.iters) {
    A=X %*% B %*% solve(t(B) %*% B)
    A=qr.Q(qr(A))
    
    A.kron=diag(rep(1,p)) %x% A
    glas=grplasso(x=A.kron, y=X.vec, index=rep(1:p,each=k), lambda=lambda, center=FALSE, 
                  standardize=FALSE, model=LinReg(), control=grpl.control(trace=0))
    B=matrix(coef(glas),p,k,byrow=TRUE)
    
    loglike=-sum((X - A %*% t(B))^2)
    penalty=sum(lambda*sqrt(k)*sqrt(rowSums(B^2)))
    loss.trace[m]=-loglike+penalty
    
    if (!quiet) 
      cat(m,"  ",zapsmall(-loglike),"   ",zapsmall(penalty),"     ",-loglike+penalty, "\n")
    
    if (m>5) {
      if ((loss.trace[m-1]-loss.trace[m])<(.0001)) 
        break
    }
  }
  zeros=sum(rowSums(abs(B))<1e-10)
#   BIC=-2*loglike+log(n)*(p+n*k+sum(abs(B)>=1e-10))
  BIC=NULL
  return(list(mu=mu,A=A,B=B,zeros=zeros,BIC=BIC,iters=m,loss.trace=loss.trace[1:m]))
}

k=2
A=matrix(svd(scale(dat,center=TRUE,scale=FALSE))$u[,1:k],nrow(dat),k)
A.kron=diag(rep(1,ncol(dat))) %x% A
lambdas <- lambdamax(x=A.kron, y=as.numeric(scale(dat,center=TRUE,scale=FALSE)), 
                    index=rep(1:ncol(dat),each=k), center=FALSE, 
                    standardize=FALSE, model=LinReg()) * 0.5^seq(10,1,-1)

zs=numeric(length(lambdas))
for (i in 1:length(lambdas)) {
  print(i)
  scp=sparse.var.pca(dat,lambda=lambdas[i],k=2,quiet=TRUE,max.iters=100)
  zs[i]=scp$zeros
}
plot(lambdas,zs)
ggmatplot(dat)
ggmatplot(outer(rep(1,n),(scp$mu)) + scp$A %*% t(scp$B))

colnames(dat)[rowSums(abs(scp$B))<1e-10]
colSums(dat[,rowSums(abs(scp$B))<1e-10])
sort(colSums(dat))

plot(udv$u[,1],A[,1])
plot(udv$u[,2],A[,2])
plot(matrix(udv$v[,1:k],p,k) %*% matrix(diag(udv$d[1:k]),k,k)[,1],B[,1])
plot(matrix(udv$v[,1:k],p,k) %*% matrix(diag(udv$d[1:k]),k,k)[,2],B[,2])
plot(-diff(loss.trace[1:m]),log='y')
