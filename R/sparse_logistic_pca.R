inv.logit.mat <- function(x, min = 0, max = 1) {
  p <- exp(x)/(1 + exp(x))
  which.large=is.na(p) & !is.na(x)
  p[which.large]=1
  p * (max - min) + min
}

sparse.logistic.pca <- function(dat,lambda=0,k=2,quiet=TRUE,max.iters=100,
                                randstart=FALSE,procrustes=TRUE) {
  # From Lee, Huang, Hu (2010)
  # Uses the uniform bound for the log likelihood
  q=as.matrix(2*dat-1)
  q[is.na(q)]<-0 # forces x to be equal to theta when data is missing
  n=nrow(dat)
  d=ncol(dat)
  
  # Initialize #
  ##################
  if (!randstart) {
    mu=colMeans(q)
    udv=svd(scale(q,center=TRUE,scale=FALSE))
    A=matrix(udv$u[,1:k],n,k)
    B=matrix(udv$v[,1:k],d,k) %*% matrix(diag(udv$d[1:k]),k,k)
  } else {
    mu=rnorm(d)
    A=matrix(runif(n*k,-1,1),n,k)
    B=matrix(runif(d*k,-1,1),d,k)
  }
  # row.names(A)=row.names(dat); row.names(B)=colnames(dat)
  loss.trace=numeric(max.iters)
  
  for (m in 1:max.iters) {
    last.mu=mu
    last.A=A
    last.B=B
    
    theta=outer(rep(1,n),mu)+A %*% t(B)
    X=as.matrix(theta+4*q*(1-inv.logit.mat(q*theta)))
    Xcross=X-A %*% t(B)
    mu=as.numeric(1/n*t(Xcross) %*% rep(1,n))
    
    theta=outer(rep(1,n),mu)+A %*% t(B)
    X=as.matrix(theta+4*q*(1-inv.logit.mat(q*theta)))
    Xstar=X-outer(rep(1,n),mu)
    if (procrustes) {
      M=svd(Xstar %*% B)
      A=M$u %*% t(M$v)
    } else {
      A=Xstar %*% B %*% solve(t(B) %*% B)
      A=qr.Q(qr(A))
    }
    
    theta=outer(rep(1,n),mu)+A %*% t(B)
    X=as.matrix(theta+4*q*(1-inv.logit.mat(q*theta)))
    Xstar=X-outer(rep(1,n),mu)
    C=t(Xstar) %*% A
    B=abs(B)/(abs(B)+4*n*lambda)*C
    
    loglike=sum(log(inv.logit.mat(q*(outer(rep(1,n),mu)+A %*% t(B))))[!is.na(dat)])
    penalty=n*lambda*sum(abs(B))
    loss.trace[m]=(-loglike+penalty)/sum(!is.na(dat))
    
    if (!quiet) 
      cat(m,"  ",zapsmall(-loglike),"   ",zapsmall(penalty),"     ",-loglike+penalty, "\n")
    
    if (m>15) {
      if ((loss.trace[m-1]-loss.trace[m])<(1e-6))
        break
    }
  }
  if (loss.trace[m-1]<loss.trace[m]) {
    mu=last.mu
    A=last.A
    B=last.B
    m=m-1
    
    loglike=sum(log(inv.logit.mat(q*(outer(rep(1,n),mu)+A %*% t(B))))[!is.na(dat)])
  }
  zeros=sum(abs(B)<1e-10)
  BIC=-2*loglike+log(n)*(d+n*k+sum(abs(B)>=1e-10))
  return(list(mu=mu,A=A,B=B,zeros=zeros,BIC=BIC,iters=m,loss.trace=loss.trace[1:m]))
}
