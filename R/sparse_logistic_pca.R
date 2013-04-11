# Andrew's first coding

inv.logit.mat <- function(x, min = 0, max = 1) {
  p <- exp(x)/(1 + exp(x))
  which.large=is.na(p) & !is.na(x)
  p[which.large]=1
  p * (max - min) + min
}

sparse.logistic.pca <- function(dat,lambda=0,k=2,quiet=TRUE,max.iters=100) {
  # From Lee, Huang, Hu (2010)
  # Does not deal with missing data
  # Uses the uniform bound for the log likelihood
  q=as.matrix(2*dat-1)
  n=nrow(dat)
  d=ncol(dat)
  
  # Initialize #
  ##################
  mu=rnorm(d)
  udv=svd(scale(dat,center=TRUE,scale=FALSE))
  A=udv$u[,1:k]
  B=udv$v[,1:k] %*% diag(udv$d[1:k])
  # row.names(A)=row.names(dat); row.names(B)=colnames(dat)
  loss.trace=numeric(max.iters)
  
  for (m in 1:max.iters) {
    theta=outer(rep(1,n),mu)+A %*% t(B)
    X=as.matrix(theta+4*q*(1-inv.logit.mat(q*theta)))
    Xcross=X-A %*% t(B)
    mu=as.numeric(1/n*t(Xcross) %*% rep(1,n))
    
    theta=outer(rep(1,n),mu)+A %*% t(B)
    X=as.matrix(theta+4*q*(1-inv.logit.mat(q*theta)))
    Xstar=X-outer(rep(1,n),mu)
    A=Xstar %*% B %*% solve(t(B) %*% B)
    A=qr.Q(qr(A))
    
    theta=outer(rep(1,n),mu)+A %*% t(B)
    X=as.matrix(theta+4*q*(1-inv.logit.mat(q*theta)))
    Xstar=X-outer(rep(1,n),mu)
    C=t(Xstar) %*% A
    B=abs(B)/(abs(B)+4*n*lambda)*C
    
    loglike=sum(log(inv.logit.mat(q*(outer(rep(1,n),mu)+A %*% t(B)))))
    penalty=n*lambda*sum(abs(B))
    loss.trace[m]=-loglike+penalty
    
    if (!quiet) 
      cat(m,"  ",zapsmall(-loglike),"   ",zapsmall(penalty),"     ",-loglike+penalty, "\n")
    
    if (m>15) {
      if ((loss.trace[m-1]-loss.trace[m])<(.0001)) 
        break
    }
    
  }
  zeros=sum(abs(B)<1e-10)
  BIC=-2*loglike+log(n)*(d+n*k+sum(abs(B)>=1e-10))
  return(list(mu=mu,A=A,B=B,zeros=zeros,BIC=BIC,iters=m,loss.trace=loss.trace[1:m]))
}
