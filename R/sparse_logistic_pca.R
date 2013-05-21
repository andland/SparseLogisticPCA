inv.logit.mat <- function(x, min = 0, max = 1) {
  p <- exp(x)/(1 + exp(x))
  which.large=is.na(p) & !is.na(x)
  p[which.large]=1
  p * (max - min) + min
}

sparse.logistic.pca <- function(dat,lambda=0,k=2,quiet=TRUE,max.iters=100,conv.crit=1e-5,
                                randstart=FALSE,procrustes=TRUE,lasso=TRUE,normalize=FALSE,
                                start.A,start.B,start.mu) {
  # From Lee, Huang, Hu (2010)
  # Uses the uniform bound for the log likelihood
  # Can only use lasso=TRUE if lambda is the same for all dimensions, 
  #   which is how this algorithm is coded
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
  if (!missing(start.B))
    B=sweep(start.B,2,sqrt(colSums(start.A^2)),"*")
  if (!missing(start.A))
    A=sweep(start.A,2,sqrt(colSums(start.A^2)),"/")
  if (!missing(start.mu))
    mu=start.mu
  
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
    if (lasso) {
      B.lse=t(Xstar) %*% A
      B=sign(B.lse)*pmax(0,abs(B.lse)-4*n*lambda)
    } else {
      C=t(Xstar) %*% A
      B=abs(B)/(abs(B)+4*n*lambda)*C
    }
    
    loglike=sum(log(inv.logit.mat(q*(outer(rep(1,n),mu)+A %*% t(B))))[!is.na(dat)])
    penalty=n*lambda*sum(abs(B))
    loss.trace[m]=(-loglike+penalty)/sum(!is.na(dat))
    
    if (!quiet) 
      cat(m,"  ",zapsmall(-loglike),"   ",zapsmall(penalty),"     ",-loglike+penalty, "\n")
    
    if (m>4) {
      if ((loss.trace[m-1]-loss.trace[m])<conv.crit)
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
  if (normalize) {
    A=sweep(A,2,sqrt(colSums(B^2)),"*")
    B=sweep(B,2,sqrt(colSums(B^2)),"/")
  }
  
  zeros=sum(abs(B)<1e-10)
  BIC=-2*loglike+log(n)*(d+n*k+sum(abs(B)>=1e-10))
  return(list(mu=mu,A=A,B=B,zeros=zeros,BIC=BIC,iters=m,loss.trace=loss.trace[1:m],lambda=lambda))
}

sparse.logistic.pca.coord <- function(dat,lambdas=0,k=2,quiet=TRUE,max.iters=100,conv.crit=1e-5,
                                randstart=FALSE,normalize=FALSE,
                                start.A,start.B,start.mu) {
  # From Lee, Huang (2013)
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
  if (!missing(start.B))
    B=sweep(start.B,2,sqrt(colSums(start.A^2)),"*")
  if (!missing(start.A))
    A=sweep(start.A,2,sqrt(colSums(start.A^2)),"/")
  if (!missing(start.mu))
    mu=start.mu
  
  # row.names(A)=row.names(dat); row.names(B)=colnames(dat)
  last.loss=1e10
  BICs=matrix(NA,length(lamabdas),k)
  
  theta=outer(rep(1,n),mu)+A %*% t(B)
  X=as.matrix(theta+4*q*(1-inv.logit.mat(q*theta)))
  Xcross=X-A %*% t(B)
  mu=as.numeric(1/n*t(Xcross) %*% rep(1,n))
  
  for (m in 1:k) {
    last.A=A
    last.B=B
    
    theta=outer(rep(1,n),mu)+A %*% t(B)
    X=as.matrix(theta+4*q*(1-inv.logit.mat(q*theta)))
    Xstar=X-outer(rep(1,n),mu)
    
    Bms=matrix(NA,d,length(lambdas))
    Ams=matrix(NA,n,length(lambdas))
    for (lambda in lambdas) {
      for (i in 1:max.iters) {
        B.lse=t(Xstar) %*% A
        B[,m]=sign(B.lse)*pmax(0,abs(B.lse)-lambda)
        
        
        loglike=sum(log(inv.logit.mat(q*(outer(rep(1,n),mu)+A %*% t(B))))[!is.na(dat)])
        penalty=0.25*lambda*sum(abs(B[,m]))
        cur.loss=(-loglike+penalty)/sum(!is.na(dat))
        
        if (!quiet) 
          cat(m,"  ",zapsmall(-loglike),"   ",zapsmall(penalty),"     ",-loglike+penalty, "\n")
        
        if (i>4) {
          if ((last.loss-cur.loss)/last.loss<conv.crit) {
            break
          }
        }
      }
      Bms[,lambda==lambdas]=B[,m]
      Ams[,lambda==lambdas]=A[,m]
      BICs[lambda==lambdas,k]=-2*loglike+log(n*d)*(sum(abs(B)>=1e-10))
    }
  }
  
  if (normalize) {
    A=sweep(A,2,sqrt(colSums(B^2)),"*")
    B=sweep(B,2,sqrt(colSums(B^2)),"/")
  }
  
  zeros=sum(abs(B)<1e-10)
  BIC=-2*loglike+log(n)*(d+n*k+sum(abs(B)>=1e-10))
  return(list(mu=mu,A=A,B=B,zeros=zeros,BIC=BIC,iters=m,loss.trace=loss.trace[1:m],lambda=lambda))
}
