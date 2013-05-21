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

sparse.logistic.pca.coord <- function(dat,lambdas=10^seq(-15,-1,len=10),k=2,quiet=TRUE,max.iters=100,conv.crit=1e-3,
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
  BICs=matrix(NA,length(lambdas),k,dimnames=list(lambdas,1:k))
  zeros.mat=matrix(NA,length(lambdas),k,dimnames=list(lambdas,1:k))
  iters=matrix(NA,length(lambdas),k,dimnames=list(lambdas,1:k))
  
  theta=outer(rep(1,n),mu)+A %*% t(B)
  X=as.matrix(theta+4*q*(1-inv.logit.mat(q*theta)))
  Xcross=X-A %*% t(B)
  mu=as.numeric(1/n*t(Xcross) %*% rep(1,n))
  
  for (m in 1:k) {
    last.A=A
    last.B=B
    
    theta=outer(rep(1,n),mu)+A %*% t(B)
    X=as.matrix(theta+4*q*(1-inv.logit.mat(q*theta)))
    Xm=X-(outer(rep(1,n),mu)+A[,-m] %*% t(B[,-m]))
    
    Bms=matrix(NA,d,length(lambdas))
    Ams=matrix(NA,n,length(lambdas))
    for (lambda in lambdas) {
      for (i in 1:max.iters) {
        if (sum(B[,m]^2)==0) {
          A[,m]=Xm %*% B[,m]
          break
        }
        A[,m]=Xm %*% B[,m]/sum(B[,m]^2)
        A[,m]=A[,m]/sqrt(sum(A[,m]^2))
        
        B.lse=t(Xm) %*% A[,m]
        B[,m]=sign(B.lse)*pmax(0,abs(B.lse)-lambda)
        
        loglike=sum(log(inv.logit.mat(q*(outer(rep(1,n),mu)+A %*% t(B))))[!is.na(dat)])
        penalty=0.25*lambda*sum(abs(B[,m]))
        cur.loss=(-loglike+penalty)/sum(!is.na(dat))
        
        if (!quiet) 
          cat(m,"  ",zapsmall(-loglike),"   ",zapsmall(penalty),"     ",-loglike+penalty, "\n")
        
#         if (!quiet & i>1 & last.loss<cur.loss) {
#           warning(paste("Loss did not decrease!!",m,which(lambda==lambdas),i))
#         }
        
        if (i>4) {
          if ((last.loss-cur.loss)/last.loss<conv.crit) {
            break
          }
        }
        last.loss=cur.loss
      }
      Bms[,lambda==lambdas]=B[,m]/ifelse(sum(B[,m]^2)==0,1,sqrt(sum(B[,m]^2)))
      Ams[,lambda==lambdas]=Xm %*% Bms[,lambda==lambdas]/ifelse(sum(Bms[,lambda==lambdas]^2)==0,1,sum(Bms[,lambda==lambdas]^2))
      
      BICs[lambda==lambdas,m]=-2*loglike+log(n*d)*(sum(abs(B)>=1e-10))
      zeros.mat[lambda==lambdas,m]=sum(abs(B[,m])<1e-10)
      iters[lambda==lambdas,m]=i
    }
    B[,m]=Bms[,which.min(BICs[,m])]
    A[,m]=Ams[,which.min(BICs[,m])]
  }
  
  if (normalize) {
    A=sweep(A,2,sqrt(colSums(B^2)),"*")
    B=sweep(B,2,sqrt(colSums(B^2)),"/")
  }
  
  zeros=sum(abs(B)<1e-10)
  BIC=-2*loglike+log(n*d)*(sum(abs(B)>=1e-10))
  return(list(mu=mu,A=A,B=B,zeros=zeros,zeros.mat=zeros.mat,BICs=BICs,BIC=BIC,lambdas=lambdas,iters=iters))
}
