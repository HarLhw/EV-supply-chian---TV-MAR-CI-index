### Helper functions

# myslice <- function(xx, K, start, end){
#   if (K==2){
#     return(xx[start:end,,,drop=FALSE])
#   } else if (K==3){
#     return(xx[start:end,,,,drop=FALSE])
#   } else {
#     stop("not support tensor mode K > 3")
#   }
# }

# mat projection
matAR.PROJ <- function(xx, dim, r, t){
  xx.mat <- matrix(xx,t,dim[1]*dim[2])
  kroneck <- t(xx.mat[2:t,]) %*% xx.mat[1:(t-1),] %*% solve(t(xx.mat[1:(t-1),]) %*% xx.mat[1:(t-1),])
  return(projection(kroneck, r, dim[1],dim[2],dim[1],dim[2]))
}

# Tensor Times List
tl <- function(x, list_mat, k = NULL){
  if (is.null(k)){
    tensor(tensor(tensor(x, list_mat[[1]], 2, 2), list_mat[[2]], 2, 2), list_mat[[3]], 2, 2)
  } else if (k == 1){
    tensor(tensor(x, list_mat[[1]], 3, 2), list_mat[[2]], 3, 2)
  } else if (k == 2){
    aperm(tensor(tensor(x, list_mat[[1]], 2, 2), list_mat[[2]], 3, 2),c(1,3,2,4))
  } else if (k == 3){
    aperm(tensor(tensor(x, list_mat[[1]], 2, 2), list_mat[[2]], 2, 2),c(1,3,4,2))
  } else {
    stop("not support tensor mode K > 3")
  }
  
}

# standard error extraction
covtosd <- function(cov, dim, R){
  K <- length(dim)
  P <- length(R)
  sd = list()
  for (p in c(1:P)){
    if (is.na(R[p])) stop("p != length(R)")
    if (R[p] == 0) next
    sd[[p]] <- lapply(1:R[p], function(j) {lapply(1:K, function(i) {list()})})
  }
  for (i in c(1:P)){
    for (j in c(1:R[i])){
      for (k in c(1:K)){
        left <- sum(dim^2)*sum(R[0:(i-1)]) + sum(dim^2)*(j-1) + sum((dim^2)[1:(k-1)])+1
        right <- sum(dim^2)*sum(R[0:(i-1)]) + sum(dim^2)*(j-1) + sum((dim^2)[1:k])
        sd[[i]][[j]][[k]] <- array(sqrt(diag(cov)[left:right]), c(dim[k], dim[k]))
      }
    }
  }
  return(sd)
}

# Permutation matrix em
em <- function(m,n,i,j){
  ## m,n,i,j set \eqn{m \times n} zero matrix with \eqn{A_{ij} = 1}
  ## return: Permutation matrix em such that \eqn{A_{ij} = 1} and other entries equals 0.
  mat <- matrix(0,m,n)
  mat[i,j] <- 1
  return(mat)
}

# Permutation matrix pm
pm <- function(m,n){
  ## m: an array of dimensions of matrices \eqn{A_1,A_2,\cdots,A_k}
  ## n: length of time
  ## return: Permutation matrix pm
  mat <- matrix(0,m*n,m*n)
  for (i in c(1:n)){
    for (j in c(1:m)){
      mat <- mat + kronecker(em(n,m,i,j),t(em(n,m,i,j)))
    }
  }
  return(mat)
}

# rearrangement operator for tensor
trearrange <- function(A,dim){
  m1 = dim[1]; m2 = dim[2]; m3 = dim[3]
  n1 = m1; n2 = m2; n3 = m3
  m <- nrow(A)
  n <- ncol(A)
  if(n!=n1*n2*n3 | m!=m1*m2*m3){
    stop("wrong dimention with your input Phi for rearrangement")
  }
  ans <- divide(A,m1,n1)
  dim <- c(m1*n1,m2*n2,m3*n3)
  t <- array(0, dim)
  for (i in c(1:m1)){
    for (j in c(1:n1)){
      t[(j-1)*m1+i,,] <- mrearrange(ans[[i]][[j]],m2,m3,n2,n3)
    }
  }
  return(t)
}


divide <- function(A,m,n){
  # the inner function of "trearrange"
  c <- dim(A)[1]/m
  l <- dim(A)[2]/n
  tmp <- lapply(1:m, function(i){
    lapply(1:n, function(j){
      A[((i-1)*c+1):(i*c),((j-1)*l+1):(j*l)]
    })
  })
  return(tmp)
}


mrearrange <- function(A,m1,m2,n1,n2){
  # the inner function of "projection"
  # A: m1m2*n1n2
  # B: m1*n1
  # C: m2*n2
  # A \approx B \otimes C
  # return RA
  m <- nrow(A)
  n <- ncol(A)
  if(n!=n1*n2 | m!=m1*m2){
    stop("error m")
  }
  ans <- matrix(NA, m1*n1, m2*n2)
  for(i in 1:m1){
    for(j in 1:n1){
      ans[(j-1)*m1+i,] <- t(as.vector(A[(i-1)*m2+1:m2,(j-1)*n2+1:n2]))
    }
  }
  return(ans)
}

projection <- function(M,r,m1,m2,n1,n2){
  # the inner function of MAR1.projection
  # M: m1m2*n1n2
  # B: m1*n1
  # C: m2*n2
  # M \approx B \otimes C
  # return B and C
  RA <- mrearrange(M,m1,m2,n1,n2)
  RA.svd <- svd(RA,nu=r,nv=r)
  A <- list()
  for (i in c(1:r)){
    A[[i]] <- list(matrix(RA.svd$v[,i] * RA.svd$d[i], m2, n2), matrix(RA.svd$u[,i], m1, n1))
  }
  for (j in c(1:r)){
    A[[j]] <- rev(A[[j]])
    a <- c()
    for (i in c(1:2)){
      m <- A[[j]][[i]]
      if (i != 2){
        a[i] <- svd(m,nu=0,nv=0)$d[1]
        A[[j]][[i]] <- m/a[i]
      } else {
        A[[j]][[i]] <- m * prod(a)
      }
    }
  }
  return(A)
}

ten.proj <- function(tt, dim, R){
  ## inner func of "TenAR.proj"
  cpd <- rTensor::cp(rTensor::as.tensor(tt), num_components = R, max_iter = 100, tol = 1e-06)
  lam <- cpd$lambdas
  A.proj <- list()
  for (j in c(1:R)){
    u1 <- cpd$U[[1]][,j]
    u2 <- cpd$U[[2]][,j]
    u3 <- cpd$U[[3]][,j]
    f1 <- sqrt(sum(cpd$U[[1]][,j]^2))
    f2 <- sqrt(sum(cpd$U[[2]][,j]^2))
    f3 <- sqrt(sum(cpd$U[[3]][,j]^2))
    a1 <- u1/f1
    a2 <- u2/f2
    a3 <- u3*f1*f2*lam[j]
    A.proj[[j]] <- list(matrix(a1,dim[1],dim[1]),
                        matrix(a2,dim[2],dim[2]),
                        matrix(a3,dim[3],dim[3]))
  }
  return(fro.order(A.proj))
}


fro.rescale <- function(A){
  r <- length(A)
  k <- length(A[[1]])
  for (j in c(1:r)){
    a <- c()
    for (i in c(1:k)){
      m <- A[[j]][[i]]
      if (i < k ){
        a[i] <- norm(m,"f")
        A[[j]][[i]] <- m/a[i]
      } else if (i == k){
        A[[j]][[i]] <- m * prod(a)
      } else {
        print("WRONG dimension")
      }
    }
  }
  return(A)
}

svd.rescale <- function(A){
  r <- length(A)
  k <- length(A[[1]])
  for (j in c(1:r)){
    a <- c()
    for (i in c(1:k)){
      m <- A[[j]][[i]]
      if (i < k ){
        a[i] <- svd(m,nu=0,nv=0)$d[1]
        A[[j]][[i]] <- m/a[i]
      } else if (i == k){
        A[[j]][[i]] <- m * prod(a)
      } else {
        print("WRONG dimension")
      }
    }
  }
  return(A)
}

eigen.rescale <- function(A){
  r <- length(A)
  k <- length(A[[1]])
  for (j in c(1:r)){
    a <- c()
    for (i in c(1:k)){
      m <- A[[j]][[i]]
      if (i < k ){
        a[i] <- eigen(m)$values[1]
        A[[j]][[i]] <- m/a[i]
      } else if (i == k){
        A[[j]][[i]] <- m * prod(a)
      } else {
        print("WRONG dimension")
      }
    }
  }
  return(A)
}

fro.order <- function(A){
  R <- length(A)
  K <- length(A[[1]])
  if (R == 1){return(A)}
  A.norm <- c()
  for (j in c(1:R)){
    A.norm[j] <- Reduce("*",lapply(c(1:K), function(k) { norm(A[[j]][[k]], 'f')}))
  }
  order.norm <- order(A.norm, decreasing=TRUE)
  A.temp <- A
  for (j in c(1:R)){
    A[[j]] <- A.temp[[order.norm[j]]]
  }
  return(A)
}

ten.dis.A <- function(A, B, R, K){
  P = length(R)
  dis <- 0
  for (p in c(1:P)){
    for (r in c(1:R[p])){
      for (k in c(1:K)){
        dis <- dis + min(sum((A[[p]][[r]][[k]] - B[[p]][[r]][[k]])^2), sum((A[[p]][[r]][[k]] + B[[p]][[r]][[k]])^2))
      }
    }
  }
  return(sqrt(dis))
}

ten.dis.phi <- function(phi.A, phi.B){
  P <- length(phi.A)
  dis <- 0
  for (i in c(1:P)){
    dis <- dis + sqrt(sum((phi.A[[i]] - phi.B[[i]])^2))
  }
  return(dis)
}


ten.res <- function(xx,A,P,R,K,t){
  L1 = 0
  for (l in c(1:P)){
    if (R[l] == 0) next
    L1 <- L1 + Reduce("+",lapply(c(1:R[l]), function(n) {rTensor::ttl(abind::asub(xx, (1+P-l):(t-l), 1, drop=FALSE), A[[l]][[n]], (c(1:K) + 1))}))
  }
  res <- abind::asub(xx, (1+P):(t), 1, drop=FALSE) - L1
  return(res)
}


M.eigen <- function(A, R, P, dim){
  phi <- list()
  PP = P
  for (i in c(1:P)){
    if (sum(R[i:length(R)]) == 0){
      PP = i-1
      break
    }
    if (R[i] == 0){
      phi[[i]] = pracma::zeros(prod(dim))
    } else {
      phi[[i]] <- Reduce("+", lapply(1:R[i], function(j) {rTensor::kronecker_list(rev(A[[i]][[j]]))}))
    }
    if (i == 1){M <- phi[[1]]} else {M <- cbind(M, phi[[i]])}
  }
  K <- dim(phi[[1]])[[1]]
  
  M <- rbind(M, cbind(diag(K*(PP-1)), array(0,c(K*(PP-1),K))))
  return(max(Mod(eigen(M, only.values = TRUE)$values)))
}

specRadius <- function(M){
  return(max(Mod(eigen(M, only.values = TRUE)$values)))
}


###################################################################################################
loglikelihood_mle <- function(SIGMA_mle,res_mle){
  ## My OWN log likelihood function for MAR(1)
  dim_mle <- dim(res_mle)
  tem_mle <- 0
  for (i in 1:dim_mle[1]) {
    immediary <- 0
    immediary <- solve(SIGMA_mle[[1]])%*%res_mle[i,,]%*%solve(SIGMA_mle[[2]])%*%t(res_mle[i,,])
    tem_mle <- tem_mle+lava::tr(immediary)
  }
  llf_mle <- -dim_mle[2]%*%dim_mle[1]%*%log(det(SIGMA_mle[[2]]))-dim_mle[3]%*%dim_mle[1]%*%log(det(SIGMA_mle[[1]]))-tem_mle
  return(llf_mle)
}
#################################################################################################


likelihood <- function(xx, A, Sigma){
  if (!(mode(xx) == "S4")) {xx <- as.tensor(xx)}
  r <- length(A[[1]])
  dd <- dim(xx)
  t <- dd[1]
  dim <- dd[-1]
  k <- length(dd[-1])
  i = 1
  res <- ten.res(xx,A,P=1,R=r,K=k,t=t)@data
  Sigma.inv <- lapply(1:k, function (i) {solve(Sigma[[i]])})
  ll <- tl(res, Sigma.inv)
  l1 <- sum(diag(tensor(ll, res, c(1:4)[-(i+1)],c(1:4)[-(i+1)])))
  l2 <- 0
  for (i in c(1:k)){
    l2 = l2 - prod(dim[-i]) * (t-1) * (log(det(Sigma[[i]])))
  }
  return((l2 - l1)/2)
}


initializer <- function(xx, k1=1, k2=1){
  PROJ = MAR1.PROJ(xx)
  if (specRadius(PROJ$A1)*specRadius(PROJ$A2) < 1){
    return(list(A1=PROJ$A1,A2=PROJ$A2))
  }
  MAR = MAR1.LS(xx)
  if (specRadius(MAR$A1)*specRadius(MAR$A2) < 1){
    return(list(A1=MAR$A1,A2=MAR$A2))
  }
  RRMAR = MAR1.RR(xx, k1, k2)
  if (specRadius(MAR1.RR$A1)*specRadius(MAR1.RR$A2) < 1){
    return(list(A1=MAR1.RR$A1,A2=MAR1.RR$A2))
  }
  stop('causality condition of initializer fails.')
}

#############################################################################################
##  row normalized to have unit variance for each AFNS indicators
rowstd <- function(data_raw_use){
  dim_data <- dim(data_raw_use)
  data_raw_use_row_norm <- array(0,dim=c(dim_data[1],dim_data[2],dim_data[3]))
  for (i in  1:dim_data[1]) {
    for (j in  1:dim_data[2]) {
      data_raw_use_row_norm[i,j,] <- data_raw_use[i,j,]/std(data_raw_use[,j,])
    }
  }
  return(data_raw_use_row_norm)
}


###############################################################################################
# My own residuals function for MAR(1) model
residuals_mar <- function(A_MLE, data_use) {
  dim_data_use <- dim(data_use)
  resid_mar <- array(0, dim = c(dim_data_use[1]-1,dim_data_use[2],dim_data_use[3]))
  
  for (i in 1:(dim_data_use[1] - 1)) {
    resid_mar[i,,] <- data_use[i + 1,,] - A_MLE[[1]][[1]][[1]] %*% data_use[i,,] %*% t(A_MLE[[1]][[1]][[2]])
  }
  return(resid_mar)
}
##############################################################################################
# initializer <- function(xx, k1=1, k2=1){
#   dim = dim(xx)[-1]
#   p = dim(xx)[2]
#   q = dim(xx)[3]
# 
#   PROJ = MAR1.PROJ(xx)
#   if (specRadius(PROJ$A1)*specRadius(PROJ$A2) < 1){
#     A1 = PROJ$A1; A2 = PROJ$A2
#     eps1 = matrix(rnorm(p^2, sd=sqrt(sum(A1^2))/(p^2)), ncol=p)
#     eps2 = matrix(rnorm(q^2, sd=sqrt(sum(A2^2))/(q^2)), ncol=q)
#     return(list(A1=PROJ$A1+eps1,A2=PROJ$A2+eps2))
#   }
#   MAR = MAR1.LS(xx)
#   if (specRadius(MAR$A1)*specRadius(MAR$A2) < 1){
#     A1 = MAR$A1; A2 = MAR$A2
#     eps1 = matrix(rnorm(p^2, sd=sqrt(sum(A1^2))/(p^2)), ncol=p)
#     eps2 = matrix(rnorm(q^2, sd=sqrt(sum(A2^2))/(q^2)), ncol=q)
#     return(list(A1=MAR$A1+eps1,A2=MAR$A2+eps2))
#   }
#   RRMAR = MAR1.RR(xx, k1, k2)
#   if (specRadius(MAR1.RR$A1)*specRadius(MAR1.RR$A2) < 1){
#     A1 = RRMAR$A1; A2 = RRMAR$A2
#     eps1 = matrix(rnorm(p^2, sd=sqrt(sum(A1^2))/(p^2)), ncol=p)
#     eps2 = matrix(rnorm(q^2, sd=sqrt(sum(A2^2))/(q^2)), ncol=q)
#     return(list(A1=MAR1.RR$A1+eps1,A2=MAR1.RR$A2+eps2))
#   }
#   stop('causality condition of initializer fails.')
# }

initializer.sig <- function(xx){
  dim = dim(xx)[-1]
  t = dim(xx)[1]
  res = tenAR.VAR(xx, P=1)$res
  SIGMA = res %*% t(res) / (t-1)
  sig = projection(SIGMA, 1, dim[1],dim[2],dim[1],dim[2])[[1]]
  if (sig[[1]][1,1] < 0){sig[[1]] = - sig[[1]]}
  if (sig[[2]][1,1] < 0){sig[[2]] = - sig[[2]]}
  return(list(Sigl.init=sig[[1]], Sigr.init=sig[[2]]))
}


likelihood.lse <- function(fres, s, d, t){
  l1 <- fres/2/s^2
  l2 <- -(t - 1)*d*log(2*pi*s^2)/2
  return(l2 - l1)
}


IC <- function(xx,res,r,t,dim){
  N <- prod(dim)
  ic <- log(sum((res)^2)/(N*t))/2 + sum(r)*log(t)/t
  return(ic)
}


#############################################################################
LMARCvalue_orig <- function(data_sim,ksiLMAR_f,numsim,LMARinterval){
  ## Calculate critical values for LMAR
  ## Input:
  ##       (1) data_sim is the generated matrix time series: 1000 x S x 5 x 11
  ##       (2) ksiLMAR is the bound for caculating critical values
  ##       (3) numsim is the number of the simulations
  ##       (4) LMARinterval is the candidate interval
  
  ## Output: 
  ##       (1) Cvalue is the critical values
  
  K <- length(LMARinterval)
  Cvalue <- matrix(1,K,1)
  
  ## 这个初始值是需要调整的，如果一开始设得很小，则可能校准出来的临界值就是这个初始值，所以一开始最好是给初始值比较大一些，然后慢慢调整。
  #CV_start <- c(35,25,20,15,15,10,10) # 2,3,4,5, 6,7,8 原来60个月
  # cv_est <- c(1000000,892,30,21,15,13,10,9,7)
  
  
  ##### 18,24,30,36,42,48,54,60, every 6 months
  #(1) CV_start <- c(40,30,20,15,15,10,10) #(underlying length is 18 months) 2,3,4,5, 6,7,8
  #(2) CV_start <- c(40,25,20,15,15,10,10) #(underlying length is 24 months) 2,3,4,5, 6,7,8
  #(3) CV_start <- c(30,20,15,15,10,10,10) #(underlying length is 30 months) 2,3,4,5, 6,7,8
  #(3) CV_start <- c(30,20,15,15,10,10,10) #(underlying length is 36 months) 2,3,4,5, 6,7,8
  #(4) CV_start <- c(30,20,15,15,10,10,10) #(underlying length is 42 months) 2,3,4,5, 6,7,8
  #(5) CV_start <- c(30,20,15,15,15,10,10) #(underlying length is 48 months) 2,3,4,5, 6,7,8
  #(6) CV_start <- c(30,20,15,15,15,10,10) #(underlying length is 54 months) 2,3,4,5, 6,7,8
  #(7) CV_start <- c(30,25,20,15,12,10,10) #(underlying length is 60 months) 2,3,4,5, 6,7,8
  
  
  
  ##### 18,22,26,31,37,45,53,60, geometrically increase
  #(1) CV_start <- c() #(underlying length is 18 months) 2,3,4,5, 6,7,8
  #(2) CV_start <- c() #(underlying length is 22 months) 2,3,4,5, 6,7,8
  #(3) CV_start <- c(25,20,20,17,17,15,13) #(underlying length is 26 months) 2,3,4,5, 6,7,8
  #(3) CV_start <- c(25,20,20,17,17,15,13) #(underlying length is 31 months) 2,3,4,5, 6,7,8
  #(4) CV_start <- c(25,20,20,17,17,15,10) #(underlying length is 37 months) 2,3,4,5, 6,7,8
  #(5) CV_start <- c(25,20,20,17,17,15,10) #(underlying length is 45 months) 2,3,4,5, 6,7,8
  #(6) CV_start <- c(25,20,20,17,17,15,10) #(underlying length is 53 months) 2,3,4,5, 6,7,8
  #(7) CV_start <- c(25,20,20,17,17,15,10) #(underlying length is 60 months) 2,3,4,5, 6,7,8
  
  
  ##### 18,24,30,36,42,48,54,60, every 6 months, without Philippine
  #(1) CV_start <- c() #(underlying length is 18 months) 2,3,4,5, 6,7,8
  #(2) CV_start <- c() #(underlying length is 24 months) 2,3,4,5, 6,7,8
  #(3) CV_start <- c(30,22,17,15,12,10,10) #(underlying length is 30 months) 2,3,4,5, 6,7,8
  #(3) CV_start <- c(30,22,17,15,12,10,10) #(underlying length is 36 months) 2,3,4,5, 6,7,8
  #(4) CV_start <- c(30,22,17,15,12,10,10) #(underlying length is 42 months) 2,3,4,5, 6,7,8
  #(5) CV_start <- c(30,22,17,15,12,10,10) #(underlying length is 48 months) 2,3,4,5, 6,7,8
  #(6) CV_start <- c(30,22,17,15,12,10,10) #(underlying length is 54 months) 2,3,4,5, 6,7,8
  #(7) CV_start <- c(30,25,20,15,12,10,10) #(underlying length is 60 months) 2,3,4,5, 6,7,8
  
  CV_start <- c(35,22,17,15,12,10,10) #(underlying length is 60 months) 2,3,4,5, 6,7,8, without Philippine
  
  
  # CV_start <- c(25,20,20,17,17,15,10) #(underlying length is 60 months) 2,3,4,5, 6,7,8， geometrically increasing
  
  Cvalue <- Cvalue%*%1000000 
  
  dist <- 0.1
  
  for (k in c(2:K)) { # We need to caculate Cvalue(k) for k=2,...,K
    # for (k in c(2:K)) { # We need to caculate Cvalue(k) for k=2,...,K  
    Cvalue[k,1] <- CV_start[k-1]
    
    LHS <- rep(0,c(K-k+1)) # We need to test every interval from k to K
    LRcompare <- matrix(0,numsim,K-k+1)
    bound <- ksiLMAR_f%*%(k-1)/(K-1)
    
    
    while(all(LHS <= rep(bound,c(K-k+1)))){ #当所有条件都满足时继续执行；当不满足其中某一条件时停止执行
      Cvalue[k,1] <- Cvalue[k,1]-dist
      LHS <- rep(0,c(K-k+1))   # we need to test every interval from k to K
      LRcompare <- matrix(0,numsim,K-k+1)
      LRcompare_tem <- matrix(0,numsim,K-k+1)
      
      if(Cvalue[k,1] < 0){ break }
      else{ 
        
        for (q in c(k:K)) {
          numk <- LMARinterval[q]
          intervalhattem <- matrix(0,numsim,K-k+1)
          
          
          for (i in c(1:numsim)) {
            yy <- data_sim[i,,,]
            
            ##这里计算adaptive mar estimator的函数
            LMARtem  <- LMARest(yy,Cvalue,LMARinterval[1:q],i)
            
            ##这里计算adaptive log likelihood function的函数
            restem <- residuals_mar(LMARtem$A_hat_1,yy[c((c(LMARinterval[S])-numk+1):c(LMARinterval[S])),,])
            LMARhat <- loglikelihood_mle(LMARtem$SIGMA_hat_1,restem)
            intervalhattem[i,q-k+1] <- LMARtem$intervalhat_1
            
            
            LLFMLE_tem <- LLFM[q,i]
            LRcompare[i,q-k+1] <- abs(sqrt(abs(LLFMLE_tem-LMARhat)))
            
            ##### Print
            print(c(i,k,q,Cvalue[k],intervalhattem[i,q-k+1],LRcompare[i,q-k+1],bound))
          }
          LHScheck <- mean(LRcompare[,q-k+1])
          if(LHScheck > bound){ break } # 第q (q=k,...,K) 个不等式不成立时，停止下来, 减少运算量
        }
        LHS <- colMeans(LRcompare)
      }
    }
    Cvalue[k,1] <- Cvalue[k,1]+dist
  }
  return(Cvalue)
}

##############################################################################
LMARest <- function(yyy,Cvalue_f,LMARinterval_ff,ii){
  ## Local MAR estimator for ii-th simulated data 
  
  kk <- length(LMARinterval_ff)
  intervalhat <- 1
  comparek <- matrix(0,kk,1)
  
  numkk <- LMARinterval_ff[1]
  # Get the MLE on 1st interval and set it as 1st optimal one
  
  item <- (ii-1)%*%S+1  ##注意：这里要乘以S，而不是kk
  
  A_hat <- A_hatM[[item]]
  SIGMA_hat <- SIGMA_hatM[[item]]
  
  for (kkk in c(2:kk)) {
    
    item <- (ii-1)%*%S+kkk ##注意：这里要乘以S，而不是kk
    
    # A_hat <- A_hatM[[item]]
    # SIGMA_hat <- SIGMA_hatM[[item]]
    
    numkk <- LMARinterval_ff[kkk] # the number of observations put into consideration
    LLFMtem <- LLFM[kkk,ii]
    
    res_tem <- residuals_mar(A_hat,yyy[c((c(LMARinterval[S])-LMARinterval_ff[kkk]+1):c(LMARinterval[S])),,])
    LLFhat <- loglikelihood_mle(SIGMA_hat,res_tem)
    
    comparek[kkk,1] <- abs(sqrt(abs(LLFMtem-LLFhat)))
    if (comparek[kkk,1] > Cvalue_f[kkk])  { break }
    else{
      A_hat <- A_hatM[[item]]
      SIGMA_hat <-SIGMA_hatM[[item]]
      intervalhat <- kkk
    }
  }
  return(list(A_hat_1 = A_hat, SIGMA_hat_1 = SIGMA_hat, intervalhat_1=intervalhat))
}



#############################################################################
LMARCvalue <- function(data_sim,ksiLMAR_f,numsim,LMARinterval){
  ## Calculate critical values for LMAR
  ## Input:
  ##       (1) data_sim is the generated matrix time series: 1000 x S x 5 x 11
  ##       (2) ksiLMAR is the bound for caculating critical values
  ##       (3) numsim is the number of the simulations
  ##       (4) LMARinterval is the candidate interval
  
  ## Output: 
  ##       (1) Cvalue is the critical values
  
  K <- length(LMARinterval)
  Cvalue <- matrix(1,K,1)
  
  ## 这个初始值是需要调整的，如果一开始设得很小，则可能校准出来的临界值就是这个初始值，所以一开始最好是给初始值比较大一些，然后慢慢调整。
  CV_start <- c(35,25,20,15,15,10,10) # 2,3,4,5, 6,7,8
  
  Cvalue <- Cvalue%*%1000000 
  
  dist <- 1 # 0.1
  
  for (k in c(2:K)) { # We need to caculate Cvalue(k) for k=2,...,K
    # for (k in c(2:K)) { # We need to caculate Cvalue(k) for k=2,...,K  
    Cvalue[k,1] <- CV_start[k-1]
    
    LHS <- rep(0,c(K-k+1)) # We need to test every interval from k to K
    LRcompare <- matrix(0,numsim,K-k+1)
    bound <- ksiLMAR_f%*%(k-1)/(K-1)
    
    
    while(all(LHS <= rep(bound,c(K-k+1)))){ #当所有条件都满足时继续执行；当不满足其中某一条件时停止执行
      Cvalue[k,1] <- Cvalue[k,1]-dist
      LHS <- rep(0,c(K-k+1))   # we need to test every interval from k to K
      LRcompare <- matrix(0,numsim,K-k+1)
      LRcompare_tem <- matrix(0,numsim,K-k+1)
      
      if(Cvalue[k,1] < 0){ break }
      else{ 
        
        for (q in c(k:K)) {
          numk <- LMARinterval[q]
          intervalhattem <- matrix(0,numsim,K-k+1)
          
          
          registerDoParallel(50)  #注册多个核心
          result_list_tem <- list()
          
          result_list_tem <- foreach(i=1:numsim)%dopar%{
            
            yy <- data_sim[i,,,]
            
            ##这里计算adaptive mar estimator的函数
            LMARtem  <- LMARest(yy,Cvalue,LMARinterval[1:q],i)
            
            ##这里计算adaptive log likelihood function的函数
            restem <- residuals_mar(LMARtem$A_hat_1,yy[c((c(LMARinterval[S])-numk+1):c(LMARinterval[S])),,])
            LMARhat <- loglikelihood_mle(LMARtem$SIGMA_hat_1,restem)
            intervalhattem[i,q-k+1] <- LMARtem$intervalhat_1
            
            
            LLFMLE_tem <- LLFM[q,i]
            LRcompare[i,q-k+1] <- abs(sqrt(abs(LLFMLE_tem-LMARhat)))
            
            
            result_tem <- LRcompare[i,q-k+1]
            return(result_tem)
            ##### Print
            #print(c(i,k,q,Cvalue[k],intervalhattem[i,q-k+1],LRcompare[i,q-k+1],bound))
          }
          
          for (jm in c(1:numsim)) {
            LRcompare_tem[jm,q-k+1] <- result_list_tem[[jm]]
          }
          
          LHScheck <- mean(LRcompare[,q-k+1])
          if(LHScheck > bound){ break } # 第q (q=k,...,K) 个不等式不成立时，停止下来, 减少运算量
        }
        
        LHS <- colMeans(LRcompare)
        
      }
    }
    
    Cvalue[k,1] <- Cvalue[k,1]+dist
    
  }
  
  return(Cvalue)
  
}


############################################################################
LMARest_cept <- function(data_real, Cvalue, LMARinterval){
  # Get the LMAR fitting model following the sequential testing procedure given the calibrated critical values
  
  ## Input:
  
  #      (1) data_real: real AFNS components, 60 by 5 by 11 matrix
  #      (2) Cvalue: the calibrated critical values
  #      (3) LMARinterval: the candidate homogeneous intervals
  
  ## Output:
  
  #      (1)  intervalhat: the index of the optimal interval
  #      (2)  A_hat: autoregressive matrix 
  #      (3)  SIGMA_hat: a list of estimated variance-covariance matrix along each dimension
  #      (3)  compare: the loglikelihood difference between MLE estimator and adaptive estimator 
  #      (4)  Sig_hat: sample covariance matrix of the residuals 
  
  K <- length(LMARinterval)
  intervalhat <- 1
  comparek <- matrix(0,K,1)
  
  numkk <- LMARinterval[1]
  
  # Get the MLE on 1st interval and set it as 1st optimal one
  est_mle_hat <- tenAR.est(data_real[(60-numkk+1):60,,],R=1,P=1,method = "MLE",niter = 10000,tol=1e-6)
  A_hat <- est_mle_hat$A
  SIGMA_hat <- est_mle_hat$SIGMA
  Sig_hat <- est_mle_hat$Sig
  res_hat <- est_mle_hat$res
  intind <- 1
  
  for (k in c(2:K)) {
    
    num <- LMARinterval[k] # the number of observations put into consideration
    
    # local ML estimator
    est_mle_hat_tem <- tenAR.est(data_real[(60-num+1):60,,],R=1,P=1,method = "MLE",niter = 10000, tol = 1e-06) 
    A_hattem <- est_mle_hat_tem$A
    SIGMA_hattem <- est_mle_hat_tem$SIGMA
    Sig_hattem <- est_mle_hat_tem$Sig
    res_hattem <- est_mle_hat_tem$res
    res_tem <- residuals_mar(A_hattem,data_real[(60-num+1):60,,])
    LLF_tem <- loglikelihood_mle(SIGMA_hattem,res_tem)
    
    # adaptive ML estimator
    res_hat <- residuals_mar(A_hat,data_real[(60-num+1):60,,])
    LLFhat <- loglikelihood_mle(SIGMA_hat,res_hat)
    
    
    comparek[k,1] <- sqrt(abs(LLF_tem-LLFhat))
    if (abs(comparek[k,1]) > Cvalue[k])
    { break }
    else{
      intind <- intind+1
      A_hat <- A_hattem
      SIGMA_hat <-SIGMA_hattem
      Sig_hat <- Sig_hattem
      res_hat <- res_hattem 
      intervalhat <- k
    }
  }
  return(list(A_hat_1 = A_hat, SIGMA_hat_1 = SIGMA_hat, Sig_hat_1 = Sig_hat, intervalhat_1=intervalhat, compare_1 = comparek, res_hat_1 = res_hat)) 
}

