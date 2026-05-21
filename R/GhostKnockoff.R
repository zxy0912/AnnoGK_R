# modified from https://github.com/biona001/ghostknockoff-gwas-reproducibility/blob/main/chen_et_al/GhostKnockoff.R

#' GhostKnockoff Preliminary Step
#'
#' @param cor.G Matrix, the correlation (LD) matrix of the genetic variants.
#' @param M Integer, the number of knockoff copies to generate. Default is 5.
#' @param method String, currently only 'sdp' is supported for solving the SDP problem.
#' @param max.size Integer, the maximum size for the strong rules. Default is 500.
#' @param corr_max Numeric, threshold of the correlation coefficient for clustering variants. Default is 0.75.
#' @param clusters Vector, optional, manually specify cluster indices. If NULL, automatic hierarchical clustering is performed.
#' @param rep.index Vector, optional, repetition indices.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{P.each}: Projection matrix for each variant.
#'   \item \code{V.left}: Cholesky decomposition of the covariance matrix for knockoff generation.
#'   \item \code{Normal_50Studies}: Matrix of randomly sampled normal variables for simulation.
#'   \item \code{M}: Number of knockoff copies.
#'   \item \code{A}: The covariance structure matrix.
#'   \item \code{A.left}: Cholesky decomposition of matrix A.
#'   \item \code{clusters}: Clustering assignments.
#'   \item \code{diag_s}: Diagonal matrix s obtained from SDP solver.
#'   \item \code{Sigma}: The input correlation matrix.
#' }
#' @importFrom Matrix Matrix forceSymmetric
#' @importFrom stats rnorm as.dist hclust cutree
#' @export
GhostKnockoff.prelim<-function(cor.G, M=5, method='sdp', max.size=500, corr_max=0.75, clusters=NULL, rep.index=NULL){
  n.G<-nrow(cor.G)
  Normal_50Studies<-matrix(rnorm(n.G*M*50),n.G*M,50)
  P.each<-matrix(0,n.G,n.G)
  
  eigen.fit<-eigen(cor.G)
  newEig <- ifelse(eigen.fit$values < 1e-5, 1e-5, eigen.fit$values)
  newMat <- eigen.fit$vectors %*% (newEig*t(eigen.fit$vectors))
  # normalize modified matrix eqn 6 from Brissette et al 2007
  newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
  cor.G<-newMat
  cor.G<-as.matrix(forceSymmetric(cor.G))
  Sigma<-cor.G
  SigmaInv<-solve(Sigma)
  
  if(length(clusters)==0){
    #clustering to identify tightly linked variants, 
    #first apply hierarchical clustering to determine number of clusters
    corr_max<-corr_max
    Sigma.distance = as.dist(1 - abs(Sigma))
    if(ncol(Sigma)>1){
      fit = hclust(Sigma.distance, method="single")
      clusters = cutree(fit, h=1-corr_max)
    }else{clusters<-1}
  }
  
  if(method != 'sdp') {
    warning("Only 'sdp' method is supported in this package. Forcing method='sdp'.")
  }

  s<-create.solve_sdp_M(Sigma,M=M)
  diag_s<-Matrix(diag(s,length(s)))
  clusters=NULL
  
  if(sum(diag_s)==0){
    V.left<-matrix(0,n.G*M,n.G*M)
  }else{
    #Sigma_k<-2*diag_s - s*t(s*SigmaInv)
    Sigma_k<-2*diag_s - diag_s%*%SigmaInv%*%diag_s
    
    V.each<-Matrix(forceSymmetric(Sigma_k-diag_s))
    
    #random part of knockoff
    V<-matrix(1,M,M)%x%V.each+diag(1,M)%x%diag_s
    V<-Matrix(forceSymmetric(V))
    
    #diag(V)<-diag(V)+rep(s,M)
    V.left<-try(t(chol(V)),silent=T)
    if(inherits(V.left, "try-error")){
      eigen.fit<-eigen(V)
      newEig <- ifelse(eigen.fit$values < 1e-5, 1e-5, eigen.fit$values)
      newMat <- eigen.fit$vectors %*% (newEig*t(eigen.fit$vectors))
      # normalize modified matrix eqn 6 from Brissette et al 2007
      newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
      V<-newMat
      V.left<-t(chol(V))
    }
  }
  #P.each<-diag(1,length(s))-s*SigmaInv
  P.each<-diag(1,n.G)-diag_s%*%SigmaInv
  Normal_50Studies<-as.matrix(V.left%*%matrix(rnorm(ncol(V.left)*50),ncol(V.left),50))
  
  #is_posdef((M+1)/M*Sigma-diag_s)
  A.each<-Sigma-diag_s
  A<-matrix(1,M+1,M+1)%x%A.each+diag(1,M+1)%x%diag_s
  #is_posdef(A)
  A<-Matrix(forceSymmetric(A))
  
  
  A.left<-try(t(chol(A)),silent=TRUE)
  if(class(A.left)[1]=="try-error"){
    #alpha<-0.05-min(eigen(A)$values)
    #A<-(1-alpha)*A+alpha*diag(nrow(A))
    eigen.fit<-eigen(A)
    newEig <- ifelse(eigen.fit$values < max(1e-5,2*min(abs(eigen.fit$values))), max(1e-5,2*min(abs(eigen.fit$values))), eigen.fit$values)
    newMat <- eigen.fit$vectors %*% (newEig*t(eigen.fit$vectors))
    # normalize modified matrix eqn 6 from Brissette et al 2007
    newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
    A<-newMat
    #A.left<-chol(newMat)
    #svd.fit<-svd.A
    #u<-svd.fit$u
    #svd.fit$d[is.na(svd.fit$d)]<-0
    #cump<-cumsum(svd.fit$d)/sum(svd.fit$d)
    #n.svd<-which(cump>=0.999)[1]
    #if(is.na(n.svd)){n.svd<-nrow(A)}
    #svd.index<-intersect(1:n.svd,which(svd.fit$d!=0))
    #A.left<-t(sqrt(svd.fit$d[svd.index])*t(u[,svd.index,drop=F]))
  }
  A.left<-as.matrix(t(chol(A)))
  #svd.A<-svd(A)
  #svd.A<-svd(A)
  
  return(list(P.each=as.matrix(P.each), V.left=V.left, Normal_50Studies=as.matrix(Normal_50Studies), M=M, A=A, A.left=A.left, clusters=clusters, diag_s=diag_s, Sigma=Sigma))
}


#' GhostKnockoff Fitting
#'
#' @param Zscore_0 Vector or Matrix, observed Z-scores from the original study.
#' @param N.effect Vector, sample sizes of the original study.
#' @param fit.prelim List, the output from GhostKnockoff.prelim.
#' @param method String, selection method for feature importance scores. Options: 'marginal', 'lasso', or 'lasso.approx.lambda'.
#' @param type String, choice of FDR control type. Options: 'fdr' or 'fwer'.
#' @param M.fwer Integer, number of copies used if type is 'fwer'. Default is 50.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{T_0}: Feature importance statistics for original variables.
#'   \item \code{T_k}: Matrix of statistics for knockoff copies.
#'   \item \code{kappa}: Calculated kappa statistic.
#'   \item \code{tau}: Calculated tau statistic.
#'   \item \code{temp.A}: Adjusted covariance matrix used in fitting.
#'   \item \code{r_all}: Concatenated vector of statistics.
#'   \item \code{lambda.seq}: Lambda sequence used in the Lasso path.
#'   \item \code{lambda}: The finally selected lambda value.
#' }
#' @importFrom stats rnorm
#' @importFrom ghostbasil ghostbasil
#' @export
GhostKnockoff.fit<-function(Zscore_0, N.effect, fit.prelim, method='susie',type='fdr',M.fwer=50){
  Zscore_0<-as.matrix(Zscore_0)
  N.effect<-as.vector(N.effect)
  Zscore_0[is.na(Zscore_0)]<-0
  N.effect[is.na(N.effect)]<-Inf
  
  M<-fit.prelim$M
  n.G<-nrow(Zscore_0)
  P.each<-fit.prelim$P.each
  Normal_50Studies<-fit.prelim$Normal_50Studies
  A<-as.matrix(fit.prelim$A)
  A.left<-fit.prelim$A.left
  V.left<-fit.prelim$V.left
  
  if(type=='fdr'){M.rep<-1}
  if(type=='fwer'){M.rep<-M.fwer}
  
  T_0<-list();T_k<-list()
  kappa<-list();tau<-list()
  for(m in 1:M.rep){
    #Normal_k<-matrix(Normal_50Studies[,m],nrow=n.G)
    Normal_k<-matrix(V.left%*%matrix(rnorm(ncol(V.left)),ncol(V.left),1),nrow=n.G)
    GK.Zscore_0<-Zscore_0
    GK.Zscore_k<-as.vector(P.each%*%GK.Zscore_0)+Normal_k
    
    if(method=='marginal'){
      T_0[[m]]<-(GK.Zscore_0)^2
      T_k[[m]]<-(GK.Zscore_k)^2
    }
    if(method=='lasso'){
      #calculate importance score
      r<-GK.Zscore_0/sqrt(N.effect)#sqrt(N.effect-1+GK.Zscore_0^2)
      r_k<-as.vector(GK.Zscore_k/sqrt(N.effect))#sqrt(N.effect-1+GK.Zscore_k^2))
      r_all<-as.matrix(c(r,r_k))
      
      nfold<-5
      nA<-N.effect*(nfold-1)/nfold;nB<-N.effect/nfold
      temp.left<-sqrt(nB/nA/N.effect)*as.matrix(A.left)
      r_all_A<-r_all+as.matrix(temp.left%*%matrix(rnorm(ncol(temp.left)),ncol(temp.left),1))
      r_all_B<-(r_all*N.effect-r_all_A*nA)/nB
      shrink=0.01#seq(0.05,1,0.05)
      beta.all<-c();parameter.set<-c()
      k<-1
      #temp.A<-(1-shrink[k])*A+diag(shrink[k],nrow(A))
      temp.A<-A+diag(shrink[k],nrow(A))
      fit.basil<-try(ghostbasil(temp.A, r_all_A, alpha=1, delta.strong.size = max(1,min(500,length(r_all_A)/20)), max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=FALSE),silent=T)
      parameter.set<-rbind(parameter.set,cbind(fit.basil$lmdas,shrink[k]))
      beta.all<-cbind(beta.all,fit.basil$betas)
      
      Get.f<-function(x){x<-as.matrix(x);return(t(x)%*%r_all_B/sqrt(t(x)%*%temp.A%*%x))}
      f.lambda<-apply(beta.all,2,Get.f)
      f.lambda[is.na(f.lambda)]<--Inf
      #beta<-beta.all[,which.max(f.lambda)]
      parameter<-parameter.set[which.max(f.lambda),]
      temp.A<-(1-parameter[2])*A+diag(parameter[2],nrow(A))
      
      lambda.all<-fit.basil$lmdas
      lambda<-fit.basil$lmdas[which.max(f.lambda)]
      lambda.seq <- lambda.all[lambda.all > lambda]
      lambda.seq <- c(lambda.seq, lambda)
      
      fit.basil<-ghostbasil(temp.A, r_all,user.lambdas=lambda.seq, alpha=1, delta.strong.size = max(1,min(500,length(r_all_A)/20)), max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=FALSE)
      beta<-fit.basil$betas[,ncol(fit.basil$betas)]
      
      T_0[[m]]<-abs(beta[1:n.G])
      T_k[[m]]<-abs(matrix(beta[-(1:n.G)],n.G,M))
    }
    if(method=='lasso.approx.lambda'){
      #calculate importance score
      r<-GK.Zscore_0/sqrt(N.effect)#sqrt(N.effect-1+GK.Zscore_0^2)
      r_k<-as.vector(GK.Zscore_k/sqrt(N.effect))#sqrt(N.effect-1+GK.Zscore_k^2))
      r_all<-as.matrix(c(r,r_k))
      
      N_eff <- as.numeric(N.effect[1])
      lambda_max<-max(abs(rnorm(length(r_all))))/sqrt(as.numeric(N_eff))
      epsilon <- .0001
      K <- 100
      lambda.all <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
                                  length.out = K)), digits = 10)
      lambda<-lambda_max*0.6
      lambda.seq <- lambda.all[lambda.all > lambda]
      lambda.seq <- c(lambda.seq, lambda)
      
      temp.A<-A+0.01*diag(1,nrow(A))
      fit.basil<-ghostbasil(temp.A, r_all,user.lambdas=lambda.seq, alpha=1, delta.strong.size = max(1,min(500,length(r_all)/20)), max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=F)
      beta<-fit.basil$betas[,ncol(fit.basil$betas)]
      
      T_0[[m]]<-abs(beta[1:n.G])
      T_k[[m]]<-abs(matrix(beta[-(1:n.G)],n.G,M))
    }
    MK.stat<-MK.statistic(T_0[[m]],T_k[[m]])
    kappa[[m]]<-MK.stat[,'kappa']
    tau[[m]]<-MK.stat[,'tau']
  }
  
  return(list(T_0=T_0,T_k=T_k,kappa=kappa,tau=tau,temp.A=temp.A,r_all=r_all,lambda.seq=lambda.seq,lambda = lambda ))
}


#' GhostKnockoff Filtering
#'
#' @description Computes knockoff statistics (\eqn{\kappa} and \eqn{\tau}) and performs 
#' FDR control to identify significant features.
#'
#' @param T_0 Vector or Matrix, feature importance statistics for the original variables.
#' @param T_k Matrix of size (p x M), feature importance statistics for the knockoff copies.
#' @param clusters Vector of size p, indicating cluster assignments for features. 
#' Defaults to 1:length(T_0).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{kappa}: Vector of signs indicating whether the original variable (0) or 
#'   which knockoff copy (1 to M) performed best.
#'   \item \code{tau}: Vector of differences between the max statistic and the 
#'   median of the other statistics.
#'   \item \code{q}: Vector of estimated q-values for FDR control.
#' }
#' @importFrom stats median
#' @export
GhostKnockoff.filter<-function (T_0,T_k,clusters=1:length(T_0)){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  M<-ncol(T_k);Rej.Bound<-10000
  
  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0
  
  which.max.alt<-function(x){
    temp.index<-which(x==max(x))
    if(length(temp.index)!=1){return(temp.index[2])}else{return(temp.index[1])}
  }
  kappa<-apply(T.temp,1,which.max.alt)-1
  
  Get.OtherMedian<-function(x){median(x[-which.max(x)])}
  tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    G.factor<-max(table(clusters[b][1:i]))
    temp_ratio<-(1/M*G.factor+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau));
  if(length(which(tau[b]>0))!=0){
    index_bound<-max(which(tau[b]>0))
    for(i in 1:length(b)){
      temp.index<-i:min(length(b),Rej.Bound,index_bound)
      if(length(temp.index)==0){next}
      q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
      if(i>Rej.Bound){break}
    }
    q[q>1]<-1
  }
  
  return(list(kappa=kappa,tau=tau,q=q))
}



#' @importFrom stats cov2cor
#' @noRd
create.solve_sdp_M <- function(Sigma, M=1, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  # Convert problem for SCS
  
  # Linear constraints
  Cl1 = rep(0,p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1,p)
  Al2 = Matrix::Diagonal(p)
  
  # Positive-definite cone
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x=d_As)
  As = As[which(Matrix::rowSums(As) > 0),]
  Cs = c((M+1)/M*G) ##change from 2 to (M+1)/M
  
  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$s=p
  K$l=2*p #not sure if it should be changed - may be not as it is the dimention of the linear part.
  
  # Objective
  b = rep(1,p)
  
  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  if(verbose) cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  if(verbose) cat("done. \n")
  
  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution. Knockoffs may lose power.')
  }
  
  # Clip solution to correct numerical errors (domain)
  s = sol$y
  s[s<0]=0
  s[s>1]=1
  
  # Compensate for numerical errors (feasibility)
  if(verbose) cat("Verifying that the solution is correct ... ")
  psd = 0
  s_eps = 1e-8
  while ((psd==0) & (s_eps<=0.1)) {
    if (is_posdef((M+1)/M*G-diag(s*(1-s_eps),length(s)),tol=1e-9)) { ##change from 2 to (M+1)/M
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  s[s<0]=0
  if(verbose) cat("done. \n")
  
  # Verify that the solution is correct
  if (all(s==0)) {
    warning('In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
  }
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}


#' @importFrom RSpectra eigs
#' @noRd
is_posdef = function(A, tol=1e-9) {
  p = nrow(matrix(A))
  
  if (p<500) {
    lambda_min = min(eigen(A,only.values=T)$values)
  }
  else {
    oldw <- getOption("warn")
    #options(warn = -1)
    lambda_min = suppressWarnings(RSpectra::eigs(A, 1, which="SM", opts=list(retvec = FALSE, maxitr=100, tol))$values)
    options(warn = oldw)
    if( length(lambda_min)==0) {
      # RSpectra::eigs did not converge. Using eigen instead."
      lambda_min = min(eigen(A,only.values=T)$values)
    }
  }
  return (lambda_min>tol*10)
}


#' @importFrom stats median
#' @noRd
MK.statistic<-function (T_0,T_k,method='median'){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0
  
  which.max.alt<-function(x){
    temp.index<-which(x==max(x))
    if(length(temp.index)!=1){return(temp.index[2])}else{return(temp.index[1])}
  }
  kappa<-apply(T.temp,1,which.max.alt)-1
  
  if(method=='max'){tau<-apply(T.temp,1,max)-apply(T.temp,1,max.nth,n=2)}
  if(method=='median'){
    Get.OtherMedian<-function(x){median(x[-which.max(x)])}
    tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  }
  return(cbind(kappa,tau))
}



