pin <-
function(object, formula=NULL,signif=NULL, corN=NULL,Rdf=NULL){
  if(is.null(signif)) signif=FALSE
  if(is.null(Rdf)) Rdf=FALSE  
  
  if(!is.null(formula)){
    transform<-formula
    #if(is.null(N)) N<-0
    pframe <- as.list(object$gammas)
    names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
    tvalue<-eval(deriv(transform[[length(transform)]], names(pframe)),pframe)
    X <- as.vector(attr(tvalue, "gradient"))
    X[object$gammas.type == 1] <- 0
    tname <- if(length(transform)==3){transform[[2]]}else ""
    n <- length(pframe)
    i <- rep(1:n, 1:n)
    j <- sequence(1:n)
    k <- 1 + (i > j)
    Vmat <- object$ai
    se <- sqrt(sum(Vmat * X[i] * X[j] * k))
    
    vv=vector() #NULL
    vv[1]=tvalue;vv[2]=se
    
    result<-data.frame(row.names=tname, Estimate=tvalue, SE=se)
    result1<-result
    result1$sig.level<-sig.level(tvalue,se)
  
    cat("\n")
    #options(digits=3)
    if(signif==TRUE){ 
      print(format(result1, nsmall=3))
      cat("---------------")
      cat("\nSig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n")    
    }else{
      if(Rdf==TRUE) print(format(vv, nsmall=3)) else print(format(result, nsmall=3))
      }
    cat("\n")
  }
  
  if(is.null(formula)){
    
    if(is.null(corN)){cat("\nAttension: since no N value, here just show fisrt one corr!!\n\n")
                    corN<-1} 
    n<-corN
    df<-summary(object)$varcomp
    tvalue<-as.vector(df[1:n,2])
    se<-as.vector(df[1:n,3])
    tname<-rownames(summary(object)$varcomp)[1:n]    
    siglevel<-sig.level(tvalue,se)
    
    result2=data.frame(row.names=tname,Estimate=tvalue, SE=se, sig.level=siglevel)
    print(format(result2, nsmall=3))
    cat("---------------")
    cat("\nSig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n\n")    
  }
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# sig.level functions

sig.level<-function(tvalue,se,...){
  n<-length(se)
  siglevel<-1:n
  for(i in 1:n){    
    sig.se<-c(se[i]*1.450,se[i]*2.326,se[i]*3.090)  
    
    if(abs(tvalue[i])>sig.se[3]) {siglevel[i]<-"***"}
     else if(abs(tvalue[i])>sig.se[2]) {siglevel[i]<-"**"}
     else if(abs(tvalue[i])>sig.se[1]) {siglevel[i]<-"*"}
     else {siglevel[i]<-"Not signif"}
  }
  siglevel
}

sig.level2=function(x){
  tt=vector() #NULL
  n=length(x)
  for(i in 1:n){
    if(abs(x[i])<0.001) tt[i]='***'
    else if(abs(x[i])<0.01) tt[i]='**'
    else if(abs(x[i])<0.05) tt[i]='*'
    else if(abs(x[i])<0.10) tt[i]='.' 
    else tt[i]=''
  }
  return(tt)
}

  
 # ----------------------------------------------------------------------------
# other functions for chapter 6
# ----------------------------------------------------------------------------

## for 6.4
pcp=function (PCA, nfactor, plot = FALSE) 
{
  m=nfactor
  W = as.matrix(PCA[[1]]^2/sum(PCA[[1]]^2))
  PCs = as.matrix(PCA$scores[, 1:m])
  PC = PCs %*% W[1:m]/sum(W[1:m])
  ans = cbind(PCs, PC = PC[, 1], rank = rank(PC[, 1]))
  if (plot) {
    plot(PCs)
    abline(h = 0, v = 0, lty = 3)
    text(PCs, label = rownames(PCs), pos = 1.1, adj = 0.5, 
         cex = 0.85)
  }
  return(ans)
}

## for 6.7.3 
discrim.dist2=function (X, G, newX = NULL, var.equal = FALSE) 
{
  disp = newX
  G = as.factor(G)
  if (is.null(newX) == TRUE) 
    newX <- X
  if (is.vector(newX) == TRUE) 
    newX <- t(as.matrix(newX))
  else if (is.matrix(newX) != TRUE) 
    newX <- as.matrix(newX)
  if (is.matrix(X) != TRUE) 
    X <- as.matrix(X)
  nx <- nrow(newX)
  newG <- matrix(rep(0, nx), nrow = 1, dimnames = list("newG", 
                                                       1:nx))
  g <- length(levels(G))
  mu <- matrix(0, nrow = g, ncol = ncol(X))
  for (i in 1:g) mu[i, ] <- colMeans(X[G == i, ])
  D <- matrix(0, nrow = g, ncol = nx)
  if (var.equal == TRUE || var.equal == T) {
    for (i in 1:g) D[i, ] <- mahalanobis(newX, mu[i, ], var(X))
  }
  else {
    for (i in 1:g) D[i, ] <- mahalanobis(newX, mu[i, ], var(X[G == i, ]))
  }
  for (j in 1:nx) {
    dmin <- Inf
    for (i in 1:g) if (D[i, j] < dmin) {
      dmin <- D[i, j]
      newG[j] <- i
    }
  }
  if (is.null(disp) == FALSE) 
    list(Dist = D, newG = newG)
  else {dd=data.frame(G = G, D = t(D), newG = t(newG))
  print(dd)
  cat("\n")
  tab=table(G=dd[,1],NewG=dd[,5])
  print(tab);cat("\nDiscriminant coincidence rate:\n")
  return(sum(diag(prop.table(tab))))
  }
}

## for 6.10
CI_CR=function (B) 
{
  RI = c(0, 0, 0.58, 0.9, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49, 
         1.51)
  Wi = weight(B)
  n = length(Wi)
  if (n > 2) {
    W = matrix(Wi, ncol = 1)
    A = matrix(B, nrow = sqrt(length(B)), ncol = sqrt(length(B)), 
               byrow = TRUE)
    AW = A %*% W
    aw = as.vector(AW)
    la_max = sum(aw/Wi)/n
    CI = (la_max - n)/(n - 1)
    CR = CI/RI[n]
    cat("\n CI=", round(CI, 4), "\n")
    cat("\n CR=", round(CR, 4), "\n")
    cat("\n la_max=", round(la_max, 4), "\n\n")
    if (CR <= 0.1) {
      cat(" Consistency test is OK！\n")
      cat("\n Wi: ", round(Wi, 4), "\n")
    }
    else {
      cat(" Please adjust the judgment matrix! \n")
      Wi = null
      break
    }
  }
  else if (n <= 2) {
    return(Wi)
  }
}

z_score=function (B, converse = FALSE) 
{
  B = as.vector(B)
  if (converse == FALSE || converse == F || converse == "") {
    min_value = min(B)
    max_value = max(B)
    z_score = (B - min_value)/(max_value - min_value) * 60 + 
      40
    z_score
  }
  else if (converse == TRUE || converse == T) {
    min_value = min(B)
    max_value = max(B)
    z_score = (max_value - B)/(max_value - min_value) * 60 + 
      40
    z_score
  }
}

z_data=function (data, converse = FALSE) 
{
  n = ncol(data)
  m = nrow(data)
  score_array = array(1:(m * n), c(m, n))
  for (i in 1:n) {
    score_array[, i] = z_score(data[, i], converse)
  }
  SCORE = as.matrix(score_array)
  dimnames(SCORE)[1] = dimnames(data)[1]
  dimnames(SCORE)[2] = dimnames(data)[2]
  round(SCORE, 4)
}

weight=function (B) 
{
  A = matrix(B, nrow = sqrt(length(B)), ncol = sqrt(length(B)), 
             byrow = TRUE)
  n = ncol(A)
  mul_collect = c(1:n)
  for (i in 1:n) mul_collect[i] = prod(A[i, ])
  weight = mul_collect^(1/n)
  weight_one = weight/sum(weight)
  round(weight_one, 4)
}

S_rank=function (data, Wi) 
{
  wight_matrix = matrix(Wi, ncol = 1, byrow = FALSE)
  score_matrix = as.matrix(data)
  Si = score_matrix %*% wight_matrix
  print(data.frame(Si = Si, ri = rank(-Si)))
  list(Si = Si)
}

# ----------------------------------------------------------------------------
# other functions for 5.7.3
# ----------------------------------------------------------------------------

relweights <- function(fit,...){                         
  R <- cor(fit$model)   
  nvar <- ncol(R)          
  rxx <- R[2:nvar, 2:nvar] 
  rxy <- R[2:nvar, 1]      
  svd <- eigen(rxx)        
  evec <- svd$vectors                           
  ev <- svd$values         
  delta <- diag(sqrt(ev))  
  lambda <- evec %*% delta %*% t(evec)        
  lambdasq <- lambda ^ 2   
  beta <- solve(lambda) %*% rxy           
  rsquare <- colSums(beta ^ 2)                   
  rawwgt <- lambdasq %*% beta ^ 2    
  import <- round((rawwgt / rsquare) * 100,2)  
  lbls <- names(fit$model[2:nvar])   
  rownames(import) <- lbls
  colnames(import) <- "Weights"
  ymax<-max(import)+20
  barplot(t(import),names.arg=lbls,
          ylab="% of R-Square",
          xlab="Predictor Variables",
          ylim=c(0,ymax),
          main=list("Relative Importance of Predictor Variables", cex = 1,
                    col = "black"),
          sub=paste("R-Square=", round(rsquare, digits=3)),
          ...) 
  return(import)
}
 
# ----------------------------------------------------------------------------
# other functions for 5.8.1.3
# ---------------------------------------------------------------------------- 
cancor2 <- function(x, y, dec = 4){
  x <- as.matrix(x); y <- as.matrix(y)
  n<-dim(x)[1]; q1<-dim(x)[2]; q2<-dim(y)[2]; q<-min(q1,q2)
  S11 <- cov(x); S12 <- cov(x, y); S21 <- t(S12); S22 <- cov(y)
  E1 <- eigen(solve(S11)%*%S12%*%solve(S22)%*%S21)
  E2 <- eigen(solve(S22)%*%S21%*%solve(S11)%*%S12)
  rsquared <-E1$values[1:q] #as.real()
  LR=pp=qq=tt <- NULL
  for ( i in 1:q ){
    LR <- c(LR, prod(1-rsquared[i:q]))
    pp <- c(pp, q1-i+1)
    qq <- c(qq, q2-i+1)
    tt <- c(tt, n-1-i+1)}
  m <- tt-0.5*(pp+qq+1); lambda <- (1/4)*(pp*qq-2); 
  s <- sqrt((pp^2*qq^2-4)/(pp^2+qq^2-5))
  F <- ((m*s-2*lambda)/(pp*qq))*((1-LR^(1/s))/LR^(1/s)); 
  df1 <- pp*qq; df2 <-(m*s-2*lambda);
  pval <- 1-pf(F, df1, df2)
  outmat<-round(cbind(sqrt(rsquared),rsquared,LR,
                      F,df1,df2, pval),dec)
  colnames(outmat) = list("R", "RSquared", "LR", "ApproxF",
                          "NumDF", "DenDF", "pvalue")
  rownames(outmat) = as.character(1:q);
  xrels<-round(cor(x,x%*%E1$vectors)[,1:q],dec)
  colnames(xrels) <- apply( cbind(rep("U", q), as.character(1:q)),
                            1, paste, collapse ="" )
  yrels <- round(cor(y, y%*%E2$vectors)[ ,1:q], dec)
  colnames(yrels) <- apply( cbind(rep("V", q), as.character(1:q)), 
                            1, paste, collapse = "")
  list( Summary = outmat, a.Coefficients = E1$vectors,
        b.Coefficients = E2$vectors,
        XUCorrelations = xrels, YVCorrelations = yrels )
}   
  
cc2 <- function (df1, df2) {
  require(CCA)
  df.ca=cc(df1, df2)
  
  # tests of canonical dimensions
  ev <- (1 - df.ca$cor^2)
  
  n <- dim(df1)[1]; p <- length(df1)
  q <- length(df2); k <- min(p, q)
  m <- n - 3/2 - (p + q)/2;  w <- rev(cumprod(rev(ev)))
  
  # initialize
  d1 <- d2 <- f <- vector("numeric", k)
  
  for (i in 1:k) {
    s <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
    si <- 1/s
    d1[i] <- p * q; d2[i] <- m * s - p * q/2 + 1
    r <- (1 - w[i]^si)/w[i]^si
    f[i] <- r * d2[i]/d1[i]
    p <- p - 1; q <- q - 1
  }
  
  pv <- pf(f, d1, d2, lower.tail = FALSE)
  dmat <- cbind(WilksL = w, F = f, df1 = d1, df2 = d2, p = pv)
  result=list(summary=df.ca,Ftest=dmat)
  return(result)
}  
