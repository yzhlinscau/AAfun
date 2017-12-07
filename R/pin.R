pin <-
function(object, formula=NULL,signif=NULL, corN=NULL,Rdf=NULL,asrV=3){
  if(is.null(signif)) signif=FALSE
  if(is.null(Rdf)) Rdf=FALSE
  
  if(!is.null(formula)){
    transform<-formula
    #if(is.null(N)) N<-0
    #ifelse(asrV==4,pframe <- as.list(object$vparameters),
           pframe <- as.list(object$gammas)#)
    names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
    tvalue<-eval(deriv(transform[[length(transform)]], names(pframe)),pframe)
    # X <- as.vector(attr(tvalue, "gradient"))
    # if(asrV==3) X[object$gammas.type == 1] <- 0
    #ifelse(asrV==4,X[object$vparameters.type == 1] <- 0,X[object$gammas.type == 1] <- 0)
    X <- matrix(as.vector(attr(tvalue, "gradient")), ncol = 1)
    
    # tname <- if(length(transform)==3){transform[[2]]}else ""
    tname <- if (length(transform) == 3) transform[[2]] else deparse(transform[[2]])
    
    #if(asrV==4) se <- as.vector(sqrt(t(X) %*% object$ai %*% X))
    #else{
      n <- length(pframe)
      i <- rep(1:n, 1:n)
      j <- sequence(1:n)
      k <- 1 + (i > j)
      Vmat <- object$ai
      se <- sqrt(sum(Vmat * X[i] * X[j] * k))
    #}

    vv=vector() #<-NULL
    vv[1]=tvalue;vv[2]=se
    
    result<-data.frame(row.names=tname, Estimate=tvalue, SE=se)
    result1<-result
    result1$sig.level<-sig.level(tvalue,se)
  
    cat("\n")
    #options(digits=3)
    if(signif==TRUE){ 
      print(result1)
      cat("---------------")
      cat("\nSig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n")    
    }else{
      if(Rdf==TRUE) print(vv) else print(result)
      }
    cat("\n")
  }
  
  if(is.null(formula)){
    
    if(is.null(corN)){cat("\nAttension: since no N value, here just show fisrt one corr!!\n\n")
                    corN<-1} 
    n<-corN
    df<-summary(object)$varcomp
    if(asrV==3){tvalue<-as.vector(df[1:n,2])
                se<-as.vector(df[1:n,3])
    }#else{
    #  tvalue<-as.vector(df[1:n,1])
    #  se<-as.vector(df[1:n,2])
    #}
    tname<-rownames(summary(object)$varcomp)[1:n]    
    siglevel<-sig.level(tvalue,se)
    

    print(data.frame(row.names=tname,Estimate=tvalue, SE=se, sig.level=siglevel))
    cat("---------------")
    cat("\nSig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n\n")    
  }
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
vpredict2=function (object, xform,signif=FALSE) 
{
  #options(digits=3)
  if (!inherits(object, "asreml")) 
    stop("Argument must be an asreml object")
  pframe <- as.list(object$vparameters)
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  tvalue <- eval(deriv(xform[[length(xform)]], names(pframe)), 
                 pframe)
  X <- matrix(as.vector(attr(tvalue, "gradient")), ncol = 1)
  tname <- if (length(xform) == 3) xform[[2]] else deparse(xform[[2]])
  se <- as.vector(sqrt(t(X) %*% object$ai %*% X))
  result<-data.frame(row.names = tname, Estimate = tvalue, SE = se)
  
  result1<-result
  result1$sig.level<-sig.level(tvalue,se)
  
  cat("\n")
  if(signif==TRUE){ 
    print(result1)
    cat("---------------")
    cat("\nSig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n")    
  }else{
    print(result)
  }
  cat("\n")

}

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
  tt=vector() #=NULL
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
