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
