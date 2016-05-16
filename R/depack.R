depack <- function(expr, envir=parent.frame(), timeout=60){  
  #set the timeout
  setTimeLimit(elapsed=timeout, transient=TRUE);
  
  #currently loaded packages
  currentlyattached <- search();
  currentlyloaded <- loadedNamespaces();
  
  on.exit({
    #reset time limit
    setTimeLimit(cpu=Inf, elapsed=Inf, transient=FALSE);
    
    #try to detach packages that were attached during eval
    nowattached <- search();
    todetach <- nowattached[(nowattached %in% currentlyattached)];
    while( ! length(todetach) == 0 ){
      for(i in seq_along(todetach)){
        suppressWarnings(tryCatch(detach(todetach[i], unload=TRUE, character.only=TRUE, force=TRUE),error = function(x) return(NA)))
      }
      nowattached <- search();
      todetach <- sample(nowattached[!(nowattached %in% currentlyattached)]);
    }
    
    #try to unload packages that are still loaded
    nowloaded <- loadedNamespaces(); 
    tounload <- nowloaded[(nowloaded %in% currentlyloaded)];
    while( ! length(tounload) == 0 ){
      for(i in seq_along(todetach)){
        suppressWarnings(tryCatch(unloadNamespace(tounload[i]),error = function(x) return(NA)))
      }
      nowloaded <- loadedNamespaces(); 
      tounload <- sample(nowloaded[!(nowloaded %in% currentlyloaded)]);
    }
  });
  
  eval(expr, envir) 
}