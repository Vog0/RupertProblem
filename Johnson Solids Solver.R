solved=rep(FALSE,92)

while(TRUE){
  for (n in which(!solved)){
    #if (any(n==c(15,37,42,44,45,57,63,64,89,92))){next} # Johnson solids we could not solve
    cat("now @",n,"\n")
    loadSolid(31+n)
    isPS=IsPolyhedronPointSymmetric(Points)
    if (isPS){solved[n]=SolveRupert(Points,NrCycles = 5)$Solution}
    if (!isPS){solved[n]=SolveRupert(Points,NrCycles = 5,NrProjections = 2^10)$Solution}
    
    cat(sum(!solved),"unsolved:",which(!solved),"\n\n\n")
  }
}




