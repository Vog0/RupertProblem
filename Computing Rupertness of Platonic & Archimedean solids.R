### This script computes the "Rupertness" of every solved point symmetric solved Platonic and Archimedean solid

Solids2Solve=c(2:5,7:11,13:17)

NrTries=rep(0,length(Solids2Solve))
NrSolved=rep(0,length(Solids2Solve))

while (TRUE){##infinite loop
  for (Nr in 1:length(Solids2Solve)){
    loadSolid(Solids2Solve[Nr])
    for (i in 1:10000){
      P=ConvexHull(Points%*%t(RandomProjection()))
      Q=ConvexHull(Points%*%t(RandomProjection()))
      if (CanFitRotation(P,Q,bothPointSymmetric = TRUE)$fits){
        NrTries[Nr]=NrTries[Nr]+1
        NrSolved[Nr]=NrSolved[Nr]+1
      } else {
        NrTries[Nr]=NrTries[Nr]+1
      }
    }
    cat("\n")
    cat("NrTries:",NrTries,"\n")
    cat("NrSolved:",NrSolved,"\n")
    cat("ratio:",NrSolved/NrTries,"\n")
    cat("\n\n")
  }
}



