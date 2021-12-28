### This script solves Platonic & Archimedean solids 
### Then it tries to find the best possible Nieuwland numbers for each solved polyhedron 
### The best parameters are stored in the variable "BestSolutions"
### The best Nieuwland numbers are stored in the variable "nieuwlands"

Solids2Solve=c(1:11,13:15,17,20:24,26:28,30)

BestSolutions=list()
nieuwlands=rep(0,length(Solids2Solve))

while (TRUE){##infinite loop
  for (Nr in 1:length(Solids2Solve)){
    loadSolid(Solids2Solve[Nr])
    sol=SolveRupert(Points,wantPlot = F)
    if (sol$Solution){
      sol=ImproveNieuwland(Points,sol,NrTries=3*10^5,wantPrint=T,NrTimeSinceLastImprovment = 30000)
      if (NieuwlandNumber(Points,sol)>nieuwlands[Nr]){
        nieuwlands[Nr]=NieuwlandNumber(Points,sol)
        BestSolutions[[Nr]]=sol
        print(nieuwlands)
      }
    }
  }
}

dat = data.frame(SolidsNames = Solids2Solve)

dat$x = unlist(lapply(BestSolutions,function(x) x$vector[1]))
dat$y = unlist(lapply(BestSolutions,function(x) x$vector[2]))
dat$alpha = unlist(lapply(BestSolutions,function(x) x$alpha))
dat$theta1 = unlist(lapply(BestSolutions,function(x) x$ThetaPhiP[1]))
dat$phi1 = unlist(lapply(BestSolutions,function(x) x$ThetaPhiP[2]))
dat$theta2 = unlist(lapply(BestSolutions,function(x) x$ThetaPhiQ[1]))
dat$phi2 = unlist(lapply(BestSolutions,function(x) x$ThetaPhiQ[2]))
dat$N = nieuwlands

dat
