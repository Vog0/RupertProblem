#######################
### Solving the TID ###
#######################

## This script solves the TID (Truncated icosidodecahedron) and optimizes its Nieuwland number

loadSolid(17) #The global variables Points and Edges are created

sol=SolveRupert(Points) #a solution of the TID
if (sol$Solution){print(sol);cat("The TID was solved :-)\n")}

NieuwlandNumber(Points, sol) #the Nieuwland number of sol

sol=ImproveNieuwland(Points,sol,NrTries=4*10^4)

### computing the new Nieuwland Number:
NieuwlandNumber(Points,sol)
PlotSolution(Points,sol)
