##
## demo file for met.plot. 
##

library(AAfun)
require(grid)
data(MET)

# example 1
# variable order: genotype,yield,site,row,col
MET2<-MET[,c(1,9,2,4:5)] 
met.plot(MET2)

# example 2
MET3<-MET[,c(1,9,2,4:7)] # add variable order on MET2: Rep, Block
met.plot(MET3,"My met trials")

