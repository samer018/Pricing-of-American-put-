
setwd("C:/Users/samer/Documents/Machine learning pour la finance")
options("scipen"=100, "digits"=10)

#Sourcing functions
source("brownien.R")
for (i in seq(from=90,to=120,by=2))
set.seed(1)
z=stock(2,0.06,0.2,i,50000,200,0.25)
S=z$prix
A=z$average

s=valuation(S,A,100,200,2,0.06,50000)
#o=cash_flow(S,A,100,200,2,0.06)
p=s$Américain




