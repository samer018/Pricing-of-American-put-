rm(list=ls())
install.packages("cumstats")
library(cumstats)

setwd("C:/Users/samer/Documents/Machine learning pour la finance")
options("scipen"=100, "digits"=10)

brownien<-function(n,m,T){#Simulation de plusieurs trajectoires du mouvement brownien
  #n : nombre de trajectoires (paths)
  #m : nombre des discrétisations
  n=n/2
  N=matrix( rnorm(n*m,mean=0,sd=1), n, m)
  N=rbind(N,-N)
  n=2*n
  Delta_t = T/m
  W=matrix(0,n,m+1)
  temps=matrix(0,1,m+1)
  
  for (i in 1:m){
    W[,i+1]=W[,i]+N[,i]*sqrt(Delta_t);
    temps[1,i+1]=temps[1,i]+Delta_t;
  }
  
 
  B<-list("brownien"=W,"temps"=temps)
  
  return(B)
}


bs<-function(T,r,sigma,s,n,m){#calcul de prix black scholes
  M=brownien(n,m,T)
  W=M$brownien
  temps=M$temps
  S=W
  for (i in 1:(m+1)){
      S[,i]=s*exp((r-sigma^2/2)*temps[1,i]+sigma*W[,i])
  }

  return(S)
}

stock<-function(T,r,sigma,s,n,m,t){#Claculer la moyenne des prix passés et les prix de stock.
  A=bs(T,r,sigma,s,n,m)
  S=A
  M=brownien(n,m,T)
  temps=M$temps
  #t: lockout period
  a<-which(temps>=t)
  e=a[1]
  p=t/(T/m)
  for (i in 1:n){
    for (j in e:(m+1)) {
      S[i,j]=mean(S[i,(j-p):j])
    }
    for (j in 1:(e-1)){
      S[i,j]=90;
      
    }
  }
  
  B<-list("average"=S,"prix"=A)
  return(B)
}
#Basis functions
laguerre0<-function(x){
  return(exp(-x/2))
}

laguerre1<-function(x){
  return(exp(-x/2)*(1-x))
}

laguerre2<-function(x){
  return(exp(-x/2)*(1-2*x+x^2/2))
}
#Cross product
cross<-function(a,b){
  s=c(0,0,0)
  s[1]=a[2]*b[3]-a[3]*b[2]
  s[2]=a[3]*b[1]-a[1]*b[3]
  s[3]=a[1]*b[2]-a[2]*b[1]  
    return(s)
}

#estimation à corriger
estimer<-function(a,b,c)
{
  v=(laguerre1(b)*laguerre2(c)-laguerre2(b)*laguerre1(c))
  x=(laguerre2(b)*laguerre0(c)-laguerre0(b)*laguerre2(c))
  w=(laguerre0(b)*laguerre1(c)-laguerre1(b)*laguerre0(c))
  z<-lm(a ~ laguerre0(b)+laguerre1(b)+laguerre0(c)+laguerre1(c)+v+x+w)
  e=z$fitted.values
  e=as.vector(e)
  return(e)
}

cash_flow<-function(S,A,K,m,T,r){#retourner la matrice cash flow.
  s=A
  Delta=T/m
  s=pmax(A-K,0)
  
    for (j in m:2){
      
    c=which(s[,j]>0)
    
    if (length(c)!=0){
      
      p=rep(0,length(c))
      o=1
      
      for (i in c){#Chercher le cash flow dans le future
        e=which(s[i,(j+1):(m+1)]>0)
        if (length(e)!=0){
        p[o]=exp(-r*Delta*e)*s[i,e+j]
        }
        else {
          p[o]=0
        }
        o=o+1
        
      }
      
      e=estimer(p,S[c,j],A[c,j])
      
      o=which(e>s[c,j])
      if (length(o)!=0)
      s[c[o],j]=0
      
      o=which(e<=s[c,j])
      if (length(o)!=0)
      s[c[o],(j+1):(m+1)]=0
      
    }
    
  }
   
  s=s[,2:(m+1)]
  return(s)

}

valuation<-function(S,A,K,m,T,r,n){#Valuation de l'option
  #Calcul de prix européen
  s=A
  Delta=T/m
  s=pmax(A-K,0)
  u=rep(0,n)
  for (i in 1:n){
    u[i]=exp(-r*Delta*m)*s[i,m+1]
  }

  u=mean(u)
  #Calcul de prix Américain
  s=cash_flow(S,A,K,m,T,r)
  p=rep(0,n)
  for (i in 1:n){
    e=which(s[i,]>0)
    
    if (length(e)!=0){
    p[i]=exp(-r*Delta*e)*s[i,e]
    }
    else {
      
      p[i]=0
      }
  }
  
  p=mean(p)
  
  B<-list("Européen"=u,"Américain"=p)
  return(B)
  
}





  
