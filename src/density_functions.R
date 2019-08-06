# Gauss Hypergeometric density with two parameters
dgaushyp2<-function(u,a,b1,b2,c,y1,y2){
  f<-(gamma(c)/(gamma(a)*gamma(c-a)))*((u^(a-1))*((1-u)^(c-a-1))*((1-u*y1)^(-b1))*((1-u*y2)^-(b2)))
  return(f)
}
# Gauss Hypergeometric density with one parameter
dgaushyp<-function(u,a,b,c,y){
  f<-(1/beta(a,c-a))*((u^(a-1))*((1-u)^(c-a-1))*((1-u*y)^(-b)))
  return(f)
}
# Density of the difference of two independent beta random variables
# the sum fo alpha (as well as the sum of the betas) should be greater
# than 1
ddifbeta<-function(x,alpha_1,beta_1,alpha_2,beta_2){
  if(x>0&x<=1){
    gaus<-integrate(dgaushyp2,lower = 0,upper=1,
                    a=beta_1,
                    b1=alpha_1+beta_1+alpha_2+beta_2-2,
                    b2=1-alpha_1,
                    c=beta_1+alpha_2,
                    y1=1-x,
                    y2=1-(x^2))
    fx<-beta(alpha_2,beta_1)*x^(beta_1+beta_2-1)*((1-x)^(alpha_2+beta_1-1))*gaus$value*(1/(beta(alpha_1,beta_1)*beta(alpha_2,beta_2)))
  }
  else if(x==0){
    if((alpha_1+alpha_2)>0&(beta_1+beta_2)>0){
      fx<-beta(alpha_1+alpha_2-1,beta_1+beta_2-1)/(beta(alpha_1,beta_1)*beta(alpha_2,beta_2))
    }
    else{
      fx<-NA
    }
  }
  else if(x>=-1&x<0){
    gaus<-integrate(dgaushyp2,lower = 0,upper=1,
                    a=beta_2,
                    b1=1-alpha_2,
                    b2=alpha_1+alpha_2+beta_1+beta_2-2,
                    c=alpha_1+beta_2,
                    y1=1-(x^2),
                    y2=1+x)
    fx<-(beta(alpha_1,beta_2))*((-x)^(beta_1+beta_2-1))*((1+x)^(alpha_1+beta_2-1))*(gaus$value)*(1/(beta(alpha_1,beta_1)*beta(alpha_2,beta_2)))
  }
  else{
    stop('x should be between -1 and 1')
  }
  return(fx)
}
# Density of the ratio of two independent beta random variables 
# the sum fo alpha (as well as the sum of the betas) should be greater
# than 1
dratiobeta<-function(x,alpha_1,beta_1,alpha_2,beta_2){
  if(x<0){
    stop("x must be greater than 0")
  }
  else if(x>0&x<=1){
    gaus<-integrate(dgaushyp,lower = 0,upper=1,
                    a=alpha_1+alpha_2,
                    b=1-beta_1,
                    c=alpha_1+alpha_2+beta_2,
                    y=x) 
    fx<-(1/(beta(alpha_1,beta_1)*beta(alpha_2,beta_2)))*(beta(alpha_1+alpha_2,beta_2))*(x^(alpha_1-1))*gaus$value
  }
  else{
    gaus<-integrate(dgaushyp,lower = 0,upper=1,
                    a=alpha_1+alpha_2,
                    b=1-beta_2,
                    c=alpha_1+alpha_2+beta_1,
                    y=1/x)
    fx<-(1/(beta(alpha_1,beta_1)*beta(alpha_2,beta_2)))*beta(alpha_1+alpha_2,beta_1)*(x^(-(1+alpha_2)))*gaus$value
  }
  return(fx)
}
# Bayesian model comparison theta_1=theta_2 with independent beta
# distributions 
# methods
# "modelcomparison" returns bayes factor in favor of alternative model
# two distributions assuming that the priors are equal for both variables
# "dif" returns savage-dickey approximation of the bayes factor in favor
# of the hypothesis that theta_1???theta_2 testing the difference between 
# two independent beta random variables at 0
# "ratio" returns savage-dickey approximation of the bayes factor in favor
# of the hypothesis that theta_1???theta_2 testing the ratio between 
# two independent beta random variables at 1
# default values initialize both prior distribution with a beta(1,1) and 
# the model comparison method
# THE METHODS MODELCOMPARISON AND DIF ARE EQUIVALENT
betabf <- function(x1,n1,x2,n2,a1=1,b1=1,a2=1,b2=1,method="modelcomparison"){
  if(method=="modelcomparison"){
    a1pos <- sum(x1)+a1
    b1pos <- n1-sum(x1)+b1
    a2pos <- sum(x2)+a2
    b2pos <- n2-sum(x2)+b2
    bf <- (beta(a1pos,b1pos)*beta(a2pos,b2pos))/(beta(a1pos+a2pos-a1,b1pos+b2pos-b1))
  }
  else if(method=="dif"){
    if((a1+a2)>1&(b1+b2)>1){
      a1pos <- sum(x1)+a1
      b1pos <- n1-sum(x1)+b1
      a2pos <- sum(x2)+a2
      b2pos <- n2-sum(x2)+b2
      priordx <-ddifbeta(0,a1,b1,a2,b2)
      posteriordx <- ddifbeta(0,a1pos,b1pos,a2pos,b2pos)
      bf <- priordx/posteriordx
    }
    else{
      stop("error: density not defined at 0 for a1+a2<=1 and/or b1+b2<=1")
    }
  }
  else if(method=="ratio"){
    if((a1+a2)>1&(b1+b2)>1){
      a1pos <- sum(x1)+a1
      b1pos <- n1-sum(x1)+b1
      a2pos <- sum(x2)+a2
      b2pos <- n2-sum(x2)+b2
      priordx <-dratiobeta(1,a1,b1,a2,b2)
      posteriordx <- dratiobeta(1,a1pos,b1pos,a2pos,b2pos)
      bf <- priordx/posteriordx
    }
    else{
      stop("error: density not defined at 1 for a1+a2<=1 and/or b1+b2<=1")
    }
  }
  return(bf)
}