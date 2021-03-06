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
    if((alpha_1+alpha_2)>1&(beta_1+beta_2)>1){
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
# Nonlinear sistem of equations to find the two points in 
# the domain of the distribution of the difference between two independent 
# beta random variables for which the density is equal to some constant k
int_f<-function(ci,k,a,b){
  ci<-sort(ci)
  y <- c()
  if(ci[1]>=-1&ci[2]<=1){
    l<-ddifbeta(ci[1],alpha_1 = a[1],beta_1=b[1],alpha_2 = a[2],beta_2 = b[2])
    u<-ddifbeta(ci[2],alpha_1 = a[1],beta_1=b[1],alpha_2 = a[2],beta_2 = b[2])
    y[1] <- l-k
    y[2] <- u-k
  }
  else{
    y[1]<-999
    y[2]<-999
  }
  return(y)
}
# Approximation to the Highest Posterior density of the difference between 
# two independen beta random variables
hpdddifbeta<-function(alpha,beta,prob=0.95,tolerance=0.00001,max.iter=1000,dens=0.1){
  mu <- alpha[1]/(alpha[1]+beta[1])-alpha[2]/(alpha[2]+beta[2])
  iter <- 0
  alt <- c()
  alt[1] <- dens
  bandera <- 1
  while(bandera==1){
    iter <- iter+1
    int_prueba<- nleqslv(c(mu-.2,mu+.2),int_f,a=alpha,b=beta,k=alt[iter])$x  
    integral <- integrate(ddifbeta,int_prueba[1],int_prueba[2],alpha_1=alpha[1],alpha_2=alpha[2],beta_1=beta[1],beta_2=beta[2])$value
    if((integral-prob)==0){
      bandera=2
      status <- 0
    }
    else if(abs(integral-prob)<=tolerance){
      bandera=2
      status <- 1
    }
    else if(iter==max.iter){
      bandera <- 2
      status <- 3
    }
    else{
      if((integral-prob)>0){
        if(iter==1){
          alt[iter+1] <- runif(1,alt[iter],ddifbeta(mu,alpha_1=alpha[1],alpha_2=alpha[2],beta_1=beta[1],beta_2=beta[2]))
        }
        else if(alt[iter]>alt[iter-1]){
          alt[iter+1] <- runif(1,alt[iter],ddifbeta(mu,alpha_1=alpha[1],alpha_2=alpha[2],beta_1=beta[1],beta_2=beta[2]))
        }
        else{
          alt[iter+1] <- runif(1,alt[iter],alt[iter-1])
        }
      }
      else{
        if(iter==1){
          alt[iter+1] <- runif(1,0,alt[iter])
        }
        else if(alt[iter]<alt[iter-1]){
          alt[iter+1] <- runif(1,0,alt[iter])
        }
        else{
          alt[iter+1] <- runif(1,alt[iter-1],alt[iter])
        }
      }
    }
  }
  results <- list()
  results$density <- alt[iter]
  results$abs.error <- abs(integral-prob)
  results$prob <- c(prob,integral)
  results$HPDI <- int_prueba
  results$convergence <- status
  return(results)
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
# Density of the relative ratio of two independent random beta variables 
# the sum fo alpha (as well as the sum of the betas) should be greater
# than 1
drelratiobeta<-function(x,alpha_1,beta_1,alpha_2,beta_2){
  if(x<0){
    stop("x must be greater than 0")
  }
  else if(x>1){
    stop("x must be lower than 1")
  }
  else if(x>0&x<=0.5){
    gaus<-integrate(dgaushyp,lower = 0,upper=1,
                    a=alpha_1+alpha_2,
                    b=1-beta_1,
                    c=alpha_1+alpha_2+beta_2,
                    y=(x/(1-x))) 
    fx<-(1/(beta(alpha_1,beta_1)*beta(alpha_2,beta_2)))*(1/((1-x)^2))*((x/(1-x))^(alpha_1-1))*beta(alpha_1+alpha_2,beta_2)*gaus$value
  }
  else{
    gaus<-integrate(dgaushyp,lower = 0,upper=1,
                    a=alpha_1+alpha_2,
                    b=1-beta_2,
                    c=alpha_1+alpha_2+beta_1,
                    y=((1-x)/x))
    fx<-(1/(beta(alpha_1,beta_1)*beta(alpha_2,beta_2)))*(1/(x^2))*(((1-x)/x)^(alpha_2-1))*beta(alpha_1+alpha_2,beta_1)*gaus$value
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
# "relative_ratio" returns returns savage-dickey approximation of the bayes factor in favor
# of the hypothesis that theta_1???theta_2 testing the ratio between 
# two independent beta random variables at 0.5
# default values initialize both prior distribution with a beta(1,1) and 
# the model comparison method
# THE METHODS MODELCOMPARISON AND DIF ARE EQUIVALENT
betabf <- function(x1,n1,x2,n2,a1=1,b1=1,a2=1,b2=1,method="modelcomparison"){
  results<-list()
  if(method=="modelcomparison"){
    a1pos <- sum(x1)+a1
    b1pos <- n1-sum(x1)+b1
    a2pos <- sum(x2)+a2
    b2pos <- n2-sum(x2)+b2
    bf <- (beta(a1pos,b1pos)*beta(a2pos,b2pos))/(beta(a1pos+a2pos-a1,b1pos+b2pos-b1))
    results$BF<-matrix(c(beta(a1pos+a2pos-a1,b1pos+b2pos-b1),
                         beta(a1pos,b1pos)*beta(a2pos,b2pos),
                         bf),nrow=3,ncol=1)
    colnames(results$BF)<-c('Value')
    rownames(results$BF)<-c('h0','h1','bf')
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
      results$BF<-matrix(c(priordx,
                           posteriordx,
                           bf),nrow=3,ncol=1)
      colnames(results$BF)<-c('Value')
      rownames(results$BF)<-c('h0','h1','bf')
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
      results$BF<-matrix(c(priordx,
                           posteriordx,
                           bf),nrow=3,ncol=1)
      colnames(results$BF)<-c('Value')
      rownames(results$BF)<-c('h0','h1','bf')
    }
    else{
      stop("error: density not defined at 1 for a1+a2<=1 and/or b1+b2<=1")
    }
  }
  else if(method=="relative_ratio"){
    if((a1+a2)>1&(b1+b2)>1){
      a1pos <- sum(x1)+a1
      b1pos <- n1-sum(x1)+b1
      a2pos <- sum(x2)+a2
      b2pos <- n2-sum(x2)+b2
      priordx <-drelratiobeta(0.5,a1,b1,a2,b2)
      posteriordx <- drelratiobeta(0.5,a1pos,b1pos,a2pos,b2pos)
      bf <- priordx/posteriordx
      results$BF<-matrix(c(priordx,
                           posteriordx,
                           bf),nrow=3,ncol=1)
      colnames(results$BF)<-c('Value')
      rownames(results$BF)<-c('h0','h1','bf')
    }
    else{
      stop("error: density not defined at 0.5 for a1+a2<=1 and/or b1+b2<=1")
    }
  }
  results$method<-method
  results$parameters<-matrix(c(a1,a1pos,
                               b1,b1pos,
                               a2,a2pos,
                               b2,b2pos),ncol=2,nrow=4,byrow=T)
  colnames(results$parameters)<-c('Prior','Posterior')
  rownames(results$parameters)<-c('alpha','beta','alpha','beta')
  return(results)
}
    