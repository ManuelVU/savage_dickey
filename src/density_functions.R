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