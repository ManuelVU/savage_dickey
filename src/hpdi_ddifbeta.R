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
#hpdddifbeta<-function(alpha,beta,prob=0.95,tolerance=0.00001,max.iter=1000,dens=0.1){
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
          alt[iter+1] <- alt[iter]+(1/2)*(ddifbeta(mu,alpha_1=alpha[1],alpha_2=alpha[2],beta_1=beta[1],beta_2=beta[2])-alt[iter])
        }
        else if(alt[iter]>alt[iter-1]){
          alt[iter+1] <- alt[iter]+(1/2)*(ddifbeta(mu,alpha_1=alpha[1],alpha_2=alpha[2],beta_1=beta[1],beta_2=beta[2])-alt[iter])
        }
        else{
          alt[iter+1] <- alt[iter]-(1/2)*abs(alt[iter]-alt[iter-1])
        }
      }
      else{
        if(iter==1){
          alt[iter+1] <- alt[iter]-(1/2)*alt[iter]
        }
        else if(alt[iter]<alt[iter-1]){
          alt[iter+1] <- alt[iter]-(1/2)*alt[iter]
        }
        else{
          alt[iter+1] <- alt[iter]+(1/2)*abs(alt[iter]-alt[iter-1])
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

hpdddifbeta(c(41,41),c(31,31),max.iter = 3000)

