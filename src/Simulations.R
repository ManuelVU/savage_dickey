# Simulations of bayes factors over different values of probability 
# and sample sizes
#### n = 5 ####
# Equal thetas theta_1=theta_2=0.5
n_1 <-  5
n_2 <-  5
theta_1 <- 0.5
theta_2 <- 0.5
sim<-50000
t1 <- c()
t2 <- c()
bf_dif <- c()
bf_mod<-c()
bf_ratio <- c()
for(i in 1:sim){
  t1[i] <- rbinom(1,n_1,theta_1)
  t2[i] <- rbinom(1,n_2,theta_2)
  bf_mod[i] <- 1/betabf(x1 = t1[i],n1 = n_1,x2 = t2[i],n2 = n_2)$BF[3]
  bf_dif[i] <- 1/betabf(x1 = t1[i],n1 = n_1,x2 = t2[i],n2 = n_2,method="dif")$BF[3]
  bf_ratio[i] <- 1/betabf(x1 = t1[i],n1 = n_1,x2 = t2[i],n2 = n_2,method="relative_ratio")$BF[3]
}

plot(t1/5-t2/5,bf_dif)
