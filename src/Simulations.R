# Simulations of bayes factors over different values of probability 
# and sample sizes
#### n = 5,10,20,40,80 ####
# Equal thetas theta_1=theta_2=0.5
n_1 <- c(5,10,20,40,80,140,200)
n_2 <- c(5,10,20,40,80,140,200)
samplesize_color <- qualitative_hcl(n=7,'Dark3')
theta_1 <- 0.5
theta_2 <- 0.5
sim<-10000
t1 <-  matrix(NA, nrow=length(n_1),ncol=sim)
t2 <-  matrix(NA, nrow=length(n_1),ncol=sim)
bf_dif <- matrix(NA, nrow=length(n_1),ncol=sim)
bf_ratio <- matrix(NA, nrow=length(n_1),ncol=sim)
for(n in 1:length(n_1)){
  for(i in 1:sim){
    t1[n,i] <- rbinom(1,n_1[n],theta_1)
    t2[n,i] <- rbinom(1,n_2[n],theta_2)
    bf_dif[n,i] <- 1/betabf(x1 = t1[n,i],n1 = n_1[n],x2 = t2[n,i],n2 = n_2[n],method="dif")$BF[3]
    bf_ratio[n,i] <- 1/betabf(x1 = t1[n,i],n1 = n_1[n],x2 = t2[n,i],n2 = n_2[n],method="relative_ratio")$BF[3]
  }
}
k <- seq(1,9,0.1)
correct_dif<-matrix(NA,nrow=length(n_1),ncol=length(k))
correct_ratio<-matrix(NA,nrow=length(n_1),ncol=length(k))
for(n in 1:length(n_1)){
  for(i in 1:length(k)){
    correct_dif[n,i] <- sum(bf_dif[n,]>k[i])/sim
    correct_ratio[n,i] <- sum(bf_ratio[n,]>k[i])/sim
  }
}
plot(0,0,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(1,9))
axis(1)
axis(2)
for(n in 1:length(n_1)){
  lines(k, correct_dif[n,],col=samplesize_color[n],lty=2,lwd=2)  
  lines(k, correct_ratio[n,],col=samplesize_color[n],lwd=2)
}

plot(abs(t1[7,]/n_1[7]-t2[7,]/n_2[7]),bf_dif[7,])
abline(h=1)
plot(abs(0.5-t1[7,]/(t1[7,]+t2[7,])),bf_ratio[7,])
abline(h=1)
test<-sequential_hcl(length(sort(unique(t1[7,]+t2[7,]))),'Dark Mint')
points(rep(0,length(test)),bf_ratio[7,which()],
       col=test,pch=16)
plot(t1[7,which(t1[7,]==t2[7,])]+t2[7,which(t1[7,]==t2[7,])],
     bf_ratio[7,which(t1[7,]==t2[7,])],col=test)
