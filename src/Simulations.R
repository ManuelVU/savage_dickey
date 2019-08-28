# Simulations of bayes factors over different values of probability 
# and sample sizes
#### Equal thetas at 0.5 ####
library(colorspace)
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
    bf_ratio2[n,i] <- 1/betabf(x1 = t1[n,i],n1 = n_1[n],x2 = t2[n,i],n2 = n_2[n],method="ratio")$BF[3]
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
pdf("results/correct_bf.pdf")
plot(0,0,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(1,9))
axis(1)
axis(2)
for(n in 1:length(n_1)){
  lines(k, correct_dif[n,],col=samplesize_color[n],lty=2,lwd=2)  
  lines(k, correct_ratio[n,],col=samplesize_color[n],lwd=2)
}
dev.off()
#### Equal thetas at .1, .3, .7, .9 n=50 ####
library(colorspace)
n_1 <- 50
n_2 <- 50
samplesize_color <- qualitative_hcl(n=4,'Dark3')
theta_1 <- c(.1, .3, .7, .9)
theta_2 <- c(.1, .3, .7, .9)
sim<-10000
t1 <-  matrix(NA, nrow=length(theta_1),ncol=sim)
t2 <-  matrix(NA, nrow=length(theta_1),ncol=sim)
bf_dif <- matrix(NA, nrow=length(theta_1),ncol=sim)
bf_ratio <- matrix(NA, nrow=length(theta_1),ncol=sim)
for(n in 1:length(theta_1)){
  for(i in 1:sim){
    t1[n,i] <- rbinom(1,n_1,theta_1[n])
    t2[n,i] <- rbinom(1,n_2,theta_2[n])
    bf_dif[n,i] <- 1/betabf(x1 = t1[n,i],n1 = n_1,x2 = t2[n,i],n2 = n_2,method="dif")$BF[3]
    bf_ratio[n,i] <- 1/betabf(x1 = t1[n,i],n1 = n_1,x2 = t2[n,i],n2 = n_2,method="relative_ratio")$BF[3]
  }
}
k <- seq(1,15,0.1)
correct_dif<-matrix(NA,nrow=length(theta_1),ncol=length(k))
correct_ratio<-matrix(NA,nrow=length(theta_1),ncol=length(k))
for(n in 1:length(theta_1)){
  for(i in 1:length(k)){
    correct_dif[n,i] <- sum(bf_dif[n,]>k[i])/sim
    correct_ratio[n,i] <- sum(bf_ratio[n,]>k[i])/sim
  }
}
# biproduct of the result that bf_ratio depends on the sum of 
# the number of successes in sample
plot(0,0,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(1,max(k)))
axis(1)
axis(2)
legend('topright',legend=theta_1,lty=1,col=samplesize_color)
for(n in 1:length(theta_1)){
  lines(k, correct_dif[n,],col=samplesize_color[n],lty=2,lwd=2)  
  lines(k, correct_ratio[n,],col=samplesize_color[n],lwd=2)
}

#### Equal thetas for ratio and dif ####
library(colorspace)
n_1 <- 50
n_2 <- 50
theta_1 <- seq(.05,.95,.05)
theta_2 <- seq(.05,.95,.05)
samplesize_color <- sequential_hcl(n=length(theta_1),"Mint")
sim<-10000
t1 <-  matrix(NA, nrow=length(theta_1),ncol=sim)
t2 <-  matrix(NA, nrow=length(theta_1),ncol=sim)
bf_dif <- matrix(NA, nrow=length(theta_1),ncol=sim)
bf_ratio <- matrix(NA, nrow=length(theta_1),ncol=sim)
for(n in 1:length(theta_1)){
  for(i in 1:sim){
    t1[n,i] <- rbinom(1,n_1,theta_1[n])
    t2[n,i] <- rbinom(1,n_2,theta_2[n])
    bf_dif[n,i] <- 1/betabf(x1 = t1[n,i],n1 = n_1,x2 = t2[n,i],n2 = n_2,method="dif")$BF[3]
    bf_ratio[n,i] <- 1/betabf(x1 = t1[n,i],n1 = n_1,x2 = t2[n,i],n2 = n_2,method="relative_ratio")$BF[3]
  }
}
k <- seq(1,15,0.1)
correct_dif<-matrix(NA,nrow=length(theta_1),ncol=length(k))
correct_ratio<-matrix(NA,nrow=length(theta_1),ncol=length(k))
for(n in 1:length(theta_1)){
  for(i in 1:length(k)){
    correct_dif[n,i] <- sum(bf_dif[n,]>k[i])/sim
    correct_ratio[n,i] <- sum(bf_ratio[n,]>k[i])/sim
  }
}

plot(0,0,axes=F, ann=F,type="n",ylim=c(0,1),xlim=c(1,max(k)))
axis(1)
axis(2, las=2)
box()
for(n in 1:length(theta_1)){
  lines(k,correct_ratio[n,],col=samplesize_color[n],lwd=1.5)
}

plot(0,0,axes=F, ann=F,type="n",ylim=c(0,1),xlim=c(1,max(k)))
axis(1)
axis(2, las=2)
box()
for(n in 1:length(theta_1)){
  lines(k,correct_dif[n,],col=samplesize_color[n],lwd=1.5,lty=2)
}

#### Unequal thetas delta=0.15 theta=0.5 ####
library(colorspace)
n_1 <- c(5,10,20,40,80,140,200)
n_2 <- c(5,10,20,40,80,140,200)
samplesize_color <- qualitative_hcl(n=7,'Dark3')
sim<-10000
theta_1 <- 0.5
theta_2 <- sort(rep(c(-1,1),sim/2))*.15+theta_1
t1 <-  matrix(NA, nrow=length(n_1),ncol=sim)
t2 <-  matrix(NA, nrow=length(n_1),ncol=sim)
bf_dif <- matrix(NA, nrow=length(n_1),ncol=sim)
bf_ratio <- matrix(NA, nrow=length(n_1),ncol=sim)
for(n in 1:length(n_1)){
  for(i in 1:sim){
    t1[n,i] <- rbinom(1,n_1[n],theta_1)
    t2[n,i] <- rbinom(1,n_2[n],theta_2[i])
    bf_dif[n,i] <- betabf(x1 = t1[n,i],n1 = n_1[n],x2 = t2[n,i],n2 = n_2[n],method="dif")$BF[3]
    bf_ratio[n,i] <- betabf(x1 = t1[n,i],n1 = n_1[n],x2 = t2[n,i],n2 = n_2[n],method="relative_ratio")$BF[3]
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

k <- seq(1,9,0.1)
correct_dif_n<-matrix(NA,nrow=length(n_1),ncol=length(k))
correct_ratio_n<-matrix(NA,nrow=length(n_1),ncol=length(k))
for(n in 1:length(n_1)){
  for(i in 1:length(k)){
    correct_dif_n[n,i] <- sum(bf_dif[n,1:(sim/2)]>k[i])/(sim*0.5)
    correct_ratio_n[n,i] <- sum(bf_ratio[n,1:(sim/2)]>k[i])/(sim*0.5)
  }
}
plot(0,0,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(1,9))
axis(1)
axis(2)
for(n in 1:length(n_1)){
  lines(k, correct_dif_n[n,],col=samplesize_color[n],lty=2,lwd=2)  
  lines(k, correct_ratio_n[n,],col=samplesize_color[n],lwd=2)
}


k <- seq(1,9,0.1)
correct_dif_p<-matrix(NA,nrow=length(n_1),ncol=length(k))
correct_ratio_p<-matrix(NA,nrow=length(n_1),ncol=length(k))
for(n in 1:length(n_1)){
  for(i in 1:length(k)){
    correct_dif_p[n,i] <- sum(bf_dif[n,((sim/2)+1):sim]>k[i])/(sim*0.5)
    correct_ratio_p[n,i] <- sum(bf_ratio[n,((sim/2)+1):sim]>k[i])/(sim*0.5)
  }
}
plot(0,0,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(1,9))
axis(1)
axis(2)
for(n in 1:length(n_1)){
  lines(k, correct_dif_p[n,],col=samplesize_color[n],lty=2,lwd=2)  
  lines(k, correct_ratio_p[n,],col=samplesize_color[n],lwd=2)
}

#### Unequal thetas delta=0.15 theta=0.2,0.4,0.6,0.8 n=50####
library(colorspace)
n_1 <- 50
n_2 <- 50
samplesize_color <- qualitative_hcl(n=4,'Dark3')
sim<-10000
theta_1 <- c(.2, .4, .6, .8)
theta_2 <- matrix(NA,nrow=length(theta_1),ncol=sim)
theta_2[1,] <- sort(rep(c(-1,1),sim/2))*.15+theta_1[1]
theta_2[2,] <- sort(rep(c(-1,1),sim/2))*.15+theta_1[2]
theta_2[3,] <- sort(rep(c(-1,1),sim/2))*.15+theta_1[3]
theta_2[4,] <- sort(rep(c(-1,1),sim/2))*.15+theta_1[4]
t1 <-  matrix(NA, nrow=length(theta_1),ncol=sim)
t2 <-  matrix(NA, nrow=length(theta_1),ncol=sim)
bf_dif <- matrix(NA, nrow=length(theta_1),ncol=sim)
bf_ratio <- matrix(NA, nrow=length(theta_1),ncol=sim)
for(n in 1:length(theta_1)){
  for(i in 1:sim){
    t1[n,i] <- rbinom(1,n_1,theta_1[n])
    t2[n,i] <- rbinom(1,n_2,theta_2[n,i])
    bf_dif[n,i] <- betabf(x1 = t1[n,i],n1 = n_1,x2 = t2[n,i],n2 = n_2,method="dif")$BF[3]
    bf_ratio[n,i] <- betabf(x1 = t1[n,i],n1 = n_1,x2 = t2[n,i],n2 = n_2,method="relative_ratio")$BF[3]
  }
}
k <- seq(1,15,0.1)
correct_dif<-matrix(NA,nrow=length(theta_1),ncol=length(k))
correct_ratio<-matrix(NA,nrow=length(theta_1),ncol=length(k))
for(n in 1:length(theta_1)){
  for(i in 1:length(k)){
    correct_dif[n,i] <- sum(bf_dif[n,]>k[i])/sim
    correct_ratio[n,i] <- sum(bf_ratio[n,]>k[i])/sim
  }
}
# biproduct of the result that bf_ratio depends on the sum of 
# the number of successes in sample
plot(0,0,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(1,max(k)))
axis(1)
axis(2)
legend('topright',legend=theta_1,lty=1,col=samplesize_color)
for(n in 1:length(theta_1)){
  lines(k, correct_dif[n,],col=samplesize_color[n],lty=2,lwd=2)  
  lines(k, correct_ratio[n,],col=samplesize_color[n],lwd=2)
}

k <- seq(1,9,0.1)
correct_dif_n<-matrix(NA,nrow=length(theta_1),ncol=length(k))
correct_ratio_n<-matrix(NA,nrow=length(theta_1),ncol=length(k))
for(n in 1:length(theta_1)){
  for(i in 1:length(k)){
    correct_dif_n[n,i] <- sum(bf_dif[n,1:(sim/2)]>k[i])/(sim*0.5)
    correct_ratio_n[n,i] <- sum(bf_ratio[n,1:(sim/2)]>k[i])/(sim*0.5)
  }
}
plot(0,0,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(1,9))
axis(1)
axis(2)
for(n in 1:length(theta_1)){
  lines(k, correct_dif_n[n,],col=samplesize_color[n],lty=2,lwd=2)  
  lines(k, correct_ratio_n[n,],col=samplesize_color[n],lwd=2)
}


k <- seq(1,9,0.1)
correct_dif_p<-matrix(NA,nrow=length(theta_1),ncol=length(k))
correct_ratio_p<-matrix(NA,nrow=length(theta_1),ncol=length(k))
for(n in 1:length(theta_1)){
  for(i in 1:length(k)){
    correct_dif_p[n,i] <- sum(bf_dif[n,((sim/2)+1):sim]>k[i])/(sim*0.5)
    correct_ratio_p[n,i] <- sum(bf_ratio[n,((sim/2)+1):sim]>k[i])/(sim*0.5)
  }
}
plot(0,0,type='n',axes=F,ann=F,ylim=c(0,1),xlim=c(1,9))
axis(1)
axis(2)
for(n in 1:length(theta_1)){
  lines(k, correct_dif_p[n,],col=samplesize_color[n],lty=2,lwd=2)  
  lines(k, correct_ratio_p[n,],col=samplesize_color[n],lwd=2)
}

