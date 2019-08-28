load("results/karina_2.Rdata")



negative <- sum(base_completah$biggerchosen[which(base_completah$pairP%%2==0)])
n_neg <- length(base_completah$biggerchosen[which(base_completah$pairP%%2==0)])
positive <- sum(base_completah$biggerchosen[which(base_completah$pairP%%2!=0)])
n_pos <- length(base_completah$biggerchosen[which(base_completah$pairP%%2!=0)])

k_bf <- betabf(positive,n_pos,negative,n_neg)

x <- seq(0,1,0.01)
dx <- c()
for(i in 1:length(x)){
  dx[i] <- dbeta(x[i],k_bf$parameters[1,2],k_bf$parameters[2,2])
}
curve(dbeta(x,357,95))
curve(dbeta(x,92,360),col='blue',add=T)

ind_bf<-c()
for(i in 1:10){
  negative <- sum(base_completah$biggerchosen[which(base_completah$pairP%%2==0&base_completah$participant==i)])
  n_neg <- length(base_completah$biggerchosen[which(base_completah$pairP%%2==0&base_completah$participant==i)])
  positive <- sum(base_completah$biggerchosen[which(base_completah$pairP%%2!=0&base_completah$participant==i)])
  n_pos <- length(base_completah$biggerchosen[which(base_completah$pairP%%2!=0&base_completah$participant==i)])
  ind_bf[i] <- betabf(positive,n_pos,negative,n_neg)$BF[3]
}

#Oyentes
negative <- sum(base_completao$biggerchosen[which(base_completao$pairP%%2==0)])
n_neg <- length(base_completao$biggerchosen[which(base_completao$pairP%%2==0)])
positive <- sum(base_completao$biggerchosen[which(base_completao$pairP%%2!=0)])
n_pos <- length(base_completao$biggerchosen[which(base_completao$pairP%%2!=0)])

k_bf <- betabf(positive,n_pos,negative,n_neg)

curve(dbeta(x,k_bf$parameters[1,2],k_bf$parameters[2,2]))
curve(dbeta(x,k_bf$parameters[3,2],k_bf$parameters[4,2]),col='blue',add=T)

#Oyentes Individual
ind_bf<-c()
for(i in 1:10){
  negative <- sum(base_completao$biggerchosen[which(base_completao$pairP%%2==0&base_completao$participant==i)])
  n_neg <- length(base_completao$biggerchosen[which(base_completao$pairP%%2==0&base_completao$participant==i)])
  positive <- sum(base_completao$biggerchosen[which(base_completao$pairP%%2!=0&base_completao$participant==i)])
  n_pos <- length(base_completao$biggerchosen[which(base_completao$pairP%%2!=0&base_completao$participant==i)])
  print(c(negative, positive))
  ind_bf[i] <- betabf(positive,n_pos,negative,n_neg)$BF[3]
}
