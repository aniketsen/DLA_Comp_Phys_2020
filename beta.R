beta <- array(data=0, dim=c(5,10,4))       #Empty array to store calculated beta values
for (trial in 1:10) {                      #Since we had 10 clusters for each size
csvfilename <- paste("Data_50000-",trial,".csv",sep="")
dat <- read.csv(file=csvfilename, header=TRUE, sep=",", colClasses=c("NULL", NA, NA))
N0 <- 50000
x0 <- dat[1:N0,1]
y0 <- dat[1:N0,2]
Nn <- 1
for(N in c(10000,15000,20000,30000,50000)) {  #For different size of the cluster
x <- x0[1:N]
y <- y0[1:N]
pn <- 1
for (p in c(0.5,0.75,0.9,0.95)) {             #For 50%, 75%, 90% and 95% last added points
m <- N*p
m0 <- N - m                                   #Gives the starting point for the corresponding percentage
Rg <- c()
for (i in seq(m0,N,100)){                     #From starting point to ending point for every 100 points added
  R <- sum(sqrt(x[1:i]^2+y[1:i]^2))/i         #Calculate Rg for cluster of i points
  Rg <- c(Rg,R)
}
n <- log(seq(m0,N,100))
lrg <- log(Rg)
fit <- lm(lrg ~ n)                            #Fit a curve between log(Rg) and log(N) to get the slope which is beta
beta[Nn,trial,pn] <- coefficients(fit)[2]     #Add the beta value for the given cluster size, given trial and given percentage
pn <- pn+1
}
Nn <- Nn+1
}
}
saveRDS(beta, file="beta_values.rds")         #Store all calculated beta values for further calculations
#require(ggplot2)
#ggdat <- data.frame(first=x, second=y)
#print(qplot(x, y, data=ggdat))
#plot(x,y,pch = 16,cex = 0.75)