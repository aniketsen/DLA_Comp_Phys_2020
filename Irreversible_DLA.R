
set.seed(100)                                                #set seed for random number generation

N <<- 10000                                                  #set number of particles in cluster

n <<- 400                                                    #set length and width of lattice

g <<- logical(n*n)                                           #create grid for positioning particles, False for lattice points with no particles

xd <- c()
yd <- c()                                                    #these will store x and y positions of the particles

counter <- 0

while(counter < N) {                                         #randomly position N particles in the lattice
  x0 <- sample(c(0:(n-1)),1)
  y0 <- sample(c(0:(n-1)),1)                                 #sample 2 points between 0 and n-1 for x and y
  if (!g[x0*n+y0+1]){                                        #if the position doesn't  contain a particle already, place a particle
    g[x0*n+y0+1] <- TRUE
    xd <- c(xd,x0)
    yd <- c(yd,y0)
    counter <- counter + 1
  }
}
clusternum <- function(m) {                                  #Function for finding the different clusters and returning the number of clusters
  counter1 <- 0
  while(length(which(m))>1) {
    i <- sample(which(m),1) - 1                              #randomly select a point in the lattice that has a particle
    xpos <- i%/%n                                            #these will store the x and y points of this cluster
    ypos <- i%%n
    m[i+1] <- FALSE                                          #since particle counted, set it to false
    counter1 <- counter1+1
    k <- 1
    l <- 1
    while(k<=l) {                                            #check the neighbouring points. if a particle exists, add it to x and y positions of the
      if(m[((xpos[k]+1)%%n)*n+ypos[k]+1]) {                     #cluster. set the grid location to false. Repeat the process again for these newly added
        m[((xpos[k]+1)%%n)*n+ypos[k]+1] <- FALSE                   #particles, untill there are no more neighbouring points.
        xpos <- c(xpos,(xpos[k]+1)%%n)
        ypos <- c(ypos,ypos[k])
        l <- l+1
      }
      if(m[((xpos[k]-1)%%n)*n+ypos[k]+1]) {
        m[((xpos[k]-1)%%n)*n+ypos[k]+1] <- FALSE
        xpos <- c(xpos,(xpos[k]-1)%%n)
        ypos <- c(ypos,ypos[k])
        l <- l+1
      }
      if(m[xpos[k]*n+((ypos[k]+1)%%n)+1]) {
        m[xpos[k]*n+((ypos[k]+1)%%n)+1] <- FALSE
        xpos <- c(xpos,xpos[k])
        ypos <- c(ypos,(ypos[k]+1)%%n)
      }
      if(m[xpos[k]*n+((ypos[k]-1)%%n)+1]) {
        m[xpos[k]*n+((ypos[k]-1)%%n)+1] <- FALSE
        xpos <- c(xpos,xpos[k])
        ypos <- c(ypos,(ypos[k]-1)%%n)
      }
      k <- k+1
    }
    assign(paste("clusx_",counter1,sep=""),xpos,envir=.GlobalEnv)    #store the x and y positions of the cluster in "clusx_(cluster number)" and
    assign(paste("clusy_",counter1,sep=""),ypos,envir=.GlobalEnv)         #"clusy_(cluster number)" variables
      
  }
  i <- which(m)
  xpos <- i%/%n
  ypos <- i%%n
  counter1 <- counter1+1
  assign(paste("clusx_",counter1,sep=""),xpos,envir=.GlobalEnv)
  assign(paste("clusy_",counter1,sep=""),ypos,envir=.GlobalEnv)
  return(counter1)
}

num <- clusternum(g)                                     #get number of clusters
clus <- seq(1,num)                                       #set of cluster numbers
st <- 1
while(num>1) {                                          #loop untill only 1 cluster remains
  sel <- sample(clus,1)                                 #select 1 cluster at random from the existing clusters
  xsel <- get(paste("clusx_",sel,sep=""))               #get the x and y positions of the cluster selected
  ysel <- get(paste("clusy_",sel,sep=""))
  u <- runif(1)                                         #generate a random number and move the cluster by 1 step in a random direction using that
  if(u <= 0.25) xsel <- (xsel+st)%%n
  else if(u<=0.5) xsel <- (xsel-st)%%n
  else if(u<=0.75) ysel <- (ysel+st)%%n
  else ysel <- (ysel-st)%%n
  assign(paste("clusx_",sel,sep=""),xsel)               #store the new position of the cluster
  assign(paste("clusy_",sel,sep=""),ysel)
  for(i in clus[clus != sel]) {                         #loop through all the remaining clusters except the selected one
    xpos <- get(paste("clusx_",i,sep=""))               #get x and y positions of the loop cluster
    ypos <- get(paste("clusy_",i,sep=""))
    flag <- 0
    for(k in 1:length(xpos)) {                          #check if any point on the loop cluster is neighbouring the selected cluster after the random walk
      ii <- which(xsel==xpos[k]+1)
      if(length(ii)!=0) {
        for(l in ii) {
          if (ysel[l]==ypos[k]) {
            flag <- 1
            break
          }
        }
      }
      ij <- which(xsel==xpos[k]-1)
      if(length(ij)!=0) {
        for(l in ij) {
          if (ysel[l]==ypos[k]) {
            flag <- 1
            break
          }
        }
      }
      ik <- which(ysel==ypos[k]+1)
      if(length(ik)!=0) {
        for(l in ik) {
          if (xsel[l]==xpos[k]) {
            flag <- 1
            break
          }
        }
      }
      il <- which(ysel==ypos[k]-1)
      if(length(il)!=0) {
        for(l in il) {
          if (xsel[l]==xpos[k]) {
            flag <- 1
            break
          }
        }
      }
      if(flag == 1) break
    }
    if(flag==1) {           #if any cluster is found at the neighbouring points
      xn <- c(xsel,xpos)
      yn <- c(ysel,ypos)
      assign(paste("clusx_",sel,sep=""),xn)         #x and y positions of that cluster is added to the selected cluster
      assign(paste("clusy_",sel,sep=""),yn)
      rm(list=paste("clusx_",i,sep=""))             #x and y positions of the 2nd cluster is deleted
      rm(list=paste("clusy_",i,sep=""))
      clus <- clus[clus != i]                       #cluster number of the one that was merged with the selected cluster is removed from the set of cluster numbers
      num <- num-1                                  #total number of clusters decreases by 1
      break
    }
      
  }
}
xf <- get(paste("clusx_",clus[1],sep=""))                        #get x and y positions of the final cluster remaining
yf <- get(paste("clusy_",clus[1],sep=""))

plot(xf,yf,pch=16,cex=0.5,xlab="x",ylab="y",xlim=c(0,n),ylim=c(0,n),xaxs="i",yaxs="i")       #plot the generated cluster

Data1 <- data.frame(x=xf,y=yf)
csvFileName1 <- paste("Data_",N,"_",n,".csv",sep="")
write.csv(Data1,file=csvFileName1)                               #x and y positions are stored in a csv file

