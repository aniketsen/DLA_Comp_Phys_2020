
set.seed(100)    #set seed for random number generation

t <- 10000      #set number of particles in cluster

ros <- ceiling(sqrt(100*t))    #set length and width of lattice
m4 <- logical(ros*ros)         #create grid for positioning particles, False for lattice points with no particles
ce <- ceiling(ros/2)           #set the origin (0,0) at the center
m4[ce*ros+ce] <- TRUE          #position seed particle at origin
xpos <- 0                      #these will store x and y positions of the particles, first one is set at (0,0)
ypos <- 0
RandWalk <- function(x,y){     #Function for random walk
  u <- runif(1)
  if (u<=0.25) x <- x+1
  else if (u<=0.5) y <- y+1
  else if (u<=0.75) x <- x-1
  else y <- y-1
  return(c(x,y))
}
nbrchk <- function(x,y) {      #Function for checking neighbouring lattice sites. Will  return 1 if at least 1 particle found, else 0.
  if (m4[(x+ce)*ros+y+ce]) return(0)
  else if(m4[(x+ce+1)*ros+y+ce] || m4[(x+ce-1)*ros+y+ce] || m4[(x+ce)*ros+y+ce+1] || m4[(x+ce)*ros+y+ce-1]) return(1)
  else return(0)
}
counter <-  1                              #counts number of particles added
while(counter<t) {                         #loop to repeat algorith untill total number of particles are placed
  rm <- max(sqrt(xpos^2+ypos^2))           #find the maximum radius of existing cluster
  th <- runif(1)*2*pi                      #generate random angle to generate a new particle at rmax*cos(theta) and rmax*sin(theta)
  xn <- ceiling((rm)*cos(th)) 
  yn <- ceiling((rm)*sin(th))
  npos <- c(xn,yn)
  flag <- 0
  while(flag!=1) {                         #do random walk untill the particle is added to cluster or deleted
    npos <- RandWalk(npos[1],npos[2])      #do one step of random walk
    dis <- sqrt(sum(npos^2))               #calculate distance of the new position from center
    if (dis<rm+2) {
      if (nbrchk(npos[1],npos[2])==1) {    #if particle is situated at less than rmax+2, check neighbouring sites
        xpos <- c(xpos,npos[1])
        ypos <- c(ypos,npos[2])
        m4[(npos[1]+ce)*ros+npos[2]+ce] <- TRUE       #if there is a particle at neighbouring site add the x and y positions and set the grid position to TRUE
        counter <- counter+1
        flag <- 1
      }
    }
   
    if (rm != 0 && dis>2*rm) break             #if particle goes to further than 2*rmax, it is discarded
  }
}


plot(xpos,ypos,pch = 16,cex = 0.75)           #plot the generated cluster

Data1 <- data.frame(x=xpos,y=ypos)
csvFileName1 <- paste("Data_",t,".csv",sep="")
write.csv(Data1,file=csvFileName1)                         #the x and y positions can be saved in a csv file
csvFileName3 <- paste("m4_",t,".rds",sep="")
saveRDS(m4, file=csvFileName3)                             #the lattice grid can be stored in a rds file

