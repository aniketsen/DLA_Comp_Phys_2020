

set.seed(100)               #set seed for random number generation

t <- 10000                  #set number of particles in cluster

ros <- ceiling(4*sqrt(t))         #set length and width of lattice
m4 <- logical(ros*ros*ros)        #create grid for positioning particles, False for lattice points with no particles
ce <- ceiling(ros/2)              #set the origin (0,0,0) at the center
m4[(ce*ros+ce)*ros+ce] <- TRUE    #position seed particle at origin
xpos <- 0
ypos <- 0                         #these will store x, y and z positions of the particles, first one is set at (0,0,0)
zpos <- 0
RandWalk <- function(x,y,z){      #Function for random walk
  u <- runif(1)
  if (u<=0.167) x <- x+1
  else if (u<=0.333) y <- y+1
  else if (u<=0.5) z <- z+1
  else if (u<=0.667) x <- x-1
  else if (u<=0.833) y <- y-1
  else z <- z-1
  return(c(x,y,z))
}
nbrchk <- function(x,y,z) {       #Function for checking neighbouring lattice sites. Will  return 1 if at least 1 particle found, else 0.
  if (m4[((x+ce)*ros+y+ce)*ros+z+ce]) return(0)
  else if(m4[((x+ce+1)*ros+y+ce)*ros+z+ce] || m4[((x+ce-1)*ros+y+ce)*ros+z+ce] || m4[((x+ce)*ros+y+ce+1)*ros+z+ce] || m4[((x+ce)*ros+y+ce-1)*ros+z+ce] || m4[((x+ce)*ros+y+ce)*ros+z+ce+1] || m4[((x+ce)*ros+y+ce)*ros+z+ce-1]) return(1)
  else return(0)
}
counter <-  1                     #counts number of particles added
while(counter<t) {                                        #loop to repeat algorith untill total number of particles are placed
  rm <- max(sqrt(xpos^2+ypos^2+zpos^2))                   #find the maximum radius of existing cluster
  nn <- rnorm(3)
  if (nn[1]==0 && nn[2]==0 && nn[3]==0) flag <- 1
  else {
  nn <- nn/sqrt(sum(nn^2))
  npos <- ceiling(nn*rm)                                  #randomly position the new particle on the sphere with radius rmax
  flag <- 0 }
  while(flag!=1) {
    npos <- RandWalk(npos[1],npos[2],npos[3])             #random walk of 1 step
    dis <- sqrt(sum(npos^2))                              #calculate distance of the new position from center
    if (dis<rm+2) {
      if (nbrchk(npos[1],npos[2],npos[3])==1) {
        xpos <- c(xpos,npos[1])
        ypos <- c(ypos,npos[2])
        zpos <- c(zpos,npos[3])
        m4[((npos[1]+ce)*ros+npos[2]+ce)*ros+npos[3]+ce] <- TRUE       #if there is a particle at neighbouring site add the x, y and z positions and set the grid position to TRUE
        counter <- counter+1
        flag <- 1
      }
    }
    
    if (rm != 0 && dis>2*rm) break                              #if particle goes to further than 2*rmax, it is discarded
  }
}

library(car)
scatter3d(x=xpos,y=ypos,z=zpos, point.col = "red", surface=FALSE)     #plot the generated cluster

Data1 <- data.frame(x=xpos,y=ypos,z=zpos)
csvFileName1 <- paste("Data_3D_",t,".csv",sep="")
write.csv(Data1,file=csvFileName1)                                    #the x, y and z positions can be saved in a csv file
csvFileName3 <- paste("m4_3D_",t,".rds",sep="")
saveRDS(m4, file=csvFileName3)                                        #the lattice grid can be stored in a rds file

