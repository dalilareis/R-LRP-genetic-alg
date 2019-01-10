library(GA)
library (stats)

#options(error=recover) #options(error=NULL)

#-------------Define initial variables (clients, depots, cars, capacities)---------------------
posClients <- matrix(c(33,1,0,0,2,30,31,29,31,32,30,31,25,31,29,31,32,0,1,0,0,2,30,31), 
                    nrow = 12, ncol = 2) # Positions of clients
colnames(posClients) <- c("X", "Y")
distClients <- as.matrix(dist(posClients)) # Distances between clients

depot <- c(0,0) # Position of depot
allPositions <- rbind(posClients, depot)
allDistances <- as.matrix(dist(allPositions))
distancesDepot <- as.matrix(allDistances[, "depot"]) # Distances between Depot and each Client

sizeChromossome <- 20 # 4 cars with 5 clients capacity
maxCli <- 5 # Capacity of each car
numCar <- 4 # Number of cars available
numCli <- 12 # Number of demands (clients requests)
unitCost <- 1 # Travelling cost /km

#tour <- c(sample.int(20, 5))

#----------------------Calculate distance for each individual gene (car)---------------------------
distanceTour <- function(tour, distMatrix, distDepot) {
  dist <- 0
  if(length(tour) != 1) {
    route <- embed(tour, 2)[, 2:1]
    dist <- sum(distMatrix[route])
  }
  totalDist <- dist + distDepot[tour[1]] + distDepot[tour[length(tour)]]
  return (totalDist * unitCost)
}

#-----------Given a chromossome, calculate the distance for all cars-------------------------------
totalDistance <- function(chromossome, distMatrix, distDepot){
  
  #Convert the chromossome into matrix where each column is a car:
  routes <- matrix(chromossome, nrow = maxCli, ncol = numCar)
  totalDist <- 0
  
  #For each car remove locations > number of demands and calculate distance
  for(car in 1:numCar){
    t <- routes[, car]
    remove <- numCli + 1:sizeChromossome #Remove numbers > numCli from tour
    tour <- setdiff(t, remove)
    if(length(tour) != 0) {
      totalDist <- totalDist + distanceTour(tour, distMatrix, distDepot)
    }
  }
  return (totalDist)
}

#-------------------Fitness Function (inverse total distance)-------------------------------------
fitness <- function(chromossome, distMatrix, distDepot){
  1 / totalDistance (chromossome, distMatrix, distDepot)
}

#----------------------------Run GA model-------------------------------------------------------
GA.fit <- ga(type = 'permutation', fitness = fitness, distMatrix = distClients, 
             distDepot = distancesDepot, min = 1, max = sizeChromossome, popSize = 30,
             maxiter = 1500, run = 200, pmutation = 0.2, monitor = NULL)

# Show distance travelled:
apply(GA.fit@solution, 1, totalDistance, distClients, distancesDepot)

# See solutions (routes)
summary(GA.fit)






