library(rstudioapi)
library(GA)
library (stats)

current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))
set.seed(123)
#browser() usar para debugging

# ---------------Aux Functions (read data from files and organize it accordingly)--------------------

getFiles <- function(folder) {
  path <- paste ("Instances", folder, sep = "/")
  files <- sort(list.files(path, pattern="dat$"), decreasing = T)
  return (files)
}

#Single values: number CLients, number Depots, car Capacity, car cost
getNumerics <- function(file, folder) { 
  path <- paste ("Instances", folder, sep = "/")
  lines <- readLines(paste(path, file, sep="/"), n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  carCapacity <- as.numeric(read.table(paste(path, file, sep="/"), header=F,nrows = 1,
                                       blank.lines.skip = TRUE, skip = 3 + numDep + 1 + numCli + 1))
  costCar <- 0
  if (folder != "Barreto") {
    costCar <- as.numeric(read.table(paste(path, file, sep="/"), header = F, nrows = 1, 
                                    blank.lines.skip = TRUE, skip = 2*3 + 3*numDep + 2*numCli + 3*1))
  }
  
  listing <- data.frame(numCli, numDep, carCapacity, costCar)
  return(listing)
}

getDistClients <- function(file, folder) {
  path <- paste ("Instances", folder, sep = "/")
  lines <- readLines(paste(path, file, sep="/"), n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  posClients <- as.matrix(read.table(paste(path, file, sep="/"), header = F, nrows = numCli, 
                                     blank.lines.skip = TRUE, skip = 3 + numDep))
  if (ncol(posClients) > 2) {
    posClients <- posClients [, 1:2]
  }
  distClients <- as.matrix(dist(posClients))
  return (distClients)
}

getDistDepots <- function(file, folder) {
  path <- paste ("Instances", folder, sep = "/")
  lines <- readLines(paste(path, file, sep="/"), n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  posDepot <- as.matrix(read.table(paste(path, file, sep="/"), header = F, nrows = numDep, 
                                   skip = 2, blank.lines.skip = TRUE))
  row.names(posDepot) <- c(paste0("Depot_", 1:nrow(posDepot)))
  posClients <- as.matrix(read.table(paste(path, file, sep="/"), header = F, nrows = numCli, 
                                     blank.lines.skip = TRUE, skip = 3 + numDep))
  if (ncol(posClients) > 2 || ncol(posDepot) > 2) {
    posClients <- posClients [, 1:2]
    posDepot <- posDepot [, 1:2]
  }
  allPositions <- rbind(posClients, posDepot)
  allDistances <- as.matrix(dist(allPositions))
  distDepots <- as.matrix(allDistances[1:nrow(posClients), (nrow(posClients)+1):nrow(allDistances)]) 
  return (distDepots)
}

getDepotCap <- function(file, folder) {
  path <- paste ("Instances", folder, sep = "/")
  lines <- readLines(paste(path, file, sep="/"), n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  if (folder != "Tuzun") {
    depotCapacity <- as.matrix(read.table(paste(path, file, sep="/"), header = F, nrows = numDep, 
                                          blank.lines.skip = TRUE, skip = 3 + numDep + numCli + 3))
    row.names(depotCapacity) <- c(paste0("Depot_", 1:nrow(depotCapacity)))
  }
  if (folder == "Tuzun") {
    demands <- as.matrix(read.table(paste(path, file, sep="/"), header = F, nrows = numDep, 
                                    blank.lines.skip = TRUE, skip = 3 + numDep + numCli + 3 + numDep + 1))
    totalOrders <- colSums(demands)
    depotCapacity <- matrix(totalOrders, nrow = numDep, ncol = 1)
    row.names(depotCapacity) <- c(paste0("Depot_", 1:nrow(depotCapacity)))
  }
  return (depotCapacity)
}

getDemands <- function(file, folder) {
  path <- paste ("Instances", folder, sep = "/")
  lines <- readLines(paste(path, file, sep="/"), n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  demands <- as.matrix(read.table(paste(path, file, sep="/"), header = F, nrows = numCli, 
                                        blank.lines.skip = TRUE, skip = 3 + numDep + numCli + 3 + numDep + 1))
  return (demands)
}

getCostDepot <- function(file, folder) {
  path <- paste ("Instances", folder, sep = "/")
  lines <- readLines(paste(path, file, sep="/"), n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  costDepot <- as.matrix(read.table(paste(path, file, sep="/"), header = F, nrows = numDep, 
                                    blank.lines.skip = TRUE, skip = 2*3 + 2*numDep + 2*numCli + 2*1))
  row.names(costDepot) <- c(paste0("Depot_", 1:nrow(costDepot)))
  return (costDepot)
}

# ---------------------Aux Functions to calculate Fitness for GA model------------------------------

#Calculate distance for each individual gene (routes for each depot)
distanceCost <- function(depot, tour, distClients, distDepot) {
  dist <- 0
  if(length(tour) != 1) {
    route <- embed(tour, 2)[, 2:1]
    dist <- sum(distClients[route])
  }
  totalDist <- dist + distDepot[tour[1], depot] + distDepot[tour[length(tour)], depot]
  return (totalDist * unitCost)
}

#Given a chromossome, calculate the distance+cost for all depots TODO: add costCar
totalCost <- function(chromossome, distClients, distDepot, costDepot, demands){
  
  #Convert the chromossome into matrix where each column is a depot:
  routes <- matrix(chromossome, nrow = avgClients, ncol = maxRotas)
  travelCost <- 0
  operationsCost <- 0
  penalty <- 0
  listDepots <- c()
  
  #init <- seq(from = 1, to = numDep*maxCarDepot, by= maxCarDepot )
  #end <- seq(from = 2, to = numDep*maxCarDepot, by= maxCarDepot )
  #depots <- rbind(routes[, init], routes[, end]) #combine routes for each depot
  
  #For each depot remove locations > number of Clients and add Depot Cost and calculate distance
    # for(depot in seq_along(depots[1, ])){
    # selectDepot <- depots[, depot]
    # remove <- numCli + 1:sizeChromossome #Remove numbers > numCli from tour
    # tour <- setdiff(selectDepot, remove)
    # if(length(tour) != 0) {
    #   operationsCost <- operationsCost + costDepot[depot, ]
    # } 
  for(car in seq_along(routes[1, ])){
    selectCar <- routes[, car]
    remove <- numCli + 1:sizeChromossome #Remove numbers > numCli from tour
    tour <- setdiff(selectCar, remove)
    
    if(length(tour) != 0) {
      depot <- ceiling(car / maxCarDepot)
      listDepots <- unique(append(listDepots, depot))
      if(!depot %in% listDepots){
        operationsCost <- operationsCost + costDepot[depot, ]
      }
      orders <- demands[c(tour), ]
      totalOrder <- sum(orders)
      if (totalOrder > carCapacity) {
        penalty <- penalty + (10*costDepot[depot, ])
      }
      operationsCost <- operationsCost + penalty + costCar
      travelCost <- travelCost + distanceCost(depot, tour, distClients, distDepot)
    }
  }
  totalCost <- operationsCost + travelCost
  return (totalCost)
}

#Fitness Function (minimize costs)
fitness <- function(chromossome, distClients, distDepot, costDepot, demands){
  1 / totalCost (chromossome, distClients, distDepot, costDepot, demands)
}

#----------------------------Run GA model & Initializations------------------------------------------


#Chromossome with X genes (X = number Depots) and each gene has Y alelles 
#(Y = max Cars per Depot * min routes necessary to satisfy total Demands of Clients)
# 1 Car equivalent to 1 route

#chromossome <- c(sample.int(sizeChromossome, (minRotas*maxCarDepot)*numDep))

#TODO: List everything, do loop for all files (add computing times to data frame)-------

unitCost <- 1
filesBar <- getFiles("Barreto")
listNumerics <- sapply(filesBar, getNumerics, "Barreto")
numCli <- as.numeric(listNumerics[1, ])
numDep <- as.numeric(listNumerics[2, ])
carCapacity <- as.numeric(listNumerics[3, ])
costCar <- as.numeric(listNumerics[4, ])
listDistClients <- sapply(filesBar, getDistClients, "Barreto")
listDistDepots <- sapply(filesBar, getDistDepots, "Barreto")
listDepotCapacity <- sapply(filesBar, getDepotCap, "Barreto")
listDemands <- sapply(filesBar, getDemands, "Barreto")
listCostDepots <- sapply(filesBar, getCostDepot, "Barreto")

#-----------------------------TESTING---------------------------------------------------------
unitCost <- 1
filesBar <- getFiles("Barreto")
distClients <- getDistClients(filesBar[4], "Barreto")
distDepot <- getDistDepots(filesBar[4], "Barreto")
costDepot <- getCostDepot(filesBar[4], "Barreto")
depotCapacity <- getDepotCap(filesBar[4], "Barreto")
listaNumerics <- getNumerics(filesBar[4], "Barreto")
numCli <- as.numeric(listaNumerics[, 1])
numDep <- as.numeric(listaNumerics[, 2])
carCapacity <- as.numeric(listaNumerics[, 3])
costCar <- as.numeric(listaNumerics[, 4])
demands <- getDemands(filesBar[4], "Barreto")
encTotal <- as.numeric(colSums(demands))
maxCarDepot <- floor(depotCapacity[1] / carCapacity)
minRotas <- ceiling(encTotal / carCapacity)
maxRotas <- numDep * maxCarDepot
maxClients <- ceiling(carCapacity / (encTotal / numCli))
minClientMaxRota <- ceiling(maxClients / (maxRotas / minRotas))
avgClients <- ceiling((maxClients + minClientMaxRota) / 2)
sizeChromossome <- maxRotas * avgClients

#GA Model
GA.fit <- ga(type = 'permutation', fitness = fitness, distClients = distClients, distDepot = distDepot, 
             costDepot = costDepot, demands = demands, min = 1, max = sizeChromossome, maxiter = 1500, 
             popSize = sizeChromossome, run = 200, pmutation = 0.2, monitor = NULL)
            

# Show best totalCost found
as.matrix(apply(GA.fit@solution, 1, totalCost, distClients, distDepot, costDepot, demands))[1,]


# See model details
summary(GA.fit)

# See solutions (routes)
GA.fit@solution[1,]

filesTuzun <- getFiles("Tuzun")
filesTuzun[1]
posCli <- getDistClients(filesTuzun[1])
fit <- kmeans(posCli, 20)




