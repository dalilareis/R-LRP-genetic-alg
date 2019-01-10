library(rstudioapi)
library(GA)
library (stats)
library(plyr)

current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))
set.seed(123)

# ---------------Aux Functions (read data from files and organize it accordingly)--------------------

getFiles <- function(folder) {
  path <- paste ("Instances", folder, sep = "/")
  files <- sort(list.files(path, pattern="dat$", full.names = T), decreasing = T)
  return (files)
}

#Single values: number CLients, number Depots, car Capacity, car cost
getNumerics <- function(file) { 
  folder <- dirname(file)
  lines <- readLines(file, n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  carCapacity <- as.numeric(read.table(file, header=F,nrows = 1, dec = ".",
                                       blank.lines.skip = TRUE, skip = 3 + numDep + 1 + numCli + 1))
  costCar <- 0
  if (folder != "Instances/Barreto") {
    costCar <- as.numeric(read.table(file, header = F, nrows = 1, 
                                     blank.lines.skip = TRUE, skip = 2*3 + 3*numDep + 2*numCli + 3*1))
  }
  
  listing <- data.frame(numCli, numDep, carCapacity, costCar)
  return(listing)
}

getDistClients <- function(file) {
  lines <- readLines(file, n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  posClients <- as.matrix(read.table(file, header = F, nrows = numCli, blank.lines.skip = TRUE, 
                                     skip = 3 + numDep))
  if (ncol(posClients) > 2) {
    posClients <- posClients [, 1:2]
  }
  distClients <- as.matrix(dist(posClients))
  return (distClients)
}

getDistDepots <- function(file) {
  lines <- readLines(file, n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  posDepot <- as.matrix(read.table(file, header = F, nrows = numDep, 
                                   skip = 2, blank.lines.skip = TRUE))
  row.names(posDepot) <- c(paste0("Depot_", 1:nrow(posDepot)))
  posClients <- as.matrix(read.table(file, header = F, nrows = numCli, 
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

getDepotCap <- function(file) {
  lines <- readLines(file, n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  depotCapacity <- as.matrix(read.table(file, header = F, nrows = numDep, blank.lines.skip = TRUE,
                                          skip = 3 + numDep + numCli + 3))
  row.names(depotCapacity) <- c(paste0("Depot_", 1:nrow(depotCapacity)))
  return (depotCapacity)
}

getDemands <- function(file) {
  lines <- readLines(file, n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  demands <- as.matrix(read.table(file, header = F, nrows = numCli, blank.lines.skip = TRUE,
                                  skip = 3 + numDep + numCli + 3 + numDep + 1))
  return (demands)
}

getCostDepot <- function(file) {
  lines <- readLines(file, n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  costDepot <- as.matrix(read.table(file, header = F, nrows = numDep, blank.lines.skip = TRUE,
                                    skip = 2*3 + 2*numDep + 2*numCli + 2*1))
  row.names(costDepot) <- c(paste0("Depot_", 1:nrow(costDepot)))
  return (costDepot)
}

sliceChrom <- function(input, by){ 
  starts <- seq(1,length(input),by)
  tt <- lapply(starts, function(y) input[y:(y+(by-1))])
  llply(tt, function(x) x[!is.na(x)])
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
  routes <- matrix(chromossome, nrow = avgClients, ncol = (numDep*maxCarDepot))
  travelCost <- 0
  operationsCost <- 0
  penalty <- 0
  listDepots <- c()
  
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

startTime <- proc.time()

instanceFolder <- c("Barreto", "Prodhon")
allFiles <- lapply(instanceFolder, getFiles)
listNumerics <- as.data.frame(t(as.data.frame(lapply(allFiles, sapply, getNumerics))))
unitCost <- 1

listDistClients <- lapply(allFiles, lapply, getDistClients)
listDistDepots <- lapply(allFiles, lapply, getDistDepots)
listDepotCapacity <- lapply(allFiles, lapply, getDepotCap)
listDemands <- lapply(allFiles, lapply, getDemands)
listCostDepots <- lapply(allFiles, lapply, getCostDepot)
header <- c("Number Clients", "Available Depots", "Used Depots", "Chromossome Size", "Iterations",
            "Execution Time", "Lowest Cost found")
results <- data.frame()


#Start loop to run model for each instance file 
pbar <- create_progress_bar('text')
files <- nrow(listNumerics)
pbar$init(files)

for (i in 1:files ) {
  numCli <- as.numeric(listNumerics[i, 1])
  numDep <- as.numeric(listNumerics[i, 2])
  carCapacity <- as.numeric(listNumerics[i, 3])
  costCar <- as.numeric(listNumerics[i, 4])
  
  if (i <= length(allFiles[[1]])){
    distClients <- listDistClients[[1]][[i]]
    distDepot <- listDistDepots[[1]][[i]]
    depotCapacity <- listDepotCapacity[[1]][[i]]
    demands <- listDemands[[1]][[i]]
    costDepot <- listCostDepots[[1]][[i]]
  }
  if (i > length(allFiles[[1]])){
    distClients <- listDistClients[[2]][[i]]
    distDepot <- listDistDepots[[2]][[i]]
    depotCapacity <- listDepotCapacity[[2]][[i]]
    demands <- listDemands[[2]][[i]]
    costDepot <- listCostDepots[[2]][[i]]
  }
  maxCarDepot <- floor(depotCapacity[1] / carCapacity)
  totalDemands <- as.numeric(colSums(demands))
  minRotas <- ceiling(totalDemands / carCapacity)
  maxRotas <- numDep * maxCarDepot
  maxClientsCar <- ceiling(carCapacity / (totalDemands / numCli)) #Trans
  minClientMaxRota <- ceiling(maxClientsCar / (maxRotas / minRotas))
  avgClients <- ceiling((maxClientsCar + minClientMaxRota) / 2)
  
  #Chromossome with X genes (X = Depot + Car associated) and each gene has Y alelles 
  #(Y = average number of Clients for each car, according to total Demands)
  sizeChromossome <- maxRotas * avgClients
  
  #popSize = 1.5*numCli

  
  GA.fit <- ga(type = 'permutation', fitness = fitness, distClients = distClients, distDepot = distDepot, 
               costDepot = costDepot, demands = demands, min = 1, max = sizeChromossome, popSize = (2*numCli),
               maxiter = 1500, run = 200, pmutation = 0.2, monitor = NULL)
  
  pbar$step()
  
  #Lowest totalCost found
  solution <- as.matrix(apply(GA.fit@solution, 1, totalCost, distClients, distDepot, costDepot, demands))[1, ]
  iterations <- GA.fit@iter
  
  #Get Number of used Depots from solution found
  chromossome <- GA.fit@solution[1,]
  carDepot <- sliceChrom(chromossome, avgClients)
  usedDep <- 0
  length(carDepot)
  for (j in 1:length(carDepot)){
    if (any(carDepot[[j]] <= numCli)){
      usedDep <- usedDep + 1
    }
  }
  endTime <- proc.time() - startTime
  
  #Save results in Data Frame
  results <- rbind(results, data.frame(numCli, numDep, usedDep, sizeChromossome, iterations, 
                                       endTime["elapsed"], solution, row.names = rownames(listNumerics[1, ])))
  
}
colnames(results) <- header
#--------------------------------Save Results to Excel File-----------------------------------------



