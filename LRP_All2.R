library(rstudioapi)
library(GA)
library (stats)
library(plyr)
library(xlsx)

set.seed(123)
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))

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
  
  if (folder != "Instances/Tuzun") {
    carCapacity <- as.numeric(read.table(file, header=F,nrows = 1, dec = ".",
                                       blank.lines.skip = TRUE, skip = 3 + numDep + 1 + numCli + 1))
  }
  
  if (folder == "Instances/Tuzun") {
    carCapacity <- as.numeric(read.table(file, header=F,nrows = 1, dec = ".",
                                       blank.lines.skip = TRUE, skip = 3 + numDep + numCli + 1))
  }
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
  distClients <- as.matrix(dist(posClients))
  return (distClients)
}

getDistDepots <- function(file) {
  lines <- readLines(file, n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  
  posDepot <- as.matrix(read.table(file, header = F, nrows = numDep, skip = 2, blank.lines.skip = TRUE))
  row.names(posDepot) <- c(paste0("Depot_", 1:nrow(posDepot)))
  
  posClients <- as.matrix(read.table(file, header = F, nrows = numCli, blank.lines.skip = TRUE, 
                                     skip = 3 + numDep))
  
  allPositions <- rbind(posClients, posDepot)
  allDistances <- as.matrix(dist(allPositions))
  distDepots <- as.matrix(allDistances[1:nrow(posClients), (nrow(posClients)+1):nrow(allDistances)]) 
  return (distDepots)
}

getDepotCap <- function(file) {
  folder <- dirname(file)
  lines <- readLines(file, n=2)
  numCli <- as.numeric(lines[1])
  numDep <- as.numeric(lines[2])
  
  depotCapacity <- as.matrix(read.table(file, header = F, nrows = numDep, blank.lines.skip = TRUE,
                                          skip = 3 + numDep + numCli + 3))
  return (depotCapacity)
}

getDemands <- function(file) {
  folder <- dirname(file)
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

#Divide chromossome in genes
sliceChrom <- function(input, by){ 
  starts <- seq(1,length(input),by)
  tt <- lapply(starts, function(y) input[y:(y+(by-1))])
  llply(tt, function(x) x[!is.na(x)])
}

# ---------------------Aux Functions to calculate Fitness for GA model------------------------------

#Calculate distance for each individual gene (routes for each car associated to a depot)
distanceCost <- function(depotIndex, tour, distClients, distDepot) {
  dist <- 0
  if(length(tour) != 1) {
    route <- embed(tour, 2)[, 2:1]
    dist <- sum(distClients[route])
  }
  totalDist <- dist + distDepot[tour[1], depotIndex] + distDepot[tour[length(tour)], depotIndex]
  return (totalDist * unitCost)
}

#Given a chromossome, calculate the distance+cost for all routes
totalCost <- function(chromossome, distClients, distDepot, costDepot, demands, depotCapacity){
  
  travelCost <- 0
  operationsCost <- 0
  penalty <- 0
  
  #---------------------------------------DEPOT----------------------------------------------
  
  #convert the chromossome into matrix where each column is a depot
  depots <- matrix(chromossome, nrow = (avgClients * maxCarDepot), ncol = numDep)
  
  #For each depot remove locations > number of Clients and add Depot Cost
  for(depot in seq_along(depots[1, ])){
    
    selectDepot <- depots[, depot]
    removeCli <- numCli + 1:sizeChromossome 
    bigTour <- setdiff(selectDepot, removeCli)
    
    if(length(bigTour) != 0) {  #Real tours only (with possible clients)
      
      #Restriction 1: Clients Demands cannot exceed Depot capacity (considering that max capacity 
      #was selected when there were different depot capacities)
      depotCap <- as.numeric(depotCapacity[depot, ])
      orders <- demands[c(bigTour), ]
      totalOrder <- sum(orders)
      if (totalOrder > depotCap) {
        penalty <- penalty + (10*costDepot[depot, ])
      }
      
      #Add Depot Cost (only when it exists tour)
      operationsCost <- operationsCost + costDepot[depot, ]
    }
  }
  
  #------------------------------------------CAR------------------------------------------------
  
  #Convert the chromossome into matrix where each column is a car
  routes <- matrix(chromossome, nrow = avgClients, ncol = maxRotas)
  
  #For each car remove numbers > Client numbers (get tours for available clients)
  for(car in seq_along(routes[1, ])){
    selectCar <- routes[, car]
    remove <- numCli + 1:sizeChromossome 
    tour <- setdiff(selectCar, remove)
    
    if(length(tour) != 0) {
      
      #Calculate Depot to which car belongs (to calculate distance and penalty)
      depotIndex <- ceiling(car / maxCarDepot)

      #Restriction 2: Clients Demands cannot exceed Car capacity
      orders <- demands[c(tour), ]
      totalOrder <- sum(orders)
      if (totalOrder > carCapacity) {
        penalty <- penalty + (10*costDepot[depotIndex, ])
      }
      
      #Increment cost (add tours)
      operationsCost <- operationsCost + penalty + costCar #Add cost of 1 car to each tour
      travelCost <- travelCost + distanceCost(depotIndex, tour, distClients, distDepot)
    }
  }
  #--------------------------------------------TOTAL--------------------------------------------
  totalCost <- operationsCost + travelCost
  return (totalCost)
}

#Fitness Function (minimize costs)
fitness <- function(chromossome, distClients, distDepot, costDepot, demands, depotCapacity){
  1 / totalCost (chromossome, distClients, distDepot, costDepot, demands, depotCapacity)
}

#----------------------------Run GA model & Initializations------------------------------------------

#Load files 
filesBar <- sapply("Barreto", getFiles)
filesProd <- sapply("Prodhon", getFiles)
allFiles <- rbind(filesBar, filesProd, sapply("Tuzun", getFiles))
 
allFiles <- sapply("Tuzun", getFiles)

#Read all numeric values from files to data frame
listNumerics <- as.data.frame(t(as.data.frame(lapply(allFiles, sapply, getNumerics))))

#Set progress bar for file loop 
pbar <- create_progress_bar('text')
files <- nrow(allFiles)
pbar$init(files)

#Set gas cost and initialize data frame to store final results
unitCost <- 1
results <- data.frame()

#Start loop to run model for each instance file 
for (i in 1:files ) {

  #------------------Read data from files and assign to variables-----------------------------
  numCli <- as.numeric(listNumerics[i, 1])
  numDep <- as.numeric(listNumerics[i, 2])
  carCapacity <- as.numeric(listNumerics[i, 3])
  costCar <- as.numeric(listNumerics[i, 4])
  distClients <- getDistClients(allFiles[i, ])
  distDepot <- getDistDepots(allFiles[i, ])
  depotCapacity <- getDepotCap(allFiles[i, ])
  demands <- getDemands(allFiles[i, ])
  costDepot <- getCostDepot(allFiles[i, ])

  #------------------Calculate variables to be used in chromossome structure-----------------------
  
  #Extract max Depot Capacity (for instances with variable capacity)
  maxDepCapacity <- as.numeric(depotCapacity[which.max(depotCapacity), ])
  
  #Max number of cars that a Depot can supply 
  maxCarDepot <- floor(maxDepCapacity / carCapacity)
  
  #Max number of cars (routes) that can exist in the system (will be the gene in chromossome)
  maxRotas <- numDep * maxCarDepot
  
  #Get total Demands and determine min cars (routes) necessary to satisfy them
  totalDemands <- as.numeric(colSums(demands))
  minRotas <- floor(totalDemands / carCapacity)

  #Max number of clients a Car can supply (based on average demand of each client)
  maxClientsCar <- ceiling(carCapacity / (totalDemands / numCli))
  
  #Min number of clients served by each car(route) if max number of routes are used
  minClientMaxRota <- ceiling(maxClientsCar / (maxRotas / minRotas))
  
  #Average Client number that a car can supply
  avgClients <- ceiling((maxClientsCar + minClientMaxRota) / 2)
  
  #Limit size of chromossome (if each depot has a lot of cars (Tuzun: 1 depot can supply all clients))
  if (maxCarDepot >= 20) {
    maxCarDepot <- floor(maxCarDepot/2)
    maxRotas <- numDep * maxCarDepot
  }
  #Check if it is still too big
  if (maxCarDepot >= 20) {
    maxCarDepot <- floor(maxCarDepot/3)
    maxRotas <- numDep * maxCarDepot
  }
  #Third check for cases with a lot of depots
  if (numDep > 10 && maxCarDepot > 10) {
    maxCarDepot <- floor(maxCarDepot/3)
    maxRotas <- numDep * maxCarDepot
  }
  #If sizeChromossome is still too big, then it is not the maxCarDepot contribution, but the clients
  # Usually, Tuzun instances still have too large chromossomes
  if (maxRotas * avgClients > 600) { #Clients are never more than 200
    avgClients <- floor(avgClients/3)
  }
  
  #Chromossome with X genes (X = Depot + Car associated) and each gene has Y alelles 
  #(Y = average number of Clients for each route, according to Demands, Capacity & Clients)
  sizeChromossome <- maxRotas * avgClients
  
  #PopSize too large with big chromossomes does not improve cost and has more CPU usage
  if (sizeChromossome >= 200 || numCli >= 200) {
    population <- 30
  }
  if(sizeChromossome < 200) {
    population <- 2 * numCli
  }
  
  #-------------------------------Run GA Model and get execution time--------------------------------
  startTime <- proc.time()
  
  GA.fit <- ga(type = 'permutation', fitness = fitness, distClients = distClients, distDepot = distDepot, 
               costDepot = costDepot, demands = demands, depotCapacity = depotCapacity, min = 1, 
               max = sizeChromossome, popSize = population, maxiter = 1500, run = 200, pmutation = 0.2, monitor = NULL)
  
  endTime <- proc.time() - startTime

  #Lowest totalCost found
  solution <- as.matrix(apply(GA.fit@solution, 1, totalCost, distClients, distDepot, costDepot, 
                              demands, depotCapacity))[1, ]
  iterations <- GA.fit@iter
  
  chromossome <- GA.fit@solution[1,]
  
  #Get Number of used Depots from solution found (combine every Z genes, Z = max car for Depot)
  carDepot <- sliceChrom(chromossome, (maxCarDepot * avgClients))
  usedDep <- 0
  for (j in 1:length(carDepot)){
    if (any(carDepot[[j]] <= numCli)){
      usedDep <- usedDep + 1
    }
  }
  
  #Get Number of used Cars from solution found (slice chromossome to get lenght=clients)
  car <- sliceChrom(chromossome, avgClients)
  usedCars <- 0
  for (k in 1:length(car)){
    if (any(car[[k]] <= numCli)){
      usedCars <- usedCars + 1
    }
  }
  
  #Save results to Data Frame
  results <- rbind(results, data.frame(numCli, numDep, usedDep, usedCars, sizeChromossome, population,
                                       iterations, endTime["elapsed"], solution, 
                                       row.names = allFiles[i, ]))
  pbar$step()
  
}
#------------------------------------Save Results--------------------------------------------------

colnames(results) <- c("Number Clients", "Available Depots", "Used Depots", "Used Cars", 
                       "Chromossome Size", "Population Size", "Iterations", "Execution Time(s)", 
                       "Lowest Cost found")

write.xlsx(results, file = "Results.xlsx", sheetName = "Tuzun", col.names = T, 
           row.names = T, append = T)

