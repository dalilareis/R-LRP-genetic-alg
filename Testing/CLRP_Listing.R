library(rstudioapi)
library(GA)
library (stats)
library(plyr)
library(xlsx)

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

#Divide chromossome in genes
sliceChrom <- function(input, by){ 
  starts <- seq(1,length(input),by)
  tt <- lapply(starts, function(y) input[y:(y+(by-1))])
  llply(tt, function(x) x[!is.na(x)])
}

# ---------------------Aux Functions to calculate Fitness for GA model------------------------------

#Calculate distance for each individual gene (routes for each car associated to a depot)
distanceCost <- function(depot, tour, distClients, distDepot) {
  dist <- 0
  if(length(tour) != 1) {
    route <- embed(tour, 2)[, 2:1]
    dist <- sum(distClients[route])
  }
  totalDist <- dist + distDepot[tour[1], depot] + distDepot[tour[length(tour)], depot]
  return (totalDist * unitCost)
}

#Given a chromossome, calculate the distance+cost for all routes
totalCost <- function(chromossome, distClients, distDepot, costDepot, demands, depotCapacity){
  
  #Convert the chromossome into matrix where each column is a car:
  routes <- matrix(chromossome, nrow = avgClients, ncol = (numDep*maxCarDepot))
  travelCost <- 0
  operationsCost <- 0
  penalty <- 0
  listDepots <- c()
  
  #For each gene (car) remove numbers > Client numbers (get tours for available clients)
  for(car in seq_along(routes[1, ])){
    selectCar <- routes[, car]
    remove <- numCli + 1:sizeChromossome 
    tour <- setdiff(selectCar, remove)
    
    if(length(tour) != 0) { #Real tours only (with possible clients)
      
      #Check which depot the car belongs to: if already listed (accounted), do not add Depot Cost
      depot <- ceiling(car / maxCarDepot)
      listDepots <- unique(append(listDepots, depot))
      if(!depot %in% listDepots){
        operationsCost <- operationsCost + costDepot[depot, ]
      }
      
      #Restriction 1: Clients Demands cannot exceed Car capacity
      orders <- demands[c(tour), ]
      totalOrder <- sum(orders)
      if (totalOrder > carCapacity) {
        penalty <- penalty + (100*costDepot[depot, ])
      }
      
      #Restriction 2: Clients Demands cannot exceed Depot capacity
      depotCap <- depotCapacity[depot, ]
      if (totalOrder > depotCap) {
        penalty <- penalty + (100*costDepot[depot, ])
      }
      
      #Increment cost (add tours)
      operationsCost <- operationsCost + penalty + costCar #Add cost of 1 car to each tour (gene=car)
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
#Load files 
filesBar <- sapply("Barreto", getFiles)
allFiles <- rbind(filesBar, sapply("Prodhon", getFiles))

#Read all numeric values from files to data frame
listNumerics <- t(as.data.frame(lapply(allFiles, sapply, getNumerics)))

#Set progress bar for file loop 
pbar <- create_progress_bar('text')
files <- nrow(allFiles)
pbar$init(files)

#Set gas cost and initialize data frame to store final results
unitCost <- 1
results <- data.frame()

#Start loop to run model for each instance file 
for (i in 1:files ) {
  
  #Read data from files and assign to variables
  numCli <- as.numeric(listNumerics[i, 1])
  numDep <- as.numeric(listNumerics[i, 2])
  carCapacity <- as.numeric(listNumerics[i, 3])
  costCar <- as.numeric(listNumerics[i, 4])
  distClients <- getDistClients(allFiles[i, ])
  distDepot <- getDistDepots(allFiles[i, ])
  depotCapacity <- getDepotCap(allFiles[i, ])
  demands <- getDemands(allFiles[i, ])
  costDepot <- getCostDepot(allFiles[i, ])
  
  #Calculate variables to be used in chromossome structure
  maxCarDepot <- floor(depotCapacity[1] / carCapacity)
  totalDemands <- as.numeric(colSums(demands))
  minRotas <- ceiling(totalDemands / carCapacity)
  maxRotas <- numDep * maxCarDepot
  maxClientsCar <- ceiling(carCapacity / (totalDemands / numCli)) #Trans
  minClientMaxRota <- ceiling(maxClientsCar / (maxRotas / minRotas))
  avgClients <- ceiling((maxClientsCar + minClientMaxRota) / 2)
  
  #Chromossome with X genes (X = Depot + Car associated) and each gene has Y alelles 
  #(Y = average number of Clients needed for each route, according to Demands, Capacity & Clients)
  sizeChromossome <- maxRotas * avgClients
  population <- 30
  
#Run GA Model and get execution time for running model
#--------------------------------------------------------------------------------------------------
  startTime <- proc.time()

  GA.fit <- ga(type = 'permutation', fitness = fitness, distClients = distClients, distDepot = distDepot, 
               costDepot = costDepot, demands = demands, min = 1, 
               max = sizeChromossome, popSize = population, maxiter = 1500, run = 200, pmutation = 0.5, monitor = NULL)
  
  endTime <- proc.time() - startTime
#---------------------------------------------------------------------------------------------------
  
  #Lowest totalCost found
  solution <- as.matrix(apply(GA.fit@solution, 1, totalCost, distClients, distDepot, costDepot, demands))[1, ]
  iterations <- GA.fit@iter
  
  #Get Number of used Depots from solution found (combine every X genes, X = max car for Depot)
  chromossome <- GA.fit@solution[1,]
  carDepot <- sliceChrom(chromossome, (maxCarDepot * avgClients))
  usedDep <- 0
  for (j in 1:length(carDepot)){
    if (any(carDepot[[j]] <= numCli)){
      usedDep <- usedDep + 1
    }
  }
  
  #Save results to Data Frame
  results <- rbind(results, data.frame(numCli, numDep, usedDep, sizeChromossome, population, iterations, 
                                       endTime["elapsed"], solution, 
                                       row.names = allFiles[i, ]))
  pbar$step()
  
}

colnames(results) <- c("Number Clients", "Available Depots", "Used Depots", "Chromossome Size", 
                       "Population Size", "Iterations", "Execution Time (s)", "Lowest Cost found")

write.xlsx(results, file = "Results.xlsx", sheetName = "Folha8", col.names = T, 
           row.names = T, append = T)
