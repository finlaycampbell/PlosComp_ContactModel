##===== Load functions and libraries =====##
source("internals.R")
load.libs()

##===== Use the data from the manuscript ======##
load("../data/store.RData")

##===== Or run the simulations yourself =====##
## This will take a long time on a single desktop computer;
## reduce the number of runs if necessary

## Specify the number of non-infectious contacts per person (psi)
psi <- seq(0, 10, 2)

## Specify the contact reporting coverage (epsilon)
eps <- seq(0, 1, 0.2)

## Create a grid of values for eps and psi
param <- expand.grid(psi = psi, eps = eps)

## Specify the number of runs per grid point
param$runs <- 100

## Run the analysis for every grid point and save raw results
for(i in 1:nrow(param)) {
  r <- run.analysis(param)
  save(r, file = paste0("../output/res_", i, ".RData"))
}

## Extract summary results from the raw results
store <- create.store("../data/")


##===== Create figures =====##
## These will be saved to the /figs folder - see R/internals for code
create.figs(store, dir = "../figs/")
