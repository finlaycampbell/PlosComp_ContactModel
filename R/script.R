##===== Load functions and libraries =====##
library(here)
source(here("R/internals.R"))

## Load libraries (warning - this will install missing packages)
load.libs()

##===== Use the data from the manuscript ======##
load(here("analysis/sim_store.rds"))
load(here("analysis/sars_store.rds"))
load(here("analysis/eg_store.rds"))
dat <- get.sars.dat()

##===== Or run the analysis yourself =====##
## This will take a long time on a single desktop computer;
## reduce the number of runs if necessary

## 1) Simulate and reconstruct outbreaks
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
  save(r, file = here(paste0("output/res_", i, ".RData")))
}

## Extract summary results from the raw results
store <- create.store(here("output/"))


## 2) Analyse the SARS outbreak data
## Specify which outbreaker settings to use
sett <- c('noctd', 'nolambda', 'fix_4', 'fix_5', 'lambda_10')

## Run the analysis
sars_store <- sapply(sett, run.sars, dat, simplify = FALSE, USE.NAMES = TRUE)


##===== Create figures =====##
## These will be saved to the /figs folder - see R/internals for code
## The visNetwork figures will open in the browser
create.figs(sim_store, sars_store, eg_store, dat, dir = here("figs/"))
