library(here)
run_dir <- here("mods","m1/")
nh <- "..//..//src//nh" 
arg <- " -steepness 0.7"
setwd(run_dir)
getwd()
run <- paste0(nh, arg)
system(run)
run
nh
