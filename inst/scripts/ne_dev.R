library(sf)
library(ncdf4)
library(dplyr)
COLORS <- palette.colors()
fvcom <- "/mnt/ecocast/corecode/R/fvcom"
#devtools::load_all(fvcom)
library(fvcom)
necofs <- "/mnt/ecocast/corecode/R/necofs"
#devtools::document(necofs)
devtools::load_all(necofs)

X <- NECOFSPhysics()
P0 = X$random_points()

if (FALSE){
  p <- sf::read_sf("/mnt/ecocast/coredata/necofs/tracks/particle_track_900_43200.gpkg") %>%
    suppressMessages(sf::st_set_crs(CB$M))
}

if (FALSE){
  tstep <- 15 * 60
  tmax <- 3600 * 1
  filename = file.path("/mnt/ecocast/coredata/necofs/tracks",
                       sprintf("particle_track_%i_%i.gpkg",tstep, tmax))
  
  p <- particle_track(X, 
                      tstep = tstep,
                      tmax = tmax,
                      filename = filename)
}

