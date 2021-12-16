cascophys = "/mnt/ecocast/corecode/R/cascophys"
devtools::load_all(cascophys)
library(sf)
library(ncdf4)
library(dplyr)
library(fvcom)

on.exit({
  CB <- 0
  ncdf4::nc_close(FV)
})

CB <- CascoBayPhysics() 
Mcb  <- CB$M
Ecb <- fvcom::fvcom_elems(CB$NC, what = 'xy')
Ccb <- sf::st_centroid(Mcb) |>
  st_transform(4326)
uri_base <- "http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_GOM/2019"
uri <- file.path(uri_base, "gom4_201901.nc")
FV <- nc_open(uri)
Mfv <- fvcom::get_mesh_geometry(FV, where = 'elems', what = 'xy')
Efv <- fvcom::fvcom_elems(FV, what = 'xy')
Cfv <- sf::st_centroid(Mfv) |>
  st_transform(4326)


if (FALSE){
  library(leaflet)
  leaflet(data = Cfv) |>
    addTiles() |>
    addCircles(radius = 1) |>
    addCircles(data = Ccb, radius = 1, color = 'orange', fillColor = "orange") 
}


#projstring <- paste0("+", ncdf4::ncatt_get(CB$NC, 0)[['CoordinateProjection']])
proj="+tmerc +datum=NAD83 +lon_0=-70d10 lat_0=42d50 k=.9999666666666667 x_0=900000 y_0=0"
wkt <- sf::st_crs(projstring)
NAD83 <- 4269

elemxy <- sapply(c("xc", "yc"),
            function(var) {ncdf4::ncvar_get(CB$NC, var)},
            simplify = FALSE) |>
  dplyr::as_tibble()

nodexy <- sapply(c("x", "y"),
                   function(var) {ncdf4::ncvar_get(CB$NC, var)},
                   simplify = FALSE) |>
  dplyr::as_tibble()

#NAD83 halfway rock
hwr <- sf::st_point(c(-70.036812, 43.658589)) |>
  sf::st_sfc(crs = (NAD83))

hwrxy <- hwr |>
  sf::st_transform(wkt)
