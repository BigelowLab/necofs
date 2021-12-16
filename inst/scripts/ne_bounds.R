nbe <- ncdf4::ncvar_get(X$NC, "nbe")
edge <- apply(nbe, 1, function(x) any(x ==0))

X$M$bounds <- "internal"
m <- X$M %>% dplyr::filter(edge)

xy <- locator(type = "l")
if (FALSE){
#dput(xy)
  xy <- list(x = c(367082.503574047, 607435.949321847, 2136581.14726906, 
    1767762.92879399, 1709746.57982038, 1709746.57982038, 470683.126741202, 
    362938.478647361), y = c(-607862.807917888, -914520.652492668, 
                             -35987.3680351906, 477871.7228739, 448863.548387097, 303822.675953079, 
                             -578854.633431085, -599574.758064516))
}

p <- do.call(cbind, xy)
p <- rbind(p, p[1,])
p <- sf::st_sfc(st_polygon(list(p)), crs = X$get_crs())

ix <- sf::st_contains(p, m, sparse = F)[1,]

opn <- m$elem[which(ix)]
cls <- m$elem[which(!ix)]

saveRDS(opn, "/mnt/ecocast/corecode/R/necofs/inst/extdata/elem_open.Rds")
saveRDS(cls, "/mnt/ecocast/corecode/R/necofs/inst/extdata/elem_closed.Rds")
