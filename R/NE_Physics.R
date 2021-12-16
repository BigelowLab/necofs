#' A class for navigating  FVCOM datasets
#' 
#' @description R6 base class
#' @export
NE_Physics <- R6::R6Class("NE_Physics",
                          
  public = list(
    #' @field filename character, path or URL
    filename = NULL,
    #' @field NC ncdf4 object
    NC = NULL,
    #' @field M sf table of mesh
    M = NULL,
    #' @field t0 POSIXct timestamp identifying the start time
    t0 = NULL,
    #' @field verbose logical for helpful messaging
    verbose = NULL,
    
    #' @param filename character either a filename or OpenDAP URL
    #' @param origin POSIXct timestamp indicating the start of the experiment
    #' @param crs 
    #' @param verbose logical for helpful messaging
    initialize = function(filename,
                          origin = as.POSIXct("1858-11-17 00:00:00", tz = 'UTC'),
                          verbose = FALSE){
      self$filename <- filename[1]
      self$verbose <- verbose[1]
      self$t0 <- origin[1]
      private$message("opening NCDF4")
      self$NC <- try(ncdf4::nc_open(self$filename))
      if (inherits(self$NC, "try-error")) stop("unable to open NCDF4 resource:", filename[1])
      private$message("retrieving mesh")
      self$M <- fvcom::get_mesh_geometry(self$NC, where = 'elems', what = 'xy',
                                         crs = self$get_crs(form = 'wkt'))
      ok <- self$append_bounds()
    },
    
    #' @description Pretty print
    print = function(){
      cat("NE_Physics class:", self$filename, "\n")
      cat("  CRS        :", self$get_crs(), "\n")
      cat("  t0         :", self$t0, "\n")
      cat("  n timesteps:", self$NC$dim$time$len, "\n")
      cat("  n elements :", self$NC$dim$nele$len, "\n")
      cat("  n nodes    :", self$NC$dim$node$len, "\n")
      cat("  verbose    :", self$verbose, "\n")
      
    }, 
    
    
    #' @description Retrieve the CRS
    #' @param form character, one of 'proj', 'wkt'
    #' @return CRS in the specified form
    get_crs = function(form = c("proj", "wkt")[1]){
      a <- ncdf4::ncatt_get(self$NC, 0)
      crs <- paste0("+", a$CoordinateProjection)
      if (tolower(form[1]) == 'wkt'){
        crs <- sf::st_crs(crs)
      }
      crs
    },
    
    #' @description retrieve the time relative to some epoch/origin
    #' @param origin POSIXct the orogin of relative timestamps
    #' @param ... other arguments for \code{\link[fvcom]{fvcom_time}}
    get_time = function(origin = self$t0, ...){
      fvcom::fvcom_time(self$NC, origin = origin[1], ...)
    },
    
    #' @description Retriev one or more random points
    #' @param n numeric, the number of points to generate
    #' @param select_time character, one of 'first' (default), 'last' or 'random'
    #' @return sf object of class POINT
    random_points = function(n = 1, select_time = c("first", "last", "random")[1]){
      random_point(self, n = n, select_time = select_time)
    },
    
    #' @description find proxy for characteristic size based upon element size
    #'
    #' NOTE if we update GEOS we could get size of inscribed circle
    #'
    #' @param x element mesh as sf
    #' @param what character
    #' \itemize{
    #' \item min or square - side length of sqaure with equivalent area
    #' \item max - side of square with double the area (or 1.4*sqrt(area))
    #' \item circle - radius of circle with equivalent area
    #'  }
    #' @return numeric characteristic size
    char_size = function(what = c("min", "max", "square", "circle")[4]){
      area <- sf::st_area(self$M)
      # area <- sf::st_inscribed_circle(x) |> sf::str_area()
      switch(tolower(what[1]),
             "min" = sqrt(area),
             "max" = sqrt(2*area),
             "circle" = sqrt(area/pi),
             "square" = sqrt(area))
    },
    
    #' @description generate kinematics metrics for a given mesh
    #'
    #' @param ofile optional file to save to (as non-spatial CSV)
    #' @return mesh (invisibly) with variables added including
    #' \itemize{
    #' \item max_u, max_v, max_w maximum velocity by element over all sigma in m/s
    #' \item area element area in m^2
    #' \item char_dim characteristic dimension, for now the radius of the circle with equivalent area in m
    #' }
    mesh_metrics = function(ofile = file.path(dirname(self$filename), "mesh-elem-metrics.csv.gz")){
      self$M$area <- sf::st_area(self$M)
      self$M$char_dim <- char_size(self$M)
      get_abs_max <- function(x){max(abs(x), na.rm = T)}
      self$M$max_u <- apply(ncdf4::ncvar_get(self$NC, "u"), 1, get_abs_max)
      self$M$max_v <- apply(ncdf4::ncvar_get(self$NC, "v"), 1, get_abs_max)
      self$M$max_w <- apply(ncdf4::ncvar_get(self$NC, "ww"), 1, get_abs_max)
      if (length(ofile[1]) > 0 && !is.na(ofile[1])){
        sf::st_drop_geometry(self$M) |>
          readr::write_csv(ofile[1])
      }
      invisible(self$M)
    },
    
    #' @description append the boundary info to the mesh table
    #' 
    #' @return mesh invisibly with appended boundary variable where
    #' \itemize{
    #'  \item{open indicates the element is boundaed on at least one side by open water}
    #'  \item{closed indicates the element is bounded on at least one side by shoreline}
    #'  \item{internal indicated the element is not on the boundary}
    #' }
    append_bounds = function(){
      b <- rep("internal", nrow(self$M))
      b[self$open_bounds] <- "open"
      b[self$closed_bounds] <- "closed"
      self$M <- dplyr::mutate(self$M, bounds = b) %>%
        dplyr::relocate(dplyr::starts_with("geometry"), .after = dplyr::last_col())
      invisible(self$M)
    }
    
    
  ), # public
  
  active = list(
    #' @field open_bounds provides the indices for open boundary elements (open water)
    open_bounds = function(){
      readRDS(system.file("extdata/elem_open.Rds", package = "necofs"))
    },
    #' @field closed_bounds provides the indices for closed boundary elements (shore)
    closed_bounds = function(){
      readRDS(system.file("extdata/elem_closed.Rds", package = "necofs"))
    }
  ), # active
  
  private = list(
    finalize = function(){
      ncdf4::nc_close(self$NC)
    },
    message = function(fmt, ...){
      if (self$verbose){
        if (missing(...)){
          cat(fmt, sep = "\n")
        } else {
          cat( sprintf(fmt, ...), sep = "\n")
        }
      } 
      invisible(NULL)
    }
  ) # private
  
  ) #NE_Physics


#' Generate  listing of one or more random points
#' 
#' @export
#' @param X NE_Physics object
#' @param n numeric the number of random points to generate
#' @param select_time character, one of 'first' (default), 'last' or 'random'
#' @return sf POINT XYZ object with \code{n} features
random_point <- function(X, 
                         n = 1,
                         select_time = c("first", "last", "random")[1]) {

  # given a 2 element vector [min,max] compute n random numbers within
  rand <- function(v, n = 1) {
    if (inherits(v, "POSIXt")){
      t0 <- v[1]
      v <- as.numeric(v)
      vd <- v - v[1]
      r <- stats::runif(n, min = vd[1], max = vd[2]) + t0
    } else {
      r <- stats::runif(n, min = v[1], max = v[2])
    }
    return(r)
  }
  
  zr <- range(X$NC$dim$siglev$vals[1,])

  p <- sf::st_sample(X$M, n)
  elem <- sf::st_contains(X$M, p, sparse = FALSE) |>
    apply(2, which)
  
  times <- X$get_time()
  thistime <- switch(tolower(select_time[1]),
    "first" = times[1],
    "last" = times[length(times)],
    rand(range(times),n))
  
  p <- sf::st_coordinates(p) |>
    dplyr::as_tibble() %>%
    dplyr::mutate(Z = rand(zr, n),
                  time = thistime, 
                  elem = elem) |>
    dplyr::relocate(dplyr::contains("elem"), .before = 1) |>
    sf::st_as_sf(coords = c("X", "Y", "Z"), 
                 dim = "XYZ",
                 crs = sf::st_crs(p))
  
  # if any are NA then try again
  while(any(is.na(p$elem))){
    ix <- is.na(p$elem)
    nbad <- sum(ix)
    nkeep <- n - nbad
    p <- p |>
      dplyr::filter(!ix) %>%
      dplyr::bind_rows(random_point(X, n = nbad))
  }
  p
}

#' Instantiate a NE_Physics R6 object
#' 
#' @export
#' @param filename character either a filename or OpenDAP URL
#' @param ... other arguments passed to \code{NE_Physics$new(filename, ...)}
#' @return NE_Physics R6 Reference Object
NECOFSPhysics <- function(filename = file.path("http://www.smast.umassd.edu:8080",
                                               "thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_GOM/2019",
                                               "gom4_201901.nc"),
                            ...){
  NE_Physics$new(filename, ...)
}