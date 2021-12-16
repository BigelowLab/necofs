#' Track a particle
#' 
#' @export
#' @param X NE_Physics object
#' @param P0 starting point, see \code{NE_Physics$random_point}
#' @param tstep numeric, seconds between each iteration
#' @param tmax numeric, number of seconds to run for
#' @param show_progress logical, if TRUE then show a progress bar
#' @param filename an optional filename to save the track to or NA to skip
#' @param reverse logical, if TRUE run the track reverse in time.  It is an error
#'   to provide the seed point with a timestamp equivalent to the first timestep of the
#'   model.  See \code{\link{random_point}} to set the seed with other timestamps
#' @param drag numeric, a three element vector of drag factors in [x, y, z], like a sinking rate in z.
#'   Drag units must be the same as the units of \code{u}, \code{v}, and \code{ww}.
#' @param overwrite logical, if TRUE allow overwriting of existing files
#' @return sf object of type POINT
particle_track <- function(X, P0 = X$random_points(),
                           tstep = 60,
                           tmax = 3600 * 12,
                           reverse = FALSE,
                           drag = c(0,0,0),
                           show_progress = FALSE,
                           filename = c("particle_track.gpkg", NA)[2],
                           overwrite = TRUE){
  
  if (FALSE){
    P0 = X$random_points()
    tstep = 60
    tmax = 600
    reverse = FALSE
    drag = c(0,0,0)
    show_progress = FALSE
    filename = c("particle_track.gpkg", NA)[1]
  }
  
  TIMES <- X$get_time() # X$NC$dim$time$vals
  if (reverse && any(P0$time <= TIMES[1])){
    stop("if operating in reverse, the seed point(s) must be later than the first model time slice")
  }
  
  NMAX = tmax/tstep
  #preallocate the points list
  P <- vector(length = NMAX, mode = "list")

  #P0 <- P0 %>%
  #  dplyr::mutate(time = TIMES[1])
  
  N <- 1
  P[[N]] <- P0
  
  if (show_progress[1]){
    pb <- utils::txtProgressBar(min = 0, max = NMAX, initial = 0, style = 3)
  }
  
  if (reverse){
    tstep <- -1 * tstep
  }
  
  while (N <= NMAX){
    if (show_progress[1]) utils::setTxtProgressBar(pb, N)
    itime <- findInterval(P[[N]]$time, TIMES)
    # get the current point's u,v and w
    uvw <- fvcom::get_elem_var(X$NC, var = c("u", "v", "ww"),
                               elem = P[[N]]$elem,
                               time = itime) |>
      as.matrix()
    # and it's coordinates
    pn <- sf::st_coordinates(P[[N]])
    # translate (aka affine shift)
    d <- pn + ((uvw[,2:4] + drag) * tstep)
    # convert to sfc 
    g <- sf::st_sfc(sf::st_point(d), crs = sf::st_crs(X$M))
    # determine which element the pount belongs to
    ix <- lengths(sf::st_contains(X$M, g)) > 0
    elem <- which(ix)
    if (length(elem) == 0){
      breakplot
    }
    # add a new point
    P[[N+1]] <- sf::st_sf(dplyr::tibble(elem = elem, 
                                        time = P[[N]]$time[1] + tstep, 
                                        geometry = g))
    N <- N + 1
  } 
  
  if (show_progress) close(pb)
      
  P <- dplyr::bind_rows(P[seq_len(N)])
  
  if (nrow(P) > 0 && !is.na(filename[1])){
    if (file.exists(filename[1]) && overwrite) ok <- file.remove(filename[1])
    tf <- tempfile(fileext = ".gpkg")
    ok <- sf::write_sf(P, tf)
    dummy <- file.copy(tf, filename[1])
    dummy <- file.remove(tf)
  }
  P
}

#' Plot a track or series of tracks
#' 
#' @export
#' @param p sf POINT tibble.  If it has a 'track' variable (column) then each is plotted
#' @param X NE_Physics object
#' @param title character plot title
#' @param filename character or NA, optional output file as PNG
#' @param ext object that defines plot extent, by default \code{p} See \code{\link[sf]{plot}}
#' @param ... other arguments for \code{\link[grDevices]{png}} 
plot_track <- function(p, X = NULL, 
                       title = "Particle Track",
                       filename = c(NA,"particle_track.png")[1],
                       ext = p,
                       ...){
  
  cols <- grDevices::palette.colors()
  if (!is.na(filename[1])){
    grDevices::png(filename[1], ...)
  }
  
  plot(sf::st_geometry(ext),
       xlab = 'Easting (m)', 
       ylab  = 'Northing (m)',
       axes = TRUE,
       main = title,
       extent = ext,
       col = "#FFFFFF")

  if (!is.null(X)){
    plot(sf::st_geometry(X$M),
         border = cols[['gray']],
         add = TRUE)
  }
  
  if (!("track" %in% colnames(p))){
    p <- p |>
      dplyr::mutate(track = 1)
  }
  
  n <- length(unique(p$track))
  colors <- cols[p$track]
  p <- p |>
    dplyr::group_by(.data$track) %>%
    dplyr::group_map(
      function(x, key){
        i <- x$track[1]
        plot(sf::st_geometry(x),
             col = cols[i+1],
             pch = ".",
             type = "l",
             lwd = 2,
             add = TRUE)
        plot(sf::st_geometry(x |> dplyr::slice(c(1,dplyr::n()))),
             col = cols[i+1],
             pch = c(1, 19),
             cex = 1.5,
             add = TRUE)
        p
      }, .keep = TRUE ) 
  
  
  if (!is.na(filename[1])){
    ok <- grDevices::dev.off()
  }
  invisible(NULL)
}


#' Load an example track
#' 
#' @export
#' @return sf POINT tibble. 
example_track <- function(){
  filename <- system.file(file.path("extdata", "particle_track_1200_604800.geojson"),
                          package = "necofs")
  suppressWarnings(sf::read_sf(filename) |>
    sf::st_set_crs("+init=nad83:1802"))
} 