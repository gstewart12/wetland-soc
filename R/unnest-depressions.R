
# This is an adaptation of Nate Jones' work...
# - https://github.com/FloodHydrology/DMV_Spatial_Analysis
# ...which is an adaptation of Wu et al. (2019)
# - https://github.com/giswqs/lidar/blob/master/lidar/slicing.py
# - recommended slicing interval = 2-3x (no smaller than) DEM vert. accuracy


unnest_depressions <- function(x, dem, merge_elev = NA, n_slices = 10, 
                               slice_interval = NA, by_quantiles = FALSE, 
                               min_area = 250, min_depth = 0.1, max_level = 3) {
  
  # Coerce to sf if necessary
  if (!inherits(x, "sf")) x <- sf::st_as_sf(x, sf_column_name = "geom")
  
  dem_dep <- x %>%
    # Only need SpatVector format for this step
    terra::vect() %>%
    terra::mask(terra::crop(dem, .), ., touches = FALSE)
  
  elev <- na.omit(terra::values(dem_dep, mat = FALSE))
  
  # If merge elevation is known then don't need to estimate
  # - Note: setting the merge elevation will return the original polygon if
  #   other criteria are not met
  if (!all(is.na(merge_elev))) {
    slice_elevs <- c(min(elev), merge_elev, max(elev))
    n_slices <- length(slice_elevs) - 1
  } else {
    if (!is.na(slice_interval)) {
      slice_elevs <- seq(min(elev), max(elev), by = slice_interval)
    } else {
      if (by_quantiles) {
        # Quantiles help pinpoint the exact merge elevation
        slice_elevs <- quantile(elev, seq(0, 1, length.out = n_slices + 1))
      } else {
        slice_elevs <- seq(min(elev), max(elev), length.out = n_slices + 1)
      }
    }
  }
  
  slices <- dem_dep %>%
    terra::classify(slice_elevs, include.lowest = TRUE) %>%
    as.factor()
  # Start at 1 instead of 0
  slices <- slices + 1
  
  nodes <- x %>%
    dplyr::transmute(
      id = as.character(id),
      has_children = FALSE,
      level = 0,
      area = as.numeric(sf::st_area(geometry)),
      depth = diff(range(elev)),
      merge_elev = NA_real_
    ) %>%
    list()
  
  for (i in (n_slices - 1):1) {
    
    objects <- slices %>%
      terra::clamp(upper = i, values = FALSE) %>%
      terra::patches(directions = 4) %>%
      terra::as.polygons() %>%
      sf::st_as_sf() %>%
      dplyr::mutate(
        area = as.numeric(sf::st_area(geometry)),
        elev = exactextractr::exact_extract(
          dem_dep, geometry, c("min", "max"), progress = FALSE
        ),
        depth = elev$max - elev$min
      ) %>%
      dplyr::select(-elev) %>%
      dplyr::filter(area >= min_area, depth >= min_depth)
    
    # Skip if no new children; quit when pieces are all crumbs
    if (!nrow(objects)) {
      break
    } else if (nrow(objects) < 2) {
      next
    }
    
    new_nodes <- list()
    for (j in 1:length(nodes)) {
      
      parent <- nodes[[j]]
      
      # Skip if already has children or maxed out levels
      if (parent$has_children | parent$level >= max_level) {
        next
      }
      
      # Get children within node
      node_slice <- dplyr::filter(
        objects, sf::st_within(geometry, parent, sparse = FALSE)
      )
      
      # Add any new children to list
      if (nrow(node_slice) > 1) {
        
        children <- node_slice %>%
          dplyr::mutate(
            # ID derived from parent's ID
            id = paste0(parent$id, "-", dplyr::row_number()),
            has_children = FALSE,
            level = stringr::str_count(id, "-"),
            merge_elev = slice_elevs[i + 1]
          )
        
        nodes[[j]]$has_children <- TRUE
        new_nodes <- c(new_nodes, dplyr::group_split(children, id))
        
      } else {
        next
      }
    }
    
    # Add all new nodes at this elevation slice
    nodes <- c(nodes, new_nodes)
  }
  
  # Return data frame of unnested depressions w/ basin info
  nodes %>% 
    dplyr::bind_rows() %>% 
    # Add "has children" indicator to ID
    dplyr::mutate(id = if_else(has_children, paste0(id, "+"), id)) %>%
    dplyr::select(id, level, area, depth, merge_elev)
}
