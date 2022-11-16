
# Allows Whitebox ins/outs to be passed directly to/from R environment
# - handles writing objects to temp files for Whitebox to access
# - handles reading and returning output file(s) from Whitebox processing

run_wbt <- function(tool, ..., sf_out = TRUE) {
  
  # Check output type in WhiteBox tool database
  data(wbttoolparameters, envir = environment())
  if (!stringr::str_starts(tool, "wbt_")) {
    tool <- paste0("wbt_", tool)
  }
  outputs <- wbttoolparameters %>%
    dplyr::filter(function_name == tool, parameter_class == "NewFile")
  
  # Find arguments that must be written as files to be accessed
  args <- list(...)
  
  rst_args <- purrr::keep(args, ~ inherits(.x, c("RasterLayer", "SpatRaster")))
  if (length(rst_args)) {
    temp_rst <- do.call(tempfile, list(fileext = rep(".tif", length(rst_args))))
    # Write rasters to file
    purrr::walk2(rst_args, temp_rst, terra::writeRaster)
    # Replace raster arguments with temp file names
    args <- purrr::list_modify(
      args, !!!as.list(rlang::set_names(temp_rst, names(rst_args)))
    )
  }
  
  vec_driver <- ".shp"
  vec_args <- purrr::keep(args, ~ inherits(.x, c("sf", "SpatVector")))
  if (length(vec_args)) {
    temp_vec <- do.call(tempfile, list(fileext = rep(vec_driver, length(vec_args))))
    is_sf <- purrr::map_lgl(vec_args, inherits, "sf")
    # Write vectors to file
    if (all(is_sf)) {
      purrr::walk2(vec_args, temp_vec, sf::write_sf)
    } else if (all(!is_sf)) {
      purrr::walk2(vec_args, temp_vec, terra::writeVector)
    } else {
      stop("Vectors must be either all sf or all SpatVector.")
    }
    
    # Replace vector arguments with temp file names
    args <- purrr::list_modify(
      args, !!!as.list(rlang::set_names(temp_vec, names(vec_args)))
    )
  }
  
  # Prepare for output
  if (!nrow(outputs)) {
    # If no listed output type, the function must update the file in place
    # - this will fail in the unlikely case that the output is a raster
    # --> multi-output functions also have no listed outputs...
    read_files <- args$input
    read_funs <- if (sf_out) list(sf::read_sf) else list(terra::vect)
  } else {
    output_names <- dplyr::pull(outputs, argument_name)
    output_types <- dplyr::pull(outputs, parameter_detail)
    vector_types <- "Polygon|Line|Point"
    read_files <- output_types %>%
      stringr::str_replace(
        "Raster", purrr::partial(tempfile, fileext = ".tif")
      ) %>%
      stringr::str_replace(
        vector_types, purrr::partial(tempfile, fileext = vec_driver)
      ) %>%
      as.list() %>%
      rlang::set_names(output_names)
    read_funs <- vector("list", length(read_files))
    read_funs <- replace(
      read_funs, stringr::str_which(output_types, "Raster"), list(terra::rast)
    )
    read_funs <- replace(
      read_funs, stringr::str_which(output_types, vector_types),
      if (sf_out) list(sf::read_sf) else list(terra::vect)
    )
    args <- c(args, read_files)
  }
  
  rlang::exec(tool, !!!args)
  
  out <- purrr::map2(read_funs, read_files, ~ rlang::exec(.x, .y))
  if (nrow(outputs)) {
    out <- rlang::set_names(
      out, stringr::str_remove(outputs$argument_name, "out_")
    )
  }
  
  # Unlist if single output
  if (length(out) == 1) out <- out[[1]]
  out
}
