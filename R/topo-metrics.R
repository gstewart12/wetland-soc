
library("readxl")
library("sf")
library("stars")
library("terra")
library("whitebox")
wbt_init("~/WBT/whitebox_tools")
library("qgisprocess")
library("lubridate", warn.conflicts = FALSE)
library("tidyverse", warn.conflicts = FALSE)

# Packages used but not loaded: exactextractr, fasterize, lwgeom, smoothr

# Functions stored elsewhere
source("R/run-wbt.R")
source("R/unnest-depressions.R")


# 1. Read data -----------------------------------------------------------------

dir_in <- "data/spatial"
dir_out <- "output/spatial"

bc_dem <- rast(file.path(dir_in, "bc-dem.tif"))
jl_dem <- rast(file.path(dir_in, "jl-dem.tif"))
dem <- merge(bc_dem, jl_dem)

bc_aoi <- read_sf(file.path(dir_in, "jones-rd-tnc"))
jl_aoi <- read_sf(file.path(dir_in, "jackson-ln-tnc"))
aoi <- bind_rows(bc_aoi, jl_aoi) %>% 
  mutate(geometry = st_snap(geometry, geometry, tolerance = 0.5)) %>% 
  summarize()


# 2. Smooth DEM errors ---------------------------------------------------------

# Bizarrely, the SAGA Lee filter produces different results in R vs QGIS
# Running it in QGIS instead
# Info about Lee filter: 
# https://sourceforge.net/p/saga-gis/discussion/790705/thread/7fa509ec/

# 1. remove buildings, dense vegetation, dirt piles, etc.
#   a) WBT remove off terrain objects (filter = 9, slope = 1.5)
#   b) WBT remove off terrain objects (filter = 7, slope = 1.0)
#   c) WBT remove off terrain objects (filter = 5, slope = 0.5)
#   d) WBT remove off terrain objects (filter = 3, slope = 0.25)
# 2. smooth terrain
#   a) SAGA multi direction lee filter (default parameters)


# 3. Delineate wetlands --------------------------------------------------------

# Read data
bc_dem_cond <- rast(file.path(dir_out, "bc-dem-roto-lee1.tif"))
jl_dem_cond <- rast(file.path(dir_out, "jl-dem-roto-lee1.tif"))
dem_cond <- merge(jl_dem_cond, bc_dem_cond)

# Helper function for finding weirdly shaped depressions
st_narrowness <- function(x) {
  
  ic <- sf::st_inscribed_circle(x)
  ic <- ic[!sf::st_is_empty(ic)]
  
  as.numeric(sf::st_area(x) / sf::st_area(ic))
}

# Conduct stochastic depression analysis
# - unfortunately this cannot be exactly reproduced; no way to set seed
bc_pdep <- run_wbt(
  "stochastic_depression_analysis",
  dem = jones_dem_cond, rmse = 0.15, range = 10, iterations = 500
)
writeRaster(
  bc_pdep, file.path(dir_out, "spatial/bc-pdep.tif"), overwrite = TRUE
)

jl_pdep <- run_wbt(
  "stochastic_depression_analysis",
  dem = jl_dem_cond, rmse = 0.15, range = 10, iterations = 500
)
writeRaster(
  bc_pdep, file.path(dir_out, "spatial/jl-pdep.tif"), overwrite = TRUE
)

# Mask depression areas with standard pdep threshold (0.8)
pdep <- merge(bc_pdep, jl_pdep)
pdep_thr <- pdep %>%
  classify(
    rbind(c(0, 0.8, 0), c(0.8, 1, 1)), include.lowest = TRUE
  ) %>%
  # Clean up edges of depression areas
  run_wbt("majority_filter", input = ., filterx = 3, filtery = 3) %>%
  run_wbt("majority_filter", input = ., filterx = 3, filtery = 3) %>%
  run_wbt("majority_filter", input = ., filterx = 3, filtery = 3)

# Extract depressions as polygons
dep <- pdep_thr %>%
  classify(cbind(0, NA)) %>%
  st_as_stars() %>%
  st_as_sf(merge = TRUE) %>%
  transmute(id = row_number()) %>%
  # Remove depressions off property
  filter(st_intersects(geometry, st_buffer(aoi, 20), sparse = FALSE)) %>%
  # Fill small holes in polygons
  qgis_run_algorithm("native:deleteholes", INPUT = ., MIN_AREA = 50) %>%
  st_as_sf() %>%
  mutate(
    area = as.numeric(st_area(geom)),
    elev = exactextractr::exact_extract(dem_cond, geom, c("min", "max")),
    depth = elev$max - elev$min
  ) %>%
  # Remove small, irregularly-shaped, and shallow depressions
  # - depth & area thresholds are sensu [1]
  # - narrow "depressions" are likely ditch fragments
  filter(area > 50, st_narrowness(geom) < 5, depth > 0.10) %>%
  select(-elev)

# Identify separate depressions that are merged at high water levels
dep_nested_auto <- dep %>%
  st_set_geometry("geometry") %>%
  group_split(id) %>%
  # Takes ~1.5 mins to run
  map(
    unnest_depressions, dem_cond, n_slices = 10, min_area = 50, min_depth = 0.10
  ) %>%
  bind_rows()

dep_nested <- dep_nested_auto %>%
  # Add unique identifier for each depression
  mutate(
    dep_id = row_number(), 
    root_id = as.integer(str_extract(id, "\\d+")),
    root = !str_detect(id, "-"),
    leaf = !str_detect(id, "\\+") & str_detect(id, "-"),
    # Label properties
    site = if_else(
      st_intersects(geometry, st_buffer(bc_aoi, 200), sparse = FALSE)[, 1], 
      "BC", "JL"
    ),
    .before = 1
  ) %>%
  select(-id)

write_sf(dep_nested, file.path(dir_out, "deps.gpkg"))


# 4. Delineate catchments ------------------------------------------------------

# Read data
dep_nested <- read_sf(file.path(dir_out, "deps.gpkg"))
bc_dem_cond <- rast(file.path(dir_out, "bc-dem-roto-lee1.tif"))
jl_dem_cond <- rast(file.path(dir_out, "jl-dem-roto-lee1.tif"))
dem_cond <- merge(bc_dem_cond, jl_dem_cond)

# Flow pointer for catchment delineation
dem_flow_catch <- dem_cond %>%
  # Fix no-flow cells on slopes
  run_wbt("fill_single_cell_pits", dem = .) %>%
  # Drain puddles (depth < 10 cm; [2]) into surrounding basins
  run_wbt(
    "fill_depressions", dem = ., flat_increment = 1e-5, max_depth = 0.10
  ) %>%
  run_wbt("d8_pointer", dem = .)

catch <- dep_nested %>%
  # Need to separate levels to avoid overlapping
  group_split(level) %>%
  map(
    fasterize::fasterize, raster::raster(dem_cond), 
    field = "dep_id", background = 0
  ) %>%
  # Determine contributing areas with wetland cells as target points
  map(~ run_wbt("watershed", d8_pntr = dem_flow_catch, pour_pts = .x)) %>%
  map(st_as_stars) %>%
  map(st_as_sf, merge = TRUE) %>%
  map(rename, dep_id = 1) %>%
  bind_rows()

dep_catch <- catch %>%
  # Remove 8-connected fragments
  with_groups(dep_id, ~ slice_max(.x, st_area(geometry))) %>%
  left_join(select(st_drop_geometry(dep_nested), dep_id:level)) %>%
  arrange(dep_id)

write_sf(dep_catch, file.path(dir_out, "dep-catch.gpkg"))


# 5. Subset study wetlands -----------------------------------------------------

dep_nested <- read_sf(file.path(dir_out, "deps.gpkg"))
wetland_info <- read_csv("output/wetland-data-ak.csv")

# Subset study wetlands
wetland_pts <- wetland_info %>%
  select(wetland, site, X, Y) %>%
  # Add indicator for specific area within site
  mutate(
    site_area = case_when(
      wetland %in% c("AI", "AJ", "QB") ~ "south",
      wetland %in% c("AB", "DK", "GN", "ND") ~ "east",
      site == "BC" ~ "north", site == "JL" ~ "west"
    )
  ) %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
  st_transform(crs = 26918)
wetlands <- dep_nested %>%
  select(-site) %>%
  st_join(wetland_pts) %>%
  select(wetland:site_area, dep_id:level, geom) %>%
  drop_na(wetland) %>%
  # If more than one nested, use the bottom-level one
  group_by(wetland) %>%
  slice_max(level, n = 1) %>%
  ungroup() %>%
  st_as_sf()

write_sf(wetlands, file.path(dir_out, "wetlands.gpkg"))


# 6. Calculate metrics ---------------------------------------------------------

# Read data
wetlands <- read_sf(file.path(dir_out, "wetlands.gpkg"))
dep_nested <- read_sf(file.path(dir_out, "deps.gpkg"))
dep_catch <- read_sf(file.path(dir_out, "dep-catch.gpkg"))
bc_dem_cond <- rast(file.path(dir_out, "bc-dem-roto-lee1.tif"))
jl_dem_cond <- rast(file.path(dir_out, "jl-dem-roto-lee1.tif"))
dem_cond <- merge(jl_dem_cond, bc_dem_cond)

## 6.1. 2D shape metrics -------------------------------------------------------

shape_metrics <- dep_nested %>%
  mutate(
    area = as.numeric(st_area(geom)),
    # Smooth rough raster edges before measuring perimeter
    perimeter = geom %>%
      smoothr::smooth(method = "ksmooth", bandwidth = 3) %>%
      lwgeom::st_perimeter() %>%
      as.numeric() %>%
      round(),
    p_a_ratio = perimeter/area,
    shape_index = perimeter / (2 * sqrt(pi * area)),
    a_axis = as.numeric(st_length(run_wbt("polygon_long_axis", input = .))),
    b_axis = as.numeric(st_length(run_wbt("polygon_short_axis", input = .))),
    orientation = run_wbt("patch_orientation", input = .)$ORIENT
  ) %>%
  st_drop_geometry()

## 6.2. Terrain metrics --------------------------------------------------------

# Slope gradient
dem_slope <- run_wbt("slope", dem = dem_cond)
# Difference from mean elevation
dem_diff <- run_wbt(
  "diff_from_mean_elev", dem = dem_cond, filterx = 500, filtery = 500
)
# Deviation from mean elevation
dem_dev_500 <- run_wbt(
  "dev_from_mean_elev", dem = dem_cond, filterx = 500, filtery = 500
)

# Get elev info from merged wetlands
dep_merge <- dep_nested %>% 
  as_tibble() %>%
  filter(!root) %>%
  left_join(select(filter(dep_nested, root), root_id, merge = geom)) %>%
  st_as_sf(sf_column_name = "merge") %>%
  rowwise() %>%
  mutate(
    max_merge_elev = exactextractr::exact_extract(dem_cond, merge, "max")
  ) %>%
  st_drop_geometry() %>%
  select(dep_id:level, max_merge_elev) %>%
  ungroup()

terrain_metrics <- dep_nested %>%
  left_join(dep_merge) %>%
  mutate(
    elev = map(exactextractr::exact_extract(dem_cond, geom), pull, value),
    mean_elev = map_dbl(elev, mean),
    min_elev = map_dbl(elev, min),
    max_elev = map_dbl(elev, max),
    mean_depth = map_dbl(elev, ~ mean(max(.x) - .x)),
    max_depth = max_elev - min_elev,
    max_merge_depth = max_merge_elev - min_elev,
    volume = map_dbl(elev, ~ sum(max(.x) - .x)),
    mean_diff = exactextractr::exact_extract(dem_diff, geom, "mean"),
    mean_dev = exactextractr::exact_extract(dem_dev_500, geom, "mean"),
    mean_slope = exactextractr::exact_extract(dem_slope, geom, "mean")
  ) %>%
  select(-elev, -max_elev, -max_merge_elev) %>%
  st_drop_geometry()

## 6.3. Profile coefficients ---------------------------------------------------

fit_profile_coef <- function(x, dem, depth_res = 0.01, min_depth = 0,
                             use_max_area = FALSE,
                             type = c("area", "volume")) {
  
  type <- match.arg(type)
  if (!inherits(x, "SpatVector")) x <- terra::vect(x)
  
  elev <- x %>%
    terra::mask(terra::crop(dem, .), ., touches = FALSE) %>%
    terra::values(mat = FALSE) %>%
    na.omit()
  rel_elev <- elev - min(elev) # the real wetland depths
  
  depth <- seq(min_depth, max(rel_elev), by = depth_res) # depths to evaluate
  depth_bins <- cut(rel_elev, c(-Inf, depth), right = FALSE)
  
  area <- cumsum(unname(table(depth_bins)))
  
  # Brooks & Hayashi method ([3]; default: Hayashi & van der Kamp [4])
  if (use_max_area) s <- max(area)
  
  if (type == "area") {
    # Original equation does not divide by max depth, but I've found it helps
    # convergence and does not affect estimation of p
    fit <- nls(
      area ~ s * (depth/max(depth))^(2/p), 
      start = if (use_max_area) list(p = 1) else list(s = 1000, p = 1)
    )
  } else if (type == "volume") {
    volume <- vector("numeric", length(depth))
    for (i in 1:length(depth)) {
      volume[i] <- sum(depth[i] - (rel_elev[rel_elev <= depth[i]]))
    }
    fit <- nls(
      volume ~ (s/(1 + 2/p)) * (depth^(1 + 2/p)/max(depth)^(2/p)), 
      start = if (use_max_area) list(p = 1) else list(s = 1000, p = 1)
    )
  }
  
  fit
}

# Plot a single fit
wetlands %>%
  filter(wetland == "AJ") %>%
  fit_profile_coef(
    dem_cond, depth_res = 0.05, min_depth = 0, use_max_area = FALSE
  ) %>% 
  broom::augment() %>%
  ggplot(aes(y = depth)) +
  geom_point(aes(x = area)) +
  geom_line(aes(x = .fitted), color = "blue")

# Fit all wetlands
profiles <- dep_nested %>% 
  nest_by(dep_id) %>% 
  mutate(
    # Minke et al. use 5 cm increments; start at 0.1 m depth [5]
    # - I'm starting at 0 since DEM is conditioned and some wetlands are shallow
    # Technically I already know max area, but since that is an estimate based
    # on the DEM I let it vary here
    fit = list(possibly(fit_profile_coef, NULL)(
      data, dem_cond, depth_res = 0.05, min_depth = 0, use_max_area = FALSE
    ))
  ) %>%
  filter(!is.null(fit)) %>%
  mutate(
    fit_data = list(broom::augment(fit)),
    s = coef(fit)[1],
    p = coef(fit)[2],
    # Hayashi & van der Kamp use RMSE as goodness of fit [4]
    rmse = sqrt(mean((fit_data$area - fit_data$.fitted)^2)),
    # RMSE as percent of max area
    rmse_pct = rmse / max(fit_data$area) * 100,
    # RMSE as percent of area at 80% max depth [5]
    depth80 = pull(slice_min(fit_data, abs(0.8 - percent_rank(depth))), depth),
    fit_data80 = list(filter(fit_data, depth <= depth80)),
    rmse80 = sqrt(mean((fit_data80$area - fit_data80$.fitted)^2)),
    rmse80_pct = rmse / max(fit_data80$area) * 100,
    n = nrow(fit_data)
  ) %>%
  select(-data, -fit, -fit_data, -depth80, -fit_data80) %>%
  ungroup()
profiles
# Results are similar using area vs. volume (as they should be)
# - sticking with area - simpler

write_csv(profiles, file.path(dir_out, "dep-profile-coefs.csv"))

## 6.4. Catchment metrics ------------------------------------------------------

catch_metrics <- dep_nested %>%
  select(dep_id) %>%
  st_drop_geometry() %>%
  left_join(dep_catch) %>%
  st_as_sf() %>%
  mutate(
    # Area in ha
    area = as.numeric(st_area(geom)) * 1e-4,
    elev = map(exactextractr::exact_extract(dem_cond, geom), pull, value),
    catch_min_elev = map_dbl(elev, min),
    catch_max_elev = map_dbl(elev, max),
    catch_mean_elev = map_dbl(elev, mean),
    catch_sd_elev = map_dbl(elev, sd),
    catch_mean_slope = exactextractr::exact_extract(dem_slope, geom, "mean")
  ) %>%
  select(-elev) %>%
  rename(catch_area = area) %>%
  st_drop_geometry()

# Gather all metrics

dep_metrics <- dep_nested %>%
  st_drop_geometry() %>%
  left_join(shape_metrics) %>%
  left_join(terrain_metrics) %>%
  left_join(catch_metrics) %>%
  left_join(select(profiles, dep_id, p))

write_csv(dep_metrics, file.path(dir_out, "dep-topo-metrics.csv"))

# Subset of metrics for analysis
topo_vars <- c(
  "area", 
  "depth" = "max_depth", 
  "catchment" = "rel_ca_area", 
  "profile" = "p",
  "rtp" = "mean_dev", 
  "shape" = "shape_index"
)

wetland_metrics <- wetlands %>%
  st_drop_geometry() %>%
  left_join(dep_metrics) %>%
  mutate(rel_ca_area = catch_area * 1e4 / area) %>%
  select(wetland, all_of(topo_vars))

write_csv(wetland_metrics, file.path(dir_out, "wetland-topo-metrics.csv"))


# References -------------------------------------------------------------------

# [1]  Kiss J and Bedard-Haughn A 2021 Predictive mapping of solute‐rich 
#      wetlands in the Canadian Prairie Pothole Region through high‐resolution 
#      digital elevation model analyses Wetlands 41 38
# [2]  Li S, MacMillan R A, Lobb D A, McConkey B G, Moulin A and Fraser W R 2011 
#      Lidar DEM error analyses and topographic depression identification in a 
#      hummocky landscape in the prairie region of Canada Geomorphology 129 
#      263–75
# [3]  Brooks R T and Hayashi M 2002 Depth-area-volume and hydroperiod 
#      relationships of ephemeral (vernal) forest pools in southern New England 
#      Wetlands 22 247–55
# [4]  Hayashi M and van der Kamp G 2000 Simple equations to represent the 
#      volume–area–depth relations of shallow wetlands in small topographic 
#      depressions Journal of Hydrology 237 74–85
# [5]  Minke A G, Westbrook C J and van der Kamp G 2010 Simplified 
#      volume-area-depth method for estimating water storage of Prairie Potholes 
#      Wetlands 30 541–51



