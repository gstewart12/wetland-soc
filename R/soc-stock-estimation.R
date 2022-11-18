
library("mpspline2")
library("readxl")
library("lubridate")
library("tidyverse")

# Packages used but not loaded: broom, rsample


# 1. Read data -----------------------------------------------------------------

dir_in <- "data/soil"
dir_out <- "output/soil"

# Forested wetland LOI, CN, and bulk density data
loi_data <- read_csv(file.path(dir_in, "loi_all.csv"))
cn_data <- read_csv(file.path(dir_in, "cn.csv"))
bd_data <- read_csv(file.path(dir_in, "bd_all.csv"))
loi_emergent <- read_csv("output/loi-data-gs.csv")


# 2. Predict SOC from LOI ------------------------------------------------------

loi_c <- left_join(loi_data, cn_data)

# Emergent wetland LOI samples
loi_samples <- loi_emergent %>%
  separate(
    sample, c("wetland", "zonecore", "top_cm", "bottom_cm"), 
    remove = FALSE, convert = TRUE
  ) %>%
  # Convert to %
  mutate(loi_mean = loi_mean * 100)

# Check linear relationship
loi_c %>%
  ggplot(aes(loi_mean, carbon_p)) +
  geom_point(na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, na.rm = TRUE)

# Define linear fit between LOI & C
loi_c_fit <- lm(carbon_p ~ loi_mean, data = loi_c)
summary(loi_c_fit) # R2 = 0.987

# Predict C in emergent samples using forested fit
loi_c_preds <- loi_samples %>% 
  broom::augment(loi_c_fit, newdata = .) %>%
  # Set negative values to 0
  # - introduces small bias but prevents large problems later on
  mutate(c_pred = pmax(.fitted, 0), .keep = "unused")


# 3. Predict BD from LOI -------------------------------------------------------

loi_bd <- loi_data %>%
  left_join(bd_data) %>%
  # Only use data from "core" or "mccauley" ("auger" is from SSURGO)
  filter(collectionmethod != "auger") %>%
  # Average duplicate cores
  group_by(wetland, horizon_class, fixed_top_cm, fixed_bottom_cm) %>% 
  summarize(
    across(c(loi_mean, bd_gpercm3), mean, na.rm = TRUE), .groups = "drop"
  )

# Check nonlinear relationship
loi_bd %>%
  mutate(horizon_class = factor(horizon_class)) %>%
  ggplot(aes(loi_mean, bd_gpercm3)) +
  geom_point(aes(color = horizon_class)) +
  geom_smooth(
    method = nls, formula = y ~ 1 / (x/a + (1 - x)/b),
    method.args = list(start = c(a = 1, b = 1)), se = FALSE
  )

# BD prediction method is mostly following [1]

# "Ideal mixing model"
# - technically should use nonlinear quantile regression because of
#   heteroskedasticity, but checked and hardly makes a difference for parameter
#   estimates
loi_bd_fit <- nls(
  bd_gpercm3 ~ 1 / (loi_mean/k1 + (1 - loi_mean)/k2), 
  data = loi_bd, start = c(k1 = 1, k2 = 1)
)
summary(loi_bd_fit)

# Bootstrap resampling for confidence intervals on predictions
set.seed(56660)
loi_bd_fit_boot <- loi_bd %>%
  rsample::bootstraps(times = 1000) %>%
  mutate(
    # Fit nonlinear model on all bootstrap samples
    fit = map(
      splits, 
      ~ nls(
        bd_gpercm3 ~ 1 / (loi_mean/k1 + (1 - loi_mean)/k2), 
        data = rsample::analysis(.x), start = c(k1 = 1, k2 = 1)
      )
    ),
    k1 = map_dbl(fit, ~ coef(.x)[1]),
    k2 = map_dbl(fit, ~ coef(.x)[2]),
    # Predict on emergent data
    pred = map(fit, broom::augment, newdata = loi_samples)
  ) 
loi_bd_fit_boot %>% summarize(across(c(k1, k2), quantile, c(0.025, 0.975)))

# Calculate confidence intervals
loi_bd_pred_int <- loi_bd_fit_boot %>%
  unnest(pred) %>%
  group_by(sample, loi_mean) %>% 
  summarize(
    # Percentile method
    bd_pred_lower = quantile(.fitted, 0.025), 
    bd_pred_upper = quantile(.fitted, 0.975),
    .groups = "drop"
  )

# Predict BD in emergent samples using forested fit
loi_bd_preds <- loi_samples %>% 
  broom::augment(loi_bd_fit, newdata = .) %>%
  rename(bd_pred = .fitted) %>%
  left_join(loi_bd_pred_int)


# 4. Estimate 1-m SOC stocks ---------------------------------------------------

# Estimate C density in emergent samples
c_dens_pred <- loi_c_preds %>%
  left_join(loi_bd_preds) %>%
  # Predicted C conc. in g g-1, C density in g cm-3
  mutate(
    c_dens_pred = c_pred/100 * bd_pred,
    c_dens_pred_lower = c_pred/100 * bd_pred_lower,
    c_dens_pred_upper = c_pred/100 * bd_pred_upper
  )

write_csv(c_dens_pred, file.path(dir_out, "soc-dens-pred-emergent.csv"))

# Convert horizons to fixed intervals
c_dens_fixed <- c_dens_pred %>%
  # Group samples by core 
  # - mpspline requires columns 1, 2, and 3 to be sample, upper, lower
  select(-sample) %>%
  separate(zonecore, c("zone", "core"), 1) %>%
  unite("sample", c(wetland, zone, core)) %>%
  # Lumping 50-75 and 75-100 intervals so all cores get estimates down to 1 m
  mpspline_tidy(
    var_name = "c_dens_pred", d = c(0, 10, 30, 50, 100), lam = 1
  ) %>%
  map(as_tibble) %>%
  map(separate, sample, c("wetland", "zone", "core")) %>%
  map_at(1:3, rename, top_cm = UD, bottom_cm = LD) %>%
  map(rename_with, tolower)

# mpspline2 output:
#   1. predicted values over input depth ranges
#   2. predicted values over output depth ranges
#   3. 1-cm increment predictions
#   4. RMSE & IQR-scaled RMSE values

# Sum over cores for total 1-m SOC stocks
c_stock_est <- c_dens_fixed %>%
  pluck("est_dcm") %>%
  rename(c_dens_pred = splined_value) %>%
  mutate(depth_cm = bottom_cm - top_cm) %>%
  group_by(wetland, zone, core) %>%
  summarize(
    c_stock_0_100_est = sum(c_dens_pred * depth_cm * 10), .groups = "drop"
  )

write_csv(c_stock_est, file.path(dir_out, "soc-stocks-est-core-emergent.csv"))

# C stock estimates for each wetland/zone, with CIs from sampling variability
c_stock_est_ci <- c_stock_est %>% 
  group_by(wetland, zone) %>% 
  summarize(
    sd_core = sd(c_stock_0_100_est), 
    c_stock_0_100_est = mean(c_stock_0_100_est), 
    n = n(),
    .groups = "drop"
  ) %>% 
  # Parametric CIs because only 3 reps
  mutate(
    c_stock_0_100_lci = c_stock_0_100_est - sd_core / sqrt(n) * qnorm(0.975), 
    c_stock_0_100_uci = c_stock_0_100_est + sd_core / sqrt(n) * qnorm(0.975)
  ) %>%
  select(wetland, zone, n_cores = n, everything())

write_csv(c_stock_est_ci, file.path(dir_out, "soc-stocks-est-emergent.csv"))


# References -------------------------------------------------------------------

# [1]  Holmquist J R, Windham-Myers L, Bliss N, Crooks S, Morris J T, Megonigal 
#      J P, Troxler T, Weller D, Callaway J, Drexler J, Ferner M C, Gonneea M E, 
#      Kroeger K D, Schile-Beers L, Woo I, Buffington K, Breithaupt J, Boyd B M, 
#      Brown L N, Dix N, Hice L, Horton B P, MacDonald G M, Moyer R P, Reay W, 
#      Shaw T, Smith E, Smoak J M, Sommerfield C, Thorne K, Velinsky D, Watson 
#      E, Grimes K W and Woodrey M 2018 Accuracy and precision of tidal wetland 
#      soil carbon mapping in the conterminous United States Sci Rep 8 9478

