
library("rdacca.hp")
library("car")
library("readxl")
library("patchwork")
library("lubridate")
library("tidyverse")

# Packages used but not loaded: corrr, gstat, khroma, sf, stargazer


# 1. Read data -----------------------------------------------------------------

core_ak <- read_csv("output/core-data-ak.csv")

wetlands <- sf::read_sf("data/spatial/wetlands.gpkg")
wetland_pts <- sf::st_centroid(wetlands)

soc_data <- core_ak %>%
  group_by(wetland, site) %>%
  summarize(across(-rep, mean), .groups = "drop") %>%
  select(wetland, site, c_stock = C_stock_100cm) %>%
  bind_cols(sf::st_coordinates(wetland_pts))

water_level_daily <- read_csv("output/water-level-daily-gapfilled.csv")
hydro_metrics <- read_csv("output/wetland-hydro-metrics.csv") %>%
  rename(duration = days_wet, exposure = date_dry, recession = rec_rate)
hydro_vars <- names(hydro_metrics)[-1]

topo_metrics <- read_csv("output/wetland-topo-metrics.csv") %>%
  rename(catchment = rel_ca_area, profile = p, shape = shape_index)
topo_vars <- names(topo_metrics)[-1]

climate <- read_csv("output/climate-vars-daily-fntower-2019.csv")

format_metric_names <- function(x, add_units = FALSE) {
  levels <- set_names(c(
    "duration", "exposure", "iqr", "max", "median", "recession", "area", 
    "catchment", "depth", "profile", "rtp", "shape"
  ))
  names(levels)[c(3, 11)] <- c("IQR", "RTP")
  if (add_units) {
    units <- c(
      "days", "d.o.y.", "m", "m", "m", "cm~day^{-1}", "m^2", "", "m", "", "", ""
    )
    units <- if_else(units != "", paste0("~(", units, ")"), units)
    names(levels) <- paste0(names(levels), units)
  }
  fct_recode(x, !!!levels)
}


# 2. Summary of all variables --------------------------------------------------

## 2.1. Statistics -------------------------------------------------------------

# SOC stocks
ggplot() + 
  geom_qq(aes(sample = scale(soc_data$c_stock))) + 
  geom_abline(color = "red")
shapiro.test(soc_data$c_stock)
plot(gstat::variogram(c_stock ~ 1, locations = ~ X + Y, data = soc_data))

# Hydro metrics
hydro_metrics %>%
  rename(maximum = max) %>%
  summarize(across(-1, list(mean = mean, min = min, max = max, sd = sd))) %>%
  pivot_longer(everything(), names_to = c("metric", "st"), names_sep = "_") %>%
  pivot_wider(names_from = st, values_from = value) %>%
  mutate(cv = sd/abs(mean) * 100)
hydro_metrics %>%
  pivot_longer(-wetland) %>%
  ggplot(aes(wetland, value)) +
  facet_wrap(~ name, scales = "free") +
  geom_point()

# Topo metrics
topo_metrics %>%
  rename(ca = catchment) %>%
  summarize(across(-1, list(mean = mean, min = min, max = max, sd = sd))) %>%
  pivot_longer(everything(), names_to = c("metric", "st"), names_sep = "_") %>%
  pivot_wider(names_from = st, values_from = value) %>%
  mutate(cv = sd/abs(mean) * 100)
topo_metrics %>%
  pivot_longer(-wetland) %>%
  ggplot(aes(wetland, value)) +
  facet_wrap(~ name, scales = "free") +
  geom_point()
# catchment has an extreme value - should be transformed
# - take the inverse since it's already a ratio
# - take negative to preserve direction of relationships
topo_metrics %>% ggplot(aes(wetland, -1/catchment)) + geom_point()

## 2.2. Distributions ----------------------------------------------------------

topo_metrics %>%
  left_join(hydro_metrics) %>%
  pivot_longer(-wetland) %>%
  ggplot(aes(value)) +
  facet_wrap(~ name, scales = "free") +
  geom_density()
summarize(hydro_metrics, across(-wetland, ~ shapiro.test(.x)$p.value))
summarize(topo_metrics, across(-wetland, ~ shapiro.test(.x)$p.value))

topo_metrics %>%
  left_join(hydro_metrics) %>%
  pivot_longer(-wetland) %>%
  ggplot() +
  facet_wrap(~ name, scales = "free") +
  geom_qq(aes(sample = scale(value))) +
  geom_abline(color = "red")
# might benefit from transformation: max -> log, profile -> log, area -> sqrt

# Compare all metrics between sites
soc_data %>%
  select(wetland, site) %>%
  left_join(topo_metrics) %>%
  mutate(catchment = -1/catchment) %>%
  left_join(hydro_metrics) %>%
  mutate(across(c(-wetland, -site), scale)) %>%
  pivot_longer(c(-wetland, -site)) %>%
  mutate(group = if_else(name %in% hydro_vars, "hydro", "topo")) %>%
  ggplot(aes(name, value, fill = site)) +
  facet_wrap(~ group, scales = "free_x") +
  geom_boxplot()


# 3. Correlations --------------------------------------------------------------

# Pairwise correlations between all metrics
corr_data <- hydro_metrics %>%
  full_join(topo_metrics) %>%
  mutate(catchment = -1/catchment)

corr_plot <- corr_data %>%
  rename(IQR = iqr, RTP = rtp) %>%
  relocate(duration, exposure, IQR, max, median, .after = 1) %>%
  relocate(catchment, depth, profile, RTP, shape, .after = last_col()) %>%
  corrr::correlate(use = "pairwise.complete.obs", method = "pearson") %>%
  corrr::shave(upper = TRUE) %>%
  corrr::stretch(na.rm = FALSE, remove.dups = FALSE) %>%
  mutate(
    across(
      c(x, y), fct_relevel, "shape", "RTP", "profile", "depth", "catchment",
      "area", "recession", "median", "max", "IQR", "exposure", "duration"
    ),
    x = fct_rev(x)
  ) %>%
  filter(y != "duration", x != "shape") %>%
  ggplot(aes(x, y, fill = r)) +
  khroma::scale_fill_BuRd(na.value = "#FFFFFF") +
  coord_equal() +
  geom_tile(height = 0.9, width = 0.9, alpha = 0.75) +
  # geom_text(aes(label = round(r, 2)), size = 3) +
  geom_text(
    aes(label = na_if(format(round(r, digits = 2), nsmall = 1), "   NA")), 
    size = 3.30
  ) +
  annotate(
    "rect", xmin = 0.5, xmax = 6.5, ymin = 0.5, ymax = 6.5, 
    fill = "transparent", color = "black", linewidth = 0.5
  ) +
  labs(x = NULL, y = NULL, fill = expression(italic("r"))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.text = element_text(face = "italic", size = 10),
    axis.ticks = element_blank(), 
    legend.key.height = unit(0.04, "npc"),
    legend.position = c(0.9, 0.75),
    panel.border = element_blank(), 
    panel.grid = element_blank()
  )
corr_plot

ggsave(
  "reports/plots/hydro-topo-corrplot.pdf", corr_plot, device = cairo_pdf,
  width = 6.00, height = 6.00, units = "in"
)

# Pairwise scatterplots between all topo/hydro metric combos
scattermat <- hydro_metrics %>%
  left_join(topo_metrics) %>%
  select(-wetland) %>%
  mutate(catchment = -1/catchment) %>%
  mutate(across(.fns = scale)) %>%
  pivot_longer(
    all_of(hydro_vars), names_to = "hydro_var", values_to = "hydro_value"
  ) %>%
  pivot_longer(
    all_of(topo_vars), names_to = "topo_var", values_to = "topo_value"
  ) %>%
  mutate(
    hydro_var = fct_recode(hydro_var, IQR = "iqr"),
    topo_var = fct_recode(topo_var, RTP = "rtp")
  ) %>%
  ggplot(aes(hydro_value, topo_value)) +
  facet_grid(
    rows = vars(topo_var), cols = vars(hydro_var), 
    scales = "free", switch = "both"
  ) +
  geom_point(shape = 21) +
  geom_smooth(se = FALSE, size = 0.4, span = 1.5) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  scale_x_continuous(expand = c(0.05, 0.05)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside"
  )
scattermat

ggsave(
  "reports/plots/hydro-topo-scattermat-supp.pdf", scattermat, 
  device = cairo_pdf, width = 6.00, height = 6.00, units = "in"
)


# 4. Single regression ---------------------------------------------------------

# Pairwise relationships between metrics & SOC stocks

# Commonality metrics
metricplot_data <- soc_data %>%
  left_join(topo_metrics) %>%
  left_join(hydro_metrics) %>%
  select(
    wetland, site, c_stock, recession, max, duration, profile, depth, shape
  ) %>%
  pivot_longer(-(1:3), names_to = "metric", values_to = "value") %>% 
  mutate(
    group = if_else(metric %in% hydro_vars, "hydro", "topo"),
    metric = fct_inorder(fct_recode(
      metric,
      `Inundation~duration~(days)` = "duration",
      `Maximum~water~level~(m)` = "max",
      `Recession~rate~(cm~day^{-1})` = "recession",
      `Depth~(m)` = "depth",
      `Basin~profile~coefficient` = "profile",
      `Shape~index` = "shape"
    )),
    well = wetland %in% hydro_metrics$wetland
  ) 
metricplot <- metricplot_data %>%
  ggplot(aes(value, c_stock, fill = site, linetype = site)) +
  facet_wrap(
    vars(metric), scales = "free_x", strip.position = "bottom", 
    labeller = "label_parsed"
  ) +
  geom_point(aes(shape = well), size = 2.35, alpha = 0.7, color = "grey10") +
  geom_smooth(
    method = "lm", formula = y ~ x, se = FALSE, size = 0.4,
    alpha = 0.5, color = "grey20"
  ) +
  # Below chunk is for adding r^2 values to plots
  # geom_text(
  #   data = metricplot_data %>% 
  #     group_by(metric, site) %>% 
  #     summarize(
  #       r2 = round(cor(value, c_stock, use = "complete.obs")^2, 2), 
  #       .groups = "drop"
  #     ) %>%
  #     mutate(
  #       x = c(-0.95, -1.1, 0.85, 0.98, rep(NA, 8)),
  #       y = c(45, 29, 46, 34, rep(NA, 8))
  #       #label = paste0("R^2 = ", r2)
#     ),
#   aes(x = x, y = y, label = paste0("italic(r^{2}==", r2, ")")), 
#   parse = TRUE, size = 2.15, color = "grey20"
# ) +
  geom_text(
    data = metricplot_data %>% 
      group_by(metric) %>% 
      summarize(
        x = max(value, na.rm = TRUE) - diff(range(value, na.rm = TRUE)) * 0.02, 
        y = max(c_stock) - diff(range(c_stock)) * 0
      ) %>% 
      mutate(label = c("a)", "b)", "c)", "d)", "e)", "f)")),
    aes(x, y, label = label), size = 3.75, nudge_y = 0.5, inherit.aes = FALSE
  ) +
  scale_y_continuous(expand = expansion(c(0.05, 0.075))) +
  scale_fill_manual(values = c("#44aa99", "#332288"), name = "Site") +
  scale_color_manual(values = c("#44aa99", "#332288"), name = "Site") +
  scale_shape_manual(values = c(23, 21)) +
  scale_linetype_manual(values = c("12", "64"), name = "Site") +
  guides(shape = "none", fill = guide_legend(override.aes = list(shape = 21))) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 0.6) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 0.6) +
  labs(x = NULL, y = bquote("SOC stock"~(kg~C~m^-2))) +
  theme_bw(base_size = 12) +
  theme(
    legend.background = element_blank(),
    legend.box.margin = margin(),
    legend.key.height = unit(0.04, "npc"),
    legend.position = c(0.06, 0.92),
    legend.spacing.x = unit(0.004, "npc"), 
    legend.spacing.y = unit(0.005, "npc"),
    legend.title = element_text(size = 11),
    panel.grid = element_blank(), 
    panel.border = element_blank(),
    strip.background = element_blank(), 
    strip.placement = "outside", 
    strip.switch.pad.wrap = unit(0.004, "npc"),
    strip.text = element_text(margin = margin(3, 3, 3, 3))
  )
metricplot

ggsave(
  "reports/plots/soc-metrics-pairwise-plot.pdf", metricplot, 
  device = cairo_pdf, width = 6.5, height = 4.5, units = "in"
)

# Other metrics
metricplot_supp <- soc_data %>% 
  select(-X, -Y) %>%
  full_join(hydro_metrics) %>% 
  full_join(topo_metrics) %>% 
  select(-duration, -max, -recession, -depth, -profile, -shape) %>%
  mutate(catchment = -1/catchment) %>%
  pivot_longer(-(1:3)) %>% 
  mutate(
    well = wetland %in% hydro_metrics$wetland,
    group = if_else(name %in% hydro_vars, "hydro", "topo"),
    name = format_metric_names(name, add_units = TRUE),
    name = fct_reorder2(name, name, group, .fun = last2, .desc = FALSE)
  ) %>%
  ggplot(aes(value, c_stock, fill = site)) + 
  facet_wrap(
    vars(name), ncol = 3, scales = "free_x", labeller = "label_parsed", 
    dir = "h", strip.position = "bottom"
  ) + 
  geom_point(aes(shape = well), size = 2.35, alpha = 0.7, color = "grey10") + 
  # geom_smooth(method = "lm", se = FALSE, size = 0.3) +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
  scale_fill_manual(values = c("#44aa99", "#332288")) +
  scale_shape_manual(values = c(23, 21)) +
  guides(
    shape = "none", 
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  labs(x = NULL, y = bquote("SOC stock"~(kg~C~m^-2))) +
  theme_bw(base_size = 12) +
  theme(
    legend.background = element_blank(),
    legend.margin = margin(), 
    legend.position = "right", 
    legend.spacing.x = unit(0.004, "npc"),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside", 
    strip.text = element_text(margin = margin(3, 3, 3, 3))
  )
metricplot_supp

ggsave(
  "reports/plots/soc-metrics-pairwise-plot-supp.pdf", metricplot_supp, 
  device = cairo_pdf, width = 6.5, height = 4.0, units = "in"
)


# 5. Multiple regression -------------------------------------------------------

## 5.1. Hydrology model --------------------------------------------------------

soc_hydro <- soc_data %>% 
  inner_join(hydro_metrics) %>%
  mutate(site = factor(site)) %>%
  as.data.frame()

# (Multi)collinearity among metrics

# PCA biplot
# biplot(prcomp(soc_hydro[, 6:11], scale. = TRUE))
# Correlations
cor(select(soc_hydro, all_of(hydro_vars)))
# Variance inflation factors
vif(lm(
  c_stock ~ exposure + duration + max + median + iqr + recession, 
  data = soc_hydro, na.action = na.fail
))
# High multicollinearity: max-median-exposure, duration-median-exposure, max-iqr

hydro_fit <- lm(
  c_stock ~ duration + max + recession + site, 
  data = soc_hydro, na.action = na.fail
)

### 5.1.1. Diagnostics ---------------------------------------------------------

# Residual normality
hydro_qq <- ggplot(hydro_fit) +
  geom_abline(linetype = "dashed") +
  stat_qq(aes(sample = .stdresid)) +
  labs(x = "Theoretical quantiles", y = "Std. residuals")
hydro_qq
shapiro.test(rstandard(hydro_fit)) # normality test OK-ish
# residuals are somewhat skewed, likely due to c stock distribution

# Residual linearity & heteroskedasticity 
hydro_resfitted <- ggplot(hydro_fit, aes(.fitted, .resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  geom_smooth(size = 0.5, alpha = 0.25) +
  labs(x = "Fitted values", y = "Residuals")
hydro_resfitted
ncvTest(hydro_fit) # heteroskedasticity test OK

hydro_fit %>% 
  broom::augment() %>% 
  mutate(site = as.integer(site)) %>% 
  pivot_longer(duration:site) %>%
  ggplot(aes(value, .resid)) +
  facet_wrap(~ name, scales = "free_x", labeller = label_parsed) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(data = ~ filter(.x, name == "site"), aes(group = value)) +
  geom_point() +
  geom_smooth(data = ~ filter(.x, name != "site"), span = 1.25)

# Residual independence
hydro_resmap <- wetland_pts %>%
  group_by(site) %>%
  mutate(xmin = sf::st_bbox(geom)[1], ymin = sf::st_bbox(geom)[2]) %>%
  bind_cols(sf::st_coordinates(.)) %>%
  sf::st_drop_geometry() %>%
  filter(wetland %in% soc_hydro$wetland) %>%
  mutate(
    `Resid.` = resid(hydro_fit),
    east = X - xmin, north = Y - ymin
  ) %>%
  ggplot(aes(east, north)) +
  facet_wrap(~ site, ncol = 2) +
  geom_point(aes(fill = `Resid.`, size = `Resid.`), shape = 21) +
  khroma::scale_fill_BuRd(breaks = seq(-10, 10, 4), limits = c(-10, 10)) +
  scale_size_area(max_size = 6, breaks = seq(-10, 10, 4), limits = c(-10, 10)) +
  guides(
    fill = guide_legend(reverse = TRUE), 
    size = guide_legend(reverse = TRUE)
  ) +
  labs(x = "East (m)", y = "North (m)") +
  coord_fixed()
hydro_resmap

hydro_diags <- hydro_qq + hydro_resfitted + hydro_resmap + 
  plot_spacer() + 
  plot_layout(design = "11112222\n33333334", heights = c(1.5, 2)) &
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme_bw() &
  theme(
    panel.grid.minor = element_blank(), 
    plot.tag = element_text(size = 11)
  )
hydro_diags

ggsave(
  "reports/plots/hydro-model-diag-plot-supp.pdf", hydro_diags, 
  device = cairo_pdf, width = 5.5, height = 4.25, units = "in"
)

### 5.1.2. Summary -------------------------------------------------------------

# Model summary & beta coefficients
stargazer::stargazer(
  lm(mutate(hydro_fit$model, across(-site, scale))), 
  type = "text", report = "vcsp", single.row = TRUE, 
  out = "reports/tables/hydro-fit.txt"
)
# Variance inflation factors
vif(hydro_fit)
# Other predictor metrics
yhat::calc.yhat(hydro_fit)[[1]]
# r^2
map_dbl(select(soc_hydro, 6:11), ~ cor(.x, soc_hydro$c_stock))^2
cor(as.integer(soc_hydro$site), soc_hydro$c_stock)^2
# r^2 separately for each site
map_depth(
  group_split(select(soc_hydro, -1), site, .keep = FALSE), 1, 
  ~ cor(.x, .x$c_stock)^2
)
map_depth(
  group_split(select(soc_hydro, -1), site, .keep = FALSE), 1, 
  ~ lm(.x$c_stock ~ .x)
)

### 5.1.3. Variance partitioning -----------------------------------------------

hydro_site_part <- rdacca.hp(
  soc_hydro$c_stock, 
  select(soc_hydro, duration, max, recession, site), 
  method = "RDA", type = "R2", var.part = TRUE
)
stargazer::stargazer(
  hydro_site_part$Hier.part, type = "text", 
  out = "reports/tables/hydro-hier-part.txt"
)
stargazer::stargazer(
  hydro_site_part$Var.part, type = "text", 
  out = "reports/tables/hydro-var-part.txt"
)

## 5.2. Topography model -------------------------------------------------------

soc_topo <- soc_data %>%
  left_join(topo_metrics) %>%
  mutate(
    catchment = -1/catchment,
    site = as.factor(site)
  ) %>%
  as.data.frame()

# (Multi)collinearity among metrics
# Correlations
cor(select(soc_topo, all_of(topo_vars)))
map_df(select(soc_topo, 6:11), ~ broom::tidy(t.test(.x ~ soc_topo$site)))
# Variance inflation factors
vif(lm(
  c_stock ~ area + depth + profile + catchment + rtp + shape, 
  data = soc_topo, na.action = na.fail
))
# PCA biplot
# biplot(prcomp(soc_topo[, 6:11], scale. = TRUE))
# Low multicollinearity, but catchment & depth are correlated
topo_fit <- lm(
  c_stock ~ area + depth + profile + rtp + shape + site, 
  data = soc_topo, na.action = na.fail
)

### 5.2.1. Diagnostics ---------------------------------------------------------

# Residual normality
topo_qq <- ggplot(topo_fit) +
  geom_abline(linetype = "dashed") +
  stat_qq(aes(sample = .stdresid)) +
  labs(x = "Theoretical quantiles", y = "Std. residuals")
topo_qq
shapiro.test(rstandard(topo_fit)) # normality test OK

# Residual linearity & heteroskedasticity 
topo_resfitted <- ggplot(topo_fit, aes(.fitted, .resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  geom_smooth(size = 0.5, alpha = 0.25) +
  labs(x = "Fitted values", y = "Residuals")
topo_resfitted
ncvTest(topo_fit) # heteroskedasticity test OK
var.test(.resid ~ site, data = broom::augment(topo_fit)) # equal variance OK

# Residual independence
topo_resmap <- wetland_pts %>%
  group_by(site) %>%
  mutate(xmin = sf::st_bbox(geom)[1], ymin = sf::st_bbox(geom)[2]) %>%
  bind_cols(sf::st_coordinates(.)) %>%
  sf::st_drop_geometry() %>%
  mutate(
    `Resid.` = resid(topo_fit),
    east = X - xmin, north = Y - ymin
  ) %>%
  ggplot(aes(east, north)) +
  facet_wrap(~ site) +
  geom_point(aes(fill = `Resid.`, size = `Resid.`), shape = 21) +
  khroma::scale_fill_BuRd(
    na.value = "#FFFFFF", breaks = seq(-15, 15, 6), limits = c(-15, 15)
  ) +
  scale_size_area(max_size = 6, breaks = seq(-15, 15, 6), limits = c(-15, 15)) +
  guides(
    fill = guide_legend(reverse = TRUE), 
    size = guide_legend(reverse = TRUE)
  ) +
  labs(x = "East (m)", y = "North (m)") +
  coord_fixed()
topo_resmap

topo_diags <- topo_qq + topo_resfitted + topo_resmap + 
  plot_spacer() + 
  plot_layout(design = "11112222\n33333334", heights = c(1.5, 2)) &
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme_bw() &
  theme(
    panel.grid.minor = element_blank(), 
    plot.tag = element_text(size = 11)
  )
topo_diags

ggsave(
  "reports/plots/topo-model-diag-plot-supp.pdf", topo_diags, 
  device = cairo_pdf, width = 5.5, height = 4.25, units = "in"
)

# Spatial dependence tests (not currently used)
# nb <- spdep::dnearneigh(wetland_pts, d1 = 0, d2 = 400)
# wt <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
# spdep::moran.plot(resid(topo_fit), wt)
# spdep::lm.morantest(topo_fit, wt)
# plot(gstat::variogram(
#   formula(topo_fit), locations = ~ X + Y, data = soc_topo, 
#   cutoff = 1500, width = 150
# ))
# 
# plot(nlme::Variogram(
#   nlme::gls(formula(topo_fit), data = soc_topo),
#   form = ~ X + Y, maxDist = 1500, collapse = "none"
# ))


### 5.2.2. Summary -------------------------------------------------------------

# Model summary & beta coefficients
stargazer::stargazer(
  lm(mutate(topo_fit$model, across(-site, scale))), 
  type = "text", report = "vcsp", single.row = TRUE, 
  out = "reports/tables/topo-fit.txt"
)
# Variance inflation factors
vif(topo_fit)
# Other predictor metrics
yhat::calc.yhat(topo_fit)[[1]]
# r^2
map_dbl(select(soc_topo, 6:11), ~ cor(.x, soc_topo$c_stock))^2
cor(as.integer(soc_topo$site), soc_topo$c_stock)^2
# r^2 separately for each site
map_depth(
  group_split(select(soc_topo, -1), site, .keep = FALSE), 1, 
  ~ cor(.x, .x$c_stock)^2
)

### 5.2.3. Variance partitioning -----------------------------------------------

topo_site_part <- rdacca.hp(
  soc_topo$c_stock, 
  select(soc_topo, area, depth, profile, rtp, shape, site), 
  method = "RDA", type = "R2", var.part = TRUE
)
stargazer::stargazer(
  topo_site_part$Hier.part, type = "text", 
  out = "reports/tables/topo-hier-part.txt"
)
stargazer::stargazer(
  topo_site_part$Var.part, type = "text", 
  out = "reports/tables/topo-var-part.txt"
)
apply(topo_site_part$Hier.part, 2, sum)


# 6. Commonality analysis ------------------------------------------------------

soc_metrics <- soc_data %>%
  inner_join(hydro_metrics) %>%
  left_join(topo_metrics) %>%
  mutate(catchment = -1/catchment)

# Choose (3?) most important metrics from each set
# - Hydro
rdacca.hp(
  soc_hydro$c_stock, 
  select(soc_hydro, duration, recession, max), 
  method = "RDA", type = "R2", var.part = FALSE
)
comm_hydro_vars <- c("recession", "duration", "max")
# - Topo
rdacca.hp(
  soc_topo$c_stock, 
  select(soc_topo, area, depth, profile, rtp, shape), 
  method = "RDA", type = "R2", var.part = FALSE
)
comm_topo_vars <- c("profile", "depth", "shape")

## 6.1. Variance components ----------------------------------------------------

# For each group
group_part <- rdacca.hp(
  soc_metrics$c_stock, 
  list(
    Hydro = select(soc_metrics, recession, max, duration), 
    Topo = select(soc_metrics, profile, depth, shape),
    Site = select(soc_metrics, site)
  ), 
  method = "RDA", type = "R2", var.part = TRUE
)
stargazer::stargazer(
  group_part$Hier.part, type = "text", 
  out = "reports/tables/group-hier-part.txt"
)
stargazer::stargazer(
  group_part$Var.part, type = "text", 
  out = "reports/tables/group-var-part.txt"
)
# small negative values reflect sampling error (Reichwein & Thompson 2006)
group_part$Hier.part

# Peres-Neto et al. (2006) say that using adj R^2 is necessary to avoid bias--
# but focused on multivariate analysis
# - regardless, it's useful to see adjusted individual contributions
rdacca.hp(
  soc_metrics$c_stock, 
  list(
    Hydro = select(soc_metrics, recession, max, duration), 
    Topo = select(soc_metrics, profile, depth, shape),
    Site = select(soc_metrics, site)
  ), 
  method = "RDA", type = "adjR2", var.part = FALSE
)

# For each metric (and site)
fit_part <- lm(
  c_stock ~ recession + duration + max + profile + depth + shape + site,
  data = soc_metrics, na.action = na.fail
)
summary(fit_part)
metric_part <- rdacca.hp(
  soc_metrics$c_stock, 
  select(soc_metrics, all_of(comm_hydro_vars), all_of(comm_topo_vars), site), 
  method = "RDA", type = "R2", var.part = TRUE
)
metric_part$Hier.part
metric_part$Var.part %>% 
  as_tibble(rownames = "vars") %>% 
  filter(str_detect(vars, "Total", negate = TRUE)) %>%
  mutate(n_vars = str_count(vars, ",") + 1, .after = 1) %>%
  arrange(desc(Fractions))

## 6.2. Components plot --------------------------------------------------------

# Prepare commonality analysis results
group_part_data <- group_part %>%
  pluck("Var.part") %>%
  as_tibble(rownames = "part") %>%
  slice(-n()) %>%
  mutate(
    part = part %>% 
      str_trim() %>% 
      str_remove_all("Unique to |Common to | and")
  ) %>% 
  left_join(as_tibble(group_part$Hier.part, rownames = "part")) %>%
  mutate(part = part %>% str_wrap(14) %>% fct_inorder() %>% fct_rev())

# Variance components bar plot

unique_ind_panel <- group_part_data %>%
  drop_na(Unique) %>%
  ggplot(aes(x = part)) +
  geom_col(
    aes(y = `I.perc(%)`), fill = NA, color = "black", width = 0.48, 
    linewidth = 0.3, linetype = "13"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.line.y = element_line(linewidth = 0.25), 
    axis.text.x = element_blank(), 
    axis.ticks = element_blank(), 
    plot.margin = margin(5, 5, 2, 5)
  )

unique_panel <- group_part_data %>%
  drop_na(Unique) %>%
  ggplot(aes(x = part)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.line.y = element_line(linewidth = 0.25), 
    axis.text.x = element_blank(), 
    axis.ticks = element_blank(), 
    plot.margin = margin(5, 5, 2, 5)
  )

common_panel <- group_part_data %>%
  filter(is.na(Unique)) %>%
  mutate(
    flag = if_else(between(`  % Total`, -10, 0), 1, NA_real_),
    Fractions = if_else(between(`  % Total`, -10, 0), 0, Fractions),
    `  % Total` = if_else(between(`  % Total`, -10, 0), 0, `  % Total`)
  ) %>%
  ggplot(aes(x = part)) +
  geom_text(aes(y = flag, label = "*"), size = 5, na.rm = TRUE) +
  labs(x = NULL, y = "Prop. total explained variance (%)") +
  theme_bw() +
  theme(
    axis.line = element_line(linewidth = 0.25), 
    axis.ticks.y = element_blank(),
    plot.margin = margin(2, 5, 5, 5)
  )

# Version with independent & commonality components
contrib_plot <- (unique_ind_panel / common_panel) +
  plot_layout(heights = c(0.44, 0.56)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  geom_col(
    aes(y = `  % Total`), color = "black", fill = "#f7bea4",
    width = 0.48, linewidth = 0.3
  ) &
  scale_y_continuous(expand = expansion(c(0, 0.05))) &
  coord_flip() &
  theme(
    legend.position = "none", 
    legend.title = element_blank(), 
    panel.border = element_blank(),
    panel.grid = element_blank(), 
    plot.tag = element_text(size = 11)
  )
contrib_plot
ggsave(
  "reports/plots/commonality-imp-plot.pdf", contrib_plot, 
  device = cairo_pdf, width = 4.00, height = 3.50, units = "in"
)

# Version with just commonality components
contrib_plot_noind <- (unique_panel / common_panel) +
  plot_layout(heights = c(0.44, 0.56)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  geom_col(
    aes(y = `  % Total`), color = "black", fill = "#f7bea4",
    width = 0.48, linewidth = 0.3
  ) &
  coord_flip(ylim = c(0, 42.93)) &
  scale_y_continuous(expand = expansion(c(0, 0.05))) &
  theme(
    legend.position = "none", 
    legend.title = element_blank(), 
    panel.border = element_blank(),
    panel.grid = element_blank(), 
    plot.tag = element_text(size = 11)
  )
contrib_plot_noind
ggsave(
  "reports/plots/commonality-plot.pdf", contrib_plot_noind, 
  device = cairo_pdf, width = 4.00, height = 3.50, units = "in"
)

# Venn diagram (not currently used)
# rad_seq <- seq(0, 2 * pi, length.out = 200)
# offset <- 0.58
# circles <- tibble(
#   x = c(cos(rad_seq) - offset, cos(rad_seq) + offset, cos(rad_seq)), 
#   y = c(sin(rad_seq) + offset, sin(rad_seq) + offset, sin(rad_seq) - offset),
#   group = rep(c("Hydro", "Topo", "Site"), each = length(rad_seq))
# ) 
# 
# part_locs <- tibble(
#   x = c(-0.78, 0.78, 0, 0, -0.48, 0.48, 0),
#   y = c(0.70, 0.70, -0.75, 0.75, -0.08, -0.08, 0.15),
#   part = c(1, 7, 3, 6, 2, 4, 5),
#   value = group_part_data$part %>% str_wrap(8) %>% fct_inorder() %>% fct_rev()
# )
# parts_plot <- circles %>%
#   ggplot(aes(x, y)) +
#   geom_polygon(
#     aes(group = group), alpha = 0.5, color = NA, fill = "#ee7733",
#     show.legend = FALSE
#   ) +
#   geom_polygon(aes(group = group), fill = NA, color = "black", size = 0.3) +
#   geom_text(
#     data = part_locs, aes(label = value), size = 3, lineheight = 1.0
#   ) +
#   coord_fixed(expand = FALSE, xlim = c(-1.62, 1.62), ylim = c(-1.62, 1.62)) +
#   labs(x = NULL, y = NULL) +
#   theme_void()
# parts_plot
# 
# common_plot <- (parts_plot | contrib_plot) +
#   plot_annotation(tag_levels = "a", tag_suffix = ")") & 
#   theme(plot.tag.position = c(0.03, 1))
# common_plot
# ggsave(
#   "reports/commonality-plot-venn.pdf", common_plot, device = cairo_pdf,
#   width = 6.50, height = 3.25, units = "in"
# )


# References -------------------------------------------------------------------



