
library("readxl")
library("lubridate")
library("tidyverse")


# 1. Read data -----------------------------------------------------------------

dir_in <- "data/water-level"
dir_out <- "output/water-level"

# Read water level data
water_level <- read_csv(file.path(dir_in, "water-level-daily.csv"))
gapfill <- read_csv(file.path(dir_in, "water-level-gapfilled.csv"))

# Subset of wetlands for analysis
wetlands <- c(
  "BB", "DB", "DK", "GN", "JA", "JB", "JC", "NB", "ND", "QB", "TA", "TB"
)


# 2. Apply gap-filling ---------------------------------------------------------

# Function for flagging times when water level bottomed out
is_well_dry <- function(x, min_tol = 0.03, delta_tol = 0.005) {
  
  minx <- if (all(is.na(x))) NA else min(x, na.rm = TRUE)
  dx <- abs(x - lag(x))
  
  is_min <- near(x, minx, tol = min_tol)
  is_flat <- near(dx, 0, tol = delta_tol)
  
  if_else(is_min & is_flat, 1L, 0L)
}

water_level_gf <- water_level %>%
  filter(name %in% paste(wetlands, "Wetland Well Shallow")) %>%
  left_join(gapfill) %>%
  mutate(
    wetland = word(name),
    year = year(date),
    wateryear = year(date %m+% months(3))
  ) %>%
  # Make sure chronological order before time series analysis
  group_by(wetland) %>%
  arrange(date, .by_group = TRUE) %>%
  # Add derivatives
  mutate(
    # Inundated?
    inund = if_else(water_level >= 0, 1, 0),
    inund_gf = coalesce(inund, inund_pred),
    water_level_gf = coalesce(water_level_gf, water_level),
    # Daily water level change
    delta_1d = water_level_gf - lag(water_level_gf),
    # Direction of water level change
    dir_1d = sign(delta_1d)
  ) %>%
  # Add year grouping since different years have different min water level
  group_by(year, .add = TRUE) %>%
  # Flag periods when well is dry (criteria = min + not changing)
  mutate(well_dry = is_well_dry(water_level_gf)) %>%
  ungroup() %>%
  select(-name, -inund_pred, -gapfill_well)

write_csv(water_level_gf, file.path(dir_out, "water-level-daily-gapfilled.csv"))

# Sanity check
water_level_gf %>%
  filter(year == 2019) %>%
  ggplot(aes(date, water_level_gf)) +
  facet_wrap(~ wetland) +
  geom_line()


# 3. Calculate yearly metrics --------------------------------------------------

water_level_gf <- read_csv(
  file.path(dir_out, "water-level-daily-gapfilled.csv")
)

hydro_summary <- water_level_gf %>%
  filter(year == 2019) %>%
  group_by(wetland) %>%
  # Need robust metrics since water level bottoms out
  summarize(
    # Number of gap-filled days
    n_gf = sum(is.na(water_level)),
    max = max(water_level_gf),
    median = median(water_level_gf),
    iqr = IQR(water_level_gf),
    mad = mad(water_level_gf),
    # Number of days inundated
    days_wet = sum(inund_gf),
    # First dry date
    date_dry = yday(first(date[which(inund_gf == 0)])),
    # First date of last yearly inundation
    date_wet = yday(last(date[which(inund_gf == 0)])) + 1,
    # Number of times wetland dried out
    n_dry = sum(abs(diff(inund_gf)))/2,
    # Robust coefficient of variation
    rcv = mad/median * 100,
    # Median absolute water level change
    med_delta = median(abs(delta_1d)),
    .groups = "drop"
  ) %>%
  arrange(wetland)

write_csv(hydro_summary, file.path(dir_out, "hydro-summary-2018-2019.csv"))


# 4. Calculate recession rates -------------------------------------------------

water_level_gf <- read_csv(
  file.path(dir_out, "water-level-daily-gapfilled.csv")
)

# Apply criteria for 'recession days'
rec_days <- water_level_gf %>%
  filter(year == 2018) %>%
  group_by(wetland) %>%
  mutate(delta_lag = lag(delta_1d)) %>%
  group_by(date) %>%
  filter(
    # All wetlands have (non gap-filled) data
    sum(is.na(water_level)) == 0,
    # Growing season (May-Oct)
    between(month(date), 5, 10),
    # No dry wells
    sum(well_dry) == 0,
    # Recession observed in all wetlands
    all(delta_1d < 0),
    # No (presumed) heavy rain the previous day
    all(delta_lag < 0.01),
    # No very dry wetlands (difference in specific yield)
    all(water_level > -0.25)
  ) %>%
  ungroup()

# How many days?
distinct(rec_days, date)

# How are daily recessions distributed?
rec_days %>%
  ggplot(aes(delta_1d)) +
  facet_wrap(~ wetland, scales = "free") +
  geom_density()

# Yearly rate is median of daily values (cm d-1)
rec_rates <- rec_days %>%
  group_by(wetland) %>% 
  summarize(
    rec_rate = median(delta_1d) * 100, 
    .groups = "drop"
  )

write_csv(rec_rates, file.path(dir_out, "recession-rates-2018.csv"))


# 5. Select metrics for analysis -----------------------------------------------

hydro_summary <- read_csv(file.path(dir_out, "hydro-summary-2018-2019.csv"))
rec_rates <- read_csv(file.path(dir_out, "recession-rates.csv"))

# Rename if needed
hydro_vars <- c(
  "median", "max", "iqr", 
  "duration" = "days_wet", 
  "exposure" = "date_dry",
  "recession" = "rec_rate"
)

hydro_metrics <- hydro_summary %>%
  left_join(rec_rates) %>%
  select(wetland, all_of(hydro_vars))
hydro_metrics

write_csv(hydro_metrics, file.path(dir_out, "wetland-hydro-metrics.csv"))


# 6. Climate & hydrograph plot -------------------------------------------------

library("patchwork")

climate <- read_csv(file.path(dir_in, "climate-vars-daily-fntower-2019.csv"))
water_level_gf <- read_csv(
  file.path(dir_out, "water-level-daily-gapfilled.csv")
)

# Climate panel: daily temp. & precip.
ppt_adj <- 20
ppt_mult <- 2
climate_plot <- climate %>%
  filter(year(date) == 2019) %>%
  ggplot(aes(date)) +
  geom_segment(
    aes(
      y = (ppt_mm / ppt_mult) - ppt_adj, xend = date, yend = -ppt_adj, 
        color = "Precip."
    ), 
    lineend = "butt", size = 0.75
  ) +
  geom_line(aes(y = tmean_c, color = "Temp."), size = 0.6) + 
  scale_y_continuous(
    name = expression("Mean temperature (" *degree*"C)"),
    sec.axis = sec_axis(
      trans = (~ (. + ppt_adj) * ppt_mult), name = "Total precip. (mm)"
    ), 
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_color_manual(
    name = NULL, values = c("Precip." = "#0077bb", "Temp." = "#cc3311")
  ) +
  labs(x = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    legend.position = c(0.115, 0.885)
  )
climate_plot

# Hydrograph panel: daily water level
water_level_plot <- water_level_gf %>%
  left_join(select(soc_data, wetland, site)) %>%
  filter(year == 2019) %>%
  ggplot(aes(date)) +
  geom_hline(aes(yintercept = 0), size = 0.4, alpha = 0.5) +
  # For means & individual lines:
  geom_line(
    aes(y = water_level_gf, group = wetland), size = 0.25, alpha = 0.2
  ) +
  stat_summary(
    aes(y = water_level_gf, color = site), geom = "line", fun = mean, size = 0.6
  ) +
  scale_color_manual(values = c("#44aa99", "#332288"), name = "Site mean") +
  # For mean & std. dev. ribbons:
  # stat_summary(
  #   aes(y = water_level_gf, fill = site), geom = "ribbon",
  #   fun.data = mean_sdl, fun.args = list(mult = 1), alpha = 0.4
  # ) +
  # scale_color_manual(values = c("#44aa99", "#332288"), name = "Site mean") +
  # scale_fill_manual(values = c("#44aa99", "#332288"), name = "Site mean") +
  labs(x = NULL, y = "Wetland water level (m)") +
  theme_bw() +
  theme(
    legend.position = c(0.925, 0.875), 
    legend.title = element_text(size = 10.5), 
    plot.margin = margin(t = 0)
  )
water_level_plot

# Stack climate & hydrograph panels
climate_water_plot <- climate_plot / water_level_plot + 
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  scale_x_date(date_labels = "%b", expand = expansion(mult = 0.025)) &
  theme(
    axis.line = element_line(color = "black", size = 0.2),
    axis.title = element_text(size = 10.5),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.03, "npc"),
    legend.spacing.x = unit(0.006, "npc"),
    panel.border = element_blank(),
    panel.grid = element_blank(), 
    plot.tag = element_text(size = 11)
  )
climate_water_plot

ggsave(
  "reports/plots/climate-waterlevel-plot.pdf", climate_water_plot, 
  device = cairo_pdf, width = 5.75, height = 5.00, units = "in"
)

