
library("tidyverse")


# 1. Read data -----------------------------------------------------------------

emergent <- read_csv("output/soil/soc-stocks-est-emergent.csv") %>%
  mutate(type = if_else(wetland %in% c("DF", "FN", "JU"), "natu", "rest"))

forested <- read_csv("output/core-data-ak.csv")
forested_mean <- forested %>%
  group_by(wetland, site) %>%
  summarize(c_stock_0_100 = mean(C_stock_100cm), .groups = "drop")

# 2. Emergent wetland variability ----------------------------------------------

# Natural vs. restored
emergent %>%
  group_by(type) %>%
  summarize(across(
    c_stock_0_100_est, list(mean = mean, min = min, max = max, sd = sd), 
    .names = "{.fn}"
  ))

# Each zone in natural vs. restored
emergent %>%
  group_by(type, zone) %>%
  summarize(across(
    c_stock_0_100_est, list(mean = mean, min = min, max = max, sd = sd), 
    .names = "{.fn}"
  )) %>%
  mutate(cv = sd / mean * 100)

# Among-wetland variability 

# Mean of CVs for each zone
emergent %>% 
  group_by(type, zone) %>% 
  summarize(cv = sd(c_stock_0_100_est) / mean(c_stock_0_100_est) * 100) %>% 
  summarize(mean_cv = mean(cv))

# Among-zone (within-wetland) variability 

# Mean of CVs for each wetland
emergent %>% 
  group_by(type, wetland) %>% 
  summarize(cv = sd(c_stock_0_100_est) / mean(c_stock_0_100_est) * 100) %>% 
  summarize(mean_cv = mean(cv))

# Statistically distinguishing between nat & rest
emergent %>%
  group_by(type, zone) %>%
  summarize(mean = mean(c_stock_0_100_est), .groups = "drop") %>% 
  t.test(mean ~ type, data = .)


# 3. Forested wetland variability ----------------------------------------------

# Mean difference between duplicate cores
forested %>%
  group_by(wetland) %>%
  summarize(diff = abs(diff(C_stock_100cm))) %>%
  summarize(mean_diff = mean(diff))

# Variability among all wetlands
forested_mean %>%
  summarize(cv = sd(c_stock_0_100) / mean(c_stock_0_100) * 100)

# Variability for each site
forested_mean %>%
  group_by(site) %>%
  summarize(cv = sd(c_stock_0_100) / mean(c_stock_0_100) * 100)


# 4. Range plot for all types/zones --------------------------------------------

# Make range data
emergent_range <- emergent %>%
  group_by(type, zone) %>%
  summarize(
    min = min(c_stock_0_100_est), 
    max = max(c_stock_0_100_est), 
    .groups = "drop"
  ) %>%
  # Hack to get everything to show up in the right place
  mutate(
    typepos = as.integer(factor(type)) + 1,
    width = 0.05,
    zonepos = as.integer(factor(zone)) * (width * 3) - (2 * width * 3),
    xpos = typepos + zonepos
  )

# Build plot
soc_varplot <- forested_mean %>%
  mutate(type = "forested", zone = "None (open water)") %>%
  ggplot(aes(type)) +
  geom_rect(
    data = emergent_range, 
    aes(
      xmin = xpos - width, xmax = xpos + width, ymin = min, ymax = max, 
      fill = zone
    )
  ) +
  geom_violin(
    # aes(y = c_stock_0_100), width = 0.3, fill = "#88ccee", color = NA
    aes(y = c_stock_0_100), width = 0.3, fill = "grey80", color = NA
  ) +
  geom_dotplot(
    aes(y = c_stock_0_100), binwidth = 1, binaxis = "y", stackdir = "center", 
    fill = "white", dotsize = 1.15, stackratio = 1, stroke = 1.4
  ) +
  geom_point(
    data = emergent, aes(y = c_stock_0_100_est, group = zone), 
    position = position_dodge(width = 0.45), 
    shape = 21, stroke = 0.75, fill = "white", size = 1.75
  ) +
  scale_x_discrete(
    limits = c("forested", "natu", "rest"), 
    labels = c("Forested\nnatural", "Emergent\nnatural*", "Emergent\nrestored*")
  ) +
  labs(
    x = NULL, 
    y = bquote("1-m SOC stock"~(kg~C~m^-2)),
    caption = "*Estimated"
  ) +
  scale_fill_manual(
    values = c("#88ccee", "#44aa99", "#ddcc77"), 
    labels = c("Open water", "Herbaceous", "Woody")
  ) +
  guides(
    fill = guide_legend(title = "Vegetation zone")
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 11), 
    legend.margin = margin(), 
    legend.position = c(0.82, 0.77),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10.5),
    panel.grid = element_blank()
  )
soc_varplot

ggsave(
  "reports/plots/soc-varplot.pdf", soc_varplot, device = cairo_pdf,
  width = 4.25, height = 3.25, units = "in"
)
