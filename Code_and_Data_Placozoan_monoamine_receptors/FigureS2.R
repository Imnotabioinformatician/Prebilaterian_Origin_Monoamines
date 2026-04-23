# Clear environment
rm(list = ls(all.names = TRUE))
gc()

# Load required packages
library(tidyverse)
library(cowplot)
library(png)
library(knitr)
library(vroom)
library(ggsci)
library(igraph)
library(viridis)
library(plotly)
library(igraph)
library(ggplot2)
library (patchwork)

# Set working directory
#setwd("~/Documents/GitHub/MonoamineDeorphanisation/Trichoplax_GPCRs")
Sys.setlocale("LC_MESSAGES", "en_GB.UTF-8")
Sys.setenv(LANG = "en_US.UTF-8")
set.seed(666)

# Read data
mixtures_raw <- vroom("data/AllMixReceptorstogether.csv")

# Preprocess and normalize
mixtures_norm <- mixtures_raw %>%
  group_by(Replicate, Receptor) %>%
  mutate(Negative10 = Negative * 6) %>%
  pivot_longer(starts_with(c("Nega", "M")), 
               names_to = "Peptide_mix", values_to = "luminescence") %>%
  filter(!is.na(luminescence)) %>%
  group_by(Replicate, Receptor) %>%
  mutate(norm_luminescence = 100 * (luminescence - first(luminescence)) / 
           (max(luminescence) - first(luminescence))) %>%
  ungroup()

# Split data
first_half <- mixtures_norm %>%
  filter(Receptor %in% c(paste0("Tadh", 111:140), "pcDNA3.1"))

second_half <- mixtures_norm %>%
  filter(Receptor %in% paste0("Tadh", 141:173))

# Common color scale
color_breaks <- c(0, 50, 100)
color_limits <- range(mixtures_norm$norm_luminescence, na.rm = TRUE)
selected_mixes <- c("M1", "M2", "M3", "M4", "M5")

# Heatmap function
plot_heatmap <- function(data, show_legend = TRUE) {
  ggplot(
    data %>% filter(Peptide_mix %in% selected_mixes),
    aes(x = Peptide_mix, y = Receptor)
  ) +
    geom_tile(aes(fill = norm_luminescence)) +
    scale_fill_gradient2(
      low = "black",
      mid = "darkgray",
      high = "#CCFF00",
      midpoint = 35,
      breaks = color_breaks,
      limits = color_limits,
      na.value = "black",
      name = NULL
    ) +
    theme_half_open() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
      axis.text.y = element_text(size = 10),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 9)
    )
}

p_legend <- plot_heatmap(first_half, show_legend = TRUE)
legend <- cowplot::get_legend(p_legend)

p1 <- plot_heatmap(first_half, show_legend = FALSE)
p2 <- plot_heatmap(second_half, show_legend = FALSE)

heatmaps <- cowplot::plot_grid(
  p1, p2,
  ncol = 2,
  align = "h"
)

final_plot <- cowplot::plot_grid(
  heatmaps,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.08)
)

print(final_plot)

# ---------------------------------------------------------
# Save figure as half-page portrait A4
# ---------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

output_base <- "figures/Trichoplax_receptor_mix_heatmap"

ggsave(
  filename = paste0(output_base, ".png"),
  plot = final_plot,
  width = 12,
  height = 14.85,
  units = "cm",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = paste0(output_base, ".pdf"),
  plot = final_plot,
  width = 12,
  height = 14.85,
  units = "cm",
  device = cairo_pdf,
  bg = "white"
)
