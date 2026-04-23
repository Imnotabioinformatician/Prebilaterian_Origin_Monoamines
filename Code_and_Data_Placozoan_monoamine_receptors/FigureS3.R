# Load necessary libraries
library(tidyverse)
library(cowplot)
library(png)
library(patchwork)
library(drc)
library(usethis)
library(gitcreds)
library(devtools)
library(knitr)
library(vroom)
library(credentials)
library(purrr)
library(ggsci)
library(igraph)
library(viridis)
library(plotly)
library(igraph)
library(ggplot2)
# Set working directory
setwd("~/Documents/GitHub/MonoamineDeorphanisation/Trichoplax_GPCRs")
Sys.setlocale("LC_MESSAGES", "en_GB.UTF-8")
Sys.setenv(LANG = "en_US.UTF-8")
set.seed(666)
# Read data
IndividualHeatplot <- vroom("data/Placozoan_GPCRAlldataIndividualCompounds.csv")

# Clean column names: Replace spaces and special characters
colnames(IndividualHeatplot) <- make.names(colnames(IndividualHeatplot))

# Pivot longer using actual column names
IndividualHeatplot <- IndividualHeatplot %>%
  pivot_longer(cols = -c(Transfection, Replicate, Receptor, Negative), 
               names_to = "Compounds", values_to = "luminescence") %>%
  filter(!is.na(luminescence)) %>%
  group_by(Replicate, Receptor) %>%
  mutate(norm_luminescence = 100 * (luminescence - first(luminescence)) / (max(luminescence) - first(luminescence)))

# Define color scale
marcascolor <- c(0, 50, 100)

# Plot
plot2 <- IndividualHeatplot %>%
  ggplot(aes(x = Compounds, y = Receptor)) +
  geom_tile(aes(fill = norm_luminescence)) +
  scale_fill_gradient2(low = "black", mid = "gray30", high = "#CCFF00", 
                       midpoint = 20, breaks = marcascolor, na.value = "black") +
  theme_half_open() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Rotate more and increase size
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank()) #+
  #guides(fill = "none")# Increase y-axis text size



plot2


#Define plot dimensions (adjust as needed, units can be "in", "cm", "mm")
plot_width <- 14
plot_height <- 10
plot_units <- "in"

# Save as PNG (good for web, presentations; raster format)
# DPI (dots per inch) controls resolution, 300 is good for general use, 600 for high-res print
ggsave(
  filename = "output_figures/Supplementary_Fig_1.png",
  plot = plot2,
  width = plot_width,
  height = plot_height,
  units = plot_units,
  dpi = 300,
  bg = "white" # Specify background color if needed (e.g., for transparency issues)
)
message("Saved plot as PNG.")



# Save as PDF (good for publications, scaling without quality loss; vector format)
# Use device = cairo_pdf for potentially better handling of fonts and special characters
ggsave(
  filename = "output_figures/Supplementary_Fig_1.pdf",
  plot = plot2,
  width = plot_width,
  height = plot_height,
  units = plot_units,
  device = cairo_pdf # Requires the Cairo package: install.packages("Cairo")
  # device = "pdf" # Alternative standard PDF device
)
message("Saved plot as PDF.")

# Save as SVG (Scalable Vector Graphics - good for web, editing in vector software)
ggsave(
  filename = "output_figures/Supplementary_Fig_1.svg",
  plot = plot2,
  width = plot_width,
  height = plot_height,
  units = plot_units
  # device = "svg" # device argument usually not needed for svg
)
message("Saved plot as SVG.")





