# Clear workspace
rm(list = ls(all.names = TRUE))
gc()

# Load required packages
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
# Set working directory
setwd("~/Documents/GitHub/MonoamineDeorphanisation/Trichoplax_GPCRs")

# Read dose-response data
DoseResCurv <- vroom::vroom("data/DeorphanisationDataPlacozoans.csv")

# Convert to long format and remove NA
AllGPCRtoplot <- DoseResCurv %>%
  pivot_longer(starts_with(c("Nega", "1e", "0", "1")),
               names_to = "concentration", values_to = "luminescence") %>%
  mutate(concentration = as.double(concentration)) %>%
  filter(!is.na(luminescence))

# Normalization function
normalize_to_ctr <- function(x) {
  100 * (x - x[1]) / (max(x) - x[1])
}

# Normalize luminescence values
AllGPCRtoplotTR1 <- AllGPCRtoplot %>%
  group_by(Replicate, Receptor, Compound) %>%
  mutate(norm_luminescence = normalize_to_ctr(luminescence))

# Function to add gridlines
grids <- function(axis = c("xy", "x", "y"), color = "slategray1", size = NULL, linetype = NULL) {
  axis <- match.arg(axis)
  grid.major <- element_line(color = color, size = size, linetype = linetype)
  grid.minor <- element_line(color = color, size = 0.25, linetype = linetype)
  switch(axis,
         xy = theme(panel.grid.major = grid.major, panel.grid.minor = grid.minor),
         x = theme(panel.grid.major.x = grid.major, panel.grid.minor.x = grid.minor),
         y = theme(panel.grid.major.y = grid.major, panel.grid.minor.y = grid.minor))
}

# Plot dose-response curves
AllGPCRtoplotTR1 %>%
  ggplot(aes(x = concentration, 
             y = norm_luminescence, 
             colour = Compound, 
             shape = Compound,  # Add this line
             group = Compound)) +
  scale_color_aaas() +
  scale_shape_manual(values = 1:length(unique(AllGPCRtoplotTR1$Compound))) + # Map shapes manually
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE, size = 1.2) +
  scale_x_log10(
    breaks = c(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2), 
    limits = c(1e-13, 1e-2),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_half_open() +
  theme(
    axis.text = element_text(size = 15), 
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16),
    axis.title = element_text(size = 18), 
    axis.title.x = element_text(margin = margin(t = 3.5), face = "bold"),
    axis.title.y = element_text(margin = margin(t = 2), face = "bold", size = 18),
    panel.background = element_blank(),
    legend.position = "left"
  ) +
  labs(x = "Concentration (M)", y = "% Luminescence") +  # Custom axis labels
  stat_summary(fun = mean, geom = "point", size = 3) + # Remove `shape` here since it's in aes()
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.3) +
  facet_wrap(vars(Receptor)) + 
  grids(linetype = "longdash")




# EC50 calculation --------------------------------------------------------

ListREceptors <- unique(AllGPCRtoplotTR1$Receptor)
ListPeptides <- unique(AllGPCRtoplotTR1$Compound)
df <- data.frame(Receptor = character(), Peptide = character(), EC50 = numeric(), stringsAsFactors = FALSE)

for (i in ListREceptors) {
  RecSubset <- filter(AllGPCRtoplotTR1, Receptor == i)
  for (k in ListPeptides) {
    SubsetRecPep <- filter(RecSubset, Compound == k)
    if (nrow(SubsetRecPep) > 0) {
      write.csv(SubsetRecPep, file = paste0("Subs_", i, "_", k, ".csv"), row.names = FALSE)
      model <- drm(
        norm_luminescence ~ concentration, 
        data = SubsetRecPep, 
        fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"))
      )
      EC50 <- model$coefficients[4]
      df <- rbind(df, data.frame(Receptor = i, Peptide = k, EC50 = EC50))
    }
  }
}

write.csv(df, file = "supplements/EC50table.csv", row.names = FALSE)
