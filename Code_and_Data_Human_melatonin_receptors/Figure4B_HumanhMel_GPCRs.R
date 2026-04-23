# Cleaned-up R script split into two figures:
# 1. Dose-response curve plot
# 2. EC50 calculations using drc's
# Clear environment
rm(list = ls(all.names = TRUE))
gc()

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(patchwork, drc, knitr, vroom, ggsci, plotly, readr, dplyr, tidyr,
       gt, webshot2, cowplot, scales, ggplot2)

# Set locale and seed
Sys.setlocale("LC_MESSAGES", 'en_GB.UTF-8')
Sys.setenv(LANG = "en_US.UTF-8")
set.seed(666)

#----------------------------
# 1. Dose-response curve plot
#----------------------------

# Load and prepare data
curve <- read.csv("Data/HumanMelatoninReceptors.csv", check.names = FALSE)
colnames(curve) <- sapply(colnames(curve), function(x) {
  if (grepl("^[0-9eE.-]+$", x)) format(as.numeric(x), scientific = TRUE) else x
})

AllGPCRtoplot <- curve %>%
  pivot_longer(starts_with(c("0", "1e", "Neg", "1")), 
               names_to = "concentration", values_to = "luminescence") %>%
  mutate(concentration = as.double(concentration)) %>%
  filter(!is.na(luminescence))

normalize_to_ctr <- function(x) {
  100 * (x - x[1]) / (max(x) - x[1])
}

grids <- function(axis = c("xy", "x", "y"), color = "slategray1", size = NULL, linetype = NULL) {
  axis <- match.arg(axis)
  grid.major <- element_line(color = color, size = 0.2, linetype = linetype)
  grid.minor <- element_line(color = color, size = 0.2, linetype = linetype)
  switch(axis,
         xy = theme(panel.grid.major = grid.major, panel.grid.minor = grid.minor),
         x = theme(panel.grid.major.x = grid.major, panel.grid.minor.x = grid.minor),
         y = theme(panel.grid.major.y = grid.major, panel.grid.minor.y = grid.minor))
}

AllGPCRtoplot <- AllGPCRtoplot %>%
  group_by(Replicate, Receptor, Peptide) %>%
  mutate(norm_luminescence = normalize_to_ctr(luminescence))

# Plot dose-response curve
AllGPCRtoplot %>% 
  ggplot(aes(x = concentration, 
             y = norm_luminescence, 
             colour = Peptide, 
             shape = Peptide,  # Add this line
             group = Peptide)) +
  scale_color_aaas() +
  scale_shape_manual(values = 1:length(unique(AllGPCRtoplot$Peptide))) + # Map shapes manually
  # Uncomment to add a boxplot, but it may not fit with multiple Compounds in the graph
  # geom_boxplot(aes(group = concentration, colour = Compound, group = Compound), 
  #              outlier.shape = NA, size = 0.15, width = 0.2) +
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
  stat_summary(fun.y = mean, geom = "point", size = 3) + # Remove `shape` here since it's in aes()
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.3) +
  facet_wrap(vars(Receptor)) + 
  grids(linetype = "longdash")


# EC50 table generation
ListREceptors <- unique(AllGPCRtoplot$Receptor)
ListPeptides <- unique(AllGPCRtoplot$Peptide)
df <- data.frame()

for (i in ListREceptors) {
  for (k in ListPeptides) {
    Subset <- filter(AllGPCRtoplot, Receptor == i, Peptide == k)
    if (nrow(Subset) > 0) {
      model <- drm(norm_luminescence ~ concentration, data = Subset,
                   fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
      EC50 <- coef(model)[4]
      df <- rbind(df, data.frame(Receptor = i, Peptide = k, EC50 = EC50))
    }
  }
}

write.csv(df, file = "supplements/EC50table2.csv", row.names = FALSE)

# GT table formatting
ec50_data <- read_csv("supplements/EC50table2.csv", col_types = cols(
  Receptor = col_character(),
  Peptide = col_character(),
  EC50 = col_double()
))

gt_table <- ec50_data %>%
  gt() %>%
  tab_header(title = "EC50 Values") %>%
  cols_label(
    Receptor = "Receptor",
    Peptide = "Compound",
    EC50 = "EC50 (M)"
  ) %>%
  fmt_scientific(columns = EC50, decimals = 2) %>%
  tab_options(
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.title.font.size = 16,
    heading.title.font.weight = "bold"
  )

gtsave(gt_table, "EC50_human2.pdf", path = "pictures/")
gtsave(gt_table, "EC50_TableHumanRec.pdf", path = "pictures/", expand = 5)


