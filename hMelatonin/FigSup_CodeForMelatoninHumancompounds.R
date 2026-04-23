rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# load packages -----------------------------------------------------------
#install.packages("pacman")
library(pacman)
library(patchwork)
library(drc)
library(knitr)
library(vroom)
library(ggsci)
library(plotly)
library(readr)
library(dplyr)
library(tidyr)
library(gt)
library(webshot2)
library(cowplot)
library(scales)
#check working directory- should be the project dir, all directories will be defined relative to this dir
getwd()
Sys.setlocale("LC_MESSAGES", 'en_GB.UTF-8')
Sys.setenv(LANG = "en_US.UTF-8")
set.seed(666)
# read data ---------------------------------------------------------------
# Load data with headers preserved
curve <- read.csv("Data/HumanMelatoninReceptors.csv", check.names = FALSE)

# View original headers
colnames(curve)
# Standardize headers to consistent scientific notation
colnames(curve) <- sapply(colnames(curve), function(x) {
  if (grepl("^[0-9eE.-]+$", x)) {
    format(as.numeric(x), scientific = TRUE)
  } else {
    x
  }
})

# Pivoting the data
AllGPCRtoplot <- curve %>%
  pivot_longer(starts_with(c("0", "1e", "Neg", "1")), 
               names_to = "concentration", values_to = "luminescence")

#convert conc valus to double
AllGPCRtoplot$concentration <- as.double(AllGPCRtoplot$concentration)

# Delete the N.A ----------------------------------------------------------
AllGPCRtoplot <- AllGPCRtoplot[!is.na(AllGPCRtoplot$luminescence), ]                 # Omit NA by column via is.na

# Functions ---------------------------------------------------------------
#FUNCTIONS to normalise to reference (zero ligand control)
normalize_to_ctr <- function(x) {
  return (100*(x - x[1]) / (max(x) - x[1]))
}


#Function to add grids to the curves 
grids <- function(axis = c("xy", "x", "y"), color = "slategray1", size = NULL, linetype = NULL)
{
  axis <- match.arg(axis)
  grid.major <- element_line(color = color, size = 0.2,
                             linetype = linetype)
  grid.minor <- element_line(color = color, size = 0.2,
                             linetype = linetype)
  
  switch(axis,
         xy = theme(panel.grid.major = grid.major, panel.grid.minor = grid.minor),
         x = theme(panel.grid.major.x = grid.major, panel.grid.minor.x = grid.minor),
         y = theme(panel.grid.major.y = grid.major, panel.grid.minor.y = grid.minor)
  )
}

#This is normalization with all the different receptors
AllGPCRtoplot <- AllGPCRtoplot %>%
  group_by(Replicate, Receptor, Peptide)%>%
  mutate('norm_luminescence'=normalize_to_ctr(luminescence))

# Plot with boxplot -------------------------------------------------------
# NEwplotting -------------------------------------------------------------


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

# Ec50 calculation --------------------------------------------------------

ListREceptors <- unique(AllGPCRtoplot$Receptor)
ListPeptides <- unique(AllGPCRtoplot$Peptide)
df <- ""
for (i in ListREceptors) {
  RecSubset <- filter(AllGPCRtoplot, Receptor == i)
  for (k in ListPeptides) {
    SubsetRecPep <- filter(RecSubset, Peptide == k)
    if (nrow(SubsetRecPep) > 0) {
      write.csv(SubsetRecPep, file = paste("Subs_", i, k), sep = ",", row.names = FALSE);
      model <- drm(
        norm_luminescence~concentration, 
        data = SubsetRecPep, 
        fct=LL.4(names =c("Slope", "Lower Limit", "Upper Limit", "ED50" ))
      )
      EC50 <- model$coefficients[4]
      dfrow <- c(Receptor = i, Peptide = k, EC50 = EC50)
      df <- rbind(df, dfrow)
      write.csv(df, file = "supplements/EC50table2.csv", sep = ",", row.names = FALSE);
    }
  }
}

# Load required libraries
# Load the gt library
library(gt)

# Create a styled table
ec50_data <- read_csv(
  "supplements/EC50table2.csv",
  col_types = cols(
    `Receptor` = col_character(),
    `Peptide` = col_character(),
    `EC50.ED50:(Intercept)` = col_double() # Ensure this is numeric
  )
)

# Step 2: Check the structure of the data to confirm types
str(ec50_data)

# Step 3: Create a styled table
gt_table <- ec50_data %>%
  gt() %>%
  tab_header(title = "EC50 Values") %>%
  cols_label(
    Receptor = "Receptor",
    Peptide = "Compound",
    `EC50.ED50:(Intercept)` = "EC50 (M)"
  ) %>%
  fmt_scientific(columns = c(`EC50.ED50:(Intercept)`), decimals = 2) %>%
  tab_options(
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    table.border.bottom.width = 2,
    table.border.top.width = 2,
    table.border.left.color  = "black",
    heading.title.font.size = 16,
    heading.title.font.weight = "bold"
  )
print(gt_table)

gtsave(gt_table, "EC50_human2.pdf", path = "pictures/")
gtsave(gt_table, "EC50_TableHumanRec22222.pdf", path = "pictures/", expand = 5)




# individualheatplot ------------------------------------------------------
# Load the data
data <- read.csv("Data/hMel_compounds.csv")

# Create an extra column with the luminescence of 10X
IndividualHeatplot <- data %>%
  group_by(Replicate, Receptor) %>%
  mutate("Negative5" = Negative * 5)

# Pivot longer for all columns except Transfection, Replicate, and Receptor
IndividualHeatplot <- IndividualHeatplot %>%
  pivot_longer(
    cols = -c(Transfection, Replicate, Receptor),  # Exclude these columns
    names_to = "Peptide_mix",
    values_to = "luminescence"
  )

IndividualHeatplot <- IndividualHeatplot[!is.na(IndividualHeatplot$luminescence), ]                 # Omit NA by column via is.na



#FUNCTIONS to normalise to 10x negative asreference (zero ligand control)
#To make normalization fairer we have to normalize, either to the highest number obtained
#in the experiment or to the 10X negative. Whatever is higher. 

normalize_to_ctr <- function(x) {
  return (100*(x - x[1]) / (max(x) - x[1]))
}

# Normalize luminescence per replicate and receptor
IndividualHeatplot2 <- IndividualHeatplot %>%
  group_by(Replicate, Receptor) %>%
  mutate(norm_luminescence = normalize_to_ctr(luminescence))

# Define color scale
marcascolor <- c(0, 50, 100)
ordeny2 <- IndividualHeatplot2 %>% pull(Receptor)

#Normalization of the data for heatplots
heatmap_plot <- IndividualHeatplot2 %>%
  filter(!Peptide_mix %in% c("Negative", "Negative5")) %>%
  ggplot(aes(x = Peptide_mix, y = Receptor)) +
  geom_tile(aes(fill = norm_luminescence)) + 
  scale_fill_gradient2(
    low = "black",
    mid = "gray30",
    na.value = "black",
    high = "#CCFF00",
    midpoint = 25,
    breaks = marcascolor
  ) +
  theme_half_open() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate labels
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank()  # Remove x-axis label
  )


# Print to viewer
print(heatmap_plot)

# Save in multiple formats
ggsave("pictures/heatmap_norm.pdf", plot = heatmap_plot, width = 8, height = 6, dpi = 300)
ggsave("pictures/heatmap_norm.jpg", plot = heatmap_plot, width = 8, height = 6, dpi = 300)
ggsave("pictures/heatmap_norm.svg", plot = heatmap_plot, width = 8, height = 6, dpi = 300)






