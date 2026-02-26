library(EcoRSI)
library(readxl)
library(ggplot2)

df <- read_excel("EcoRSI_trial_dataset_1.xlsx")

result <- computeRSI(
  data = df,
  umd_label = "UMD",
  re_label = "RE",
  positive_indicators = c("N","P","K","OC","floral_richness","faunal_richness"),
  negative_indicators = c("Mn","Fe","pH","Pb","Cr"),
  index_col = "Index",
  site_col = "site"
)
View(result)

# Ensure correct site order
result$Site <- factor(result$Site, 
                      levels = c("UMD", "PMD", "AMD", "RE"))

# Boxplot with colors
ggplot(result, aes(x = Site, y = RSI, fill = Site)) +
  geom_boxplot(alpha = 0.8) +
  theme_grey() +
  labs(
    title = "Restoration State Index (RSI) Across Sites",
    x = "Site Category",
    y = "RSI Value"
  ) +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

