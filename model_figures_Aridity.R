#Created by: ANNA DOMÉNECH-PASCUAL
# April 2024

rm(list = ls())
dev.off()

#Load packages ----
library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)


results <- read_excel("model_results_Aridity_v2.xlsx")

results$Sign <- ifelse(results$Standardized > 0, "Positive", "Negative")

results$Std_perc <- abs(results$Std_perc)
results <- subset(results, Ranks <= 5)

# results <- subset(results, Model != "BB")
# results <- subset(results, Model != "FB")
# results <- subset(results, Model != "MB")
# results <- subset(results, Model != "alphabeta")
# results <- subset(results, Model != "xylcbh")
# results <- subset(results, Model != "Cenz")

results$Model <- str_replace(results$Model,'MB','Microbial Biomass')
results$Model <- str_replace(results$Model,'FOS','PHOS')
results$Model <- str_replace(results$Model,'GLA','NAG')
results$Variables <- str_replace(results$Variables,'C_N','C/N')
results$Variables <- str_replace(results$Variables,'Soil_Temp','STemp')
results$Variables <- str_replace(results$Variables,'Peak_B','Peak B')
results$Variables <- str_replace(results$Variables, 'L_TN','LTN')
results$Variables <- str_replace(results$Variables,'LTC_LTN','LTC/LTN')

#Reorder rows for ggarrange plots in desired order:
new_order <- c(34:38, 29:33, 16:23, 11:15, 24:28, 1:10)
results <- results[new_order, ]


unique_models <- unique(results$Model)
plots <- list()

for (model in unique_models) {
  model_data <- subset(results, Model == model)
  model_data <- arrange(model_data, Std_perc)
  
  plot <- ggplot(model_data, aes(x = reorder(Variables, Std_perc), y = Std_perc, fill = Sign)) +
    geom_bar(stat = "identity") +
    labs(title = paste(model), x = "Variables") +
    coord_flip()+
    theme_bw()+
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())+
    scale_fill_manual(values = c(Positive = "#5F8249", Negative = "#CE5C17"))+
    ylim(0,60)
  
  plots[[model]] <- plot
}


combined_plots <- wrap_plots(plots, nrow = 2) +
  theme(legend.position = "bottom",
        legend.justification = c("right", "bottom"))+
  guides(fill = guide_legend(reverse = TRUE))

print(combined_plots)

# ggsave("Figures/1 GRADIENT/Drivers_models_Aridity_v2.png", width = 12, height = 6, dpi = 300)



# POSTER:

unique_models <- unique(results$Model)
plots <- list()

for (model in unique_models) {
  model_data <- subset(results, Model == model)
  model_data <- arrange(model_data, Std_perc)
  
  plot <- ggplot(model_data, aes(x = reorder(Variables, Std_perc), y = Std_perc, fill = Sign)) +
    geom_bar(stat = "identity") +
    labs(title = paste(model), x = "Variables") +
    coord_flip()+
    theme_bw()+
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 15),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12))+
    scale_fill_manual(values = c(Positive = "#5F8249", Negative = "#E97132"))+
    ylim(0,60)
  
  plots[[model]] <- plot
}


combined_plots <- wrap_plots(plots, nrow = 2) +
  theme(legend.position = "bottom",
        legend.justification = c("right", "bottom"))+
  guides(fill = guide_legend(reverse = TRUE))

print(combined_plots)


ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/- Congressos/4. ISME (2024)/POSTER/Figures and Images/high resolution figures",
       "Drivers_models_Aridity_v2.png",
       width = 24, height = 12, dpi = 1000)

ggsave(path = "C:/Users/Anna/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
       "Drivers_models_Aridity_v2.png",
       width = 12, height = 6, dpi = 1000)


