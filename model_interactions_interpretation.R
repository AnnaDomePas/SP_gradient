#INTERACTIONS INTERPRETATION 
# https://cran.r-project.org/web/packages/interactions/vignettes/interactions.html


#Load packages ----
library(dplyr)
library(lme4)
library(interactions)
library(ggpubr)
library(ggplot2)

# Load and edit data ----
my_data         = read.csv("SP_metadata_2021.csv", sep=";")

my_data$LTC_LTN <- my_data$L_TC/my_data$L_TN
my_data$xylcbh <- my_data$xyl + my_data$cbh
my_data$alphabeta <- my_data$alpha + my_data$beta

#NEW VARIABLE: ARIDITY
my_data$Aridity <- (1 - my_data$AI)

#To replace NA values with a mean of the other values of the Site:
for (i in which(sapply(my_data, is.numeric))) {
  for (j in which(is.na(my_data[, i]))) {
    my_data[j, i] <- mean(my_data[my_data[, "Site"] == my_data[j, "Site"], i],  na.rm = TRUE)
  }
}

max_values      = apply(my_data[,c(8:34,40,41,76,81:93)],2,max)
my_data.1       = as.data.frame(cbind(my_data[,c(2,96,46:53,64,94,95)],
                                      my_data[,c(8:34,40,41,76,81:93)] / as.list(max_values)))
a.1             = as.data.frame(colnames(my_data.1))

corr            = as.data.frame(cor(my_data.1[,14:56]))

my_data.1 <-my_data.1 %>%
  mutate(Cenz = select(., 3:6) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(MB = select(.,41:42) %>% rowSums(na.rm = TRUE))





#RESPIRATION ----
my_data.resp = my_data.1[-c(9,21,35,39,48),]

respiration.full.2r = lmer(Respiration ~ 1 + Water_content + Clay + BIX + Aridity + BIX:Aridity + 
                             Clay:Aridity + (1 | Site),
                           data = my_data.resp)

# EXEMPLE:
# El plot s'alimenta del model sencer, de manera que té en compte totes
# les altres variables del model. Es basa en efectes marginals.
# Són plots que ajuden a veure molt ràpid el sentit/direcció de la interacció.


# resp.bix <- interact_plot(respiration.full.2r,
#                            pred = BIX,
#                            modx = Aridity,
#                            modx.values = 'terciles',
#                            modx.labels=c('Low','Medium','High'),
#                            legend.main='Aridity level')
# 
# resp.bix

# En els llocs més humits (low aridity), a més BIX clarament menys
# respiració, mentres que si ens movem a llocs més àrids l'efecte disminueix
# i fins i tot sembla que s'inverteix



resp.clay <- interact_plot(respiration.full.2r,
              pred = Clay,
              modx = Aridity,
              modx.values = 'terciles',
              modx.labels=c('Low','Medium','High'),
              legend.main='Aridity level',
              colors = c("deepskyblue4", "green", "#C8404A"))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.3, 0.2, 0, 0),"cm"))+
  labs(y = "Respiration")

resp.clay

# L'efecte de clay no és massa important als llocs humits (low aridity),
# però que a mesura que ens desplacem a llocs més àrids el
# contingut de clay té un efecte negatiu en la respiració.



#ALPHABETA ----
my_data.alphabeta = my_data.1[-c(6,9,37,45),]

alphabeta.full.2r = lmer(log(alphabeta) ~ 1 + Silt + L_TN + HIX + Soil_Temp + C_N + SR +  
                           FB + PO43 + Water_content + SO42 + Aridity +
                           E4.E6 + BIX + LTC_LTN + Litter + BB + SOM + NH4 +
                           C_N:Aridity + PO43:Aridity + L_TN:Aridity + HIX:Aridity+
                           SO42:Aridity + Aridity:Litter + Aridity:BIX + SR:Aridity + (1|Site), data = my_data.alphabeta)                                         


alphabeta.hix <- interact_plot(alphabeta.full.2r,
                               pred = HIX,
                               modx = Aridity,
                               modx.values = 'terciles',
                               modx.labels=c('Low','Medium','High'),
                               legend.main='Aridity level',
                               colors = c("deepskyblue4", "green", "#C8404A"))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.3, 0.1, 0, 0.2),"cm"))+
  labs(y = "log(ALPHABETA)")


alphabeta.CN <- interact_plot(alphabeta.full.2r,
                               pred = C_N,
                               modx = Aridity,
                               modx.values = 'terciles',
                               modx.labels=c('Low','Medium','High'),
                               legend.main='Aridity level',
                              colors = c("deepskyblue4", "green", "#C8404A"))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.3, 0, 0, 0.5),"cm"))+
  labs(x = "C/N", y = "log(ALPHABETA)")


#Trying to put them in the same grid:
# alphabeta.hix <- interact_plot(alphabeta.full.2r,
#                                pred = HIX,
#                                modx = Aridity,
#                                modx.values = 'terciles',
#                                modx.labels=c('Low','Medium','High'),
#                                legend.main= "Aridity level",
#                                colors = c("deepskyblue4", "green", "#C8404A"))+
#   theme(legend.position = "none")+
#   labs(y = "log(ALPHABETA)")
# 
# alphabeta.hix
# 
# 
# alphabeta.CN <- interact_plot(alphabeta.full.2r,
#                               pred = C_N,
#                               modx = Aridity,
#                               modx.values = 'terciles',
#                               modx.labels=c('Low','Medium','High'),
#                               legend.main= "Aridity level",
#                               colors = c("deepskyblue4", "green", "#C8404A")) +
#   theme(legend.position = "none")+
#   scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
#   labs(x = "C/N", y = NULL)
# 
# 
# plot_list = list(alphabeta.hix,alphabeta.CN) 
# 
# alphabeta <- ggarrange(plotlist=plot_list, labels = c("a", "b"),
#           common.legend = TRUE, 
#           legend = "right")






#XYLCBH ----

xylcbh.full.2 = lmer(xylcbh ~ 1 + L_TN + HIX + Aridity + L_TN:Aridity + C_N + Soil_Temp + E4.E6 +  
                       Silt + FB + TC + Silt:Aridity + (1|Site), data = my_data.1)   

xylcbh.silt <- interact_plot(xylcbh.full.2,
                              pred = Silt,
                              modx = Aridity,
                              modx.values = 'terciles',
                              modx.labels=c('Low','Medium','High'),
                              legend.main='Aridity level',
                             colors = c("deepskyblue4", "green", "#C8404A"))+
  theme(legend.position = "none")+
  labs(y = "XYLCBH")

xylcbh.silt




#Saving all plots together ----
plot_list_all = list(resp.clay, xylcbh.silt, alphabeta.hix, alphabeta.CN) 

all <- ggarrange(plotlist=plot_list_all,
                 labels = c("A", "B", "C", "D"),
                 common.legend = TRUE, 
                 legend = "right")
all <- all + theme(plot.background = element_rect(fill = "white"))


interaction_plots_models <- ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/ANALISIS/R/SP_gradient/Figures/1 GRADIENT",
                                   "interactions_models.png", width = 9, height = 9, dpi = 400, device = "png")








# # ___________________________________________________________________________________________________________________
# #Load packages ----
# library(dplyr)
# library(lme4)
# 
# # Load and edit data ----
# my_data         = read.csv("SP_metadata_2021.csv", sep=";")
# 
# my_data$LTC_LTN <- my_data$L_TC/my_data$L_TN
# my_data$xylcbh <- my_data$xyl + my_data$cbh
# my_data$alphabeta <- my_data$alpha + my_data$beta
# 
# #NEW VARIABLE: ARIDITY
# my_data$Aridity <- (1 - my_data$AI)
# 
# #To replace NA values with a mean of the other values of the Site:
# for (i in which(sapply(my_data, is.numeric))) {
#   for (j in which(is.na(my_data[, i]))) {
#     my_data[j, i] <- mean(my_data[my_data[, "Site"] == my_data[j, "Site"], i],  na.rm = TRUE)
#   }
# }
# 
# max_values      = apply(my_data[,c(8:34,40,41,76,81:93)],2,max)
# my_data.1       = as.data.frame(cbind(my_data[,c(2,96,46:53,64,94,95)],
#                                       my_data[,c(8:34,40,41,76,81:93)] / as.list(max_values)))
# a.1             = as.data.frame(colnames(my_data.1))
# 
# corr            = as.data.frame(cor(my_data.1[,14:56]))
# 
# my_data.1 <-my_data.1 %>%
#   mutate(Cenz = select(., 3:6) %>% rowSums(na.rm = TRUE)) %>% 
#   mutate(MB = select(.,41:42) %>% rowSums(na.rm = TRUE))
# 
# 
# #RESPIRATION ----
# # > Clay:Aridity ----
# # interaction.plot(
# #   x.factor = my_data$Clay,
# #   trace.factor = my_data$Respiration,
# #   response = my_data$Aridity,
# #   fun = mean,
# #   ylab = "Aridity",
# #   xlab = "Clay",
# #   trace.label = "Respiration"
# # )
# 
# my_data.resp = my_data.1[-c(9,21,35,39,48),]
# 
# respiration.full.2r = lmer(Respiration ~ 1 + Water_content + Clay + BIX + Aridity + BIX:Aridity + 
#                              Clay:Aridity + (1 | Site),
#                            data = my_data.resp)  
# 
# 
# 
# 
# # Create a grid of values for Clay, Aridity, Water_content, and BIX
# clay_values <- seq(min(my_data.resp$Clay), max(my_data.resp$Clay), length.out = 50)
# aridity_values <- seq(min(my_data.resp$Aridity), max(my_data.resp$Aridity), length.out = 50)
# water_content_values <- seq(min(my_data.resp$Water_content), max(my_data.resp$Water_content), length.out = 50)
# bix_values <- seq(min(my_data.resp$BIX), max(my_data.resp$BIX), length.out = 50)
# 
# grid <- expand.grid(Clay = clay_values, Aridity = aridity_values, Water_content = water_content_values, BIX = bix_values)
# 
# # Predict Respiration values for the grid
# pred <- predict(respiration.full.2r, newdata = grid, re.form = NA)
# 
# # Combine predicted values with the grid
# grid$respiration <- pred
# 
# # Plot the interaction
# interaction_plot <- plot(grid$Clay, grid$Aridity, type = "n", xlab = "Clay", ylab = "Aridity", main = "Interaction Plot: Respiration vs. Clay and Aridity")
# points(grid$Clay, grid$Aridity, col = "blue", pch = 20, cex = 0.5)
# 
# # Reshape the data for contour plotting
# contour_data <- with(grid, contour(x = sort(unique(Clay)), y = sort(unique(Aridity)), z = matrix(respiration, nrow = length(unique(Aridity)), byrow = TRUE)))
# 
# # Plot the contour
# contour(x = contour_data$x, y = contour_data$y, z = contour_data$z, add = TRUE, labels = contour_data$level)# Add legend
# legend("topright", legend = "Respiration", col = "blue", pch = 20)
# 
# # Display the plot
# interaction_plot
