
# Reference paper: https://bsssjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/ejss.13419

## SP gradient-----------
rm(list = ls())

# LOAD LIBRARY ----
library("ggplot2")
library("scales")
library("tidyr")
library("dplyr")
library("ggpmisc")
library("ggpubr")
library("plotly")
library(plyr)
library(RColorBrewer)
library(devtools)
# install_github("vqv/ggbiplot")
library(ggbiplot)
library(vegan)
library(FSA)
# install.packages("DescTools")
library(DescTools)
library(rcompanion)

# New Liner mixed model ####
library(lme4)
library(car)
# remotes::install_version("MuMIn", "1.46.0")
library(MuMIn)
library(cAIC4)
library(domir)

my_data         = read.csv("SP_metadata_2021.csv", sep=";")

my_data$LTC_LTN <- my_data$L_TC/my_data$L_TN
my_data$xylcbh <- my_data$xyl + my_data$cbh
my_data$alphabeta <- my_data$alpha + my_data$beta

#To replace NA values with a mean of the other values of the Site:
for (i in which(sapply(my_data, is.numeric))) {
  for (j in which(is.na(my_data[, i]))) {
    my_data[j, i] <- mean(my_data[my_data[, "Site"] == my_data[j, "Site"], i],  na.rm = TRUE)
  }
}

# Checking data outliers ####

ggplot(my_data, aes(x=as.factor(AI), y=Soil_Temp, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Soil_Temp)

ggplot(my_data, aes(x=as.factor(AI), y=Water_activity, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Water_activity)

ggplot(my_data, aes(x=as.factor(AI), y=SOM, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$SOM)

ggplot(my_data, aes(x=as.factor(AI), y=pH, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$pH)

ggplot(my_data, aes(x=as.factor(AI), y=TOC, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$TOC)

ggplot(my_data, aes(x=as.factor(AI), y=C_N, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$C_N)

ggplot(my_data, aes(x=as.factor(AI), y=TC, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$TC)

ggplot(my_data, aes(x=as.factor(AI), y=TN, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$TN)

ggplot(my_data, aes(x=as.factor(AI), y=NH4, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$NH4)

ggplot(my_data, aes(x=as.factor(AI), y=PO43, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$PO43)

ggplot(my_data, aes(x=as.factor(AI), y=SO42, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$SO42)

ggplot(my_data, aes(x=as.factor(Site), y=Sand, fill=as.factor(Site))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Sand)

ggplot(my_data, aes(x=as.factor(Site), y=Silt, fill=as.factor(Site))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Silt)

ggplot(my_data, aes(x=as.factor(Site), y=Clay, fill=as.factor(Site))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Clay)

ggplot(my_data, aes(x=as.factor(AI), y=Litter, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Litter)

ggplot(my_data, aes(x=as.factor(AI), y=L_TC, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$L_TC)

ggplot(my_data, aes(x=as.factor(AI), y=L_TN, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$L_TN)

ggplot(my_data, aes(x=as.factor(AI), y=BB, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$BB)

ggplot(my_data, aes(x=as.factor(AI), y=FB, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$FB)

ggplot(my_data, aes(x=as.factor(AI), y=ShannonEEA, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$ShannonEEA)

ggplot(my_data, aes(x=as.factor(AI), y=SR, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$SR)

ggplot(my_data, aes(x=as.factor(AI), y=E2.E3, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$E2.E3)

ggplot(my_data, aes(x=as.factor(AI), y=E3.E4, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$E3.E4)

ggplot(my_data, aes(x=as.factor(AI), y=E4.E6, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$E4.E6)

ggplot(my_data, aes(x=as.factor(AI), y=FI, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$FI)

ggplot(my_data, aes(x=as.factor(AI), y=BIX, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$BIX)

ggplot(my_data, aes(x=as.factor(AI), y=Peak_A, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Peak_A)

ggplot(my_data, aes(x=as.factor(AI), y=Peak_C, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Peak_C)

ggplot(my_data, aes(x=as.factor(AI), y=Peak_M, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Peak_M)

ggplot(my_data, aes(x=as.factor(AI), y=Peak_T, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Peak_T)

ggplot(my_data, aes(x=as.factor(AI), y=Peak_B, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$Peak_B)

ggplot(my_data, aes(x=as.factor(AI), y=HIX, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$HIX)

# Data normalization ####
# We are normalizing all the variables to use the coefficients to assess the
# importance of the variables for the model

a               = as.data.frame(colnames(my_data))
max_values      = apply(my_data[,c(8:34,40,41,76,81:93)],2,max)
my_data.1       = as.data.frame(cbind(my_data[,c(7,46:53,64,94,95)],
                                      my_data[,c(8:34,40,41,76,81:93)] / as.list(max_values)))
a.1             = as.data.frame(colnames(my_data.1))

corr            = as.data.frame(cor(my_data.1[,13:55]))

# Parameters correlated: Water activity-water content; TOC-TN,FB,BB,TC;
# Sand-silt,clay; Peak_A-Peak_M, Peak_C; Peak_T-Peak_B; E2.E3-E3.E4

my_data.1 <-my_data.1 %>%
  mutate(Cenz = select(., 3:6) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(MB = select(.,21:22) %>% rowSums(na.rm = TRUE))

#1rst. Full model
#2nd. Remove variables with Estimates <0.1 --> full model II
#3rd. Remove variables with p-val no sig. --> full model III
#(Repeat 3 until everything is sig.)
#(Keep checking AIC and R2)


# Respiration ####

# Full model interactions ####

respiration.full.I = lmer(Respiration ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                                    L_TN+BB+FB+BIX+Soil_Temp+
                                                    Water_content+pH+PO43+Litter+
                                                    SR+E2.E3+FI+HIX+Peak_A+
                                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         
   
summary(respiration.full.I)
Anova(respiration.full.I)
cAIC(respiration.full.I)
AIC(respiration.full.I)
r.squaredGLMM(respiration.full.I)

qqnorm(residuals(respiration.full.I))
scatter.smooth(residuals(respiration.full.I) ~ fitted(respiration.full.I))

# Full model II  interactions ####

respiration.full.II = lmer(Respiration ~ altitude+(Silt+Clay+L_TC+
                                                     BB+TOC+L_TN+
                                                     Water_content+
                                                     SR+FI)*AI + (1|Site), data = my_data.1)                                         

summary(respiration.full.II)
Anova(respiration.full.II)
cAIC(respiration.full.II)
AIC(respiration.full.II)
r.squaredGLMM(respiration.full.II)

qqnorm(residuals(respiration.full.II))
scatter.smooth(residuals(respiration.full.II) ~ fitted(respiration.full.II))

# Full model III  interactions ####

respiration.full.III = lmer(Respiration ~ altitude+Water_content+Silt+Clay+Silt:AI+
                              Clay:AI+SR:AI+(1|Site), data = my_data.1)                                         

summary(respiration.full.III)
Anova(respiration.full.III)
cAIC(respiration.full.III)
AIC(respiration.full.III)
r.squaredGLMM(respiration.full.III)

qqnorm(residuals(respiration.full.III))
scatter.smooth(residuals(respiration.full.III) ~ fitted(respiration.full.III))

# Importance assessment ####

domin(Respiration ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("altitude","Water_content","Silt","Silt:AI","Clay","Clay:AI","SR:AI"), 
      consmodel = "(1|Site)")



# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(Respiration ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("altitude","Water_content","Silt","Silt:AI","Clay","Clay:AI","SR:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(respiration.full.III)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("altitude","WC","Silt","Silt:AI","Clay","Clay:AI","SR:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)

# 
# Respi <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
#   geom_bar(stat = "identity",
#            color = "black") +
#   labs(
#     y = "Predictor Variables",
#     x = "Standardized dominance"
#   ) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
#             position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
#             color = "black", hjust = -0.3) +  # Adjust label position
#   ggtitle("Respiration") +
#   theme(legend.position = "none")+
#   theme(legend.title=element_blank())+
#   scale_fill_manual(values = c("#B54545", "#5F8249"), name = "Standardized",
#                     labels = c("Negative", "Positive"))+
#   theme(axis.title.x =element_blank())+
#   xlim(0, 1)
# Respi

Respi <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle("Respiration") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)


Respi

# Enzymes ####

ggplot(my_data, aes(x=as.factor(AI), y=alpha, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$alpha)

ggplot(my_data, aes(x=as.factor(AI), y=beta, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$beta)

ggplot(my_data, aes(x=as.factor(AI), y=xyl, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$xyl)

ggplot(my_data, aes(x=as.factor(AI), y=cbh, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$cbh)

ggplot(my_data, aes(x=as.factor(AI), y=gla, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$gla)

ggplot(my_data, aes(x=as.factor(AI), y=fos, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$fos)

ggplot(my_data, aes(x=as.factor(AI), y=leu, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$leu)

ggplot(my_data, aes(x=as.factor(AI), y=phe, fill=as.factor(AI))) +
  geom_boxplot(alpha=0.7)
plot(my_data$AI,my_data$phe)

# alpha ####

# Test alpha####

plot(my_data.1$Soil_Temp,my_data.1$alpha)
plot(my_data.1$Peak_T,my_data.1$alpha)

library(mgcv)

gam_mod = gam(alpha ~ s(BB) + s(as.factor(my_data.1$Site), bs = "re"),data = my_data.1, method = "REML")
gam.check(gam_mod)
summary(gam_mod)
plot(resid(gam_mod)~fitted(gam_mod))

# Full model interactions ####

alpha.full.I = lmer(alpha ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                        L_TN+BB+FB+BIX+Soil_Temp+
                                        Water_content+pH+PO43+Litter+
                                        SR+E2.E3+FI+HIX+Peak_A+
                                        Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(alpha.full.I, tol = 1e-4)


summary(alpha.full.I)
Anova(alpha.full.I)
cAIC(alpha.full.I)
AIC(alpha.full.I)
r.squaredGLMM(alpha.full.I)

qqnorm(residuals(alpha.full.I))
scatter.smooth(residuals(alpha.full.I) ~ fitted(alpha.full.I))

# Full model interactions II ####

alpha.full.II = lmer(alpha ~ (Clay+L_TC+L_TN+BB+Water_content+pH+SR+FI)*AI + 
                       (1|Site), data = my_data.1)                                         

summary(alpha.full.II)
Anova(alpha.full.II)
cAIC(alpha.full.II)
AIC(alpha.full.II)
r.squaredGLMM(alpha.full.II)

qqnorm(residuals(alpha.full.II))
scatter.smooth(residuals(alpha.full.II) ~ fitted(alpha.full.II))

# Full model interactions III ####

alpha.full.III = lmer(alpha ~ (Clay+L_TC+L_TN+BB+Water_content)*AI + 
                        (1|Site), data = my_data.1) 

summary(alpha.full.III)
Anova(alpha.full.III)
cAIC(alpha.full.III)
AIC(alpha.full.III)
r.squaredGLMM(alpha.full.III)

qqnorm(residuals(alpha.full.III))
scatter.smooth(residuals(alpha.full.III) ~ fitted(alpha.full.III))

# Full model interactions 4 ####

alpha.full.4 = lmer(alpha ~ L_TC+L_TC:AI+BB+BB:AI+(1|Site), data = my_data.1) 

summary(alpha.full.4)
Anova(alpha.full.4)
cAIC(alpha.full.4)
AIC(alpha.full.4)
r.squaredGLMM(alpha.full.4)

qqnorm(residuals(alpha.full.4))
scatter.smooth(residuals(alpha.full.4) ~ fitted(alpha.full.4))

# Importance assessment ####

domin(alpha ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("L_TC","L_TC:AI","BB","BB:AI"), 
      consmodel = "(1|Site)")


# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(alpha ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("L_TC","L_TC:AI","BB","BB:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(alpha.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("LTC","LTC:AI","BB","BB:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


alpha <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle(expression(alpha-glucosidase)) +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)+
  theme(axis.title.y =element_blank())

alpha

# beta ####

# Full model interactions ####

beta.full.I = lmer(beta ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                      L_TN+BB+FB+BIX+Soil_Temp+
                                      Water_content+pH+PO43+Litter+
                                      SR+E2.E3+FI+HIX+Peak_A+
                                      Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(beta.full.I, tol = 1e-4)


summary(beta.full.I)
Anova(beta.full.I)
cAIC(beta.full.I)
AIC(beta.full.I)
r.squaredGLMM(beta.full.I)

qqnorm(residuals(beta.full.I))
scatter.smooth(residuals(beta.full.I) ~ fitted(beta.full.I))

# Full model II interactions ####

beta.full.II = lmer(beta ~altitude+(TOC+L_TC+BB+pH+Soil_Temp+Water_content+pH+
                                      PO43+Litter+HIX)*AI + (1|Site), data = my_data.1)                                         
isSingular(beta.full.II, tol = 1e-4)


summary(beta.full.II)
Anova(beta.full.II)
cAIC(beta.full.II)
AIC(beta.full.II)
r.squaredGLMM(beta.full.II)

qqnorm(residuals(beta.full.I))
scatter.smooth(residuals(beta.full.I) ~ fitted(beta.full.I))

# Full model III interactions ####

beta.full.III = lmer((beta) ~ HIX+Soil_Temp + (1|Site), data = my_data.1)                                         
isSingular(beta.full.III, tol = 1e-4)

summary(beta.full.III)
Anova(beta.full.III)
cAIC(beta.full.III)
AIC(beta.full.III)
r.squaredGLMM(beta.full.III)

qqnorm(residuals(beta.full.III))
scatter.smooth(residuals(beta.full.III) ~ fitted(beta.full.III))

# Importance assessment ####

domin(beta ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("HIX","Soil_Temp"), 
      consmodel = "(1|Site)")


# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(beta ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("HIX","Soil_Temp"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(beta.full.III)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("HIX","STemp"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


beta <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle(expression(beta-glucosidase)) +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)+
  theme(axis.title.y =element_blank())

beta




# xyl ####

# Full model interactions ####

xyl.full.I = lmer(xyl ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                      L_TN+BB+FB+BIX+Soil_Temp+
                                      Water_content+pH+PO43+Litter+
                                      SR+E2.E3+FI+HIX+Peak_A+
                                      Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(xyl.full.I, tol = 1e-4)


summary(xyl.full.I)
Anova(xyl.full.I)
cAIC(xyl.full.I)
AIC(xyl.full.I)
r.squaredGLMM(xyl.full.I)

qqnorm(residuals(xyl.full.I))
scatter.smooth(residuals(xyl.full.I) ~ fitted(xyl.full.I))

# Full model II interactions ####

xyl.full.II = lmer(xyl ~ (TOC+L_TC+L_TN+BB+FB+Soil_Temp)*AI + (1|Site), data = my_data.1)                                         
isSingular(xyl.full.II, tol = 1e-4)

summary(xyl.full.II)
Anova(xyl.full.II)
cAIC(xyl.full.II)
AIC(xyl.full.II)
r.squaredGLMM(xyl.full.II)

qqnorm(residuals(xyl.full.II))
scatter.smooth(residuals(xyl.full.II) ~ fitted(xyl.full.II))

# Full model III interactions #### 

xyl.full.III = lmer(xyl ~ L_TN + L_TN:AI + L_TC + L_TC:AI + (1|Site), data = my_data.1)                                         
isSingular(xyl.full.III, tol = 1e-4)

summary(xyl.full.III)
Anova(xyl.full.III)
cAIC(xyl.full.III)
AIC(xyl.full.III)
r.squaredGLMM(xyl.full.III)

qqnorm(residuals(xyl.full.III))
scatter.smooth(residuals(xyl.full.III) ~ fitted(xyl.full.III))

# Importance assessment ####

domin(xyl ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("L_TN","L_TN:AI","L_TC","L_TC:AI"), 
      consmodel = "(1|Site)")


# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(xyl ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("L_TN","L_TN:AI","L_TC","L_TC:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(xyl.full.III)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("LTN","LTN:AI","LTC","LTC:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


xyl <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle("Xylosidase") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)+
  theme(axis.title.y =element_blank())

xyl







# cbh ####

# Full model interactions ####

cbh.full.I = lmer(cbh ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                    L_TN+BB+FB+BIX+Soil_Temp+
                                    Water_content+pH+PO43+Litter+
                                    SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(cbh.full.I, tol = 1e-4)


summary(cbh.full.I)
Anova(cbh.full.I)
cAIC(cbh.full.I)
AIC(cbh.full.I)
r.squaredGLMM(cbh.full.I)

qqnorm(residuals(cbh.full.I))
scatter.smooth(residuals(cbh.full.I) ~ fitted(cbh.full.I))

# Full model II interactions ####

cbh.full.II = lmer(cbh ~ (TOC+L_TC+L_TN+BB+PO43)*AI + (1|Site), data = my_data.1)                                         
isSingular(cbh.full.II, tol = 1e-4)


summary(cbh.full.II)
Anova(cbh.full.II)
cAIC(cbh.full.II)
AIC(cbh.full.II)
r.squaredGLMM(cbh.full.II)

qqnorm(residuals(cbh.full.II))
scatter.smooth(residuals(cbh.full.II) ~ fitted(cbh.full.II))

# Full model III interactions ####

cbh.full.III = lmer(cbh ~ (L_TC+L_TN)*AI + (1|Site), data = my_data.1)                                         
isSingular(cbh.full.III, tol = 1e-4)


summary(cbh.full.III)
Anova(cbh.full.III)
cAIC(cbh.full.III)
AIC(cbh.full.III)
r.squaredGLMM(cbh.full.III)

qqnorm(residuals(cbh.full.III))
scatter.smooth(residuals(cbh.full.III) ~ fitted(cbh.full.III))

# Importance assessment ####

domin(cbh ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("L_TN","L_TN:AI","L_TC","L_TC:AI","AI"), 
      consmodel = "(1|Site)")


# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(cbh ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("L_TN","L_TN:AI","L_TC","L_TC:AI","AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(cbh.full.III)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("LTN","LTN:AI","LTC","LTC:AI","AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


cbh <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle("Cellobiohydrolase") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)





# gla ####

# Full model interactions ####

gla.full.I = lmer(gla ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                    L_TN+BB+FB+BIX+Soil_Temp+
                                    Water_content+pH+PO43+Litter+
                                    SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(gla.full.I, tol = 1e-4)


summary(gla.full.I)
Anova(gla.full.I)
cAIC(gla.full.I)
AIC(gla.full.I)
r.squaredGLMM(gla.full.I)

qqnorm(residuals(gla.full.I))
scatter.smooth(residuals(gla.full.I) ~ fitted(gla.full.I))

# Full model II interactions ####

gla.full.II = lmer(gla ~ altitude+(SO42+BB+FB+BIX+Water_content+PO43+
                                    SR+E2.E3+FI+HIX+Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(gla.full.II, tol = 1e-4)


summary(gla.full.II)
Anova(gla.full.II)
cAIC(gla.full.II)
AIC(gla.full.II)
r.squaredGLMM(gla.full.II)

qqnorm(residuals(gla.full.II))
scatter.smooth(residuals(gla.full.II) ~ fitted(gla.full.II))

# Full model III interactions ####

gla.full.III = lmer(gla ~ (BB+Water_content+FI+PO43+E2.E3+SO42+FB+SR+FI+HIX)*AI+ (1|Site), data = my_data.1)                                         
isSingular(gla.full.III, tol = 1e-4)

summary(gla.full.III)
Anova(gla.full.III)
cAIC(gla.full.III)
AIC(gla.full.III)
r.squaredGLMM(gla.full.III)

qqnorm(residuals(gla.full.III))
scatter.smooth(residuals(gla.full.III) ~ fitted(gla.full.III))

# Full model 4 interactions #### BB+Water_content+PO43+E2.E3+SO42:AI+FB:AI+Water_content:AI+SR:AI+FI:AI+HIX:AI+(1|Site)

gla.full.4 = lm(gla ~ Water_content+Water_content:AI+SR:AI, data = my_data.1)                                         
isSingular(gla.full.4, tol = 1e-4)

summary(gla.full.4)
Anova(gla.full.4)
cAIC(gla.full.4)
AIC(gla.full.III)
r.squaredGLMM(gla.full.4)

qqnorm(residuals(gla.full.4))
scatter.smooth(residuals(gla.full.4) ~ fitted(gla.full.4))

# Importance assessment ####

domin(gla ~ 1, 
      lm, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("Water_content","Water_content:AI","SR:AI"))


# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(gla ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("Water_content","Water_content:AI","SR:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(gla.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("WC","WC:AI","SR:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


gla <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle(expression(beta-glucosaminidase)) +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)+
  theme(axis.title.y =element_blank())

gla






# fos ####

# Full model interactions ####

fos.full.I = lmer(fos ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                    L_TN+BB+FB+BIX+Soil_Temp+
                                    Water_content+pH+PO43+Litter+
                                    SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(fos.full.I, tol = 1e-4)


summary(fos.full.I)
Anova(fos.full.I)
cAIC(fos.full.I)
AIC(fos.full.I)
r.squaredGLMM(fos.full.I)

qqnorm(residuals(fos.full.I))
scatter.smooth(residuals(fos.full.I) ~ fitted(fos.full.I))

# Full model II interactions ####

fos.full.II = lmer(fos ~ altitude+(Silt+Clay+NH4+
                                    L_TN+FB+BIX+Soil_Temp+
                                    Water_content+pH+
                                    SR+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(fos.full.II, tol = 1e-4)


summary(fos.full.II)
Anova(fos.full.II)
cAIC(fos.full.II)
AIC(fos.full.II)
r.squaredGLMM(fos.full.II)

qqnorm(residuals(fos.full.II))
scatter.smooth(residuals(fos.full.II) ~ fitted(fos.full.II))

# Full model III interactions ####

fos.full.III = lmer(fos ~ BIX+Water_content+SR+Peak_A+AI+Silt:AI+Clay:AI+Soil_Temp:AI+
                      pH:AI+FI:AI+Peak_A:AI+Peak_T:AI+(1|Site), data = my_data.1)                                         
isSingular(fos.full.III, tol = 1e-4)


summary(fos.full.III)
Anova(fos.full.III)
cAIC(fos.full.III)
AIC(fos.full.III)
r.squaredGLMM(fos.full.III)

qqnorm(residuals(fos.full.III))
scatter.smooth(residuals(fos.full.III) ~ fitted(fos.full.III))

# Full model 4 interactions ####

fos.full.4 = lmer(fos ~ BIX+Water_content+Peak_A+Soil_Temp:AI+
                      pH:AI+FI:AI+(1|Site), data = my_data.1)                                         
isSingular(fos.full.4, tol = 1e-4)

summary(fos.full.4)
Anova(fos.full.4)
cAIC(fos.full.4)
AIC(fos.full.4)
r.squaredGLMM(fos.full.4)

qqnorm(residuals(fos.full.4))
scatter.smooth(residuals(fos.full.4) ~ fitted(fos.full.4))

# Importance assessment ####

domin(fos ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("BIX","Water_content","Peak_A","Soil_Temp:AI","pH:AI","FI:AI"), 
      consmodel = "(1|Site)")


# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(fos ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("BIX","Water_content","Peak_A","Soil_Temp:AI","pH:AI","FI:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(fos.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("BIX","WC","PeakA","STemp:AI","pH:AI","FI:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


fos <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle("Phosphatase") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)+
  theme(axis.title.y =element_blank())

fos






# leu ####

# Full model interactions ####

leu.full.I = lmer(leu ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                    L_TN+BB+FB+BIX+Soil_Temp+
                                    Water_content+pH+PO43+Litter+
                                    SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(leu.full.I, tol = 1e-4)


summary(leu.full.I)
Anova(leu.full.I)
cAIC(leu.full.I)
AIC(leu.full.I)
r.squaredGLMM(leu.full.I)

qqnorm(residuals(leu.full.I))
scatter.smooth(residuals(leu.full.I) ~ fitted(leu.full.I))

# Full model II interactions ####

leu.full.II = lmer(leu ~ (TOC+NH4+BIX+Water_content+SR+E2.E3+Peak_A)*AI + 
                     (1|Site), data = my_data.1)                                         
isSingular(leu.full.II, tol = 1e-4)


summary(leu.full.II)
Anova(leu.full.II)
cAIC(leu.full.II)
AIC(leu.full.II)
r.squaredGLMM(leu.full.II)

qqnorm(residuals(leu.full.II))
scatter.smooth(residuals(leu.full.II) ~ fitted(leu.full.II))

# Full model III interactions ####

leu.full.III = lmer(leu ~ TOC + NH4 + Water_content + SR + Peak_A + 
                     (1|Site), data = my_data.1)                                         
isSingular(leu.full.III, tol = 1e-4)


summary(leu.full.III)
Anova(leu.full.III)
cAIC(leu.full.III)
AIC(leu.full.III)
r.squaredGLMM(leu.full.III)

qqnorm(residuals(leu.full.III))
scatter.smooth(residuals(leu.full.III) ~ fitted(leu.full.III))

# Full model 4 interactions ####

leu.full.4 = lmer(leu ~ TOC + Water_content + SR + Peak_A + 
                      (1|Site), data = my_data.1)                                         
isSingular(leu.full.4, tol = 1e-4)


summary(leu.full.4)
Anova(leu.full.4)
cAIC(leu.full.4)
AIC(leu.full.4)
r.squaredGLMM(leu.full.4)

qqnorm(residuals(leu.full.4))
scatter.smooth(residuals(leu.full.4) ~ fitted(leu.full.4))

# Importance assessment ####

domin(leu ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("TOC","Water_content","SR","Peak_A"), 
      consmodel = "(1|Site)")


# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(leu ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("TOC","Water_content","SR","Peak_A"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(leu.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("TOC","WC","SR","PeakA"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


leu <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle("Leu-aminopeptidase") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)+
  theme(axis.title.y =element_blank())

leu







# phe ####

# Full model interactions ####

phe.full.I = lmer(phe ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                    L_TN+BB+FB+BIX+Soil_Temp+
                                    Water_content+pH+PO43+Litter+
                                    SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(phe.full.I, tol = 1e-4)


summary(phe.full.I)
Anova(phe.full.I)
cAIC(phe.full.I)
AIC(phe.full.I)
r.squaredGLMM(phe.full.I)

qqnorm(residuals(phe.full.I))
scatter.smooth(residuals(phe.full.I) ~ fitted(phe.full.I))

# Full model II interactions ####

phe.full.II = lmer(phe ~ (L_TC+pH+FB+FI)*AI + (1|Site), data = my_data.1)                                         
isSingular(phe.full.II, tol = 1e-4)

summary(phe.full.II)
Anova(phe.full.II)
cAIC(phe.full.II)
AIC(phe.full.II)
r.squaredGLMM(phe.full.II)

qqnorm(residuals(phe.full.II))
scatter.smooth(residuals(phe.full.II) ~ fitted(phe.full.II))

# Full model III interactions ####

phe.full.III = lmer(phe ~ L_TC + pH + AI + L_TC:AI + (1|Site), data = my_data.1)                                         
isSingular(phe.full.III, tol = 1e-4)

summary(phe.full.III)
Anova(phe.full.III)
cAIC(phe.full.III)
AIC(phe.full.III)
r.squaredGLMM(phe.full.III)

qqnorm(residuals(phe.full.III))
scatter.smooth(residuals(phe.full.III) ~ fitted(phe.full.III))

# Importance assessment ####

domin(phe ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("L_TC","pH","AI","L_TC:AI"), 
      consmodel = "(1|Site)")



# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(phe ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("L_TC","pH","AI","L_TC:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(phe.full.III)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("LTC","pH","AI","LTC:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


phe <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle("Phenol oxidase") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)

phe








# Carbon enzymes ####

# Full model interactions ####

Cenz.full.I = lmer(Cenz ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                    L_TN+BB+FB+BIX+Soil_Temp+
                                    Water_content+pH+PO43+Litter+
                                    SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(Cenz.full.I, tol = 1e-4)


summary(Cenz.full.I)
Anova(Cenz.full.I)
cAIC(Cenz.full.I)
AIC(Cenz.full.I)
r.squaredGLMM(Cenz.full.I)

qqnorm(residuals(Cenz.full.I))
scatter.smooth(residuals(Cenz.full.I) ~ fitted(Cenz.full.I))

# Full model II interactions ####

Cenz.full.II = lmer(Cenz ~ (L_TC+BB+Soil_Temp+HIX)*AI + (1|Site), data = my_data.1)                                         
isSingular(Cenz.full.II, tol = 1e-4)

summary(Cenz.full.II)
Anova(Cenz.full.II)
cAIC(Cenz.full.II)
AIC(Cenz.full.II)
r.squaredGLMM(Cenz.full.II)

qqnorm(residuals(Cenz.full.II))
scatter.smooth(residuals(Cenz.full.II) ~ fitted(Cenz.full.II))

# Full model III interactions ####

Cenz.full.III = lmer(Cenz ~ BB+Soil_Temp+HIX+Soil_Temp:AI + (1|Site), data = my_data.1)                                         
isSingular(Cenz.full.III, tol = 1e-4)

summary(Cenz.full.III)
Anova(Cenz.full.III)
cAIC(Cenz.full.III)
AIC(Cenz.full.III)
r.squaredGLMM(Cenz.full.III)

qqnorm(residuals(Cenz.full.III))
scatter.smooth(residuals(Cenz.full.III) ~ fitted(Cenz.full.III))

# Full model 4 interactions ####

Cenz.full.4 = lmer(Cenz ~ Soil_Temp+HIX + (1|Site), data = my_data.1)                                         
isSingular(Cenz.full.4, tol = 1e-4)

summary(Cenz.full.4)
Anova(Cenz.full.4)
cAIC(Cenz.full.4)
AIC(Cenz.full.4)
r.squaredGLMM(Cenz.full.4)

qqnorm(residuals(Cenz.full.4))
scatter.smooth(residuals(Cenz.full.4) ~ fitted(Cenz.full.4))

# Importance assessment ####

domin(Cenz ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("Soil_Temp","HIX"), 
      consmodel = "(1|Site)")



# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(Cenz ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("Soil_Temp","HIX"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(Cenz.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("SoilTemp","HIX"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


Cenz <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle("C enzymes") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)







# BB ####

# Full model interactions ####

BB.full.I = lmer(BB ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                    L_TN+BIX+Soil_Temp+
                                    Water_content+pH+PO43+Litter+
                                    SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(BB.full.I, tol = 1e-4)


summary(BB.full.I)
Anova(BB.full.I)
cAIC(BB.full.I)
AIC(BB.full.I)
r.squaredGLMM(BB.full.I)

qqnorm(residuals(BB.full.I))
scatter.smooth(residuals(BB.full.I) ~ fitted(BB.full.I))

# Full model II interactions ####

BB.full.II = lmer(BB ~ (SO42+PO43+Clay+pH+E2.E3+Peak_A+Peak_T)*AI + 
                    (1|Site), data = my_data.1)                                         
isSingular(BB.full.II, tol = 1e-4)


summary(BB.full.II)
Anova(BB.full.II)
cAIC(BB.full.II)
AIC(BB.full.II)
r.squaredGLMM(BB.full.II)

qqnorm(residuals(BB.full.II))
scatter.smooth(residuals(BB.full.II) ~ fitted(BB.full.II))

# Full model III interactions ####

BB.full.III = lmer(BB ~ SO42 + PO43 + Clay + Peak_A + SO42:AI + Clay:AI + Peak_A:AI +
                    (1|Site), data = my_data.1)                                         
isSingular(BB.full.III, tol = 1e-4)


summary(BB.full.III)
Anova(BB.full.III)
cAIC(BB.full.III)
AIC(BB.full.III)
r.squaredGLMM(BB.full.III)

qqnorm(residuals(BB.full.III))
scatter.smooth(residuals(BB.full.III) ~ fitted(BB.full.III))

# Full model 4 interactions ####

BB.full.4 = lmer(BB ~ SO42 + PO43 + Clay + Peak_A + SO42:AI + 
                     (1|Site), data = my_data.1)                                         
isSingular(BB.full.4, tol = 1e-4)


summary(BB.full.4)
Anova(BB.full.4)
cAIC(BB.full.4)
AIC(BB.full.4)
r.squaredGLMM(BB.full.4)

qqnorm(residuals(BB.full.4))
scatter.smooth(residuals(BB.full.4) ~ fitted(BB.full.4))

# Importance assessment ####

domin(BB ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("SO42","PO43","Clay","Peak_A","SO42:AI"), 
      consmodel = "(1|Site)")



# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(BB ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("SO42","PO43","Clay","Peak_A","SO42:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(BB.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("SO42","PO43","Clay","PeakA","SO42:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


BB <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle("Bacteria Biomass") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)+
  theme(axis.title.y =element_blank())

BB






# FB ####

# Full model interactions ####

FB.full.I = lmer(FB ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                  L_TN+BIX+Soil_Temp+
                                  Water_content+pH+PO43+Litter+
                                  SR+E2.E3+FI+HIX+Peak_A+
                                  Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(FB.full.I, tol = 1e-4)


summary(FB.full.I)
Anova(FB.full.I)
cAIC(FB.full.I)
AIC(FB.full.I)
r.squaredGLMM(FB.full.I)

qqnorm(residuals(FB.full.I))
scatter.smooth(residuals(FB.full.I) ~ fitted(FB.full.I))

# Full model II interactions ####

FB.full.II = lmer(FB ~ (TOC+Silt+Clay+NH4+SO42+PO43+Water_content+SR+FI+Litter+
                          Peak_A)*AI + (1|Site), data = my_data.1)                                         
isSingular(FB.full.II, tol = 1e-4)


summary(FB.full.II)
Anova(FB.full.II)
cAIC(FB.full.II)
AIC(FB.full.II)
r.squaredGLMM(FB.full.II)

qqnorm(residuals(FB.full.II))
scatter.smooth(residuals(FB.full.II) ~ fitted(FB.full.II))

# Full model III interactions ####

FB.full.III = lmer(FB ~ TOC+Clay+NH4+Water_content+SR+Litter+Peak_A+FI+
                     AI+TOC:AI+SO42:AI+PO43:AI+Peak_A:AI+
                     (1|Site), data = my_data.1)                                         
isSingular(FB.full.III, tol = 1e-4)


summary(FB.full.III)
Anova(FB.full.III)
cAIC(FB.full.III)
AIC(FB.full.III)
r.squaredGLMM(FB.full.III)

qqnorm(residuals(FB.full.III))
scatter.smooth(residuals(FB.full.III) ~ fitted(FB.full.III))

# Full model 4 interactions ####

FB.full.4 = lmer(FB ~ TOC+Clay+Water_content+SR+Peak_A+
                     AI+(1|Site), data = my_data.1)                                         
isSingular(FB.full.4, tol = 1e-4)


summary(FB.full.4)
Anova(FB.full.4)
cAIC(FB.full.4)
AIC(FB.full.4)
r.squaredGLMM(FB.full.4)

qqnorm(residuals(FB.full.4))
scatter.smooth(residuals(FB.full.4) ~ fitted(FB.full.4))

# Importance assessment ####

domin(FB ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("TOC","Clay","Water_content","SR","Peak_A","AI"), 
      consmodel = "(1|Site)")


# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(FB ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("TOC","Clay","Water_content","SR","Peak_A","AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(FB.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("TOC","Clay","WC","SR","PeakA","AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


FB <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle("Fungal Biomass") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)+
  theme(axis.title.y =element_blank())

FB








# MB ####

# Full model interactions ####

MB.full.I = lmer(MB ~ altitude+(TOC+Silt+Clay+L_TC+NH4+SO42+
                                  L_TN+BIX+Soil_Temp+
                                  Water_content+pH+PO43+Litter+
                                  SR+E2.E3+FI+HIX+Peak_A+
                                  Peak_T)*AI + (1|Site), data = my_data.1)                                         
isSingular(MB.full.I, tol = 1e-4)


summary(MB.full.I)
Anova(MB.full.I)
cAIC(MB.full.I)
AIC(MB.full.I)
r.squaredGLMM(MB.full.I)

qqnorm(residuals(MB.full.I))
scatter.smooth(residuals(MB.full.I) ~ fitted(MB.full.I))

# Full model II interactions ####

MB.full.II = lmer(MB ~ (TOC+Clay+NH4+SO42+Water_content+PO43+Litter+SR+FI+
                          Peak_A)*AI + (1|Site), data = my_data.1)                                         
isSingular(MB.full.II, tol = 1e-4)


summary(MB.full.II)
Anova(MB.full.II)
cAIC(MB.full.II)
AIC(MB.full.II)
r.squaredGLMM(MB.full.II)

qqnorm(residuals(MB.full.II))
scatter.smooth(residuals(MB.full.II) ~ fitted(MB.full.II))

# Full model III interactions ####

MB.full.III = lmer(MB ~ TOC+Clay+NH4+Water_content+PO43+SR+FI+Peak_A+AI+
                     TOC:AI+SO42:AI+PO43:AI+SR:AI+Peak_A:AI+(1|Site), data = my_data.1)                                         
isSingular(MB.full.III, tol = 1e-4)


summary(MB.full.III)
Anova(MB.full.III)
cAIC(MB.full.III)
AIC(MB.full.III)
r.squaredGLMM(MB.full.III)

qqnorm(residuals(MB.full.III))
scatter.smooth(residuals(MB.full.III) ~ fitted(MB.full.III))

# Full model 4 interactions ####

MB.full.4 = lmer(MB ~ TOC+Clay+Water_content+SR+Peak_A+AI+
                     (1|Site), data = my_data.1)                                         
isSingular(MB.full.4, tol = 1e-4)


summary(MB.full.4)
Anova(MB.full.4)
cAIC(MB.full.4)
AIC(MB.full.4)
r.squaredGLMM(MB.full.4)

qqnorm(residuals(MB.full.4))
scatter.smooth(residuals(MB.full.4) ~ fitted(MB.full.4))

# Importance assessment ####

domin(MB ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("TOC","Clay","Water_content","SR","AI","Peak_A"), 
      consmodel = "(1|Site)")




# Parameter ranking plot ----
# # Change the signs of the variables
# dominance_output <- domin(MB ~ 1, 
#                           lmer, 
#                           list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
#                           data = my_data.1, 
#                           sets = list("TOC","Clay","Water_content","SR","AI","Peak_A"), 
#                           consmodel = "(1|Site)") # Replace with your actual function
# 
# # Extracting General Dominance Standardized Ranks
# general_dominance <- dominance_output$Standardized
# general_dominance_ranks <- dominance_output$Ranks
# 
# dominance_data <- data.frame(Standardized = general_dominance)
# dominance_data$Ranks <- general_dominance_ranks
# dominance_data <- dominance_data[order(dominance_data$Ranks), ]
# 
# # Extracting estimates from the summary output
# estimates <- summary(MB.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept
# 
# # Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
# aligned_dominance <- sign(estimates) * abs(dominance_data)
# 
# components <- data.frame(Variables = c("TOC","Clay","Water_content","SR","AI","Peak_A"))
# 
# aligned_dominance <- cbind(aligned_dominance, components)
# 
# 
# # Reorder the levels of the 'Variables' column based on 'Standardized' values
# aligned_dominance <- aligned_dominance %>% 
#   arrange(Standardized)  # Arrange in descending order of 'Standardized'
# 
# # Convert 'Variables' to a factor with levels based on the ordered 'Variables'
# aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
#                                       levels = aligned_dominance$Variables)
# 
# # Plotting
# ggplot(aligned_dominance, aes(x = Standardized, y = Variables)) +
#   geom_bar(stat = "identity", fill = "#5F8249", color = "black") +
#   labs(
#     y = "Predictor Variables",
#     x = "Standardized dominance"
#   ) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
#             position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
#             color = "black", hjust = ifelse(aligned_dominance$Standardized >= 0, -0.3, 1.2)) +  # Adjust label position
#   ggtitle("MB model") +
#   theme(legend.position = "none")  # Remove the legend






#Change the signs of the variables
dominance_output <- domin(MB ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("TOC","Clay","Water_content","SR","AI","Peak_A"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(MB.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("TOC","Clay","WC","SR","AI","PeakA"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)

# Plotting
# ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals))) +
#   geom_bar(stat = "identity", fill = ifelse(aligned_dominance$Standardized >= 0, "#5F8249", "#CE5C17"),
#            color = "black") +
#   labs(
#     y = "Predictor Variables",
#     x = "Standardized dominance"
#   ) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
#             position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
#             color = "black", hjust = -0.3) +  # Adjust label position
#   ggtitle("MB model") +
#   theme(legend.position = "bottom") 


MB <- ggplot(aligned_dominance, aes(x = abs_vals, y = reorder(Variables, abs_vals), fill = factor(Standardized >= 0))) +
  geom_bar(stat = "identity",
           color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), -round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = -0.3) +  # Adjust label position
  ggtitle("Microbial Biomass") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values = c("#CE5C17", "#5F8249"), name = "Standardized",
                    labels = c("Negative", "Positive"))+
  theme(axis.title.x =element_blank())+
  xlim(0, 1)+
  theme(axis.title.y =element_blank())

MB


#ALL parameter ranking together ----


all <- ggarrange(Respi, alpha, beta, xyl, cbh, gla, fos, leu, phe, BB, FB,
                 ncol = 4, nrow = 3,
                 common.legend = TRUE,
                 legend = "right")
all

ggsave(path = "Figures/1 GRADIENT", "model_parameters.png", width = 20, height = 12, dpi = 300)

all_ag <- ggarrange(Respi, gla, fos, leu, phe, MB,
                 ncol = 3, nrow = 2,
                 common.legend = TRUE,
                 legend = "right")

all_ag
# IF I USE THIS ONE --> CHECK TITLES FOR OY !!!!!!!!!!!!!!!!!!!!!!!!!!!

ggsave(path = "Figures/1 GRADIENT", "model_parameters_ag.png", width = 20, height = 12, dpi = 300)



#. --------------------------------
# TRYING THINGS -----
## Parameter plots with negative values -----

# Respiration ----
#Change the signs of the variables
dominance_output <- domin(Respiration ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("altitude","Water_content","Silt","Silt:AI","Clay","Clay:AI","SR:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(respiration.full.III)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("altitude","WC","Silt","Silt:AI","Clay","Clay:AI","SR:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)

Respi2 <- ggplot(aligned_dominance, aes(x = Standardized, y = Variables)) +
  geom_bar(stat = "identity", fill = "#5F8249", color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = ifelse(aligned_dominance$Standardized >= 0, -0.3, 1.2)) +  # Adjust label position
  ggtitle("Respiration") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  theme(axis.title.x =element_blank())+
  xlim(-1, 1)+
  theme(axis.title.y =element_blank())

Respi2


# gla ----
#Change the signs of the variables
dominance_output <- domin(gla ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("Water_content","Water_content:AI","SR:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(gla.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("WC","WC:AI","SR:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)

gla2 <- ggplot(aligned_dominance, aes(x = Standardized, y = Variables)) +
  geom_bar(stat = "identity", fill = "#5F8249", color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = ifelse(aligned_dominance$Standardized >= 0, -0.3, 1.2)) +  # Adjust label position
  ggtitle(expression(beta-glucosaminidase)) +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  theme(axis.title.x =element_blank())+
  xlim(-1, 1)+
  theme(axis.title.y =element_blank())
gla2



# fos ----
#Change the signs of the variables
dominance_output <- domin(fos ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("BIX","Water_content","Peak_A","Soil_Temp:AI","pH:AI","FI:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(fos.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("BIX","WC","PeakA","STemp:AI","pH:AI","FI:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


fos2 <- ggplot(aligned_dominance, aes(x = Standardized, y = Variables)) +
  geom_bar(stat = "identity", fill = "#5F8249", color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = ifelse(aligned_dominance$Standardized >= 0, -0.3, 1.2)) +  # Adjust label position
  ggtitle("Phosphatase") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  theme(axis.title.x =element_blank())+
  xlim(-1, 1)+
  theme(axis.title.y =element_blank())




# leu ----
#Change the signs of the variables
dominance_output <- domin(leu ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("TOC","Water_content","SR","Peak_A"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(leu.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("TOC","WC","SR","PeakA"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


leu2 <- ggplot(aligned_dominance, aes(x = Standardized, y = Variables)) +
  geom_bar(stat = "identity", fill = "#5F8249", color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = ifelse(aligned_dominance$Standardized >= 0, -0.3, 1.2)) +  # Adjust label position
  ggtitle("Leu-aminopeptidase") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  theme(axis.title.x =element_blank())+
  xlim(-1, 1)+
  theme(axis.title.y =element_blank())






# Phe ----
#Change the signs of the variables
dominance_output <- domin(phe ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("L_TC","pH","AI","L_TC:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(phe.full.III)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("LTC","pH","AI","LTC:AI"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>% 
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables, 
                                      levels = aligned_dominance$Variables)

aligned_dominance$abs_vals <- abs(aligned_dominance$Standardized)


phe2 <- ggplot(aligned_dominance, aes(x = Standardized, y = Variables)) +
  geom_bar(stat = "identity", fill = "#5F8249", color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = ifelse(aligned_dominance$Standardized >= 0, -0.3, 1.2)) +  # Adjust label position
  ggtitle("Phenol oxidase") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  theme(axis.title.x =element_blank())+
  xlim(-1, 1)+
  theme(axis.title.y =element_blank())



# MB ----
# Change the signs of the variables
dominance_output <- domin(MB ~ 1,
                          lmer,
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"),
                          data = my_data.1,
                          sets = list("TOC","Clay","Water_content","SR","AI","Peak_A"),
                          consmodel = "(1|Site)") # Replace with your actual function

# Extracting General Dominance Standardized Ranks
general_dominance <- dominance_output$Standardized
general_dominance_ranks <- dominance_output$Ranks

dominance_data <- data.frame(Standardized = general_dominance)
dominance_data$Ranks <- general_dominance_ranks
dominance_data <- dominance_data[order(dominance_data$Ranks), ]

# Extracting estimates from the summary output
estimates <- summary(MB.full.4)$coefficients[, "Estimate"][-1]  # Exclude intercept

# Automatically aligning signs of General Dominance Standardized Ranks with Estimate values
aligned_dominance <- sign(estimates) * abs(dominance_data)

components <- data.frame(Variables = c("TOC","Clay","WC","SR","AI","PeakA"))

aligned_dominance <- cbind(aligned_dominance, components)


# Reorder the levels of the 'Variables' column based on 'Standardized' values
aligned_dominance <- aligned_dominance %>%
  arrange(Standardized)  # Arrange in descending order of 'Standardized'

# Convert 'Variables' to a factor with levels based on the ordered 'Variables'
aligned_dominance$Variables <- factor(aligned_dominance$Variables,
                                      levels = aligned_dominance$Variables)

# Plotting
MB2 <- ggplot(aligned_dominance, aes(x = Standardized, y = Variables)) +
  geom_bar(stat = "identity", fill = "#5F8249", color = "black") +
  labs(
    y = "Predictor Variables",
    x = "Standardized dominance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Standardized >= 0, round(Standardized, 2), round(Standardized, 2))),
            position = position_dodge(width = 0.9), vjust = ifelse(aligned_dominance$Standardized >= 0, 0, 0),
            color = "black", hjust = ifelse(aligned_dominance$Standardized >= 0, -0.3, 1.2)) +  # Adjust label position
  ggtitle("Microbial Biomass") +
  theme(legend.position = "none")+
  theme(legend.title=element_blank())+
  theme(axis.title.x =element_blank())+
  xlim(-1, 1)+
  theme(axis.title.y =element_blank())
MB2

#All parameter ranking together ----
all_ag2 <- ggarrange(Respi2, gla2, fos2, leu2, phe2, MB2,
                    ncol = 3, nrow = 2,
                    common.legend = TRUE,
                    legend = "right")

all_ag2
# IF I USE THIS ONE --> CHECK TITLES FOR OY !!!!!!!!!!!!!!!!!!!!!!!!!!!


ggsave(path = "Figures/1 GRADIENT", "model_parameters2.png", width = 20, height = 12, dpi = 300)



