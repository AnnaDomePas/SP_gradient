
# Reference paper: https://bsssjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/ejss.13419

# New Liner mixed model ####
library(lme4)
library(car)
library(MuMIn)

my_data         = read.csv("SP_metadata_2021.csv", sep=";")

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

my_data          = my_data %>% mutate(C_N.1 = C_N/max(C_N)) %>% mutate(NH4.1 = NH4/max(NH4)) %>%
  mutate(PO43.1 = PO43/max(PO43)) %>% mutate(SO42.1 = SO42/max(SO42)) %>% 
  mutate(Litter.1 = Litter/max(Litter)) %>% mutate(HIX.1 = HIX/max(HIX)) %>% 
  mutate(E2.E3.1 = E2.E3/max(E2.E3)) %>% mutate(X16S.1 = X16S/max(X16S)) %>% 
  mutate(ITS2.1 = ITS2/max(ITS2)) %>% mutate(altitude.1 = altitude/max(altitude)) %>%
  mutate(Water_content.1 = Water_content/max(Water_content)) %>% mutate(Soil_Temp.1 = Soil_Temp/max(Soil_Temp)) %>%
  mutate(pH.1 = pH/max(pH))
a                = as.data.frame(colnames(my_data))

# Respiration ####

# Full model ####

respiration.full = lmer(Respiration ~ (AI+altitude.1+Soil_Temp.1+Water_content.1+
                                         pH.1+Silt+Clay+C_N.1+TOC+TC+TN+PO43.1+
                                         Litter.1+L_TC+L_TN+BB+FB+ShannonEEA+
                                         SR+E2.E3.1+Peak_A+HIX.1+FI) + (1|Site), data = my_data)
summary(respiration.full)
Anova(respiration.full)
AIC(respiration.full)
r.squaredGLMM(respiration.full)

qqnorm(residuals(respiration.full))
scatter.smooth(residuals(respiration.full) ~ fitted(respiration.full))

# Full model interactions ####

respiration.full.I = lmer(Respiration ~ (altitude.1+Soil_Temp.1+Water_content.1+
                                           pH.1+Silt+Clay+C_N.1+TOC+TC+TN+PO43.1+
                                           Litter.1+L_TC+L_TN+BB+FB+ShannonEEA+
                                           SR+E2.E3.1+Peak_A+HIX.1+FI)*AI + (1|Site), data = my_data)                                         
   
summary(respiration.full.I)
Anova(respiration.full.I)
AIC(respiration.full.I)
r.squaredGLMM(respiration.full.I)

qqnorm(residuals(respiration.full.I))
scatter.smooth(residuals(respiration.full.I) ~ fitted(respiration.full.I))

# Full model A  interactions ####

respiration.full.Ia = lmer(Respiration ~ (altitude.1+Soil_Temp.1+Water_content.1+
                                           pH.1+Clay+C_N.1+PO43.1+
                                           Litter.1+BB+FB+
                                           E2.E3.1)*AI + (1|Site), data = my_data)
summary(respiration.full.Ia)
Anova(respiration.full.Ia)
AIC(respiration.full.Ia)
r.squaredGLMM(respiration.full.Ia)

qqnorm(residuals(respiration.full.Ia))
scatter.smooth(residuals(respiration.full.Ia) ~ fitted(respiration.full.Ia))

# Simplified 1 interactions ####

respiration.1.I = lmer(Respiration ~ (Clay+Water_content.1+pH.1+E2.E3.1)*AI + (1|Site), data = my_data)
summary(respiration.1.I)
Anova(respiration.1.I)
AIC(respiration.1.I)
r.squaredGLMM(respiration.1.I)

qqnorm(residuals(respiration.1.I))
scatter.smooth(residuals(respiration.1.I) ~ fitted(respiration.1.I))

# Simplified 2 interactions ####

respiration.2.I = lmer(Respiration ~ (Clay+Water_content.1+pH.1)*AI + (1|Site), data = my_data)
summary(respiration.2.I)
Anova(respiration.2.I)
AIC(respiration.2.I)
r.squaredGLMM(respiration.2.I)

qqnorm(residuals(respiration.2.I))
scatter.smooth(residuals(respiration.2.I) ~ fitted(respiration.2.I))

# Simplified 3 interactions ####

respiration.3.I = lmer(Respiration ~ (Clay+Water_content.1)*AI + (1|Site), data = my_data)
summary(respiration.3.I)
Anova(respiration.3.I)
AIC(respiration.3.I)
r.squaredGLMM(respiration.3.I)

qqnorm(residuals(respiration.3.I))
scatter.smooth(residuals(respiration.3.I) ~ fitted(respiration.3.I))

# Simplified 4 interactions ####

respiration.4.I = lmer(Respiration ~ (Clay)*AI + (1|Site), data = my_data)
summary(respiration.4.I)
Anova(respiration.4.I)
AIC(respiration.4.I)
r.squaredGLMM(respiration.4.I)

qqnorm(residuals(respiration.4.I))
scatter.smooth(residuals(respiration.4.I) ~ fitted(respiration.4.I))

# New Linear Model ####

my_data.1 = my_data %>% group_by(Site) %>% summarise(across(everything(), list(mean)))

# Full model linear interactions ####

respiration.full.linear.I = lm(Respiration_1 ~ altitude.1_1:AI_1+Soil_Temp.1_1:AI_1+
                                 Water_content.1_1:AI_1+pH.1_1:AI_1+Clay_1:AI_1+
                                 C_N.1_1:AI_1+PO43.1_1:AI_1+BB_1:AI_1+FB_1:AI_1+
                                 E2.E3.1_1:AI_1, data = my_data.1)  
summary(respiration.full.linear.I)
Anova(respiration.full.linear.I)
qqnorm(residuals(respiration.full.linear.I))
scatter.smooth(residuals(respiration.full.linear.I) ~ fitted(respiration.full.linear.I))

# Full model A  linear interactions ####

respiration.full.linear.Ia = lm(Respiration_1 ~ Soil_Temp.1_1:AI_1+
                                  Water_content.1_1:AI_1+
                                  C_N.1_1:AI_1+PO43.1_1:AI_1+BB_1:AI_1+FB_1:AI_1+
                                  E2.E3.1_1:AI_1, data = my_data.1)
summary(respiration.full.linear.Ia)
Anova(respiration.full.linear.Ia)
qqnorm(residuals(respiration.full.linear.Ia))
scatter.smooth(residuals(respiration.full.linear.Ia) ~ fitted(respiration.full.linear.Ia))
