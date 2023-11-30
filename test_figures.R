
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
# We are normalizing all the variables to use the coefficients to assess the
# importance of the variables for the model

a               = as.data.frame(colnames(my_data))
max_values      = apply(my_data[,c(8,10,13,20,22,30,76,81,82,83,85,92)],2,max)
my_data.1       = as.data.frame(cbind(my_data[,c(2,64,46,47,48,49,7,9,17,18,19,27,28,29,33,34,40,41,86,87,88,89,90,91,74)],
                                      my_data[,c(8,10,13,20,22,30,76,81,82,83,85,92)] / as.list(max_values)))
a.1             = as.data.frame(colnames(my_data.1))

corr            = as.data.frame(cor(my_data.1[,8:37]))

# Parameters correlated: Water activity-water content; TOC-TN,FB,BB,TC;
# Sand-silt,clay; Peak_A-Peak_M, Peak_C; Peak_T-Peak_B; E2.E3-E3.E4

# Respiration ####

# Full model interactions ####

respiration.full.I = lmer(Respiration ~ altitude+(TOC+ShannonEEA+Silt+Clay+L_TC+
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

library(domir)
library(MuMIn)

domin(Respiration ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("altitude","Water_content","Silt","Silt:AI","Clay","Clay:AI","SR:AI"), 
      consmodel = "(1|Site)")

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

# alpha

# Full model interactions ####

alpha.full.I = lmer(alpha ~ altitude+(TOC+ShannonEEA+Silt+Clay+L_TC+
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

alpha.full.II = lmer(alpha ~ (TOC+Clay+L_TC+L_TN+BB+BIX+Soil_Temp+Water_content+pH+PO43+
                                        SR+FI)*AI + (1|Site), data = my_data.1)                                         

summary(alpha.full.II)
Anova(alpha.full.II)
cAIC(alpha.full.II)
AIC(alpha.full.II)
r.squaredGLMM(alpha.full.II)

qqnorm(residuals(alpha.full.II))
scatter.smooth(residuals(alpha.full.II) ~ fitted(alpha.full.II))

# Full model interactions III ####

alpha.full.III = lmer(alpha ~ L_TC+L_TC:AI+L_TN:AI+BB+BB:AI+FI+AI + 
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

# Full model interactions 4 ####

alpha.full.4 = lmer(alpha ~ L_TC+L_TC:AI+BB+BB:AI+(1|Site), data = my_data.1) 

summary(alpha.full.4)
Anova(alpha.full.4)
cAIC(alpha.full.4)
AIC(alpha.full.4)
r.squaredGLMM(alpha.full.4)

qqnorm(residuals(alpha.full.4))
scatter.smooth(residuals(alpha.full.4) ~ fitted(alpha.full.4))

domin(alpha ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("L_TC","L_TC:AI","BB","BB:AI"), 
      consmodel = "(1|Site)")
