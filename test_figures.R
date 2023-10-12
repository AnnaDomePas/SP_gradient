a = as.data.frame(colnames(my_data))
library(car)

# Soil characterization along the gradient ####

# Soil Temperature
test = aov(Soil_Temp ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(Soil_Temp ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "Soil_Temp",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# Water content
test = aov(Water_content ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(Water_content ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "Water_content",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# Soil organic matter (%)
test = aov(SOM_perc ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(SOM_perc ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "SOM_perc",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# pH (%)
test = aov(pH ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(pH ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "pH",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# TOC_perc (%)
test = aov(TOC_perc ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(TOC_perc ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "TOC_perc",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# TN_perc (%)
test = aov(TN_perc ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(TN_perc ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "TN_perc",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# C_N
test = aov(C_N ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(C_N ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "C_N",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# NH4
test = aov(NH4 ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(NH4 ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "NH4",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# PO43
test = aov(PO43 ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(PO43 ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "PO43",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# SO42
test = aov(SO42 ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(SO42 ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "SO42",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# Litter
test = aov(Litter ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(Litter ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "Litter",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# Litter
test = aov(Litter ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(Litter ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "Litter",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# L_TC_perc
test = aov(L_TC_perc ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(L_TC_perc ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "L_TC_perc",
          fill = "AI") + stat_compare_means(method = "anova") + 
  scale_fill_gradient(low="blue", high="red")

# L_TN_perc
test = aov(L_TN_perc ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(L_TN_perc ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "L_TN_perc",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# CO2_dark
test = aov(CO2_dark ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(CO2_dark ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "CO2_dark",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# CH4_ave
test = aov(CH4_ave ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(CH4_ave ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "CH4_ave",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# N2O
test = aov(N2O ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(N2O ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "N2O",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# BB_2
test = aov(BB_2 ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(BB_2 ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "BB_2",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# BB_2
test = aov(BB_2 ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(BB_2 ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "BB_2",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# Respiration
test = aov(Respiration ~ Site,data = my_data)
shapiro.test(test$residuals)
bartlett.test(Respiration ~ Site, data = my_data)
par(mfrow = c(2, 2))
plot(test)
par(mfrow = c(1, 1))
summary(test)
ggboxplot(my_data, x = "Site", y = "Respiration",
          fill = "AI") + stat_compare_means() + 
  scale_fill_gradient(low="blue", high="red")

# Linear model between soil aridity and physicochemical variables ####

model_1 = lm(AI~Water_content,data=my_data)
model_2 = lm(AI~Water_content+Soil_Temp,data=my_data) 
model_3 = lm(AI~Water_content+Soil_Temp+SOM_perc,data=my_data)
model_4 = lm(AI~Water_content+Soil_Temp+SOM_perc+pH,data=my_data)
model_5 = lm(AI~Water_content+Soil_Temp+SOM_perc+pH+Litter,data=my_data)
model_6 = lm(AI~Water_content+Soil_Temp+SOM_perc+pH+Litter+PO43,data=my_data)
model_7 = lm(AI~Water_content+SOM_perc+pH+PO43,data=my_data)

anova(model_1,model_2,model_3,model_4,model_5,model_6,model_7)
AIC(model_1,model_2,model_3,model_4,model_5,model_6,model_7)

shapiro.test(model_7$residuals)
plot(model_7, which = 1)
plot(model_7, which = 2)
summary(model_7)

# PCA ####

results          = prcomp(my_data[,c(8,9,10,11,13,14,16,20,21,22,23,24,25,26,30,31,32)],
                  scale = TRUE)
results$rotation = -1*results$rotation
results$rotation
biplot(results, scale = 0)

results$x        = -1*results$x
results$x

results$sdev^2 / sum(results$sdev^2)

# Manova test

y_variables = as.data.frame((my_data[,c(35,42:53,56:64)]))
corre       = as.data.frame(cor(y_variables))
mod         = manova(as.matrix(y_variables) ~ as.factor(my_data$AI))
summary(mod,summary=TRUE)
summary.aov(mod)

mod         = manova(as.matrix(y_variables) ~ as.factor(my_data$AI),
                     subset(as.factor(my_data$AI) %in% c(1,2)))
