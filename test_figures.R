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

# PCA ####

results          = prcomp(my_data[,c(8,9,10,11,13,14,16,20,21,22,23,24,25,26,30,31,32)],
                  scale = TRUE)
results$rotation = -1*results$rotation
results$rotation
biplot(results, scale = 0)

results$x        = -1*results$x
results$x

results$sdev^2 / sum(results$sdev^2)

# Manova ####

# Manova Asumptions ####

y_variables = as.data.frame((my_data[,c(2,35,42:53,56:64)]))

library(rstatix)
library(broom)

# Outliers ####
out = as.data.frame(y_variables %>% group_by(Site) %>% identify_outliers("CO2_dark"))
# There are outliers but we are keeping them for the analysis

# Multivariate outliers ####
y_variables_2       = scale(y_variables[,c(2:6,8,9,11:13,15:23)],center = FALSE)
mahal               = mahalanobis(y_variables_2, colMeans(y_variables_2), cov(y_variables_2))
p_val               = pchisq(mahal, df=3, lower.tail=FALSE)
y_variables_2       = as.data.frame(cbind(y_variables_2,mahal,p_val))
# There are outliers but we are keeping them for the analysis

# Normality ####
y_variables_2      = as.data.frame(cbind(y_variables$Site,y_variables_2))
ggqqplot(y_variables_2, "CO2_dark", facet.by = "y_variables$Site",
         ylab = "CO2_dark", ggtheme = theme_bw())
ggqqplot(y_variables_2, "chla", facet.by = "y_variables$Site",
         ylab = "chla", ggtheme = theme_bw())
ggqqplot(y_variables_2, "chlb", facet.by = "y_variables$Site",
         ylab = "chlb", ggtheme = theme_bw())
ggqqplot(y_variables_2, "carotene", facet.by = "y_variables$Site",
         ylab = "carotene", ggtheme = theme_bw())
ggqqplot(y_variables_2, "EPS", facet.by = "y_variables$Site",
         ylab = "EPS", ggtheme = theme_bw())
ggqqplot(y_variables, "alpha", facet.by = "Site",
         ylab = "alpha", ggtheme = theme_bw())
ggqqplot(y_variables_2, "beta", facet.by = "y_variables$Site",
         ylab = "beta", ggtheme = theme_bw())
ggqqplot(y_variables_2, "xyl", facet.by = "y_variables$Site",
         ylab = "xyl", ggtheme = theme_bw())
ggqqplot(y_variables, "cbh", facet.by = "Site",
         ylab = "cbh", ggtheme = theme_bw())
ggqqplot(y_variables_2, "gla", facet.by = "y_variables$Site",
         ylab = "gla", ggtheme = theme_bw())
ggqqplot(y_variables_2, "fos", facet.by = "y_variables$Site",
         ylab = "fos", ggtheme = theme_bw())
ggqqplot(y_variables_2, "leu", facet.by = "y_variables$Site",
         ylab = "leu", ggtheme = theme_bw())
ggqqplot(y_variables, "phe", facet.by = "Site",
         ylab = "phe", ggtheme = theme_bw())
ggqqplot(y_variables_2, "mcrA", facet.by = "y_variables$Site",
         ylab = "mcrA", ggtheme = theme_bw())
ggqqplot(y_variables_2, "pmoA", facet.by = "y_variables$Site",
         ylab = "pmoA", ggtheme = theme_bw())
ggqqplot(y_variables_2, "nifH", facet.by = "y_variables$Site",
         ylab = "nifH", ggtheme = theme_bw())
ggqqplot(y_variables_2, "AOA", facet.by = "y_variables$Site",
         ylab = "AOA", ggtheme = theme_bw())
ggqqplot(y_variables_2, "AOB", facet.by = "y_variables$Site",
         ylab = "AOB", ggtheme = theme_bw())
ggqqplot(y_variables_2, "qnorB", facet.by = "y_variables$Site",
         ylab = "qnorB", ggtheme = theme_bw())
ggqqplot(y_variables_2, "nosZ", facet.by = "y_variables$Site",
         ylab = "nosZ", ggtheme = theme_bw())
ggqqplot(y_variables_2, "phoD", facet.by = "y_variables$Site",
         ylab = "phoD", ggtheme = theme_bw())
ggqqplot(y_variables_2, "Respiration", facet.by = "y_variables$Site",
         ylab = "Respiration", ggtheme = theme_bw())
# I think most of the variables are normally distributed and MANOVA is robust with
# slightly violations of normality. 

# Multivariate normality ####
mshapiro_test((y_variables_2[,2:22]))
# Not really normal

# Multicollinearity ####
cor.mat <- y_variables[,2:23] %>% cor_mat()
cor.mat %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot(label = TRUE)
# New selected variables
y_variables_3 = as.data.frame((y_variables[,c(1:3,5:20,22,23)]))

# Linearity assumption ####
library(GGally)
linear_ass <- y_variables_3 %>% group_by(Site) %>% doo(~ggpairs(.) + theme_bw(), result = "plots")
linear_ass$plots

# Homogeneity of covariances ####
# Balance design and we could continue with the analysis. We will use the 
# Pillai’s multivariate statistic as a more robust metric

# Homogeneity of variance assumption ####
a = as.data.frame(colnames(y_variables_3))
y_variables_3 %>% gather(key = "variable", value = "value", a[2:21,]) %>%
  group_by(variable) %>% levene_test(sqrt(value) ~ as.factor(Site))
# Squaring the variables solves almost completely the homogeneity of variances

# Manova test ####
mod         = manova(as.matrix(y_variables_3[,2:21]) ~ (my_data$AI))
summary(mod,summary=TRUE)

# Univariate one-way ANOVA ####
grouped.data <- y_variables_3 %>% gather(key = "variable", value = "value", a[2:21,]) %>%
  group_by(variable)
grouped.data %>% kruskal_test(value ~ my_data$AI)

# What variables to use in the MANOVA as covariables
x_variables = as.data.frame((my_data[,c(2,5:8,10,11,13:16,20:26,30:32)]))

cor.mat <- x_variables[,2:21] %>% cor_mat()
cor.mat %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot(label = TRUE)
# New selected variables
x_variables_1 = as.data.frame((x_variables[,c(1,5:7,10,12:15,19:21)]))

# VIF ####
library(car)
cor.mat <- x_variables_1[,2:12] %>% cor_mat()
model <- lm(cbind(y_variables_3$CO2_dark+y_variables_3$chla+y_variables_3$carotene+
                    y_variables_3$EPS) 
            ~ x_variables_1$Soil_Temp+
              x_variables_1$Water_content+x_variables_1$SOM_perc+x_variables_1$TC_perc+
              x_variables_1$C_N+x_variables_1$NH4+x_variables_1$PO43+x_variables_1$SO42+
              x_variables_1$Litter+x_variables_1$L_TC_perc+x_variables_1$L_TN_perc+as.numeric(my_data$AI))
car::vif(model)

library(pls)
model <- plsr(y_variables_3$CO2_dark~x_variables_1$Soil_Temp+
                x_variables_1$Water_content+x_variables_1$SOM_perc+x_variables_1$TC_perc+
                x_variables_1$C_N+x_variables_1$NH4+x_variables_1$PO43+x_variables_1$SO42+
                x_variables_1$Litter+x_variables_1$L_TC_perc+x_variables_1$L_TN_perc, scale=TRUE, validation="CV")
summary(model)
coefficients = coef(model)

cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
R2(model)

# Reference paper: https://bsssjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/ejss.13419