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


## Working directory
# setwd("G:/GRADCATCH/ANALISIS/R")
# getwd()
# version


# IMPORT DATA ----
my_data <- read.csv("SP_metadata_2021.csv", sep=",")

#To replace NA values with a mean of the other values of the Site:
for (i in which(sapply(my_data, is.numeric))) {
  for (j in which(is.na(my_data[, i]))) {
    my_data[j, i] <- mean(my_data[my_data[, "Site"] == my_data[j, "Site"], i],  na.rm = TRUE)
  }
}


# IMPORT FUNCTIONS ----
# > DATA SUMMARY (SD) ----
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



#.----

# UCI STAY ANALYSES ----


# One-way MANOVA ----

#Condition for MANOVA: the dataset must have more
# observations (rows) per group in the independent
# variable than a number of the dependent variables.

manova_data <- my_data[,c(2,35,37,38,40:64)]

# 1rst we separate the dependent variables from the other ones:
dependent_vars <- as.matrix(manova_data[,-c(1)])
# dependent variables needs to be entered as a matrix, not a dataframe
independent_var <- as.factor(manova_data[,c(1)])
# independent variable needs to be a factor (with different levels)

manova_model <- manova(dependent_vars ~ independent_var, data = manova_data)
summary(manova_model)

# By default, MANOVA in R uses Pillai’s Trace test statistic.
# The P-value is practically zero, which means we can safely
# reject the null hypothesis in the favor of the alternative one:
# at least one group mean vector differs from the rest.


# we could also measure the effect size. 
# One metric often used with MANOVA is Partial Eta Squared.
# It measures the effect the independent variable has on the
#dependent variables. If the value is 0.14 or greater, we can
# say the effect size is large.

library(effectsize)

eta_squared(manova_model)

# The value is 0.83, which means the effect size is large.
# It’s a great way to double-check the summary results of a MANOVA test.



# Post hoc test: Linear Discriminant Analysis (LDA)
# It finds a linear combination of features that best separates two
# or more groups.

library(MASS)

manova_lda <- lda(independent_var ~ dependent_vars, CV = F)
manova_lda


# The snippet below uses the predict() function to get the linear
# discriminants and combines them with our independent variable:
lda_df <- data.frame(
  sites = manova_data[, "Site"],
  lda = predict(manova_lda)$x
)
lda_df


# The final step in this post-hoc test is to visualize the above
#  lda_df as a scatter plot. Ideally, we should see one or multiple
#  groups stand out:


lda_df$sites <- factor(lda_df$sites, levels = c("SP08", "SP01", "SP02",
                                                "SP07","SP06","SP03",
                                                "SP12", "SP11", "SP04",
                                                "SP09", "SP10","SP05"))

mycolors2<-c("#10449F","#51B7DF","#00FFFF","#00A01D","#064700",
             "#7A7615", "#583200", "#C24A0A", "#F5C92D","#FA0C00",
             "#7D1809", "#290500")

ggplot(lda_df) +
  geom_point(aes(x = lda.LD1, y = lda.LD2, fill=sites), colour= "black", pch=21, size = 5) +
  theme_classic()+
  scale_fill_manual(values = mycolors2)





library(rstatix)
pwc <- manova_data %>%
  gather(key = "variables", value = "value", CO2_dark,CH4_ave,N2O,BB,FB,chla,chlb,carotene,
         EPS, alpha, beta, xyl, cbh, gla, fos, leu, phe,
         X16S, ITS2, mcrA, pmoA, nifH, AOA, AOB, qnorB,
         nosZ, phoD, Respiration) %>%
  group_by(variables) %>%
  games_howell_test(value ~ Site) %>%
  select(-estimate, -conf.low, -conf.high) # Remove details
pwc





# https://www.r-bloggers.com/2022/01/manova-in-r-how-to-implement-and-interpret-one-way-manova/
#  EXEMPLE WITH IRIS DATA:
# data("iris")
# dependent_vars <- cbind(iris$Sepal.Length, iris$Sepal.Width, iris$Petal.Length, iris$Petal.Width)
# independent_var <- iris$Species
# manova_model <- manova(dependent_vars ~ independent_var, data = iris)
# summary(manova_model)
# eta_squared(manova_model)
# iris_lda <- lda(independent_var ~ dependent_vars, CV = F)
# iris_lda
# lda_df <- data.frame(
#   species = iris[, "Species"],
#   lda = predict(iris_lda)$x)
# lda_df
# ggplot(lda_df) +
#   geom_point(aes(x = lda.LD1, y = lda.LD2, color = species), size = 4) +
#   theme_classic()



# > MANOVA assumptions ----
#DO WE FOLLOW ALL THE ASSUMPTIONS FOR MANOVA?

#https://www.datanovia.com/en/lessons/one-way-manova-in-r/#data-preparation

# 1. Extreme outliers:
# (trying for CO2_dark)
library(rstatix)
a <- manova_data %>%
  group_by(Site) %>%
  identify_outliers(CO2_dark)
# THERE ARE SOME EXTREME OUTLIERS!!
# Note that, in the situation where you have extreme outliers,
# this can be due to: 1) data entry errors, measurement errors or unusual values.
# You can include the outlier in the analysis anyway if you do not believe the result
# will be substantially affected. This can be evaluated by comparing the result
# of the MANOVA with and without the outlier.
# Remember to report in your written results section any decisions
# you make regarding any outliers you find.


# 2.Multivariate outliers
manova_data %>%
  group_by(Site) %>%
  mahalanobis_distance(-id) %>%
  filter(is.outlier == TRUE) %>%
  as.data.frame()
# Does not work


manova_data %>%
  group_by(Site) %>%
  shapiro_test(CO2_dark) %>%
  arrange(variable)
# 1 site (SP11) does not have normality
# for the CO2_dark

apply(dependent_vars,2,shapiro.test)
# Any variable has normality?

# WE SHOULD TRANSFORM VARIABLES


# 3. Seeing if there is multicollineality:
# pvalue > 0.05 --> there is collineality
manova_data %>% cor_test(CO2_dark,CH4_ave,N2O,BB,FB,chla,chlb,carotene,
                         EPS, alpha, beta, xyl, cbh, gla, fos, leu, phe,
                         X16S, ITS2, mcrA, pmoA, nifH, AOA, AOB, qnorB,
                         nosZ, phoD, Respiration)

# In the situation, where you have multicollinearity, you could consider
# removing one of the outcome variables that is highly correlated.



# 4. Lineality assumption
library(GGally)
results <- manova_data %>%
  select(CO2_dark,CH4_ave,N2O,BB,FB,chla,chlb,carotene,
         EPS, alpha, beta, xyl, cbh, gla, fos, leu, phe,
         X16S, ITS2, mcrA, pmoA, nifH, AOA, AOB, qnorB,
         nosZ, phoD, Respiration, Site) %>%
  group_by(Site) %>%
  doo(~ggpairs(.) + theme_bw(), result = "plots")
results$plots


# 5. Homogenity of covariances:
box_m(manova_data[, c("CO2_dark", "CH4_ave")], manova_data$Site)

# Note that, if you have balanced design (i.e., groups with similar sizes),
# you don’t need to worry too much about violation of the homogeneity
# of variances-covariance matrices and you can continue your analysis.
# However, having an unbalanced design is problematic.
# Possible solutions include:
# 1) transforming the dependent variables;
# 2) running the test anyway, but using Pillai’s multivariate
# statistic instead of Wilks’ statistic.


# 6. Homogenity of variance:
manova_data %>% 
  gather(key = "variable", value = "value", CO2_dark,CH4_ave,N2O,BB,FB,chla,chlb,carotene,
         EPS, alpha, beta, xyl, cbh, gla, fos, leu, phe,
         X16S, ITS2, mcrA, pmoA, nifH, AOA, AOB, qnorB,
         nosZ, phoD, Respiration) %>%
  group_by(variable) %>%
  levene_test(value ~ Site)


#  Note that, if you do not have homogeneity of variances,
#  you can try to transform the outcome (dependent) variable
# to correct for the unequal variances.

# Alternatively, you can continue, but accept a lower level of
# statistical significance (alpha level) for your MANOVA result.
# Additionally, any follow-up univariate ANOVAs will need to be corrected
# for this violation (i.e., you will need to use different post-hoc tests).




#.----
# Multiple lineal regressions ----
library(car)
library(MASS)

# > Respiration ----
scatterplotMatrix(~AI+Soil_Temp+MAP+MAT+Water_activity+Water_content+SOM+pH+TOC+TC+TN+C_N+Clay+Sand+Silt, data=my_data)
#no todas las variables tienen una distribucion en campana

R <- lm(Respiration ~ AI+Soil_Temp+MAP+MAT+Water_activity+Water_content+SOM+pH+TOC+TC+TN+C_N+Clay+Sand+Silt + Site, my_data)

par(mfrow=c(2,2))
plot(R)
#Seems pretty fine
par(mfrow=c(1,1))


vif(R)

#There is an error of correlated variables

respi <- my_data[,c(5:10,12,13,17:20,27:29,64)]
cor(respi)

#When two variables have a correlation coefficient of 1
# it results in an error, as the two variables are perfectly
#correlated

#VALUES >0.8 (or 0.7 if you are more strict) are considered
# highly correlated

library(GGally)
dev.off()
ggpairs(respi)

#We could exclude variables correlated for >0.7
# (i.e., MAP and MAT are 0.767 cor)

# We eliminate variable with > 0.7 cor, trying to keep the ones
# with higher cor with the dependent variable (i.e. Respiration)


my_data2 <- my_data[,c(10,13,18,20,27,64)]
ggpairs(my_data2)

R2 <- lm(Respiration ~ Water_content+pH+TC+C_N+Sand + Site, my_data)
summary(R2)
vif(R2)
# Problem is with SITE


#Find the linearly dependent variables from your model:
ld.vars <- attributes(alias(R)$Complete)$dimnames[[1]]
ld.vars

vif_values <- vif(R)

# VIF VALUES > 5 need to be removed as
# it indicates potentially severe correlation
# between a given predictor variable and other
# predictor variables in the model.
# In this case, the coefficient estimates and
# p-values in the regression output are likely unreliable. 

#create horizontal bar chart to display each VIF value
barplot(vif_values, main = "VIF Values", horiz = FALSE, col = "steelblue")
#add vertical line at 5
abline(h = 5, lwd = 3, lty = 2)

# Variables with VIF > 5 probably have high correlation
# values in the cor matrix.

respi <- my_data[,c(5:10,12,13,17:20,27:29,64)]
cor(respi)
ggpairs(respi)

# Repeating MLR without high VIF variables
R <- lm(Respiration ~ AI+Soil_Temp+MAP+MAT+Water_activity+Water_content+SOM+pH+TOC+TC+TN+C_N, my_data)
summary(R)
vif_values <- vif(R)
barplot(vif_values, main = "VIF Values", horiz = FALSE, col = "steelblue")
abline(h = 5, lwd = 3, lty = 2)

R <- lm(Respiration ~ Soil_Temp+MAT+Water_activity+Water_content+SOM+pH+TOC+TC+TN+C_N, my_data)
summary(R)
vif_values <- vif(R)
barplot(vif_values, main = "VIF Values", horiz = FALSE, col = "steelblue")
abline(h = 5, lwd = 3, lty = 2)

R <- lm(Respiration ~ Soil_Temp+MAT+Water_activity+Water_content+SOM+pH+TC+C_N, my_data)
summary(R)
vif_values <- vif(R)
barplot(vif_values, main = "VIF Values", horiz = FALSE, col = "steelblue")
abline(h = 5, lwd = 3, lty = 2)

R <- lm(Respiration ~ Soil_Temp+MAT+Water_activity+Water_content+SOM+pH+C_N, my_data)
summary(R)
vif_values <- vif(R)
barplot(vif_values, main = "VIF Values", horiz = FALSE, col = "steelblue")
abline(h = 5, lwd = 3, lty = 2)

R <- lm(Respiration ~ Water_content+SOM+pH+C_N, my_data)
summary(R)
vif_values <- vif(R)
barplot(vif_values, main = "VIF Values", horiz = FALSE, col = "steelblue")
abline(h = 5, lwd = 3, lty = 2)


summary(R)
par(mfrow=c(2,2))
plot(R)
#Regular, could be better
par(mfrow=c(1,1))


stepb <- stepAIC(R, direction="backward")

#step backwards AIC removes SOM from our model (decreasing AIC, which is better)


#Final model:
R <- lm(Respiration ~ Water_content+pH+C_N, my_data) 
summary(R)


#To know which variable has more influence
# to the dependent variable (i.e Respiration)
# the partial standarized coeficients (Beta)
# Beta = coef * SD(X)/SD(Y)

coef_n_WC <- R$coeff["Water_content"]
coef_n_pH <- R$coeff["pH"]
coef_n_C_N <- R$coeff["C_N"]

Beta_n_WC <- coef_n_WC*sd(my_data$Water_content)/sd(my_data$Respiration)
Beta_n_pH <- coef_n_pH*sd(my_data$pH)/sd(my_data$Respiration)
Beta_n_C_N <- coef_n_C_N*sd(my_data$C_N)/sd(my_data$Respiration)

Beta_n_WC
Beta_n_pH
Beta_n_C_N
#Higher beta IN ABSOLUTE VALUE = higher
# influence on dependent variable (i.e Respi)

# In this case its pH > Water content > C/N






R <- lm(Respiration ~ AI*Soil_Temp*MAP*MAT*Water_activity*Water_content*SOM*pH*TOC*TC*TN*C_N*Clay*Sand*Silt, my_data)
summary(R)

# "Coefficients not defined because of singularities"
# indicates that two or more predictor variables in the model
# have a perfect linear relationship and thus not
# every regression coefficient in the model can be estimated



# summary(R)
# AIC(R)
# mpar <- step(R)
# 
# R1 <- lm(Respiration ~ AI*Soil_Temp*MAP*MAT*Water_activity*Water_content*SOM*pH*TOC*TC*TN*C_N*Clay*Sand*Silt + Site, my_data)
# summary(R1)
# R2 <- lm(Respiration ~ AI + Site, my_data)
# summary(R2)
# R3 <- lm(Respiration ~ AI*Soil_Temp + Site, my_data)
# summary(R3)
# R4 <- lm(Respiration ~ AI*Soil_Temp*MAP + Site, my_data)
# summary(R4)
# R5 <- lm(Respiration ~ AI*Soil_Temp*MAP*MAT + Site, my_data)
# summary(R5)
# m <- step(R5)
# R6 <- lm(Respiration ~ AI*Soil_Temp*MAP*MAT*Water_activity + Site, my_data)
# summary(R6)
# R7 <- lm(Respiration ~ AI*Soil_Temp*MAP*MAT*Water_content + Site, my_data)
# summary(R7)
# #Si pongo WC y WA juntos se buguea





#.----
# PERMANOVA ----
my_data <- read.csv("SP_metadata_2021.csv", sep=",")

#To replace NA values with a mean of the other values of the Site:
for (i in which(sapply(my_data, is.numeric))) {
  for (j in which(is.na(my_data[, i]))) {
    my_data[j, i] <- mean(my_data[my_data[, "Site"] == my_data[j, "Site"], i],  na.rm = TRUE)
  }
}

perma_data <- my_data[,c(21:23,35:38,40:64,74)]
perma_data2 <- my_data[,c(2,21:23,35:38,40:64,74)]


# https://www.youtube.com/watch?v=1ETBgbXl-BM&ab_channel=RiffomonasProject
site <- data.frame(perma_data2[,c(1)]) 
site <- setNames(site, c("Site"))

# I HAVE NOT STANDARIZED THE DATA AND I SHOULD, AS
# ALL VARIABLES HAVE DIFFERENT UNITS

#without negative values because if not i can not
# obtain matrix distance with Bray-curtis
perma_data_trans <- log(perma_data+1-min(perma_data))
perma_data2_trans <- cbind(site,perma_data_trans)

summary(perma_data_trans)
# Summary shows that Bray-curtis is not a good distance
# method, as although all variables are the same order,
# enzymes for example have lost their differences between
# them

per.dist <- as.matrix(vegdist(perma_data_trans, method="bray"))

# I wanna compare with my Site variable
# In data should i put a dataframe with my Site variable,
# it is not important if there is more variables in that dataset
# it will only take the one specified after ~
adonis2(per.dist ~ Site, data = perma_data2_trans, permutations = 10000, method="bray")





## SHIT THAT DIDNT WORK OUT:

# # DAY 1
# # https://www.youtube.com/watch?v=1ETBgbXl-BM&ab_channel=RiffomonasProject
# # # CAN NOT DO BRAY CURTIS BECAUSE OF THE NEGATIVE DATA
# per.dist <- vegdist(perma_data, method="bray")
# # per.div <- adonis2(perma_data ~ ., data = perma_data2, permutations = 999, method="bray")
# # 
# # DAY 2
# # https://www.youtube.com/watch?v=xyufizOpc5I&ab_channel=RiffomonasProject
# # library(tidyverse)
# # # # shared <- as.matrix(perma_data)
# # class(shared)
# # # # library(vegan)
# # # # set.seed(19760620)
# # CAN NOT DO BRAY CURTIS BECAUSE OF THE NEGATIVE DATA
# # dist <- vegdist(shared, method = "bray")
# # nmds <- metaMDS(dist)
# # # # nmds <- metaMDS(shared, autotransform=FALSE)
# # scores(nmds) %>% 
# #   as_tibble(rownames="Site") %>% 
# #   ggplot(aes(x=NMDS1, y= NMDS2)) +
# #   geom_point()


# . ----
# PCA with all explanatory variables ----

#Selecting variables:
pca_data <- my_data[,c(2,21:23,35:38,40:64,74)]

pca_data <- pca_data %>%
  group_by(Site) %>%
  summarise_all("mean")

site_order <- my_data[,c(2,7)] 
site_order <- site_order[!duplicated(site_order), ] #Erase duplicated lines from dataframe
pca_data <- pca_data[order(site_order$AI, decreasing = T),]
pca_data$Site <- factor(pca_data$Site, levels = pca_data$Site[order(site_order$AI)])

pcr2 <- pca_data[,c(-1)]


#Select column with levels (Site)
site <- factor(pca_data$Site, levels = pca_data$Site)
site

pc <- prcomp(na.omit(pcr2), center = TRUE,
             scale. = TRUE) 

plot(pc, type = "l")
plot(pc)
summary(pc)

mycolors2<-c("#10449F","#51B7DF","#00FFFF","#00A01D","#064700",
             "#7A7615", "#583200", "#C24A0A", "#F5C92D","#FA0C00",
             "#7D1809", "#290500")

library(ggfortify)
library(ggplot2)
autoplot(pc, data=pca_data, 
         loadings = TRUE, loadings.colour = 'brown',
         loadings.label.colour='brown', loadings.label = TRUE,
         loadings.label.size = 7,
         loadings.label.repel=TRUE)+
  theme_classic()+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 5)+
  scale_fill_manual(values = mycolors2)+
  ggtitle("All response variables")+
  theme(legend.title = element_blank(),
        legend.text=element_text(size = 12),
        title = element_text(size = 15,face="bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15, face="plain"))

# ggsave(path = "Figures","PCA_response_means.png", width = 10, height = 8, dpi = 300)







# _______________________________________ ----
# .----
# PHYSICOCHEMICAL PARAMETERS ----
#.----


# Site order depending on OX variable ----
#Preparing data:
site_order <- my_data[,c(2,5:7,76)] 
site_order <- site_order[!duplicated(site_order), ] #Erase duplicated lines from dataframe

AI_ord <- site_order[order(site_order$AI, decreasing = T),]
MAT_ord <- site_order[order(site_order$MAT, decreasing = F),]
MAP_ord <- site_order[order(site_order$MAP, decreasing = T),]
ALT_ord <- site_order[order(site_order$altitude, decreasing = F),]

AI_ord$Site <- factor(AI_ord$Site, levels = AI_ord$Site[order(AI_ord$AI)])
MAT_ord$Site <- factor(MAT_ord$Site, levels = MAT_ord$Site[order(MAT_ord$MAT, decreasing = TRUE)])
MAP_ord$Site <- factor(MAP_ord$Site, levels = MAP_ord$Site[order(MAP_ord$MAP)])
ALT_ord$Site <- factor(ALT_ord$Site, levels = ALT_ord$Site[order(ALT_ord$altitude)])


## AI plot ----
Sites_by_AI <- ggplot(AI_ord, aes(x=Site, y=AI)) +
  geom_point(size= 3, color ="gold") +
  ylab("Aridity Index") +
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(size = 20, face = "bold", colour = "black")) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0,1.4))+
  scale_y_continuous(breaks = breaks_width(0.2))

Sites_by_AI


## MAT plot ----
Sites_by_MAT <- ggplot(MAT_ord, aes(x=Site, y=MAT)) +
  geom_point(size= 3, color ="tomato") +
  ylab("MAT (ºC)") +
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(size = 20, face = "bold", colour = "black")) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(9,18))+
  scale_y_continuous(breaks = breaks_width(3))

Sites_by_MAT


## MAP plot ----
Sites_by_MAP <- ggplot(MAP_ord, aes(x=Site, y=MAP)) +
  geom_point(size= 3, color ="cornflowerblue") +
  ylab("MAP (mm)") +
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(size = 20, face = "bold", colour = "black")) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(260,1400))+
  scale_y_continuous(breaks = breaks_width(350))

Sites_by_MAP

## All together ----

allplots <- ggarrange(Sites_by_AI, Sites_by_MAT, Sites_by_MAP,
                      ncol = 3, nrow = 1)
allplots <- annotate_figure(allplots,
                            bottom = text_grob(("Sites"), color = "black", size = 20, face = "bold"))

allplots

ggsave("SP_sites_by_X_axis.png", width = 20, height = 6, dpi = 300)



## Ibuttons ----
#To see ibuttons data and plots check folder "CLIMATIC_DATA/ibuttons"
#There there are also the plots for the ibuttons and the ibuttons plus MAT.



## Altitude ----
Sites_by_ALT <- ggplot(ALT_ord, aes(x=Site, y=altitude)) +
  geom_point(size= 3, color ="darkmagenta") +
  ylab("Altitude (masl)") +
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(size = 20, face = "bold", colour = "black")) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0,1200))+
  scale_y_continuous(breaks = breaks_width(200))

Sites_by_ALT


ggsave(path = "Figures/1 GRADIENT","SP_sites_by_ALT.png", width = 7, height = 6, dpi = 300)






#. ----
#Phsco variables by different OX ----
## Prepare data for plots ----
New_data <- my_data %>% gather(Physio, values, c(9:11,13:16,75)) #With TOC/TN and without TC/TN
head(New_data)

# reorder parameters
New_data$Physio <- factor(New_data$Physio, levels = c("Water_activity",
                                                      "Water_content", "SOM_perc","pH",
                                                      "TOC_perc", "TC_perc", "TN_perc",
                                                      "TOC_TN"))

#IF I NEED PHYSICOCHEMICAL DATA AS A FACTOR FOR REPRESENTING IT:
#New_data$AI <- as.factor(New_data$AI)

# IMPORT FUNCTION DATA_SUMMARY !!!

# summarize the data
New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site", "MAT", "MAP", "AI","Physio"))




## Plot ordered by AI ----

## EXPLANATION OF N:
# To know my N in the ggplots, I tried first to do one LM for one variable
# and obtain the N
a <- lm(Water_activity ~ AI, my_data)
summary(a)

length(a$residuals)
# In this case, if you do the LM separately from the plot, N=60 (by replicates)

aaa <- ggplot(my_data, aes(x=AI, y=Water_activity))+
  geom_point(size= 3, color ="grey")+
  geom_smooth(method='lm', formula= y~x)+
  stat_cor(label.x = 0, label.y.npc="bottom",
          aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
          p.accuracy = 0.001, r.accuracy = 0.01,
          color = "blue", size = 4.5)

aaa
# So if I add this LM to the plot, I obtain the p-value and R with the N = 60
# But if I look into the plots altogether (below), I see that the
# p-value and the R2 are slightly different for the water_activity
# Then I interpret that when I plot all the variables together (the 8 plots) and
# add the geom_smooth directly to the ggplot, N = 12 (by means).


New_data$AI = as.numeric(New_data$AI)
pphysio <- ggplot(New_data, aes(x=AI, y=values, color = Physio, 
                                 fill = Physio)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd), 
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ Physio , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0, label.y.npc="bottom",
                                  aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                                  p.accuracy = 0.001, r.accuracy = 0.01,
                                  color = "blue", size = 4.5) 
  # stat_regline_equation(label.x = 0.1, label.y = 2.5, 
  #                       #aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~`,`~"))
  # )

pphysio

# save the plot
ggsave(path = "Figures/1 GRADIENT","SP_physio2_AI.png", width = 20, height = 10, dpi = 300)


### Without SP6 ----
New_data <- subset(New_data, Site != "SP06")
pphysio <- ggplot(New_data, aes(x=AI, y=values, color = Physio, 
                                fill = Physio)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd), 
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ Physio , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0, label.y.npc="bottom",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5) 


pphysio

# save the plot
ggsave("SP_physio_AI_noSP6.png", width = 20, height = 10, dpi = 300)




## Plot ordered by MAT ----

pphysio <- ggplot(New_data, aes(x=MAT, y=values, color = Physio, 
                                fill = Physio)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd), 
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Mean Anual Temperature (ºC)") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ Physio , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(9,18))+
  scale_x_continuous(breaks = c(9, 12, 15, 18),
                     labels = c(9, 12, 15, 18)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey")+
  stat_cor(label.x = 13.5, label.y.npc="bottom",  
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
# stat_regline_equation(label.x = 0.1, label.y = 2.5, 
#                       #aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~`,`~"))
# )
pphysio


# save the plot
ggsave("SP_physio_MAT.png", width = 20, height = 10, dpi = 300)




## Plot ordered by MAP ----

pphysio <- ggplot(New_data, aes(x=MAP, y=values, color = Physio, 
                                fill = Physio)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd), 
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Mean Anual Precipitation (mm)") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ Physio , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(260,1400))+
  scale_x_continuous(breaks = breaks_width(350)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 400, label.y.npc="bottom",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
# stat_regline_equation(label.x = 0.1, label.y = 2.5, 
#                       #aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~`,`~"))
# )
pphysio


# save the plot
ggsave("SP_physio_MAP.png", width = 20, height = 10, dpi = 300)



# . ----
# VARIABLES BY AI ----

## Litter ----

# gathering data
New_data <- my_data %>% gather(litter, values, c(30:32))
head(New_data)

## Without SP6 ---- !!!!!!!!!!!!!!!!!
New_data <- subset(New_data, Site != "SP06")

# reorder parameters
New_data$litter <- factor(New_data$litter, levels = c("Litter","L_TC","L_TN"))

# IMPORT FUNCTION DATA_SUMMARY !!!

# summarize the data
New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site", "AI","litter"))

# make a plot
New_data$values = as.numeric(New_data$values)
litter <- ggplot(New_data, aes(x=AI, y=values, color = litter, 
                                fill = litter)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ litter , nrow = 1, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.65, label.y.npc="bottom",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)

litter

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_litter_AI.png", width = 16, height = 6, dpi = 300)
ggsave(path = "Figures/1 GRADIENT", "SP_litter_AI_noSP6.png", width = 16, height = 6, dpi = 300)



## Gas fluxes ----

New_data <- my_data %>% gather(Flux, Gas, c(35, 37:38))

## !!!!!!!!!!!!!!!Without SP6 !!!!!!!!!!!!!!!!!
New_data <- subset(New_data, Site != "SP06")


# reorder parameters
New_data$Flux <- factor(New_data$Flux, levels = c("CO2_dark",  "CH4_ave", "N2O"))
levels(New_data$Flux) <- c("CO[2]",  "CH[4]", "N2O")

# IMPORT FUNCTION DATA_SUMMARY !!!

New_data <- data_summary(New_data, varname="Gas", 
                          groupnames=c("AI",  "Flux"))

# plot
pCO2CH4N2O_2 <- ggplot(New_data, aes(x= AI, y=Gas, color = Flux,  fill = Flux)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=Gas-sd, ymax=Gas+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ Flux , nrow = 1, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.65, label.y.npc="bottom",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)

pCO2CH4N2O_2

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_gas_AI.png", width = 16, height = 6, dpi = 300)
ggsave(path = "Figures/1 GRADIENT", "SP_gas_AI_noSP6.png", width = 16, height = 6, dpi = 300)



## Enzymes ----
# my_data <- read.csv("SP_metadata_2021.csv", sep=",")
# my_data <- read.csv("SP_metadata_2021_SD.csv", sep=",")

New_data <- my_data %>% gather(enzyme, values, c(46:53))

## !!!!!!!!!!!!!!!Without SP6 !!!!!!!!!!!!!!!!!
# New_data <- subset(New_data, Site != "SP06")

# reorder parameters
New_data$enzyme <- factor(New_data$enzyme, levels = c("alpha", "beta", 
                                                      "xyl", "cbh", 
                                                      "gla", "fos",
                                                      "leu","phe" ))
# IMPORT FUNCTION DATA_SUMMARY !!!

New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "enzyme"))
# make a plot
New_data$values = as.numeric(New_data$values)
enzyme <- ggplot(New_data, aes(x=AI, y=values, color = enzyme, 
                               fill = enzyme)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data=subset(New_data, Site != "SP06"), aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
    xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ enzyme , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.6, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
enzyme

# save the plot
# ggsave(path = "Figures/1 GRADIENT", "SP_enzymes_AI_sd.png", width = 16, height = 8, dpi = 300)
# ggsave(path = "Figures/1 GRADIENT", "SP_enzymes_AI_noSP6_sd.png", width = 16, height = 8, dpi = 300)



#With SP06 but out of the correlation and colored differently FOR ALL of them:
enzyme <- ggplot(New_data, aes(x=AI, y=values, color = enzyme, 
                               fill = enzyme)) +
  geom_point(data=subset(New_data, Site != "SP06"), size= 3, color ="grey") +
  geom_point(data=subset(New_data, Site == "SP06"), alpha=.2, size= 3, color ="grey") +
  # geom_point(size= 3, color ="grey") +
  geom_pointrange(data=subset(New_data, Site != "SP06"), aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  geom_pointrange(data=subset(New_data, Site == "SP06"), aes(ymin=values-sd, ymax=values+sd),
                  color = "grey", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ enzyme , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(data=subset(New_data, Site != "SP06"), method = "lm", formula = y~x,
              color = "blue", fill = "grey") +
  stat_cor(data=subset(New_data, Site != "SP06"), label.x = 0.6, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
enzyme

ggsave(path = "Figures/1 GRADIENT", "SP_enzymes_AI_sp06excluded.png", width = 16, height = 8, dpi = 300)





#With SP06 but out of the correlation and colored
#differently FOR LEU ONLY:
enzyme <- ggplot(New_data, aes(x=AI, y=values, color = enzyme, 
                               fill = enzyme)) +
  geom_point(data=subset(New_data, Site != "SP06" | enzyme != "leu"), size= 3, color ="grey") +
  geom_point(data=subset(New_data, Site == "SP06" & enzyme == "leu"), alpha=.2, size= 3, color ="grey") +
  # geom_point(size= 3, color ="grey") +
  geom_pointrange(data=subset(New_data, Site != "SP06" | enzyme != "leu"), aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  geom_pointrange(data=subset(New_data, Site == "SP06" & enzyme == "leu"), aes(ymin=values-sd, ymax=values+sd),
                  color = "grey", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ enzyme , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(data=subset(New_data, Site != "SP06" | enzyme != "leu"), method = "lm", formula = y~x,
              color = "blue", fill = "grey") +
  stat_cor(data=subset(New_data, Site != "SP06" | enzyme != "leu"), label.x = 0.6, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
enzyme

# ggsave(path = "Figures/1 GRADIENT", "SP_enzymes_AI_sp06_leu_excluded.png", width = 16, height = 8, dpi = 300)
# ggsave(path = "Figures/1 GRADIENT", "SP_enzymes_AI_sp06_leu_excluded_sd.png", width = 16, height = 8, dpi = 300)



# ggplotly(enzyme) #To show which site is which, add "label=Site" inside aes.







### Shannon EEA ----
New_data <- my_data %>% gather(Shannon, values, c(74))

# reorder parameters
New_data$Shannon <- factor(New_data$Shannon, levels = c("ShannonEEA" ))

# IMPORT FUNCTION DATA_SUMMARY !!!

New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "Shannon"))
# make a plot
New_data$values = as.numeric(New_data$values)
sha <- ggplot(New_data, aes(x=AI, y=values, color = Shannon, 
                               fill = Shannon)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ Shannon , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.6, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
sha

ggsave(path = "Figures/1 GRADIENT", "SP_shannon_AI.png", width = 6, height = 6, dpi = 300)


#### No AI on OX ----

# OX: SOM
shannon_om <- ggplot(my_data, aes(x=SOM, y=ShannonEEA)) +
  geom_point(size= 3, color ="darkgrey") +
  theme(axis.title.x = element_text(size = 15, face="bold", colour = "black")) +
  theme(axis.title.y = element_text(size = 15, face="bold", colour = "black")) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  xlab("Soil Organic Matter") + ylab("Shannon index EEA")+
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.1, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)

shannon_om





## qPCR ----
New_data <- my_data %>% gather(genes, Copy, c(54:63))

## !!!!!!!!!!!!!!!Without SP6 !!!!!!!!!!!!!!!!!
New_data <- subset(New_data, Site != "SP06")

# reorder the sites
New_data$genes <- factor(New_data$genes, levels = c("X16S", "ITS2", "mcrA",
                                                    "pmoA", "nifH", "AOA",
                                                    "AOB", "qnorB", "nosZ",
                                                    "phoD"))
# IMPORT FUNCTION DATA_SUMMARY !!!

New_data <- data_summary(New_data, varname="Copy", 
                         groupnames=c("Site","AI", "genes"))

# make a plot
New_data$Copy = as.numeric(New_data$Copy)
genes <- ggplot(New_data, aes(x=AI, y=Copy, color = genes, 
                               fill = genes)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=Copy-sd, ymax=Copy+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ genes , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
genes

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_qPCR_AI.png", width = 20, height = 10, dpi = 300)
ggsave(path = "Figures/1 GRADIENT", "SP_qPCR_AI_noSP6.png", width = 20, height = 10, dpi = 300)






## Microbial biomasses ----
New_data <- my_data %>% gather(mb, values, c(40:41))

## !!!!!!!!!!!!!!!Without SP6 !!!!!!!!!!!!!!!!!
New_data <- subset(New_data, Site != "SP06")


# reorder parameters
New_data$mb <- factor(New_data$mb, levels = c("BB", "FB"))

# IMPORT FUNCTION DATA_SUMMARY !!!

New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "mb"))
# make a plot
New_data$values = as.numeric(New_data$values)
biomasses <- ggplot(New_data, aes(x=AI, y=values, color = mb, 
                               fill = mb)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ mb , nrow = 1, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.1, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
biomasses

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_MB_AI.png", width = 12, height = 6, dpi = 300)
ggsave(path = "Figures/1 GRADIENT", "SP_MB_AI_noSP6.png", width = 12, height = 6, dpi = 300)



###  16s, ITS2 + BB, FB ----
New_data <- my_data %>% gather(biogene, values, c(40,41,54,55))

## !!!!!!!!!!!!!!!Without SP6 !!!!!!!!!!!!!!!!!
New_data <- subset(New_data, Site != "SP06")

New_data$biogene <- factor(New_data$biogene, levels = c("BB","FB","X16S", "ITS2"))
New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "biogene"))

# make a plot
New_data$values = as.numeric(New_data$values)
ratio1 <- ggplot(New_data, aes(x=AI, y=values, color = biogene, 
                               fill = biogene)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ biogene , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))+
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.1, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
ratio1

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_microbial-comparison_AI.png", width = 10, height = 8, dpi = 300)
ggsave(path = "Figures/1 GRADIENT", "SP_microbial-comparison_noSP6_AI.png", width = 10, height = 8, dpi = 300)








##  Pigments ----

New_data <- my_data %>% gather(pigments, values, c(42:45))

## !!!!!!!!!!!!!!!Without SP6 !!!!!!!!!!!!!!!!!
New_data <- subset(New_data, Site != "SP06")

# reorder parameters
New_data$pigments <- factor(New_data$pigments, levels = c("chla", "chlb", 
                                                          "carotene","EPS"))

# IMPORT FUNCTION DATA_SUMMARY !!!

New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "pigments"))
# make a plot
New_data$values = as.numeric(New_data$values)
pigments <- ggplot(New_data, aes(x=AI, y=values, color = pigments, 
                                  fill = pigments)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ pigments , nrow = 1, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.1, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
pigments

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_pigments_AI.png", width = 20, height = 6, dpi = 300)
ggsave(path = "Figures/1 GRADIENT", "SP_pigments_AI_noSP6.png", width = 20, height = 6, dpi = 300)




## Respiration ----
New_data <- my_data %>% gather(respi, values, c(64))

## !!!!!!!!!!!!!!!Without SP6 !!!!!!!!!!!!!!!!!
New_data <- subset(New_data, Site != "SP06")

# reorder parameters
New_data$respi <- factor(New_data$respi, levels = c("Respiration"))


# IMPORT FUNCTION DATA_SUMMARY !!!

New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "respi"))

# make a plot
New_data$values = as.numeric(New_data$values)
respiration <- ggplot(New_data, aes(x=AI, y=values, color = respi, 
                                 fill = respi)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ respi , nrow = 1, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.1, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
respiration

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_respiration_AI.png", width = 6, height = 6, dpi = 300)






## Other: Ratios ----

### BB/BF and 16S/ITS2 ----
New_data <- my_data %>% gather(bactfung, values, c(65,66))

# ## !!!!!!!!!!!!!!!Without SP6 !!!!!!!!!!!!!!!!!
# New_data <- subset(New_data, Site != "SP06")

# reorder parameters
New_data$bactfung <- factor(New_data$bactfung, levels = c("BB_BF", "X16S_ITS2"))


# IMPORT FUNCTION DATA_SUMMARY !!!

New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "bactfung"))

# make a plot
New_data$values = as.numeric(New_data$values)
ratio <- ggplot(New_data, aes(x=AI, y=values, color = bactfung, 
                                    fill = bactfung)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ bactfung , nrow = 1, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.1, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
ratio

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_ratio_bactfung_AI.png", width = 12, height = 6, dpi = 300)




### BB/16S and BF/ITS2 ----
New_data <- my_data %>% gather(biogene, values, c(67,68))

# ## !!!!!!!!!!!!!!!Without SP6 !!!!!!!!!!!!!!!!!
# New_data <- subset(New_data, Site != "SP06")

# reorder parameters
New_data$biogene <- factor(New_data$biogene, levels = c("BB_16S", "FB_ITS2"))


# IMPORT FUNCTION DATA_SUMMARY !!!

New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "biogene"))

# make a plot
New_data$values = as.numeric(New_data$values)
ratio1 <- ggplot(New_data, aes(x=AI, y=values, color = biogene, 
                              fill = biogene)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ biogene , nrow = 1, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
ratio1

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_ratio_biogene_AI.png", width = 12, height = 6, dpi = 300)



##### No AI on OX----
ratio1 <- ggplot(my_data, aes(x=X16S, y=BB)) +
  geom_point(size= 3, color ="darkgrey") +
  theme(axis.title.x = element_text(size = 15, face="bold", colour = "black")) +
  theme(axis.title.y = element_text(size = 15, face="bold", colour = "black")) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  xlab("16S (copies/gDW)") + ylab("BB (mg C bact/gDW)")+
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.1, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)

ratio1


ratio2 <- ggplot(my_data, aes(x=ITS2, y=FB)) +
  geom_point(size= 3, color ="darkgrey") +
  theme(axis.title.x = element_text(size = 15, face="bold", colour = "black")) +
  theme(axis.title.y = element_text(size = 15, face="bold", colour = "black")) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  xlab("ITS2 (copies/gDW)") + ylab("FB (mg C fung/gDW)")+
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.1, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
ratio2

ratios12 <- ggarrange(ratio1, ratio2,
                    ncol = 2, nrow = 1,
                    common.legend = FALSE)
ratios12
# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_BB16S-FBITS2_AI.png", width = 12, height = 6, dpi = 300)





### Chla/chlb, carotene/chla, EPS/chla ----
New_data <- my_data %>% gather(pig, values, c(69,70,71))
New_data$pig <- factor(New_data$pig, levels = c("chla_chlb", "carotene_chla","EPS_chla"))
New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "pig"))

New_data$values = as.numeric(New_data$values)
ratio3 <- ggplot(New_data, aes(x=AI, y=values, color = pig, 
                              fill = pig)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ pig , nrow = 1, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.1, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
ratio3

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_ratio_pigments_AI.png", width = 16, height = 6, dpi = 300)




### EPS/BB and EPS/16S ----
New_data <- my_data %>% gather(pigbac, values, c(72,73))
New_data$pigbac <- factor(New_data$pigbac, levels = c("EPS_BB", "EPS_16S"))
New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "pigbac"))

New_data$values = as.numeric(New_data$values)
ratio3 <- ggplot(New_data, aes(x=AI, y=values, color = pigbac, 
                               fill = pigbac)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data = New_data, aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ pigbac , nrow = 1, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  geom_smooth(method = "lm", color = "blue", fill = "grey") +
  stat_cor(label.x = 0.1, label.y.npc="top",
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           p.accuracy = 0.001, r.accuracy = 0.01,
           color = "blue", size = 4.5)
ratio3

# save the plot
ggsave(path = "Figures/1 GRADIENT", "SP_EPS-bact_AI.png", width = 12, height = 6, dpi = 300)





#. ----

# Enzymes PCA ----

enz <- my_data[,c(2,46:53)]

#To replace NA values with a mean of the other values of the Site:
for (i in which(sapply(enz, is.numeric))) {
  for (j in which(is.na(enz[, i]))) {
    enz[j, i] <- mean(enz[enz[, "Site"] == enz[j, "Site"], i],  na.rm = TRUE)
  }
}

enz2 <- enz[,c(2:9)]



### Op. 1, means ----
enz_mit <- ddply(enz, .(Site), summarize,
                 alpha=mean(alpha),
                 beta=mean(beta),
                 xyl=mean(xyl),
                 cbh=mean(cbh),
                 gla=mean(gla),
                 fos=mean(fos),
                 leu=mean(leu),
                 phe=mean(phe))

# pcr_mit <- pcr %>% 
#   group_by(Site) %>%
#   summarise_all("mean")

site_order <- my_data[,c(2,7)] 
site_order <- site_order[!duplicated(site_order), ] #Erase duplicated lines from dataframe
enz_mit <- enz_mit[order(site_order$AI, decreasing = T),]
enz_mit$Site <- factor(enz_mit$Site, levels = enz_mit$Site[order(site_order$AI)])

enz_mit2 <- enz_mit[,c(2:9)]


#Select column with levels (Site)
site <- paste("SP", order(enz_mit[, 1], decreasing =  TRUE), sep='')
site <- factor(site, levels = site)

pc <- prcomp(na.omit(enz_mit2), center = TRUE,
              scale. = TRUE) 

plot(pc, type = "l")
plot(pc)
summary(pc)

# colourCount = length(unique(enz_mit$Site))
# mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(colourCount)

mycolors2<-c("#10449F","#51B7DF","#00FFFF","#00A01D","#064700",
             "#7A7615", "#583200", "#C24A0A", "#F5C92D","#FA0C00",
             "#7D1809", "#290500")

eea <- ggbiplot(pc, obs.scale = 1, var.scale = 1, 
              ellipse = FALSE, fill=site,
              varname.size = 5,
              circle = TRUE, alpha=0) +
  theme_classic()+
  scale_fill_manual(values = mycolors2)+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 5)+
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=15))+
  ggtitle("Enzymes Spain (PCA)")+
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=15))
eea


# save the plot
# ggsave(path = "Figures/1 GRADIENT","SP_PCA_enzymes.png", width = 7, height = 6, dpi = 300)



### Op. 2, replicates ----
#Select column with levels (Site)
site <- enz$Site
site #To know the order of the sites and then to write them in order in the ggbiplot


pc <- prcomp(na.omit(enz2), center = TRUE,
             scale. = TRUE) 

plot(pc, type = "l")
plot(pc)
summary(pc)

mycolors2<-c("#10449F","#51B7DF","#00FFFF","#00A01D","#064700",
             "#7A7615", "#583200", "#C24A0A", "#F5C92D","#FA0C00",
             "#7D1809", "#290500")

eea <- ggbiplot(pc, obs.scale = 1, var.scale = 1, 
                ellipse = FALSE, fill=site,
                varname.size = 5,
                circle = TRUE, alpha=0) +
  theme_classic()+
  scale_fill_manual(values = mycolors2,
                    breaks=c('SP08', 'SP01', 'SP02','SP07',
                             'SP06', 'SP03','SP12','SP11',
                             'SP04','SP09','SP10','SP05'))+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 5)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))+
  ggtitle("Enzymes Spain (PCA)")+
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=15))
eea


# save the plot
# ggsave(path = "Figures/1 GRADIENT","SP_PCA_enzymes_replicates.png", width = 7, height = 6, dpi = 300)


### Op.1+2, means with SD ----
#Obtain scores:
scores <- as.data.frame(pc$x[,1:2])
scores$Site <- my_data[,2]
scores <- scores[,c(3,1,2)] #Reorder to have Site as the first column

scores_msd <- aggregate(cbind(PC1,PC2) ~ Site, data = scores, 
                                     FUN = function(scores) c(mean = mean(scores), sd = sd(scores)))

#Si no ho faig així no em deixa canviar els noms de les columnes ni graficar
#perquè interpreta que mean i sd són la mateixa columna dividida en 2 o algo raro
PC1<- as.data.frame(scores_msd[,2])
PC2<- as.data.frame(scores_msd[,3])

names(PC1)[names(PC1) == "mean"] <- "PC1_mean"
names(PC1)[names(PC1) == "sd"] <- "PC1_sd"
names(PC2)[names(PC2) == "mean"] <- "PC2_mean"
names(PC2)[names(PC2) == "sd"] <- "PC2_sd"

PC1$Site <- scores_msd$Site
PC1 <- PC1[,c(3,1,2)]
PC2$Site <- scores_msd$Site
PC2 <- PC2[,c(3,1,2)]

scores_msd <- cbind(PC1, PC2[,2:3])




mycolors2<-c("#10449F","#51B7DF","#00FFFF","#00A01D","#064700",
             "#7A7615", "#583200", "#C24A0A", "#F5C92D","#FA0C00",
             "#7D1809", "#290500")



# site_order <- my_data[,c(2,7)] 
# site_order <- site_order[!duplicated(site_order), ] #Erase duplicated lines from dataframe
# scores_msd <- scores_msd[order(site_order$AI, decreasing = T),]
# scores_msd$Site <- factor(scores_msd$Site, levels = scores_msd$Site[order(site_order$AI)])
# 
# site <- paste("SP", order(scores_msd[, 1], decreasing =  TRUE), sep='')
# site <- factor(site, levels = site)


pca_sd <- ggplot(scores_msd, aes(x=PC1_mean, y=PC2_mean)) +
  geom_point(aes(fill=Site), colour= "black", pch=21, size = 5)+
  geom_errorbarh(aes(xmax = PC1_mean + PC1_sd, xmin = PC1_mean - PC1_sd))+
  geom_errorbar(aes(ymax = PC2_mean + PC2_sd, ymin = PC2_mean - PC2_sd))+
  theme_classic()+
  scale_fill_manual(values = mycolors2,
                    breaks=c('SP08', 'SP01', 'SP02','SP07',
                             'SP06', 'SP03','SP12','SP11',
                             'SP04','SP09','SP10','SP05'))+
  ylab("PC2 (25.9% explained var.)") +
  xlab("PC1 (36.2% explained var.)") +
  theme(axis.text=element_text(size=12),
      axis.title=element_text(size=15))+
  ggtitle("Enzymes Spain (PCA)")+
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=15))
pca_sd

ggsave(path = "Figures/1 GRADIENT","SP_PCA_enzymes_mean_sd.png", width = 7, height = 6, dpi = 300)


  
  
  



#### X Op. 3 ----
library(ggfortify)
pc <- autoplot(pc,
               data = enz_mit,
               colour = 'Site')
pc

#### X Op. 4 ---- 

enz.pca <- rda(enz_mit2, scale = TRUE) 

summary(enz.pca)

barplot(enz.pca$CA$eig, ylab='Eigenvalues', las=1,
        main='Kaiser-Gutman')
abline(h=mean(enz.pca$CA$eig), lty=2, col='dark red')

par(mfrow=c(1,2))
biplot(enz.pca, scaling = 1, main = "Scaling 1")
biplot(enz.pca, main = "Scaling 2")  # Default scaling 2

#.----



# qPCR PCA ----
## Without ITS2 and 16S
pcr <- my_data[,c(2,56:63)]

#To replace NA values with a mean of the other values of the Site:
for (i in which(sapply(pcr, is.numeric))) {
  for (j in which(is.na(pcr[, i]))) {
    pcr[j, i] <- mean(pcr[pcr[, "Site"] == pcr[j, "Site"], i],  na.rm = TRUE)
  }
}

pcr2 <- pcr[,-1]


### Op. 1, means ----
pcr_mit <- ddply(pcr, .(Site), summarize,
                 mcrA=mean(mcrA),
                 pmoA=mean(pmoA),
                 nifH=mean(nifH),
                 AOA=mean(AOA),
                 AOB=mean(AOB),
                 qnorB=mean(qnorB),
                 nosZ=mean(nosZ),
                 phoD=mean(phoD))

# pcr_mit <- pcr %>% 
#   group_by(Site) %>%
#   summarise_all("mean")

pcr_mit2 <- pcr_mit[,c(2:9)]

site_order <- my_data[,c(2,7)] 
site_order <- site_order[!duplicated(site_order), ] #Erase duplicated lines from dataframe
pcr_mit <- pcr_mit[order(site_order$AI, decreasing = T),]
pcr_mit$Site <- factor(pcr_mit$Site, levels = pcr_mit$Site[order(site_order$AI)])

pcr_mit2 <- pcr_mit[,c(2:9)]




#Select column with levels (Site)
site <- paste("SP", order(pcr_mit[, 1], decreasing =  TRUE), sep='')
site <- factor(site, levels = site)

pc2 <- prcomp(na.omit(pcr_mit2), center = TRUE,
             scale. = TRUE) 

plot(pc2, type = "l")
plot(pc2)
summary(pc2)


mycolors2<-c("#10449F","#51B7DF","#00FFFF","#00A01D","#064700",
             "#7A7615", "#583200", "#C24A0A", "#F5C92D","#FA0C00",
             "#7D1809", "#290500")

qpcr <- ggbiplot(pc2, obs.scale = 1, var.scale = 1, 
                ellipse = FALSE, fill=site,
                varname.size = 5,
                circle = TRUE, alpha=0) +
  theme_classic()+
  scale_fill_manual(values = mycolors2)+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 5)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))+
  ggtitle("qPCR Spain (PCA)")+
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=15))
qpcr


# save the plot
ggsave(path = "Figures/1 GRADIENT","SP_PCA_qPCR.png", width = 7, height = 6, dpi = 300)





### Op. 2, replicates ----
#Select column with levels (Site)
site <- pcr$Site
site #To know the order of the sites and then to write them in order in the ggbiplot


pc <- prcomp(pcr2, center = TRUE,
             scale. = TRUE) 

plot(pc, type = "l")
plot(pc)
summary(pc)

mycolors2<-c("#10449F","#51B7DF","#00FFFF","#00A01D","#064700",
             "#7A7615", "#583200", "#C24A0A", "#F5C92D","#FA0C00",
             "#7D1809", "#290500")

qpcr <- ggbiplot(pc, obs.scale = 1, var.scale = 1, 
                ellipse = FALSE, fill=site,
                varname.size = 5,
                circle = TRUE, alpha=0) +
  theme_classic()+
  scale_fill_manual(values = mycolors2,
                    breaks=c('SP08', 'SP01', 'SP02','SP07',
                             'SP06', 'SP03','SP12','SP11',
                             'SP04','SP09','SP10','SP05'))+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 5)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))+
  ggtitle("qPCR Spain (PCA)")+
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=15))
qpcr




# save the plot
ggsave(path = "Figures/1 GRADIENT","SP_PCA_qPCR_replicates.png", width = 7, height = 6, dpi = 300)






#.----

# ALL data PCA ----
data <- my_data[,c(2,5:74)]

#To replace NA values with a mean of the other values of the Site:
for (i in which(sapply(data, is.numeric))) {
  for (j in which(is.na(data[, i]))) {
    data[j, i] <- mean(data[data[, "Site"] == data[j, "Site"], i],  na.rm = TRUE)
  }
}

data2 <- data[,c(-1)]

### Op. 1, means ----
data_mit <- data %>% 
  group_by(Site) %>%
  summarise_all("mean")

data_mit2 <- data_mit[,c(-1)]


#Select column with levels (Site)
site <- data_mit[,1]
site <- site$Site
site <- factor(site, levels = site)

all <- prcomp(na.omit(data_mit2), center = TRUE,
              scale. = TRUE) 

plot(all, type = "l")
plot(all)
summary(all)


mycolors2<-c("#10449F","#51B7DF","#00FFFF","#00A01D","#064700",
             "#7A7615", "#583200", "#C24A0A", "#F5C92D","#FA0C00",
             "#7D1809", "#290500")

all_plot <- ggbiplot(all, obs.scale = 1, var.scale = 1, 
                 ellipse = FALSE, fill=site,
                 varname.size = 5,
                 circle = TRUE, alpha=0) +
  theme_classic()+
  scale_fill_manual(values = mycolors2,
                    breaks=c('SP08', 'SP01', 'SP02','SP07',
                             'SP06', 'SP03','SP12','SP11',
                             'SP04','SP09','SP10','SP05'))+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 5)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))+
  ggtitle("All data Spain (PCA)")+
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=15))
all_plot


# save the plot
ggsave(path = "Figures/1 GRADIENT","SP_PCA_all.png", width = 7, height = 6, dpi = 300)







### Op. 2, replicates ----
#Select column with levels (Site)

site <- data[, 1]

data2 <- data[,-1]

pc <- prcomp(na.omit(data2), center = TRUE,
             scale. = TRUE) 

plot(pc, type = "l")
plot(pc)
summary(pc)

mycolors2<-c("#10449F","#51B7DF","#00FFFF","#00A01D","#064700",
             "#7A7615", "#583200", "#C24A0A", "#F5C92D","#FA0C00",
             "#7D1809", "#290500")

all <- ggbiplot(pc, var.scale = 1,varname.size = 5,
                obs.scale = 1,
                 ellipse = FALSE, fill=site,
                 circle = TRUE, alpha=0) +
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))+
  ggtitle("All data Spain (PCA)")+
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=15))+
  scale_fill_manual(values = mycolors2,
                    breaks=c('SP08', 'SP01', 'SP02','SP07',
                             'SP06', 'SP03','SP12','SP11',
                             'SP04','SP09','SP10','SP05'))+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 5)
all

# save the plot
ggsave(path = "Figures/1 GRADIENT","SP_PCA_all_replicate.png", width = 7, height = 6, dpi = 300)




#.----
# FSCOQA PCA ----
#BY MEANS
pca_data <- my_data[,c(2,5:10,12,13,17:20,27:29)]

pca_data <- pca_data %>%
  group_by(Site) %>%
  summarise_all("mean")

site_order <- my_data[,c(2,7)] 
site_order <- site_order[!duplicated(site_order), ] #Erase duplicated lines from dataframe
pca_data <- pca_data[order(site_order$AI, decreasing = T),]
pca_data$Site <- factor(pca_data$Site, levels = pca_data$Site[order(site_order$AI)])

pcr2 <- pca_data[,c(-1)]


#Select column with levels (Site)
site <- factor(pca_data$Site, levels = site)
site

pc <- prcomp(na.omit(pcr2), center = TRUE,
             scale. = TRUE) 

plot(pc, type = "l")
plot(pc)
summary(pc)

mycolors2<-c("#10449F","#51B7DF","#00FFFF","#00A01D","#064700",
             "#7A7615", "#583200", "#C24A0A", "#F5C92D","#FA0C00",
             "#7D1809", "#290500")

library(ggfortify)
autoplot(pc, data=pca_data, 
         loadings = TRUE, loadings.colour = 'brown',
         loadings.label.colour='brown', loadings.label = TRUE,
         loadings.label.size = 7,
         loadings.label.repel=TRUE)+
  theme_classic()+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 5)+
  scale_fill_manual(values = mycolors2)+
  ggtitle("Physicochemical variables")+
  theme(legend.title = element_blank(),
        legend.text=element_text(size = 12),
        title = element_text(size = 15,face="bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15, face="plain"))

# ggsave(path = "Figures","PCA_Fscq_means.png", width = 10, height = 8, dpi = 300)




#.----
# (All enzymes all gradients) ----
# Data: summary_out_value_20221221 (Antic folder)
library(readxl)
library(dplyr)
my_data <- read_excel("Antic/summary_out_value_20221221.xlsx")

my_data <- rename(data,
               alpha = alpha_umolMUF_gDW_h,
               beta = beta_umolMUF_gDW_h,
               xyl = xyl_umolMUF_gDW_h,
               cbh = cbh_umolMUF_gDW_h,
               gla = gla_umolMUF_gDW_h,
               fos = fos_umolMUF_gDW_h,
               leu = leu_umolAMC_gDW_h,
               phe = umolDIQC_gDW_h_phe
)
step <- my_data[my_data$Site == 'SP3',]
step$Site <- rep('EU5', 5)
my_data <- rbind(my_data,step)
my_data <- my_data[order(my_data$Site),]
my_data$Gradient <- substr(my_data$Site, 1,2)


New_data <- my_data %>% gather(enzyme, values, c(3:10))

## !!!!!!!!!!!!!!!Without SP6 !!!!!!!!!!!!!!!!!
# New_data <- subset(New_data, Site != "SP06")

# reorder parameters
New_data$enzyme <- factor(New_data$enzyme, levels = c("alpha", "beta", 
                                                      "xyl", "cbh", 
                                                      "gla", "fos",
                                                      "leu","phe" ))


# meta <- read_excel("Antic/Coordinates_alt.xlsx")
meta <- read_csv("C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/ANALISIS/CLIMATIC_DATA/Extracted_climaticData/results/extractedClimaticVariables.csv")
Ox <- data.frame(Gradient=c('SP', 'AL', 'EU', 'GL', 'SA'), Variable=c('AI', 'AI', 'AI', 'AI', 'AI')) 


# variables <- c('alpha', 'beta', 'gla', 'fos', 'leu', 'xyl', 'cbh', 'phenox')
# unitats <- c(rep(expression(paste(mu,'mol MUF DW'^-1*'h'^-1)),7),   expression(paste(mu,'mol DIQC DW'^-1*'h'^-1)))
# linia_mtext <- c(3,3,3,2,3,3,3,2)
meta <- rename(meta,
               aaa = Site,
               Site = Sites)

New_data <- merge(New_data, meta[,c(2,11)], by.x=c("Site"), by.y=c("Site"), sort=T, all.x = TRUE)


# IMPORT FUNCTION DATA_SUMMARY !!!

New_data <- data_summary(New_data, varname="values", 
                         groupnames=c("Site","AI", "enzyme"))
# make a plot
New_data$values = as.numeric(New_data$values)
enzyme <- ggplot(New_data, aes(x=AI, y=values, color = enzyme, 
                               fill = enzyme)) +
  geom_point(size= 3, color ="grey") +
  geom_pointrange(data=subset(New_data, Site != "SP06"), aes(ymin=values-sd, ymax=values+sd),
                  color = "black", stroke = 1, size = 0.6, shape = 0) +
  xlab("Aridity Index") +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black")) +
  facet_wrap( .~ enzyme , nrow = 2, scales = "free_y", labeller = label_parsed) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black", face = "bold"))+
  theme(axis.text.y = element_text(size = 18, color = "black", face = "bold"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  theme(legend.position = "none")+
  # coord_cartesian(xlim = c(0,1.4))+
  # scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
  #                    labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
   geom_smooth(color = "blue", fill = "grey") #+
  # stat_cor(label.x = 0.6, label.y.npc="top",
  #          aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
  #          p.accuracy = 0.001, r.accuracy = 0.01,
  #          color = "blue", size = 4.5)
enzyme

# save the plot
ggsave(path = "Figures/1 GRADIENT", "ALL_enzymes_loess.png", width = 16, height = 8, dpi = 300)
 
