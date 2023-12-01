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
library(rcompanion)


## Working directory
# setwd("G:/GRADCATCH/ANALISIS/R")
# getwd()
# version


# IMPORT DATA ----
my_data <- read.csv("SP_metadata_2021.csv", sep=";")

#To replace NA values with a mean of the other values of the Site:
for (i in which(sapply(my_data, is.numeric))) {
  for (j in which(is.na(my_data[, i]))) {
    my_data[j, i] <- mean(my_data[my_data[, "Site"] == my_data[j, "Site"], i],  na.rm = TRUE)
  }
}

#Data without ratios, percentages....
data <- my_data[,c(5:10,12:13,17:23,27:30,33:35,37:38,40:64,76,79:92)]

#Dummies....
dummy <- my_data[,c(5:7,76)]

#Independent variables without dummies....
indepe <- my_data[,c(8:10,12:13,17:23,27:30,33:34,79:92)]

#Independent var with dummies...
duindepe <- my_data[,c(5:7,76,8:10,12:13,17:23,27:30,33:34,79:92)]

#Dependen vars....
depe <- my_data[,c(35,37:38,40:64)]


#Functional dependent vars ....
# MICROBIAL BIOMASSES (aggregated)
# RESPIRATION
# EEA (aggregated)
func <- my_data[,c(2,40,41,46:53,64)]

func_ag <-func %>%
  mutate(Cenz = select(., 4:7) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(MB = select(.,2:3) %>% rowSums(na.rm = TRUE))
func_ag <- func_ag[,-c(2:7)]



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


data_summary2 <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# _______________________________________ ----
# UCI STAY ANALYSES **** ----
#.----

# Clustering sites ----
# https://bookdown.org/stephi_gascon/bookdown-demo-master_-_multivariant/_book/cluster-analysis.html

# > By climatic vars. ----
library(ade4) 
library(vegan)  
library(gclus) 
library(cluster)
library(RColorBrewer)  
library(labdsv)
library(leaflet)

# 1. Create dissimilarity matrix:
# https://stats.stackexchange.com/questions/80377/which-distance-to-use-e-g-manhattan-euclidean-bray-curtis-etc

## Trec les repliques per site pq totes tenen igual valors
dummy2 <- as.data.frame(unique(dummy))

boxplot(dummy2)
# Vars. with different units and scales
# Standardization its needed

# If not, MAP and altitude would have much more
# influence on the distance matrix only because of
# their magnitude and variation in absolut terms.

env.std <- scale(dummy2)
boxplot(env.std)
# Nicer

# Using Euclidean distance:
env.de <- dist(env.std)
attr(env.de, "Labels") <- unique(my_data$Site)
env.de


# 2. Analysis of hierarchical cluster
par(mfrow = c(2, 2))

# Compute single linkage agglomerative clustering
env.de.single <- hclust(env.de, method = "single")
plot(env.de.single, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Single linkage")

# Compute complete-linkage agglomerative clustering
env.de.complete <- hclust(env.de, method = "complete")
plot(env.de.complete, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Complete linkage")

# Compute UPGMA agglomerative clustering
env.de.UPGMA <- hclust(env.de, method = "average")
plot(env.de.UPGMA, 
     labels = unique(my_data$Site), 
     main = "Euclidean - UPGMA")

# Compute Ward’s Minimum Variance Clustering
env.de.ward <- hclust(env.de, method = "ward.D2")
plot(env.de.ward, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Ward")
dev.off()

# 3. Comparison dendograms by cophenetic distance
# HIGHEST CORRELATION = BEST DENDOGRAM

# Single linkage clustering
env.de.single.coph <- cophenetic(env.de.single)
cor(env.de, env.de.single.coph)

# Complete linkage clustering
env.de.comp.coph <- cophenetic(env.de.complete)
cor(env.de, env.de.comp.coph)

# Average clustering
env.de.UPGMA.coph <- cophenetic(env.de.UPGMA)
cor(env.de, env.de.UPGMA.coph)

# Ward clustering
env.de.ward.coph <- cophenetic(env.de.ward)
cor(env.de, env.de.ward.coph)

# Ward > UPGMA > Complete > Single


# 4. Comparison dendograms by Gower distance
# LOWEST VALUE = BEST DENDOGRAM
(gow.dist.single <- sum((env.de - env.de.single.coph) ^ 2))
(gow.dist.comp <- sum((env.de - env.de.comp.coph) ^ 2))
(gow.dist.UPGMA <- sum((env.de - env.de.UPGMA.coph) ^ 2))
(gow.dist.ward <- sum((env.de - env.de.ward.coph) ^ 2))

# UPGMA >> Complete > Single >>>>>> Ward

# Considering both results from points 3 and 4,
# we consider the UPGMA as the best dendogram option.

plot(env.de.UPGMA, hang=-1, labels = unique(my_data$Site),  main = "Euclidean - UPGMA")

# png(file = "Figures/1 GRADIENT/cluster_by_ENV.png", width = 700, height = 600)



# > By climatic vars. WITHOUT MAP or ALTITUDE ----
library(ade4) 
library(vegan)  
library(gclus) 
library(cluster)
library(RColorBrewer)  
library(labdsv)
library(leaflet)

# 1. Create dissimilarity matrix:
# https://stats.stackexchange.com/questions/80377/which-distance-to-use-e-g-manhattan-euclidean-bray-curtis-etc

## Trec les repliques per site pq totes tenen igual valors
dummy2 <- as.data.frame(unique(dummy[,-c(1,4)]))

boxplot(dummy2)
# Vars. with different units and scales
# Standardization its needed

# If not, MAP and altitude would have much more
# influence on the distance matrix only because of
# their magnitude and variation in absolut terms.

env.std <- scale(dummy2)
boxplot(env.std)
# Nicer

# Using Euclidean distance:
env.de <- dist(env.std)
attr(env.de, "Labels") <- unique(my_data$Site)
env.de


# 2. Analysis of hierarchical cluster
par(mfrow = c(2, 2))

# Compute single linkage agglomerative clustering
env.de.single <- hclust(env.de, method = "single")
plot(env.de.single, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Single linkage")

# Compute complete-linkage agglomerative clustering
env.de.complete <- hclust(env.de, method = "complete")
plot(env.de.complete, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Complete linkage")

# Compute UPGMA agglomerative clustering
env.de.UPGMA <- hclust(env.de, method = "average")
plot(env.de.UPGMA, 
     labels = unique(my_data$Site), 
     main = "Euclidean - UPGMA")

# Compute Ward’s Minimum Variance Clustering
env.de.ward <- hclust(env.de, method = "ward.D2")
plot(env.de.ward, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Ward")

par(mfrow = c(1,1))

# 3. Comparison dendograms by cophenetic distance
# HIGHEST CORRELATION = BEST DENDOGRAM

# Single linkage clustering
env.de.single.coph <- cophenetic(env.de.single)
cor(env.de, env.de.single.coph)

# Complete linkage clustering
env.de.comp.coph <- cophenetic(env.de.complete)
cor(env.de, env.de.comp.coph)

# Average clustering
env.de.UPGMA.coph <- cophenetic(env.de.UPGMA)
cor(env.de, env.de.UPGMA.coph)

# Ward clustering
env.de.ward.coph <- cophenetic(env.de.ward)
cor(env.de, env.de.ward.coph)

# Ward > UPGMA > Complete > Single


# 4. Comparison dendograms by Gower distance
# LOWEST VALUE = BEST DENDOGRAM
(gow.dist.single <- sum((env.de - env.de.single.coph) ^ 2))
(gow.dist.comp <- sum((env.de - env.de.comp.coph) ^ 2))
(gow.dist.UPGMA <- sum((env.de - env.de.UPGMA.coph) ^ 2))
(gow.dist.ward <- sum((env.de - env.de.ward.coph) ^ 2))

# UPGMA >> Complete > Single >>>>>> Ward

# Considering both results from points 3 and 4,
# we consider the UPGMA as the best dendogram option.

plot(env.de.UPGMA, hang=-1, labels = unique(my_data$Site),  main = "Euclidean - UPGMA")

# png(file = "Figures/1 GRADIENT/cluster_by_ENV.png", width = 700, height = 600)




# # > By climatic vars. WITHOUT ALTITUDE ----
library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(RColorBrewer)
library(labdsv)
library(leaflet)

# 1. Create dissimilarity matrix:
# https://stats.stackexchange.com/questions/80377/which-distance-to-use-e-g-manhattan-euclidean-bray-curtis-etc

## Trec les repliques per site pq totes tenen igual valors
dummy2 <- as.data.frame(unique(dummy))
dummy2 <- dummy2[,-c(4)]

boxplot(dummy2)

env.std <- scale(dummy2)
boxplot(env.std)


# Using Euclidean distance:
env.de <- dist(env.std)
attr(env.de, "Labels") <- unique(my_data$Site)
env.de

# 2. Analysis of hierarchical cluster
par(mfrow = c(2, 2))

# Compute single linkage agglomerative clustering
env.de.single <- hclust(env.de, method = "single")
plot(env.de.single, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Single linkage")

# Compute complete-linkage agglomerative clustering
env.de.complete <- hclust(env.de, method = "complete")
plot(env.de.complete, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Complete linkage")

# Compute UPGMA agglomerative clustering
env.de.UPGMA <- hclust(env.de, method = "average")
plot(env.de.UPGMA, 
     labels = unique(my_data$Site), 
     main = "Euclidean - UPGMA")

# Compute Ward’s Minimum Variance Clustering
env.de.ward <- hclust(env.de, method = "ward.D2")
plot(env.de.ward, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Ward")


# 3. Comparison dendograms by cophenetic distance
# HIGHEST CORRELATION = BEST DENDOGRAM

# Single linkage clustering
env.de.single.coph <- cophenetic(env.de.single)
cor(env.de, env.de.single.coph)

# Complete linkage clustering
env.de.comp.coph <- cophenetic(env.de.complete)
cor(env.de, env.de.comp.coph)

# Average clustering
env.de.UPGMA.coph <- cophenetic(env.de.UPGMA)
cor(env.de, env.de.UPGMA.coph)

# Ward clustering
env.de.ward.coph <- cophenetic(env.de.ward)
cor(env.de, env.de.ward.coph)

# UPGMA > Complete > Ward > Single


# 4. Comparison dendograms by Gower distance
# LOWEST VALUE = BEST DENDOGRAM
(gow.dist.single <- sum((env.de - env.de.single.coph) ^ 2))
(gow.dist.comp <- sum((env.de - env.de.comp.coph) ^ 2))
(gow.dist.UPGMA <- sum((env.de - env.de.UPGMA.coph) ^ 2))
(gow.dist.ward <- sum((env.de - env.de.ward.coph) ^ 2))

# UPGMA >>> Complete > Single > Ward


# Considering both results from points 3 and 4,
# we consider the UPGMA as the best dendogram option.

plot(env.de.UPGMA, hang=-1, labels = unique(my_data$Site),  main = "Euclidean - UPGMA")

png(file = "Figures/1 GRADIENT/cluster_by_ENV_withoutAltitude.png", width = 700, height = 600)


# # > By climatic vars. WITHOUT SP06 ----
# library(ade4) 
# library(vegan)  
# library(gclus) 
# library(cluster)
# library(RColorBrewer)  
# library(labdsv)
# library(leaflet)
# 
# # 1. Create dissimilarity matrix:
# # https://stats.stackexchange.com/questions/80377/which-distance-to-use-e-g-manhattan-euclidean-bray-curtis-etc
# 
# ## Trec les repliques per site pq totes tenen igual valors
# dummy2 <- as.data.frame(unique(dummy))
# dummy2 <- dummy2[-c(6),]
# 
# boxplot(dummy2)
# # Vars. with different units and scales
# # Standardization its needed
# 
# # If not, MAP and altitude would have much more
# # influence on the distance matrix only because of
# # their magnitude and variation in absolut terms.
# 
# env.std <- scale(dummy2)
# boxplot(env.std)
# # Nicer
# 
# # Using Euclidean distance:
# env.de <- dist(env.std)
# my_data2 <- subset(my_data,Site != "SP06")
# attr(env.de, "Labels") <- unique(my_data2$Site)
# env.de
# 
# 
# # 2. Analysis of hierarchical cluster
# par(mfrow = c(2, 2))
# 
# # Compute single linkage agglomerative clustering
# env.de.single <- hclust(env.de, method = "single")
# plot(env.de.single, 
#      labels = unique(my_data2$Site), 
#      main = "Euclidean - Single linkage")
# 
# # Compute complete-linkage agglomerative clustering
# env.de.complete <- hclust(env.de, method = "complete")
# plot(env.de.complete, 
#      labels = unique(my_data2$Site), 
#      main = "Euclidean - Complete linkage")
# 
# # Compute UPGMA agglomerative clustering
# env.de.UPGMA <- hclust(env.de, method = "average")
# plot(env.de.UPGMA, 
#      labels = unique(my_data2$Site), 
#      main = "Euclidean - UPGMA")
# 
# # Compute Ward’s Minimum Variance Clustering
# env.de.ward <- hclust(env.de, method = "ward.D2")
# plot(env.de.ward, 
#      labels = unique(my_data2$Site), 
#      main = "Euclidean - Ward")
# dev.off()
# 
# # 3. Comparison dendograms by cophenetic distance
# # HIGHEST CORRELATION = BEST DENDOGRAM
# 
# # Single linkage clustering
# env.de.single.coph <- cophenetic(env.de.single)
# cor(env.de, env.de.single.coph)
# 
# # Complete linkage clustering
# env.de.comp.coph <- cophenetic(env.de.complete)
# cor(env.de, env.de.comp.coph)
# 
# # Average clustering
# env.de.UPGMA.coph <- cophenetic(env.de.UPGMA)
# cor(env.de, env.de.UPGMA.coph)
# 
# # Ward clustering
# env.de.ward.coph <- cophenetic(env.de.ward)
# cor(env.de, env.de.ward.coph)
# 
# 
# # 4. Comparison dendograms by Gower distance
# # LOWEST VALUE = BEST DENDOGRAM
# (gow.dist.single <- sum((env.de - env.de.single.coph) ^ 2))
# (gow.dist.comp <- sum((env.de - env.de.comp.coph) ^ 2))
# (gow.dist.UPGMA <- sum((env.de - env.de.UPGMA.coph) ^ 2))
# (gow.dist.ward <- sum((env.de - env.de.ward.coph) ^ 2))






# > By Microbial community composition ----
# >> 16S ----
library(ade4) 
library(vegan)  
library(gclus) 
library(cluster)
library(RColorBrewer)  
library(labdsv)
library(leaflet)

library(readxl)
AbuPhyl <- read_excel("C:/Users/Anna/OneDrive - Universitat de Girona/GRADCATCH/ANALISIS/amplicon_sequencing_2021_JD/16S/AbuPhyl.xlsx")


# AbuPhyl_16S <- AbuPhyl_16S[,c(4, 61:63, 8, 66,67)]
# 
# AbuPhyl_16S$Altitude[AbuPhyl_16S$Site == '02'] <- '664.5'
# AbuPhyl_16S$Altitude <- as.numeric(AbuPhyl_16S$Altitude)
# AbuPhyl_16S$MAP <- as.numeric(AbuPhyl_16S$MAP)
# AbuPhyl_16S$MAT <- as.numeric(AbuPhyl_16S$MAT)
# AbuPhyl_16S$AI <- as.numeric(AbuPhyl_16S$AI)

AbuPhyl_16S <- AbuPhyl[,c(4,66,67)]

library(reshape2)
AbuPhyl_16S <- dcast(AbuPhyl_16S,Site~variable)

AbuPhyl_16S <- AbuPhyl_16S[,-c(1)]

# 1. Create dissimilarity matrix:

boxplot(AbuPhyl_16S)
# Vars. with same units but different scales
# a transformation might be needed

env.std <- log(AbuPhyl_16S+1)
boxplot(env.std)
# Nicer

# Using Bray-Curtis distance:
env.de <- vegdist(env.std, method = 'bray')
attr(env.de, "Labels") <- unique(my_data$Site)
env.de


# 2. Analysis of hierarchical cluster
par(mfrow = c(2, 2))

# Compute single linkage agglomerative clustering
env.de.single <- hclust(env.de, method = "single")
plot(env.de.single, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Single linkage")

# Compute complete-linkage agglomerative clustering
env.de.complete <- hclust(env.de, method = "complete")
plot(env.de.complete, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Complete linkage")

# Compute UPGMA agglomerative clustering
env.de.UPGMA <- hclust(env.de, method = "average")
plot(env.de.UPGMA, 
     labels = unique(my_data$Site), 
     main = "Euclidean - UPGMA")

# Compute Ward’s Minimum Variance Clustering
env.de.ward <- hclust(env.de, method = "ward.D2")
plot(env.de.ward, 
     labels = unique(my_data$Site), 
     main = "Euclidean - Ward")
par(mfrow = c(1,1))

# 3. Comparison dendograms by cophenetic distance
# HIGHEST CORRELATION = BEST DENDOGRAM

# Single linkage clustering
env.de.single.coph <- cophenetic(env.de.single)
cor(env.de, env.de.single.coph)

# Complete linkage clustering
env.de.comp.coph <- cophenetic(env.de.complete)
cor(env.de, env.de.comp.coph)

# Average clustering
env.de.UPGMA.coph <- cophenetic(env.de.UPGMA)
cor(env.de, env.de.UPGMA.coph)

# Ward clustering
env.de.ward.coph <- cophenetic(env.de.ward)
cor(env.de, env.de.ward.coph)

# UPGMA > Complete > Ward > Single


# 4. Comparison dendograms by Gower distance
# LOWEST VALUE = BEST DENDOGRAM
(gow.dist.single <- sum((env.de - env.de.single.coph) ^ 2))
(gow.dist.comp <- sum((env.de - env.de.comp.coph) ^ 2))
(gow.dist.UPGMA <- sum((env.de - env.de.UPGMA.coph) ^ 2))
(gow.dist.ward <- sum((env.de - env.de.ward.coph) ^ 2))

# UPGMA >> Complete >> Single >>> Ward

# Considering both results from points 3 and 4,
# we consider the UPGMA as the best dendogram option.

plot(env.de.UPGMA, hang=-1, labels = unique(my_data$Site),  main = "Euclidean - UPGMA")

png(file = "Figures/1 GRADIENT/cluster_by_16S.png", width = 700, height = 600)




#.----
# Feature Selection (CARET PACKAGE) ----
# https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
# https://bookdown.org/rehk/stm1001_dsm_t1_introduction_to_machine_learning_in_r/machine-learning-in-r-using-the-caret-package.html

# load the library
library(mlbench)
library(caret)

# > Remove Redundant Features ----
# ONLY FOR INDEPENDENT VARIABLES, WITH DUMMIES ***********
set.seed(7)
nearZeroVar(duindepe, saveMetrics = TRUE)
comboInfo <- findLinearCombos(duindepe)
comboInfo

duindepe_cor <- cor(duindepe)
duindepe_highlyCor <- findCorrelation(duindepe_cor, cutoff = 0.75)
print(duindepe_highlyCor)

clean_duindepe <- duindepe[,-c(7,2,5,6,27,3,23,1,10,8,12,29,17,11,32,33,34)]

clean_duindepe$Site <- my_data$Site
clean_duindepe <- clean_duindepe[, c("Site",  
                              names(clean_duindepe)[names(clean_duindepe) != "Site"])] 

rm(duindepe_cor, comboInfo)

model <- lm(cbind(depe$CO2_dark+depe$CH4_ave+depe$N2O+depe$BB+
                    depe$FB+depe$chla+depe$chlb+depe$carotene+depe$EPS+
                    depe$EPS+depe$alpha+depe$beta+depe$beta+depe$xyl+depe$cbh+
                    depe$gla+depe$fos+depe$leu+depe$phe+depe$X16S+depe$ITS2+
                    depe$mcrA+depe$pmoA+depe$nifH+depe$nifH+depe$AOA+depe$AOB+
                    depe$qnorB+depe$nosZ+depe$phoD+depe$Respiration)~clean_duindepe$altitude+
              clean_duindepe$pH+clean_duindepe$C_N+
              clean_duindepe$NH4+clean_duindepe$PO43+clean_duindepe$SO42+clean_duindepe$Silt+
              clean_duindepe$Clay+
              clean_duindepe$Litter+clean_duindepe$L_TC+clean_duindepe$L_TN+
              clean_duindepe$s350.400+clean_duindepe$SR+
              clean_duindepe$E2.E3+clean_duindepe$E4.E6+
              clean_duindepe$BIX+clean_duindepe$Peak_A+clean_duindepe$Peak_B+clean_duindepe$HIX
)

car::vif(model)
# vif_values <- vif(model)
# barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")
#CONFIRMAMOS QUE NINGUNA VARIABLES INDEPENDIENTE TIENE MULTICOLINEALIDAD (VIF>5)

model2 <- lm(cbind(depe$CO2_dark+depe$CH4_ave+depe$N2O+depe$BB+
                     depe$FB+depe$chla+depe$chlb+depe$carotene+depe$EPS+
                     depe$EPS+depe$alpha+depe$beta+depe$beta+depe$xyl+depe$cbh+
                     depe$gla+depe$fos+depe$leu+depe$phe+depe$X16S+depe$ITS2+
                     depe$mcrA+depe$pmoA+depe$nifH+depe$nifH+depe$AOA+depe$AOB+
                     depe$qnorB+depe$nosZ+depe$phoD+depe$Respiration)~clean_duindepe$altitude+
               clean_duindepe$pH+clean_duindepe$C_N+
               clean_duindepe$NH4+clean_duindepe$PO43+clean_duindepe$SO42+clean_duindepe$Silt+
               clean_duindepe$Clay+
               clean_duindepe$Litter+clean_duindepe$L_TC+clean_duindepe$L_TN+clean_duindepe$SR+
               clean_duindepe$E2.E3+clean_duindepe$E4.E6+
               clean_duindepe$BIX+clean_duindepe$Peak_A+clean_duindepe$Peak_B+clean_duindepe$HIX
)
car::vif(model2)

model3 <- lm(cbind(depe$CO2_dark+depe$CH4_ave+depe$N2O+depe$BB+
                     depe$FB+depe$chla+depe$chlb+depe$carotene+depe$EPS+
                     depe$EPS+depe$alpha+depe$beta+depe$beta+depe$xyl+depe$cbh+
                     depe$gla+depe$fos+depe$leu+depe$phe+depe$X16S+depe$ITS2+
                     depe$mcrA+depe$pmoA+depe$nifH+depe$nifH+depe$AOA+depe$AOB+
                     depe$qnorB+depe$nosZ+depe$phoD+depe$Respiration)~clean_duindepe$altitude+
               clean_duindepe$pH+clean_duindepe$C_N+
               clean_duindepe$NH4+clean_duindepe$PO43+clean_duindepe$SO42+clean_duindepe$Silt+
               clean_duindepe$Clay+
               clean_duindepe$Litter+clean_duindepe$L_TC+clean_duindepe$L_TN+clean_duindepe$SR+
               clean_duindepe$E2.E3+clean_duindepe$E4.E6+
               clean_duindepe$BIX+clean_duindepe$Peak_A+clean_duindepe$HIX
)
car::vif(model3)

model4 <- lm(cbind(depe$CO2_dark+depe$CH4_ave+depe$N2O+depe$BB+
                     depe$FB+depe$chla+depe$chlb+depe$carotene+depe$EPS+
                     depe$EPS+depe$alpha+depe$beta+depe$beta+depe$xyl+depe$cbh+
                     depe$gla+depe$fos+depe$leu+depe$phe+depe$X16S+depe$ITS2+
                     depe$mcrA+depe$pmoA+depe$nifH+depe$nifH+depe$AOA+depe$AOB+
                     depe$qnorB+depe$nosZ+depe$phoD+depe$Respiration)~clean_duindepe$altitude+
               clean_duindepe$pH+clean_duindepe$C_N+
               clean_duindepe$NH4+clean_duindepe$PO43+clean_duindepe$SO42+clean_duindepe$Silt+
               clean_duindepe$Clay+
               clean_duindepe$Litter+clean_duindepe$L_TC+clean_duindepe$L_TN+clean_duindepe$SR+
               clean_duindepe$E2.E3+clean_duindepe$E4.E6+clean_duindepe$Peak_A+clean_duindepe$HIX
)
car::vif(model4)

model5 <- lm(cbind(depe$CO2_dark+depe$CH4_ave+depe$N2O+depe$BB+
                     depe$FB+depe$chla+depe$chlb+depe$carotene+depe$EPS+
                     depe$EPS+depe$alpha+depe$beta+depe$beta+depe$xyl+depe$cbh+
                     depe$gla+depe$fos+depe$leu+depe$phe+depe$X16S+depe$ITS2+
                     depe$mcrA+depe$pmoA+depe$nifH+depe$nifH+depe$AOA+depe$AOB+
                     depe$qnorB+depe$nosZ+depe$phoD+depe$Respiration)~clean_duindepe$altitude+
               clean_duindepe$C_N+
               clean_duindepe$NH4+clean_duindepe$PO43+clean_duindepe$SO42+clean_duindepe$Silt+
               clean_duindepe$Clay+
               clean_duindepe$Litter+clean_duindepe$L_TC+clean_duindepe$L_TN+clean_duindepe$SR+
               clean_duindepe$E2.E3+clean_duindepe$E4.E6+clean_duindepe$Peak_A+clean_duindepe$HIX
)
car::vif(model5)


clean_duindepe <- clean_duindepe[,-c(3,13,17,19)]

rm(duindepe,dummy,model,model2,model3,model4,model5)


# JUST WITH THE INDEPENDENTS, WITHOUT DUMMIES *****

set.seed(7)
nearZeroVar(indepe, saveMetrics = TRUE)
comboInfo <- findLinearCombos(indepe)
comboInfo

indepe_cor <- cor(indepe)
indepe_highlyCor <- findCorrelation(indepe_cor, cutoff = 0.75)
print(indepe_highlyCor)

clean_indepe <- indepe[,-c(3,1,23,19,6,4,8,13,25,7,28,29,30)]

clean_indepe$Site <- my_data$Site
clean_indepe <- clean_indepe[, c("Site",  
                                     names(clean_indepe)[names(clean_indepe) != "Site"])] 

rm(indepe_cor, comboInfo)

model <- lm(cbind(depe$CO2_dark+depe$CH4_ave+depe$N2O+depe$BB+
                    depe$FB+depe$chla+depe$chlb+depe$carotene+depe$EPS+
                    depe$EPS+depe$alpha+depe$beta+depe$beta+depe$xyl+depe$cbh+
                    depe$gla+depe$fos+depe$leu+depe$phe+depe$X16S+depe$ITS2+
                    depe$mcrA+depe$pmoA+depe$nifH+depe$nifH+depe$AOA+depe$AOB+
                    depe$qnorB+depe$nosZ+depe$phoD+depe$Respiration)~clean_indepe$Water_activity+
              clean_indepe$pH+clean_indepe$C_N+
              clean_indepe$NH4+clean_indepe$PO43+clean_indepe$SO42+clean_indepe$Silt+
              clean_indepe$Clay+
              clean_indepe$Litter+clean_indepe$L_TC+clean_indepe$L_TN+
              clean_indepe$s350.400+clean_indepe$SR+
              clean_indepe$E2.E3+clean_indepe$E4.E6+
              clean_indepe$BIX+clean_indepe$Peak_A+clean_indepe$Peak_B+clean_indepe$HIX
)
car::vif(model)

model2 <- lm(cbind(depe$CO2_dark+depe$CH4_ave+depe$N2O+depe$BB+
                    depe$FB+depe$chla+depe$chlb+depe$carotene+depe$EPS+
                    depe$EPS+depe$alpha+depe$beta+depe$beta+depe$xyl+depe$cbh+
                    depe$gla+depe$fos+depe$leu+depe$phe+depe$X16S+depe$ITS2+
                    depe$mcrA+depe$pmoA+depe$nifH+depe$nifH+depe$AOA+depe$AOB+
                    depe$qnorB+depe$nosZ+depe$phoD+depe$Respiration)~clean_indepe$Water_activity+
              clean_indepe$pH+clean_indepe$C_N+
              clean_indepe$NH4+clean_indepe$PO43+clean_indepe$SO42+clean_indepe$Silt+
              clean_indepe$Clay+
              clean_indepe$Litter+clean_indepe$L_TC+clean_indepe$L_TN+
              clean_indepe$SR+
              clean_indepe$E2.E3+clean_indepe$E4.E6+
              clean_indepe$BIX+clean_indepe$Peak_A+clean_indepe$Peak_B+clean_indepe$HIX
)
car::vif(model2)

model3 <- lm(cbind(depe$CO2_dark+depe$CH4_ave+depe$N2O+depe$BB+
                     depe$FB+depe$chla+depe$chlb+depe$carotene+depe$EPS+
                     depe$EPS+depe$alpha+depe$beta+depe$beta+depe$xyl+depe$cbh+
                     depe$gla+depe$fos+depe$leu+depe$phe+depe$X16S+depe$ITS2+
                     depe$mcrA+depe$pmoA+depe$nifH+depe$nifH+depe$AOA+depe$AOB+
                     depe$qnorB+depe$nosZ+depe$phoD+depe$Respiration)~clean_indepe$Water_activity+
               clean_indepe$pH+clean_indepe$C_N+
               clean_indepe$NH4+clean_indepe$PO43+clean_indepe$SO42+clean_indepe$Silt+
               clean_indepe$Clay+
               clean_indepe$Litter+clean_indepe$L_TC+clean_indepe$L_TN+
               clean_indepe$SR+
               clean_indepe$E2.E3+clean_indepe$E4.E6+
               clean_indepe$Peak_A+clean_indepe$Peak_B+clean_indepe$HIX
)
car::vif(model3)

model4 <- lm(cbind(depe$CO2_dark+depe$CH4_ave+depe$N2O+depe$BB+
                     depe$FB+depe$chla+depe$chlb+depe$carotene+depe$EPS+
                     depe$EPS+depe$alpha+depe$beta+depe$beta+depe$xyl+depe$cbh+
                     depe$gla+depe$fos+depe$leu+depe$phe+depe$X16S+depe$ITS2+
                     depe$mcrA+depe$pmoA+depe$nifH+depe$nifH+depe$AOA+depe$AOB+
                     depe$qnorB+depe$nosZ+depe$phoD+depe$Respiration)~clean_indepe$Water_activity+
               clean_indepe$pH+clean_indepe$C_N+
               clean_indepe$NH4+clean_indepe$PO43+clean_indepe$SO42+clean_indepe$Silt+
               clean_indepe$Clay+
               clean_indepe$Litter+clean_indepe$L_TC+clean_indepe$L_TN+
               clean_indepe$SR+
               clean_indepe$E2.E3+clean_indepe$E4.E6+
               clean_indepe$Peak_A+clean_indepe$HIX
)
car::vif(model4)

model5 <- lm(cbind(depe$CO2_dark+depe$CH4_ave+depe$N2O+depe$BB+
                     depe$FB+depe$chla+depe$chlb+depe$carotene+depe$EPS+
                     depe$EPS+depe$alpha+depe$beta+depe$beta+depe$xyl+depe$cbh+
                     depe$gla+depe$fos+depe$leu+depe$phe+depe$X16S+depe$ITS2+
                     depe$mcrA+depe$pmoA+depe$nifH+depe$nifH+depe$AOA+depe$AOB+
                     depe$qnorB+depe$nosZ+depe$phoD+depe$Respiration)~clean_indepe$Water_activity+
               clean_indepe$C_N+
               clean_indepe$NH4+clean_indepe$PO43+clean_indepe$SO42+clean_indepe$Silt+
               clean_indepe$Clay+
               clean_indepe$Litter+clean_indepe$L_TC+clean_indepe$L_TN+
               clean_indepe$SR+
               clean_indepe$E2.E3+clean_indepe$E4.E6+
               clean_indepe$Peak_A+clean_indepe$HIX
)
car::vif(model5)


clean_indepe <- clean_indepe[,-c(3,13,17,19)]

rm(indepe,indepe_cor,comboInfo,model,model2,model3,model4,model5)





# #WITH DEPENDENT AND INDEPENDENT VARIABLES ALTOGETHER ***********
# # ensure the results are repeatable
# set.seed(7)
# 
# #1. Remove zero- and near zero-variance predictors
# nearZeroVar(data, saveMetrics = TRUE)
# # There are none in data so it returns an empty vector (integer(0)).
# # # nzv <- nearZeroVar(data)
# # # data_filtered <- data[,nzv] #not needed because there is no variables near 0 var.
# 
# #2. Remove linear dependencies
# comboInfo <- findLinearCombos(data)
# data[,-comboInfo$remove]
# 
# #3.Remove correlated predictors
# data_cor <- cor(data)
# data_highlyCor <- findCorrelation(data_cor, cutoff = 0.75)
# print(data_highlyCor)
# #It gives back the number of the columns from the 
# # database USED FOR THE CORRELATION (data_cor)
# # that should be removed, as they are highly correlated with other variables
# # CHECK THAT CUTOFF LEVEL IS ARBITRARY!
# 
# # Remember to check that you are removing the correct columns.
# # findCorrelation gives the columns to be removed from the dataset used
# # on the correlation (data_cor), not the original dataset (data)
# clean_data <- data[,-c(2,39,46,6,9,11,7,4,42,5,43,47,25,3,1,10,26,27,29,16,31)]
# clean_data <- data[,-c(6,2,39,46,42,9,4,43,11,47,7,5,3,1,25,55,51,10,26,27,57,29,16,60,61,63,32)]
# 
# clean_data$Site <- my_data$Site
# clean_data <- clean_data[ , c("Site",  
#                                  names(clean_data)[names(clean_data) != "Site"])] 
# 








# > Rank Features by Importance
# 
# # ensure results are repeatable
# set.seed(7)
# 
# clean_data$Site <- factor(clean_data$Site, levels = c("SP08", "SP01", "SP02",
#                                                     "SP07","SP06","SP03",
#                                                     "SP12", "SP11", "SP04",
#                                                     "SP09", "SP10","SP05"))
# 
# X_clean_data <- clean_data[,c(1,2:11,30)]
# Y_clean_data <- clean_data[,c(1,12:29)]
# 
# # prepare training scheme
# control <- trainControl(method="repeatedcv", number=10, repeats=3)
# 
# #I plot them separately because if not it is difficult to interpret
# #***FOR INDEPENDENT VARS.****
# # train the model
# TG = expand.grid(k=1:3,size=seq(5,20,by=5))
# model <- train(Site~., data=X_clean_data, method="lvq", preProcess="scale", trControl=control, tuneGrid=TG)
# 
# # estimate variable importance
# importance <- varImp(model, scale=FALSE)
# 
# # summarize importance
# print(importance)
# 
# # plot importance
# plot(importance)
# 
# 
# #***FOR DEPENDENT VARS.****
# # train the model
# TG = expand.grid(k=1:3,size=seq(5,20,by=5))
# model <- train(Site~., data=Y_clean_data, method="lvq", preProcess="scale", trControl=control, tuneGrid=TG)
# 
# # estimate variable importance
# importance <- varImp(model, scale=FALSE)
# 
# # summarize importance
# print(importance)
# 
# # plot importance
# plot(importance)

#.----




# Regressions ----
library(car)
library(MASS)

#CO2_dark

m1 <- lm(data$CO2_dark ~ ., clean_indepe[,-1])
summary(m1)
stepb <- stepAIC(m1, direction="backward")

m2 <- lm(data$Respiration ~ ., clean_indepe[,-1])
summary(m2)
stepb <- stepAIC(m2, direction="backward")


# clean_idepe_depe <- cbind(clean_indepe, depe)
# 
# varlist <- c(colnames(depe))
# models <- lapply(varlist, function(x) {
#   lm(substitute(i ~ ., list(i = as.name(x))), data = clean_idepe_depe[,-1])
# })
# 
# lapply(models, glance)



# .----
# PERMANOVA ----

# https://www.youtube.com/watch?v=1ETBgbXl-BM&ab_channel=RiffomonasProject
site <- data.frame(my_data[,c(2)]) 
site <- setNames(site, c("Site"))

perma_data <- depe
perma_data2 <- cbind(site, depe)

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


# .----



# STATISTICAL TESTS ----
# 1. MANOVA for dependent variables ----
dummy_variables = as.data.frame((my_data[,c(2,5:7,27:29,76)]))
dummy_y <- my_data[,c(2,5:7,27:29,76, 35:38,40:64,74)]
y_variables <- my_data[,c(2,35:38,40:64,74)]

# 1.1 MANOVA assumptions ----
library(rstatix)
library(broom)

# > Outliers ####
out = as.data.frame(y_variables %>% group_by(Site) %>% identify_outliers("CO2_dark"))
out
# There are outliers but we are keeping them for the analysis

# > Multivariate outliers ####
# Absence of multivariate outliers is checked by
# assessing Mahalanobis Distances among the dependent vars.
y_variables_2       = scale(y_variables[,-c(1)],center = FALSE)
mahal               = mahalanobis(y_variables_2, colMeans(y_variables_2), cov(y_variables_2))
p_val               = pchisq(mahal, df=3, lower.tail=FALSE)
y_variables_2       = as.data.frame(cbind(y_variables_2,mahal,p_val))
# There are outliers (see p-values) but we are keeping them for the analysis


# > Normality ####
y_variables_2      = as.data.frame(cbind(y_variables$Site,y_variables_2))

ggqqplot(y_variables_2, "CO2_dark", facet.by = "y_variables$Site",
         ylab = "CO2_dark", ggtheme = theme_bw())
ggqqplot(y_variables_2, "CO2_light", facet.by = "y_variables$Site",
         ylab = "CO2_light", ggtheme = theme_bw())
ggqqplot(y_variables_2, "CH4_ave", facet.by = "y_variables$Site",
         ylab = "CH4_ave", ggtheme = theme_bw())
ggqqplot(y_variables_2, "N2O", facet.by = "y_variables$Site",
         ylab = "N2O", ggtheme = theme_bw())
ggqqplot(y_variables_2, "BB", facet.by = "y_variables$Site",
         ylab = "BB", ggtheme = theme_bw())
ggqqplot(y_variables_2, "FB", facet.by = "y_variables$Site",
         ylab = "FB", ggtheme = theme_bw())
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
ggqqplot(y_variables_2, "X16S", facet.by = "y_variables$Site",
         ylab = "X16S", ggtheme = theme_bw())
ggqqplot(y_variables_2, "ITS2", facet.by = "y_variables$Site",
         ylab = "ITS2", ggtheme = theme_bw())
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
ggqqplot(y_variables_2, "ShannonEEA", facet.by = "y_variables$Site",
         ylab = "ShannonEEA", ggtheme = theme_bw())
# Most of the variables seem normally distributed and MANOVA is robust with
# slightly violations of normality. 


# > Multivariate normality ####
mshapiro_test((y_variables_2[,2:31]))
# Not normal


# > Multicollinearity ####
cor.mat <- y_variables[,2:31] %>% cor_mat()
cor.mat %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot(label = FALSE)
cor.mat %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot(label = TRUE)
#  A cross at each matrix position where the correlation coefficient is not significant.
# If i put TRUE on label the p-values will appear.

# New selected variables
# Deleted variables:
# alpha, CO2_light, chla, carotene, X16S, nosZ, nifH, qnorB, ShannonEEA
y_variables_3 = as.data.frame((y_variables[,-c(12, 3, 8,10, 20, 24, 28, 27, 31)]))
cor.mat <- y_variables_3[,2:22] %>% cor_mat()
cor.mat %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot(label = FALSE)
cor.mat %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot(label = TRUE)
# No more dependent variables correlated >80%


# > Linearity assumption ####
library(GGally)
linear_ass <- y_variables_3 %>% group_by(Site) %>% doo(~ggpairs(.) + theme_bw(), result = "plots")
linear_ass$plots

# > Homogeneity of covariances ####
# Balance design and we could continue with the analysis. We will use the 
# Pillai’s multivariate statistic as a more robust metric

# > Homogeneity of variance assumption ####
a = as.data.frame(colnames(y_variables_3))
y_variables_3 %>% gather(key = "variable", value = "value", a[2:22,]) %>%
  group_by(variable) %>% levene_test(sqrt(value) ~ as.factor(Site))
# Squaring the variables solves almost completely the homogeneity of variances


# 1.2 MANOVA test ####
mod <- manova(as.matrix(y_variables_3[,2:22]) ~ (my_data$AI))
summary(mod,summary=TRUE)

data <- cbind(dummy_variables,y_variables_3[,-c(1)])

attach(data)
mod <- manova(as.matrix(y_variables_3[,2:22]) ~ AI+Site)
summary(mod,summary=TRUE)
detach(data)

str(data)
data$Site <- factor(data$Site)
data$MAP <- factor(data$MAP)
data$MAT <- factor(data$MAT)
data$AI <- factor(data$AI)
data$Sand <- factor(data$Sand)
data$Silt <- factor(data$Silt)
data$Clay <- factor(data$Clay)
data$altitude <- factor(data$altitude)
str(data)

mod <- manova(as.matrix(y_variables_3[,2:22]) ~ data$Site)
summary(mod,summary=TRUE)

mod <- manova(as.matrix(y_variables_3[,2:22]) ~ data$Clay + data$Silt + data$altitude)
summary(mod,summary=TRUE)

mod <- manova(as.matrix(y_variables_3[,2:22]) ~ data$AI + data$Site + data$Sand + data$Clay + data$Silt + data$altitude)
summary(mod,summary=TRUE)





# 
# #packages
# library("car")
# library("MASS")
# #
# data <- cbind(dummy_variables,y_variables_3[,-c(1)])
# 
# str(data)
# data$Site <- factor(data$Site)
# data$MAP <- factor(data$MAP)
# data$MAT <- factor(data$MAT)
# data$AI <- factor(data$AI)
# data$Sand <- factor(data$Sand)
# data$Silt <- factor(data$Silt)
# data$Clay <- factor(data$Clay)
# data$altitude <- factor(data$altitude)
# str(data)
# 
# #
# attach(data)
# mod1 <- lm(cbind(CO2_dark, CH4_ave, N2O,BB, FB, chlb, EPS, beta, xyl, cbh, gla, fos, leu, phe,
#            ITS2, mcrA, pmoA, AOA, AOB, phoD, Respiration) ~ Site+MAP+MAT+AI+altitude+Clay+Sand+Silt)
# 
# summary(mod1)
# 
# Manova(mod1, multivariate = T, type=c("III"), test=("Wilks"))
# 
# 
# 
# 
# 
# 
# mod1 <- lm(cbind(CO2_dark, CH4_ave, N2O,BB, FB, chlb, EPS, beta, xyl, cbh, gla, fos, leu, phe,
#                  ITS2, mcrA, pmoA, AOA, AOB, phoD, Respiration) ~ Site)
# 
# summary(mod1)
# 
# Manova(mod1, multivariate = T, type=c("III"), test=("Wilks"))
# 
# detach(data)
# 
# 
# 
# 
# 
# 
# 
# #
# DV <- cbind(y_variables_3$CO2_dark, y_variables_3$CH4_ave, y_variables_3$N2O,
#             y_variables_3$BB, y_variables_3$FB, y_variables_3$chlb,
#             y_variables_3$EPS, y_variables_3$beta, y_variables_3$xyl,
#             y_variables_3$cbh, y_variables_3$gla, y_variables_3$fos,
#             y_variables_3$leu, y_variables_3$phe, y_variables_3$ITS2,
#             y_variables_3$mcrA, y_variables_3$pmoA, y_variables_3$AOA,
#             y_variables_3$AOB, y_variables_3$phoD, y_variables_3$Respiration)
# 
# # DV <- na.exclude(DV)
# 
# dummy_y$Site <- as.factor(dummy_y$Site)
# dummy_y$MAP <- as.factor(dummy_y$MAP)
# dummy_y$MAT <- as.factor(dummy_y$MAT)
# dummy_y$AI <- as.factor(dummy_y$AI)
# dummy_y$Sand <- as.factor(dummy_y$Sand)
# dummy_y$Clay <- as.factor(dummy_y$Clay)
# dummy_y$Silt <- as.factor(dummy_y$Silt)
# dummy_y$altitude <- as.factor(dummy_y$altitude)
# 
# 
# # output = lm(DV ~ Site*MAP*MAT*AI*Sand*Clay*Silt*altitude, data = dummy_y,
# #             contrast=list(Site=contr.sum, MAP=contr.sum, MAT=contr.sum,
# #                           AI=contr.sum, Sand=contr.sum, Clay=contr.sum,
# #                           Silt=contr.sum, altitude=contr.sum))
# 
# output = lm(DV ~ Site+MAP+MAT+AI+Sand+Clay+Silt+altitude, data = dummy_y,
#             contrast=list(Site=contr.sum, MAP=contr.sum, MAT=contr.sum,
#                           AI=contr.sum, Sand=contr.sum, Clay=contr.sum,
#                           Silt=contr.sum, altitude=contr.sum))
# 
# summary(output)
# library(car)
# manova_out <- Manova(output, type ="III")

# .----

# 2. Kruskal-Wallis (AI (Site) vs dependent vars.) ----
# > Enzymes ----
# >> Figure 3. Enzyme by bars and post-hoc letters ----
dummy_variables = as.data.frame((my_data[,c(2,5:7,27:29,76)]))
dummy_y <- my_data[,c(2,5:7,27:29,76, 35:38,40:64,74)]
y_variables <- my_data[,c(2,35:38,40:64,74)]

dummy_y$Site <- factor(dummy_y$Site)

require(agricolae)
library(tidyverse)

data <- my_data[,c(2,7,46:53)]

#IMPORT FUNCTION DATA SUMMARY 2

alpha <- data_summary2(data, varname="alpha", 
                       groupnames=c("Site", "AI"))
beta <- data_summary2(data, varname="beta", 
                      groupnames=c("Site", "AI"))
xyl <- data_summary2(data, varname="xyl", 
                     groupnames=c("Site", "AI"))
cbh <- data_summary2(data, varname="cbh", 
                     groupnames=c("Site", "AI"))
gla <- data_summary2(data, varname="gla", 
                     groupnames=c("Site", "AI"))
fos <- data_summary2(data, varname="fos", 
                     groupnames=c("Site", "AI"))
leu <- data_summary2(data, varname="leu", 
                     groupnames=c("Site", "AI"))
phe <- data_summary2(data, varname="phe", 
                     groupnames=c("Site", "AI"))
alpha$Site <- factor(alpha$Site, levels = c("SP08", "SP01", "SP02",
                                            "SP07","SP06","SP03",
                                            "SP12", "SP11", "SP04",
                                            "SP09", "SP10","SP05"))
beta$Site <- factor(beta$Site, levels = c("SP08", "SP01", "SP02",
                                          "SP07","SP06","SP03",
                                          "SP12", "SP11", "SP04",
                                          "SP09", "SP10","SP05"))
xyl$Site <- factor(xyl$Site, levels = c("SP08", "SP01", "SP02",
                                        "SP07","SP06","SP03",
                                        "SP12", "SP11", "SP04",
                                        "SP09", "SP10","SP05"))
cbh$Site <- factor(cbh$Site, levels = c("SP08", "SP01", "SP02",
                                        "SP07","SP06","SP03",
                                        "SP12", "SP11", "SP04",
                                        "SP09", "SP10","SP05"))
fos$Site <- factor(fos$Site, levels = c("SP08", "SP01", "SP02",
                                        "SP07","SP06","SP03",
                                        "SP12", "SP11", "SP04",
                                        "SP09", "SP10","SP05"))
gla$Site <- factor(gla$Site, levels = c("SP08", "SP01", "SP02",
                                        "SP07","SP06","SP03",
                                        "SP12", "SP11", "SP04",
                                        "SP09", "SP10","SP05"))
leu$Site <- factor(leu$Site, levels = c("SP08", "SP01", "SP02",
                                        "SP07","SP06","SP03",
                                        "SP12", "SP11", "SP04",
                                        "SP09", "SP10","SP05"))
phe$Site <- factor(phe$Site, levels = c("SP08", "SP01", "SP02",
                                        "SP07","SP06","SP03",
                                        "SP12", "SP11", "SP04",
                                        "SP09", "SP10","SP05"))
library(ggplot2)
plot.theme1 <- theme_classic() +
  theme(text=element_text(size=15),
        axis.title.x = element_text(size = rel(1.2), angle = 00, margin = margin(t=8)),
        axis.title.y = element_text(size = rel(1.2), angle = 90, margin = margin(t=8)),
        plot.title = element_text(size=22),
        legend.position = "none",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))

#alpha
k <- kruskal(dummy_y$alpha, dummy_y$Site, console = TRUE,
             p.adj=c("bonferroni"))
#p.adj = "none" is t-student
(t_comp <- k$means %>% 
    rownames_to_column(var = "Site") %>%
    rename(alpha = dummy_y.alpha) %>%
    as_tibble() %>% 
    left_join(as_tibble(k$groups), by = c("rank" = "dummy_y$alpha")))
alpha <- merge(alpha, t_comp[,c(1,3,11)], by = "Site")

palpha <-ggplot(alpha, aes(x=Site, y=alpha, fill=AI)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=alpha-sem, ymax=alpha+sem), width=.2)+
  ggtitle("alpha") +
  plot.theme1+
  scale_fill_AI(discrete = FALSE, palette = "Sites")+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1))) +
  geom_text(data = alpha, aes(x = Site, y = alpha+sem+0.001, label = groups), size = 6, color = "red")
palpha


#beta
k <- kruskal(dummy_y$beta, dummy_y$Site, console = TRUE,
             p.adj=c("bonferroni"))
#p.adj = "none" is t-student
(t_comp <- k$means %>% 
    rownames_to_column(var = "Site") %>%
    rename(beta = dummy_y.beta) %>%
    as_tibble() %>% 
    left_join(as_tibble(k$groups), by = c("rank" = "dummy_y$beta")))
beta <- merge(beta, t_comp[,c(1,3,11)], by = "Site")

pbeta <-ggplot(beta, aes(x=Site, y=beta, fill=AI)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=beta-sem, ymax=beta+sem), width=.2)+
  ggtitle("beta") +
  plot.theme1+
  scale_fill_AI(discrete = FALSE, palette = "Sites")+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))+
  geom_text(data = beta, aes(x = Site, y = beta+sem+0.007, label = groups), size = 6, color = "red")
pbeta


#cbh
k <- kruskal(dummy_y$cbh, dummy_y$Site, console = TRUE,
             p.adj=c("bonferroni"))
#p.adj = "none" is t-student
(t_comp <- k$means %>% 
    rownames_to_column(var = "Site") %>%
    rename(cbh = dummy_y.cbh) %>%
    as_tibble() %>% 
    left_join(as_tibble(k$groups), by = c("rank" = "dummy_y$cbh")))
cbh <- merge(cbh, t_comp[,c(1,3,11)], by = "Site")

pcbh <-ggplot(cbh, aes(x=Site, y=cbh, fill=AI)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=cbh-sem, ymax=cbh+sem), width=.2)+
  ggtitle("cbh") +
  plot.theme1+
  scale_fill_AI(discrete = FALSE, palette = "Sites")+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))+
  geom_text(data = cbh, aes(x = Site, y = cbh+sem+0.0005, label = groups), size = 6, color = "red")
pcbh



#xyl
k <- kruskal(dummy_y$xyl, dummy_y$Site, console = TRUE,
             p.adj=c("bonferroni"))
#p.adj = "none" is t-student
(t_comp <- k$means %>% 
    rownames_to_column(var = "Site") %>%
    rename(xyl = dummy_y.xyl) %>%
    as_tibble() %>% 
    left_join(as_tibble(k$groups), by = c("rank" = "dummy_y$xyl")))
xyl <- merge(xyl, t_comp[,c(1,3,11)], by = "Site")

pxyl <-ggplot(xyl, aes(x=Site, y=xyl, fill=AI)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=xyl-sem, ymax=xyl+sem), width=.2)+
  ggtitle("xyl") +
  plot.theme1+
  scale_fill_AI(discrete = FALSE, palette = "Sites")+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))+
  geom_text(data = xyl, aes(x = Site, y = xyl+sem+0.0014, label = groups), size = 6, color = "red")
pxyl



#gla
k <- kruskal(dummy_y$gla, dummy_y$Site, console = TRUE,
             p.adj=c("bonferroni"))
#p.adj = "none" is t-student
(t_comp <- k$means %>% 
    rownames_to_column(var = "Site") %>%
    rename(gla = dummy_y.gla) %>%
    as_tibble() %>% 
    left_join(as_tibble(k$groups), by = c("rank" = "dummy_y$gla")))
gla <- merge(gla, t_comp[,c(1,3,11)], by = "Site")

pgla <-ggplot(gla, aes(x=Site, y=gla, fill=AI)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=gla-sem, ymax=gla+sem), width=.2)+
  ggtitle("gla") +
  plot.theme1+
  scale_fill_AI(discrete = FALSE, palette = "Sites")+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))+
  geom_text(data = gla, aes(x = Site, y = gla+sem+0.008, label = groups), size = 6, color = "red")
pgla


#fos
k <- kruskal(dummy_y$fos, dummy_y$Site, console = TRUE,
             p.adj=c("bonferroni"))
#p.adj = "none" is t-student
(t_comp <- k$means %>% 
    rownames_to_column(var = "Site") %>%
    rename(fos = dummy_y.fos) %>%
    as_tibble() %>% 
    left_join(as_tibble(k$groups), by = c("rank" = "dummy_y$fos")))
fos <- merge(fos, t_comp[,c(1,3,11)], by = "Site")

pfos <-ggplot(fos, aes(x=Site, y=fos, fill=AI)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=fos-sem, ymax=fos+sem), width=.2)+
  ggtitle("fos") +
  plot.theme1+
  scale_fill_AI(discrete = FALSE, palette = "Sites")+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))+
  geom_text(data = fos, aes(x = Site, y = fos+sem+0.04, label = groups), size = 6, color = "red")
pfos

#leu
k <- kruskal(dummy_y$leu, dummy_y$Site, console = TRUE,
             p.adj=c("bonferroni"))
#p.adj = "none" is t-student
(t_comp <- k$means %>% 
    rownames_to_column(var = "Site") %>%
    rename(leu = dummy_y.leu) %>%
    as_tibble() %>% 
    left_join(as_tibble(k$groups), by = c("rank" = "dummy_y$leu")))
leu <- merge(leu, t_comp[,c(1,3,11)], by = "Site")

pleu <-ggplot(leu, aes(x=Site, y=leu, fill=AI)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=leu-sem, ymax=leu+sem), width=.2)+
  ggtitle("leu") +
  plot.theme1+
  scale_fill_AI(discrete = FALSE, palette = "Sites")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  labs(y = expression(paste("μmol AMC·" ~  g ~ DW^-1,"·" ~ h^-1)))+
  geom_text(data = leu, aes(x = Site, y = leu+sem+0.012, label = groups), size = 6, color = "red")
pleu


#phe
k <- kruskal(dummy_y$phe, dummy_y$Site, console = TRUE,
             p.adj=c("bonferroni"))
#p.adj = "none" is t-student
(t_comp <- k$means %>% 
    rownames_to_column(var = "Site") %>%
    rename(phe = dummy_y.phe) %>%
    as_tibble() %>% 
    left_join(as_tibble(k$groups), by = c("rank" = "dummy_y$phe")))
phe <- merge(phe, t_comp[,c(1,3,11)], by = "Site")

pphe <-ggplot(phe, aes(x=Site, y=phe, fill=AI)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=phe-sem, ymax=phe+sem), width=.2)+
  ggtitle("phe") +
  plot.theme1+
  scale_fill_AI(discrete = FALSE, palette = "Sites")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
  labs(y = expression(paste("μmol DIQC·" ~  g ~ DW^-1,"·" ~ h^-1)))+
  geom_text(data = phe, aes(x = Site, y = phe+sem+0.2, label = groups), size = 6, color = "red")
pphe


all <- ggarrange(palpha, pbeta, pxyl, pcbh, pfos, pgla, pleu, pphe,
                 nrow = 4, ncol = 2,
                 common.legend = TRUE,
                 legend = "right")
all
ggsave(path = "Figures/1 GRADIENT","enzymes3.png",
       width = 30, height = 20, dpi = 300)


# . ----
# IV REDUCTION & PLSR ----
# What variables to use in the MANOVA as covariables
x_variables = as.data.frame((my_data[,c(2,5:8,10,12:13,17:23,27:30)]))

??cor_mat
library(rstatix)
cor.mat <- x_variables[,2:19] %>% cor_mat()
cor.mat %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot(label = TRUE)

# < 0.8 criteria
# New selected variables:
x_variables_1 = as.data.frame((x_variables[,-c(2,3,6,8,9,10,11,16)]))
# Removed variables: MAP, pH, TOC, TC, TN, MAT, Water content
cor.mat <- x_variables_1[,2:11] %>% cor_mat()
cor.mat %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot(label = TRUE)


# < 0.70 criteria
# New selected variables:
x_variables_1 = as.data.frame((x_variables[,-c(2:3,5,6,7,8,9,10,11,16)]))
# Removed variables: MAP, MAT, Soil temp, pH,SOM, Water content, TOC, TC, TN, Sand
cor.mat <- x_variables_1[,2:9] %>% cor_mat()
cor.mat %>% cor_reorder() %>% pull_lower_triangle() %>% cor_plot(label = TRUE)


# VIF ####
# < 0.8 criteria
library(car)
cor.mat <- x_variables_1[,2:11] %>% cor_mat()
model <- lm(cbind(y_variables_3$CO2_dark+y_variables_3$chla+y_variables_3$carotene+
                    y_variables_3$EPS) 
            ~ x_variables_1$AI+x_variables_1$Soil_Temp+x_variables_1$SOM+
              x_variables_1$C_N+x_variables_1$NH4+x_variables_1$PO43+x_variables_1$SO42+x_variables_1$Silt+
              x_variables_1$Clay+x_variables_1$Litter)
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



#. ----






# .----
# FIGURES ----
# 1. Physicochemical plots ----

New_data <- my_data[,c(2,7,8,10,12:13,17:23,27:30)]

New_data <- setNames(New_data, c("Site","AI","Soiltemp","WaterContent","SOM","pH",
                                 "TOC", "TC", "TN",
                                 "CN","NH4","PO43","SO42",
                                 "Sand","Silt","Clay","Litter"))
New_data1 <- New_data %>%
  gather(Physio, values, c(3:6))
New_data2 <- New_data %>%
  gather(Physio, values, c(7:10))
New_data3 <- New_data %>%
  gather(Physio, values, c(11:13,17))
# New_data4 <- New_data %>%
#   gather(Physio, values, c(14:17))


New_data1 <- data_summary(New_data1, varname="values", 
                         groupnames=c("Site","AI","Physio"))
New_data2 <- data_summary(New_data2, varname="values", 
                          groupnames=c("Site","AI","Physio"))
New_data3 <- data_summary(New_data3, varname="values", 
                          groupnames=c("Site","AI","Physio"))
# New_data4 <- data_summary(New_data4, varname="values", 
                          # groupnames=c("Site","AI","Physio"))

New_data1$AI = as.numeric(New_data1$AI)
New_data2$AI = as.numeric(New_data2$AI)
New_data3$AI = as.numeric(New_data3$AI)
# New_data4$AI = as.numeric(New_data4$AI)

New_data1$Site <- factor(New_data1$Site, levels = c("SP08", "SP01", "SP02",
                                                "SP07","SP06","SP03",
                                                "SP12", "SP11", "SP04",
                                                "SP09", "SP10","SP05"))
New_data2$Site <- factor(New_data2$Site, levels = c("SP08", "SP01", "SP02",
                                                    "SP07","SP06","SP03",
                                                    "SP12", "SP11", "SP04",
                                                    "SP09", "SP10","SP05"))
New_data3$Site <- factor(New_data3$Site, levels = c("SP08", "SP01", "SP02",
                                                    "SP07","SP06","SP03",
                                                    "SP12", "SP11", "SP04",
                                                    "SP09", "SP10","SP05"))
# New_data4$Site <- factor(New_data4$Site, levels = c("SP08", "SP01", "SP02",
#                                                     "SP07","SP06","SP03",
#                                                     "SP12", "SP11", "SP04",
#                                                     "SP09", "SP10","SP05"))

mycolors2<-c("#427681","#3891A6","#9BBC79","#CCD263","#E5DD58",
             "#FDE74C", "#EC9E57", "#E3655B", "#DB5461","#D84652",
             "#7D1809", "#290500")


p1 <- ggplot(New_data1, aes(x=AI, y=values, color = Site)) +
  geom_point(size= 5, shape=16) +
  geom_pointrange(data = New_data1, aes(ymin=values-sd, ymax=values+sd), 
                  color = "black", stroke = 1, size = 1, shape = 1) +
  xlab("Aridity Index") +
  scale_color_manual(values = mycolors2)+
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, colour = "black")) +
  facet_wrap(.~ Physio , nrow = 1, scales = "free_y",
             labeller = labeller(Physio = c("pH" = "pH", 
                                            "Soiltemp" = "Soil Temp.",
                                            "SOM" = "SOM",
                                            "WaterContent" = "Water content"))) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black"))+
  theme(axis.text.y = element_text(size = 18, color = "black"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  # theme(legend.position = "none")+
  coord_cartesian(xlim = c(0.10,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))

p1

ggsave(path = "Figures/1 GRADIENT","SP_fig1_1.png", width = 20, height = 5, dpi = 300)


p2 <- ggplot(New_data2, aes(x=AI, y=values, color = Site)) +
  geom_point(size= 5, shape=16) +
  geom_pointrange(data = New_data2, aes(ymin=values-sd, ymax=values+sd), 
                  color = "black", stroke = 1, size = 1, shape = 1) +
  xlab("Aridity Index") +
  scale_color_manual(values = mycolors2)+
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, colour = "black")) +
  facet_wrap(.~ Physio , nrow = 1, scales = "free_y",
             labeller = labeller(Physio = c("CN" = "C/N", 
                                            "TC" = "TC",
                                            "TN" = "TN",
                                            "TOC" = "TOC"))) +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black"))+
  theme(axis.text.y = element_text(size = 18, color = "black"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  # theme(legend.position = "none")+
  coord_cartesian(xlim = c(0.10,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
p2
ggsave(path = "Figures/1 GRADIENT","SP_fig1_2.png", width = 20, height = 5, dpi = 300)


p3 <- ggplot(New_data3, aes(x=AI, y=values, color = Site)) +
  geom_point(size= 5, shape=16) +
  geom_pointrange(data = New_data3, aes(ymin=values-sd, ymax=values+sd), 
                  color = "black", stroke = 1, size = 1, shape = 1) +
  xlab("Aridity Index") +
  scale_color_manual(values = mycolors2)+
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x = element_text(size = 20, colour = "black")) +
  facet_wrap(.~ Physio , nrow = 1, scales = "free_y") +
  theme(strip.background =element_rect(fill="light grey")) +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(size = 18, angle = 0, color = "black"))+
  theme(axis.text.y = element_text(size = 18, color = "black"))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
  # theme(legend.position = "none")+
  coord_cartesian(xlim = c(0.10,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
p3
ggsave(path = "Figures/1 GRADIENT","SP_fig1_3.png", width = 20, height = 5, dpi = 300)



# p4 <- ggplot(New_data4, aes(x=AI, y=values, color = Site)) +
#   geom_point(size= 5, shape=16) +
#   geom_pointrange(data = New_data4, aes(ymin=values-sd, ymax=values+sd), 
#                   color = "black", stroke = 1, size = 1, shape = 1) +
#   xlab("Aridity Index") +
#   scale_color_manual(values = mycolors2)+
#   theme(axis.title.y=element_blank()) +
#   theme(axis.title.x = element_text(size = 20, colour = "black")) +
#   facet_wrap(.~ Physio , nrow = 1, scales = "free_y")+
#   theme(strip.background =element_rect(fill="light grey")) +
#   theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0)) +
#   theme(panel.background = element_blank()) +
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
#   theme(axis.text.x = element_text(size = 18, angle = 0, color = "black"))+
#   theme(axis.text.y = element_text(size = 18, color = "black"))+
#   theme(plot.margin = unit(c(0.5, 0.5, 0.3, 0.5), "cm")) + #top, right, bottom, left
#   # theme(legend.position = "none")+
#   coord_cartesian(xlim = c(0.10,1.4))+
#   scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
#                      labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
# p4
# ggsave(path = "Figures/1 GRADIENT","SP_fig1_4.png", width = 20, height = 5, dpi = 300)


# Granulometry data:
granu <- my_data[,c(2,24:26)]

colnames(granu) <- c('Site','Sand','Silt','Clay')

granu <- granu %>% 
  gather(key="Granulometry", value="Percentage", 2:4) %>% 
  unique()

granu$Site <- factor(granu$Site, levels = c("SP05", "SP10", "SP09",
                                            "SP04","SP11","SP12",
                                            "SP03", "SP06", "SP07",
                                            "SP02", "SP01","SP08"))

plot.theme2 <- theme_classic() +
  theme(text=element_text(size=15),
        axis.title.x = element_text(size = rel(1.2), angle = 00, margin = margin(t=8)),
        axis.title.y = element_text(size = rel(1.2), angle = 90, margin = margin(t=8)),
        plot.title = element_text(size=22),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))

ggplot(granu, aes(x = granu$Site, y = granu$Percentage, fill=granu$Granulometry)) + 
  coord_flip() + 
  plot.theme2 +
  geom_bar(stat = "identity") + 
  labs(y = "Percentage", x = "Sites")+
  scale_fill_manual(values = c("Silt" = "#fec44f",
                               "Sand" = "#ec7014",
                               "Clay"="#993404"))

ggsave(path = "Figures/1 GRADIENT","fig1_5.png",
       width = 7, height = 6, dpi = 300)


# **************************************************************************************

# 2. Physicochemical PCA (without AI) ----

data <- my_data[,c(2,8,10,12,13,17:23,27:30,85,87,90,92)]
for (i in which(sapply(data, is.numeric))) {
  for (j in which(is.na(data[, i]))) {
    data[j, i] <- mean(data[data[, "Site"] == data[j, "Site"], i],  na.rm = TRUE)
  }
}
library(dplyr)
data <- rename(data,
                  Temperature = Soil_Temp,
                  WC = Water_content,
                  CN = C_N)

data2 <- data[,c(-1)]

#By means:
data_mit <- data %>% 
  group_by(Site) %>%
  summarise_all("mean")
data_mit2 <- data_mit[,c(-1)]
site <- data_mit[,1]
site <- site$Site
site <- factor(site, levels = site)
all <- prcomp(na.omit(data_mit2), center = TRUE,
              scale. = TRUE) 
plot(all, type = "l")
plot(all)
summary(all)


library(factoextra)
# Eigenvalues
eig.val <- get_eigenvalue(all)
eig.val

# Results for Variables
res.var <- get_pca_var(all)
res.var$coord          # Coordinates
res.var$contrib # Contributions to the PCs
# # To add a dataframe with a table of OTU
contributions <- res.var$contrib[,1:2]
contributionsdf = as.data.frame(contributions[,1:2])
write.csv(contributionsdf, "Figures/1 GRADIENT/PCA_fsca_contributions.csv", row.names=TRUE)

res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(all)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

a <- fviz_contrib(all, choice = "var", axes = 1, top = 10)
b <- fviz_contrib(all, choice = "var", axes = 2, top = 10)
ggarrange(a,b, ncol=2) + theme_classic()

ggsave(path = "Figures/1 GRADIENT","PCA_fsca_contributions.png", width = 10, height = 6, dpi = 300)


mycolors2<-c("#427681","#3891A6","#9BBC79","#CCD263","#E5DD58",
             "#FDE74C", "#EC9E57", "#E3655B", "#DB5461","#D84652",
             "#7D1809", "#290500")
library(ggbiplot)
all_plot <- ggbiplot(all, obs.scale = 1, var.scale = 1, 
                     ellipse = TRUE, fill=site,
                     varname.adjust = 2.5,
                     circle = TRUE, alpha=0,
                     loadings.label.repel=TRUE) +
  theme_classic()+
  scale_fill_manual(values = mycolors2,
                    breaks=c('SP08', 'SP01', 'SP02','SP07',
                             'SP06', 'SP03','SP12','SP11',
                             'SP04','SP09','SP10','SP05'))+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 6)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))+
  ggtitle("Physicochemical variables")+
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=15))
all_plot

#Place the arrows in the forefront of the points
all_plot$layers <- c(all_plot$layers, all_plot$layers[[2]])

#The options for styling the plot within the function itself are somewhat limited, but since it produces a 
#ggplot object, we can re-specify the necessary layers. The following code should work on any object 
#output from ggbiplot. First we find the geom segment and geom text layers:

seg <- which(sapply(all_plot$layers, function(x) class(x$geom)[1] == 'GeomSegment'))
txt <- which(sapply(all_plot$layers, function(x) class(x$geom)[1] == 'GeomText'))

#We can change the colour and width of the segments by doing
all_plot$layers[[seg[1]]]$aes_params$colour <- 'darkred' 
all_plot$layers[[seg[2]]]$aes_params$colour <- 'darkred'

#Labels
# Extract loadings of the variables
PCAloadings <- data.frame(Variables = rownames(all$rotation), all$rotation)

#To change the labels to have a gray background, we need to overwrite the geom_text layer with a geom_label layer:
all_plot$layers[[txt]] <- geom_label(aes(x = xvar, y = yvar, label = PCAloadings$Variables,
                                         angle = 0.45, hjust = 0.5, fontface = "bold"), 
                                     label.size = NA,
                                     color= 'darkred',
                                     data = all_plot$layers[[txt]]$data,
                                     fill = '#dddddd80')
all_plot

ggsave(path = "Figures/1 GRADIENT","SP_fig2_1.png", width = 7, height = 6, dpi = 300)



## >> Obtaining scores PC1 and PC2 ----
scores <- as.data.frame(all$x[,1:2])
scores$Site <- unique(data[,1])
scores <- scores[,c(3,1,2)] #Reorder to have Site as the first column




# 2.2. Physicochemical PCA (with AI, MAT & MAP) ----

data <- my_data[,c(2,5:8,10,12,13,17:23,27:30)]
for (i in which(sapply(data, is.numeric))) {
  for (j in which(is.na(data[, i]))) {
    data[j, i] <- mean(data[data[, "Site"] == data[j, "Site"], i],  na.rm = TRUE)
  }
}
library(dplyr)
data <- rename(data,
               Temperature = Soil_Temp,
               WC = Water_content,
               CN = C_N)

data2 <- data[,c(-1)]

#By means:
data_mit <- data %>% 
  group_by(Site) %>%
  summarise_all("mean")
data_mit2 <- data_mit[,c(-1)]
site <- data_mit[,1]
site <- site$Site
site <- factor(site, levels = site)
all <- prcomp(na.omit(data_mit2), center = TRUE,
              scale. = TRUE) 
plot(all, type = "l")
plot(all)
summary(all)


library(factoextra)
# Eigenvalues
eig.val <- get_eigenvalue(all)
eig.val

# Results for Variables
res.var <- get_pca_var(all)
res.var$coord          # Coordinates
res.var$contrib # Contributions to the PCs
# # To add a dataframe with a table of OTU
contributions <- res.var$contrib[,1:2]
contributionsdf = as.data.frame(contributions[,1:2])
# write.csv(contributionsdf, "Figures/1 GRADIENT/PCA_fsca_contributions_AIMATMAP.csv", row.names=TRUE)

res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(all)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

a <- fviz_contrib(all, choice = "var", axes = 1, top = 10)
b <- fviz_contrib(all, choice = "var", axes = 2, top = 10)
ggarrange(a,b, ncol=2) + theme_classic()

# ggsave(path = "Figures/1 GRADIENT","PCA_fsca_contributions_AIMATMAP.png", width = 10, height = 6, dpi = 300)


mycolors2<-c("#427681","#3891A6","#9BBC79","#CCD263","#E5DD58",
             "#FDE74C", "#EC9E57", "#E3655B", "#DB5461","#D84652",
             "#7D1809", "#290500")
library(ggbiplot)
all_plot <- ggbiplot(all, obs.scale = 1, var.scale = 1, 
                     ellipse = TRUE, fill=site,
                     varname.adjust = 2.5,
                     circle = TRUE, alpha=0,
                     loadings.label.repel=TRUE) +
  theme_classic()+
  scale_fill_manual(values = mycolors2,
                    breaks=c('SP08', 'SP01', 'SP02','SP07',
                             'SP06', 'SP03','SP12','SP11',
                             'SP04','SP09','SP10','SP05'))+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 6)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))+
  ggtitle("Physicochemical variables")+
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=15))
all_plot

#Place the arrows in the forefront of the points
all_plot$layers <- c(all_plot$layers, all_plot$layers[[2]])

#The options for styling the plot within the function itself are somewhat limited, but since it produces a 
#ggplot object, we can re-specify the necessary layers. The following code should work on any object 
#output from ggbiplot. First we find the geom segment and geom text layers:

seg <- which(sapply(all_plot$layers, function(x) class(x$geom)[1] == 'GeomSegment'))
txt <- which(sapply(all_plot$layers, function(x) class(x$geom)[1] == 'GeomText'))

#We can change the colour and width of the segments by doing
all_plot$layers[[seg[1]]]$aes_params$colour <- 'darkred' 
all_plot$layers[[seg[2]]]$aes_params$colour <- 'darkred'

#Labels
# Extract loadings of the variables
PCAloadings <- data.frame(Variables = rownames(all$rotation), all$rotation)

#To change the labels to have a gray background, we need to overwrite the geom_text layer with a geom_label layer:
all_plot$layers[[txt]] <- geom_label(aes(x = xvar, y = yvar, label = PCAloadings$Variables,
                                         angle = 0.45, hjust = 0.5, fontface = "bold"), 
                                     label.size = NA,
                                     color= 'darkred',
                                     data = all_plot$layers[[txt]]$data,
                                     fill = '#dddddd80')
all_plot


# ggsave(path = "Figures/1 GRADIENT","SP_fig2_3.png", width = 7, height = 6, dpi = 300)




# (3. See 2./Enzymes statistical tests)
# data <- my_data[,c(2,7,46:53)]
# 
# data_summary2 <- function(data, varname, groupnames){
#   require(plyr)
#   summary_func <- function(x, col){
#     c(mean = mean(x[[col]], na.rm=TRUE),
#       sd = sd(x[[col]], na.rm=TRUE),
#       sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
#   }
#   data_sum<-ddply(data, groupnames, .fun=summary_func,
#                   varname)
#   data_sum <- rename(data_sum, c("mean" = varname))
#   return(data_sum)
# }
# 
# alpha <- data_summary2(data, varname="alpha", 
#                     groupnames=c("Site", "AI"))
# beta <- data_summary2(data, varname="beta", 
#                      groupnames=c("Site", "AI"))
# xyl <- data_summary2(data, varname="xyl", 
#                     groupnames=c("Site", "AI"))
# cbh <- data_summary2(data, varname="cbh", 
#                     groupnames=c("Site", "AI"))
# gla <- data_summary2(data, varname="gla", 
#                     groupnames=c("Site", "AI"))
# fos <- data_summary2(data, varname="fos", 
#                     groupnames=c("Site", "AI"))
# leu <- data_summary2(data, varname="leu", 
#                     groupnames=c("Site", "AI"))
# phe <- data_summary2(data, varname="phe", 
#                     groupnames=c("Site", "AI"))
# 
# alpha$Site <- factor(alpha$Site, levels = c("SP08", "SP01", "SP02",
#                                                     "SP07","SP06","SP03",
#                                                     "SP12", "SP11", "SP04",
#                                                     "SP09", "SP10","SP05"))
# beta$Site <- factor(beta$Site, levels = c("SP08", "SP01", "SP02",
#                                             "SP07","SP06","SP03",
#                                             "SP12", "SP11", "SP04",
#                                             "SP09", "SP10","SP05"))
# xyl$Site <- factor(xyl$Site, levels = c("SP08", "SP01", "SP02",
#                                             "SP07","SP06","SP03",
#                                             "SP12", "SP11", "SP04",
#                                             "SP09", "SP10","SP05"))
# cbh$Site <- factor(cbh$Site, levels = c("SP08", "SP01", "SP02",
#                                             "SP07","SP06","SP03",
#                                             "SP12", "SP11", "SP04",
#                                             "SP09", "SP10","SP05"))
# fos$Site <- factor(fos$Site, levels = c("SP08", "SP01", "SP02",
#                                             "SP07","SP06","SP03",
#                                             "SP12", "SP11", "SP04",
#                                             "SP09", "SP10","SP05"))
# gla$Site <- factor(gla$Site, levels = c("SP08", "SP01", "SP02",
#                                             "SP07","SP06","SP03",
#                                             "SP12", "SP11", "SP04",
#                                             "SP09", "SP10","SP05"))
# leu$Site <- factor(leu$Site, levels = c("SP08", "SP01", "SP02",
#                                             "SP07","SP06","SP03",
#                                             "SP12", "SP11", "SP04",
#                                             "SP09", "SP10","SP05"))
# phe$Site <- factor(phe$Site, levels = c("SP08", "SP01", "SP02",
#                                             "SP07","SP06","SP03",
#                                             "SP12", "SP11", "SP04",
#                                             "SP09", "SP10","SP05"))
#                          
# 
# library(ggplot2)
# plot.theme1 <- theme_classic() +
#   theme(text=element_text(size=15),
#         axis.title.x = element_text(size = rel(1.2), angle = 00, margin = margin(t=8)),
#         axis.title.y = element_text(size = rel(1.2), angle = 90, margin = margin(t=8)),
#         plot.title = element_text(size=22),
#         legend.position = "none",
#         legend.title = element_text(size = 18),
#         legend.text = element_text(size = 15),
#         axis.text.x = element_text(size=15),
#         axis.text.y = element_text(size=15))
# 
# palpha <-ggplot(alpha, aes(x=Site, y=alpha, fill=AI)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=alpha-sem, ymax=alpha+sem), width=.2)+
#   ggtitle("alpha") +
#   plot.theme1+
#   scale_fill_AI(discrete = FALSE, palette = "Sites")+
#   theme(axis.title.x = element_blank())+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))
# 
# pbeta <-ggplot(beta, aes(x=Site, y=beta, fill=AI)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=beta-sem, ymax=beta+sem), width=.2)+
#   ggtitle("beta") +
#   plot.theme1+
#   scale_fill_AI(discrete = FALSE, palette = "Sites")+
#   theme(axis.title.x = element_blank())+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))
# 
# pxyl <-ggplot(xyl, aes(x=Site, y=xyl, fill=AI)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=xyl-sem, ymax=xyl+sem), width=.2)+
#   ggtitle("xyl") +
#   plot.theme1+
#   scale_fill_AI(discrete = FALSE, palette = "Sites")+
#   theme(axis.title.x = element_blank())+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))
# 
# pcbh <-ggplot(cbh, aes(x=Site, y=cbh, fill=AI)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=cbh-sem, ymax=cbh+sem), width=.2)+
#   ggtitle("cbh") +
#   plot.theme1+
#   scale_fill_AI(discrete = FALSE, palette = "Sites")+
#   theme(axis.title.x = element_blank())+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))
# 
# pgla <-ggplot(gla, aes(x=Site, y=gla, fill=AI)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=gla-sem, ymax=gla+sem), width=.2)+
#   ggtitle("gla") +
#   plot.theme1+
#   scale_fill_AI(discrete = FALSE, palette = "Sites")+
#   theme(axis.title.x = element_blank())+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))
# 
# pfos <-ggplot(fos, aes(x=Site, y=fos, fill=AI)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=fos-sem, ymax=fos+sem), width=.2)+
#   ggtitle("fos") +
#   plot.theme1+
#   scale_fill_AI(discrete = FALSE, palette = "Sites")+
#   theme(axis.title.x = element_blank())+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   labs(y = expression(paste("μmol MUF·" ~  g ~ DW^-1,"·" ~ h^-1)))
# 
# pleu <-ggplot(leu, aes(x=Site, y=leu, fill=AI)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=leu-sem, ymax=leu+sem), width=.2)+
#   ggtitle("leu") +
#   plot.theme1+
#   scale_fill_AI(discrete = FALSE, palette = "Sites")+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   labs(y = expression(paste("μmol AMC·" ~  g ~ DW^-1,"·" ~ h^-1)))
# 
# pphe <-ggplot(phe, aes(x=Site, y=phe, fill=AI)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_errorbar(aes(ymin=phe-sem, ymax=phe+sem), width=.2)+
#   ggtitle("phe") +
#   plot.theme1+
#   scale_fill_AI(discrete = FALSE, palette = "Sites")+
#   scale_y_continuous(labels = scales::number_format(accuracy = 0.001))+
#   labs(y = expression(paste("μmol DIQC·" ~  g ~ DW^-1,"·" ~ h^-1)))
# 
# # all <- ggarrange(palpha, pbeta, pxyl, pcbh, pfos, pgla, pleu, pphe,
# #                   nrow = 8, ncol = 1,
# #                   common.legend = TRUE,
# #                   legend = "right")
# # 
# # ggsave(path = "Figures/1 GRADIENT","enzymes.png",
# #        width = 10, height = 30, dpi = 300)
# 
# all <- ggarrange(palpha, pbeta, pxyl, pcbh, pfos, pgla, pleu, pphe,
#                  nrow = 4, ncol = 2,
#                  common.legend = TRUE,
#                  legend = "right")
# 
# ggsave(path = "Figures/1 GRADIENT","enzymes2.png",
#        width = 30, height = 20, dpi = 300)


#.----
# One-way MANOVA ----

#Condition for MANOVA: the dataset must have more
# observations (rows) per group in the independent
# variable than a number of the dependent variables.

manova_data <- my_data[,c(2,35,37,38,40:64)]

# 1rst we separate the dependent variables from the other ones:
dependent_vars <- as.matrix(manova_data[,-c(1)])
# dependent variables needs to be entered as a matrix, not as a dataframe
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
# PCA with FUNCTIONAL explan. vars. ----

pca_data <- func %>%
  group_by(Site) %>%
  summarise_all("mean")

site_order <- my_data[,c(2,7)] 
site_order <- site_order[!duplicated(site_order), ] #Erase duplicated lines from dataframe
pca_data <- pca_data[order(site_order$AI, decreasing = T),]
pca_data$Site <- factor(pca_data$Site, levels = c("SP08","SP01","SP02","SP07", "SP06" ,"SP03" ,"SP12", "SP11" ,"SP04", "SP09", "SP10" ,"SP05"))

pcr2 <- pca_data[,c(-1)]


#Select column with levels (Site)
site <- factor(pca_data$Site, levels = pca_data$Site)
site

pc <- prcomp(na.omit(pcr2), center = TRUE,
             scale. = TRUE) 

# 
# loadings <- pc$rotation
# 
scores = as.data.frame(pc$x)
scores$AI <- site_order$AI
scores$Site <- site_order$Site
# 
# 
# 
# p5 <- ggplot()+
#   geom_point(data = scores, aes(x = PC1, y = PC2, size=4, color=AI))+
#   geom_text(aes(scores$PC1, scores$PC2,label = scores$Site), size = 3, hjust = 1.7) +
#   scale_color_AI(discrete = FALSE, palette = "Sites")+
#   plot.theme1+
#   ggtitle("PCoA") +
#   guides(size = "none")
# 
# grid.arrange(p5,ncol=1)




pca_data$AI <- site_order$AI

scores = as.data.frame(scale(pc$x))
scores$AI <- site_order$AI
scores$Site <- site_order$Site

autoplot(pc, data=pca_data, 
         loadings = TRUE, loadings.colour = 'brown',
         loadings.label.colour='brown', loadings.label = TRUE,
         loadings.label.size = 7,
         loadings.label.repel=TRUE,
         color = pca_data$AI)+
  plot.theme1+
  geom_point(aes(fill=AI), colour= "black", pch=21, size = 5)+
  scale_fill_AI(discrete = FALSE, palette = "Sites")+
  ggtitle("Functional response variables")+
  geom_text(aes(label = scores$Site), size = 4, hjust = 1.5) +
  theme(legend.title = element_blank(),
        legend.text=element_text(size = 12),
        title = element_text(size = 15,face="bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15, face="plain"))+
  scale_x_continuous(expand = c(0.1, 0.1))

# ggsave(path = "Figures/1 GRADIENT","PCA_func_response_means.png", width = 10, height = 8, dpi = 300)





# PCA with FUNCTIONAL AGGREGATED explan. vars. ----

func <-func %>%
  mutate(Cenz = select(., 4:7) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(MB = select(.,2:3) %>% rowSums(na.rm = TRUE))
func <- func[,-c(2:7)]

pca_data <- func %>%
  group_by(Site) %>%
  summarise_all("mean")

site_order <- my_data[,c(2,7)] 
site_order <- site_order[!duplicated(site_order), ] #Erase duplicated lines from dataframe
pca_data <- pca_data[order(site_order$AI, decreasing = T),]
pca_data$Site <- factor(pca_data$Site, levels = c("SP08","SP01","SP02","SP07", "SP06" ,"SP03" ,"SP12", "SP11" ,"SP04", "SP09", "SP10" ,"SP05"))

pcr2 <- pca_data[,c(-1)]


#Select column with levels (Site)
site <- factor(pca_data$Site, levels = pca_data$Site)
site

pc <- prcomp(na.omit(pcr2), center = TRUE,
             scale. = TRUE) 

# 
# loadings <- pc$rotation
# 
scores = as.data.frame(pc$x)
scores$AI <- site_order$AI
scores$Site <- site_order$Site
# 
# 
# 
# p5 <- ggplot()+
#   geom_point(data = scores, aes(x = PC1, y = PC2, size=4, color=AI))+
#   geom_text(aes(scores$PC1, scores$PC2,label = scores$Site), size = 3, hjust = 1.7) +
#   scale_color_AI(discrete = FALSE, palette = "Sites")+
#   plot.theme1+
#   ggtitle("PCoA") +
#   guides(size = "none")
# 
# grid.arrange(p5,ncol=1)




pca_data$AI <- site_order$AI

scores = as.data.frame(scale(pc$x))
scores$AI <- site_order$AI
scores$Site <- site_order$Site

autoplot(pc, data=pca_data, 
         loadings = TRUE, loadings.colour = 'brown',
         loadings.label.colour='brown', loadings.label = TRUE,
         loadings.label.size = 7,
         loadings.label.repel=TRUE,
         color = pca_data$AI)+
  plot.theme1+
  geom_point(aes(fill=AI), colour= "black", pch=21, size = 5)+
  scale_fill_AI(discrete = FALSE, palette = "Sites")+
  ggtitle("Functional response variables")+
  geom_text(aes(label = scores$Site), size = 4, hjust = 1.5) +
  theme(legend.title = element_blank(),
        legend.text=element_text(size = 12),
        title = element_text(size = 15,face="bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15, face="plain"))+
  scale_x_continuous(expand = c(0.1, 0.1))

# ggsave(path = "Figures/1 GRADIENT","PCA_func_response_means.png", width = 10, height = 8, dpi = 300)




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

# New_data <- data_summary(New_data, varname="values", 
#                          groupnames=c("Site","AI", "enzyme"))
New_data <- data_summary2(New_data, varname="values", 
                          groupnames=c("Site","AI", "enzyme"))
# make a plot
New_data$values = as.numeric(New_data$values)
enzyme <- ggplot(New_data, aes(x=AI, y=values, color = AI)) +
  geom_point(size= 3, shape = 15) +
  geom_pointrange(data=subset(New_data), aes(ymin=values-sem, ymax=values+sem),
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
  theme(legend.position = "right")+
  coord_cartesian(xlim = c(0,1.4))+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))+
  scale_color_AI(discrete = FALSE, palette = "Sites")
  # geom_smooth(method = "lm", color = "blue", fill = "grey") +
  # stat_cor(label.x = 0.6, label.y.npc="top",
  #          aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
  #          p.accuracy = 0.001, r.accuracy = 0.01,
  #          color = "blue", size = 4.5)
enzyme

# save the plot
# ggsave(path = "Figures/1 GRADIENT", "SP_enzymes_AI_se.png", width = 16, height = 8, dpi = 300)
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

# ggsave(path = "Figures/1 GRADIENT","SP_PCA_enzymes_mean_sd.png", width = 7, height = 6, dpi = 300)


  
  
  



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
# WEOM PCA ----
#BY MEANS
pca_data <- my_data[,c(2,79:92)]

pca_data <- pca_data %>%
  group_by(Site) %>%
  summarise_all("mean")

site_order <- my_data[,c(2,7)] 
site_order <- site_order[!duplicated(site_order), ] #Erase duplicated lines from dataframe
pca_data <- pca_data[order(site_order$AI, decreasing = T),]
pca_data$Site <- factor(pca_data$Site, levels = pca_data$Site[order(site_order$AI)])

pcr2 <- pca_data[,c(-1)]

#Select column with levels (Site)
site <- factor(pca_data$Site, levels = c("SP08","SP01","SP02",
                                         "SP07","SP06","SP03",
                                         "SP12","SP11","SP04",
                                         "SP09","SP10","SP05"))
site

pc <- prcomp(na.omit(pcr2), center = TRUE,
             scale. = TRUE) 

plot(pc, type = "l")
plot(pc)
summary(pc)

mycolors2<-c("#427681","#3891A6","#9BBC79","#CCD263","#E5DD58",
             "#FDE74C", "#EC9E57", "#E3655B", "#DB5461","#D84652",
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

ggsave(path = "Figures/1 GRADIENT","PCA_weom.png", width = 10, height = 8, dpi = 300)



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
 
#__________________ -----
# SHIT .-----

# New_data <- my_data[,c(2,7,8,10,12:13,17:23,27:30)]
# 
# New_data <- setNames(New_data, c("Site","AI","Soil temp.","Water Content","SOM","pH",
#                                  "TOC", "TC", "TN",
#                                  "C/N","NH4","PO43","SO42",
#                                  "Sand","Silt","Clay","Litter"))
# New_data1 <- New_data %>% 
#   gather(Physio, values, c(3:6))
# New_data2 <- New_data %>% 
#   gather(Physio, values, c(7:10))
# New_data3 <- New_data %>% 
#   gather(Physio, values, c(11:13))
# New_data4 <- New_data %>% 
#   gather(Physio, values, c(14:17))
# 
# 
# library(formattable)
# New_data1$AI <-formattable(New_data1$AI ,format="f",digits=2)
# New_data2$AI <-formattable(New_data2$AI ,format="f",digits=2)
# New_data3$AI <-formattable(New_data3$AI ,format="f",digits=2)
# New_data4$AI <-formattable(New_data4$AI ,format="f",digits=2)
# 
# New_data1$AI <- as.factor(New_data1$AI)
# New_data2$AI <- as.factor(New_data2$AI)
# New_data3$AI <- as.factor(New_data3$AI)
# New_data4$AI <- as.factor(New_data4$AI)
# 
# 
# New_data1$Site <- factor(New_data1$Site, levels = c("SP08", "SP01", "SP02",
#                                                     "SP07","SP06","SP03",
#                                                     "SP12", "SP11", "SP04",
#                                                     "SP09", "SP10","SP05"))
# New_data2$Site <- factor(New_data2$Site, levels = c("SP08", "SP01", "SP02",
#                                                     "SP07","SP06","SP03",
#                                                     "SP12", "SP11", "SP04",
#                                                     "SP09", "SP10","SP05"))
# New_data3$Site <- factor(New_data3$Site, levels = c("SP08", "SP01", "SP02",
#                                                     "SP07","SP06","SP03",
#                                                     "SP12", "SP11", "SP04",
#                                                     "SP09", "SP10","SP05"))
# New_data4$Site <- factor(New_data4$Site, levels = c("SP08", "SP01", "SP02",
#                                                     "SP07","SP06","SP03",
#                                                     "SP12", "SP11", "SP04",
#                                                     "SP09", "SP10","SP05"))
# 
# mycolors2<-c("#427681","#3891A6","#9BBC79","#CCD263","#E5DD58",
#              "#FDE74C", "#EC9E57", "#E3655B", "#DB5461","#D84652",
#              "#7D1809", "#290500")
# 
# 
# 
# p1 <- ggplot(New_data1, aes(x=AI, y=values, fill=Site))+
#   scale_fill_manual(values = mycolors2)+
#   geom_boxplot(width = 0.1)+
#   facet_wrap(.~Physio, scales = "free", nrow = 2)+
#   xlab("Aridity Index")+
#   theme(axis.title.y=element_blank(),
#         strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold"),
#         panel.background = element_blank(),
#         panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1),
#         axis.text.x = element_text(angle= 45, hjust = 1, size = 18, color = "black"),
#         axis.text.y = element_text(size = 18, color = "black"))+
#   coord_cartesian(xlim = c(0,1.4))+
#   scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
#                      labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
# p1
# 
# 
# 
# 
# 
# New_data <- my_data[,c(2,7,8,10,12:13,17:23,27:30)]
# 
# New_data <- setNames(New_data, c("Site","AI","Soil temp.","Water Content","SOM","pH",
#                                  "TOC", "TC", "TN",
#                                  "C/N","NH4","PO43","SO42",
#                                  "Sand","Silt","Clay","Litter"))
# 
# New_data <- New_data %>%
#   gather(Physio, values, c(3:17))
# 
# 
# New_data$AI <- as.factor(New_data$AI)
# # New_data$AI = as.numeric(New_data$AI)
# 
# New_data$Site <- factor(New_data$Site, levels = c("SP08", "SP01", "SP02",
#                                                   "SP07","SP06","SP03",
#                                                   "SP12", "SP11", "SP04",
#                                                   "SP09", "SP10","SP05"))
# 
# 
mycolors2<-c("#427681","#3891A6","#9BBC79","#CCD263","#E5DD58",
             "#FDE74C", "#EC9E57", "#E3655B", "#DB5461","#D84652",
             "#7D1809", "#290500")
# 
# 
# 
# ggplot(New_data, aes(x=AI, y=values, fill=Site))+
#   scale_fill_manual(values = mycolors2)+
#   geom_boxplot()+
#   facet_wrap(.~Physio, scales = "free")+
#   xlab("Aridity Index")+
#   theme(axis.title.y=element_blank(),
#         strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold"),
#         panel.background = element_blank(),
#         panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1),
#         axis.text.x = element_text(angle= 45, hjust = 1, size = 18, color = "black"),
#         axis.text.y = element_text(size = 18, color = "black"))+
#   scale_x_discrete(breaks = c(0, 0.5, 1, 1.5))
# 
# 
# 
# # save the plot
# # ggsave(path = "Figures/1 GRADIENT","SP_physio2_AI.png", width = 20, height = 10, dpi = 300)

