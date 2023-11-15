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

my_data$Köppen <- my_data$Site
my_data$Köppen <- gsub("SP05", "Dry Arid", my_data$Köppen) #SP05
my_data$Köppen <- gsub("SP10", "Dry Arid", my_data$Köppen) #SP10
my_data$Köppen <- gsub("SP09", "Dry Semi-arid", my_data$Köppen) #SP09
my_data$Köppen <- gsub("SP04", "Dry Semi-arid", my_data$Köppen) #SP04
my_data$Köppen <- gsub("SP11", "Temperate (Dry Hot Summer)", my_data$Köppen) #SP11
my_data$Köppen <- gsub("SP12", "Temperate (Dry Hot Summer)", my_data$Köppen) #SP12
my_data$Köppen <- gsub("SP03", "Temperate (Dry Hot Summer)", my_data$Köppen) #SP03
my_data$Köppen <- gsub("SP06", "Temperate (No dry season)", my_data$Köppen) #SP06
my_data$Köppen <- gsub("SP07", "Temperate (No dry season)", my_data$Köppen) #SP07
my_data$Köppen <- gsub("SP02", "Temperate (Dry Warm Summer)", my_data$Köppen) #SP02
my_data$Köppen <- gsub("SP01", "Temperate (No dry season)", my_data$Köppen) #SP01
my_data$Köppen <- gsub("SP08", "Temperate (Dry Warm Summer)", my_data$Köppen) #SP08


#Data without ratios, percentages....
data <- my_data[,c(93,5:10,12:13,17:23,27:30,33:35,37:38,40:64,76,79:92)]

#Dummies....
dummy <- my_data[,c(5:7,76)]

#Independent variables without dummies....
indepe <- my_data[,c(8:10,12:13,17:23,27:30,33:34,86,92)]

#Independent var with dummies...
duindepe <- my_data[,c(5:7,76,8:10,12:13,17:23,27:30,33:34,86,92)]

#Dependen vars....
depe <- my_data[,c(35,37:38,40:64)]



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

# # Feature Selection (CARET PACKAGE) ----
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

clean_duindepe <- duindepe[,-c(2,7,5,6,12,8,10,3,1,19,17)]

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
              clean_duindepe$pH+clean_duindepe$TC+clean_duindepe$C_N+
              clean_duindepe$NH4+clean_duindepe$PO43+clean_duindepe$SO42+clean_duindepe$Silt+
              clean_duindepe$Litter+clean_duindepe$L_TC+clean_duindepe$L_TN+
              clean_duindepe$BIX+clean_duindepe$HIX
)

car::vif(model)
# vif_values <- vif(model)
# barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")
#CONFIRMAMOS QUE NINGUNA VARIABLES INDEPENDIENTE TIENE MULTICOLINEALIDAD (VIF>5)

rm(duindepe,dummy,model)



# FIGURES ----
# Physicochemical PCA (without AI) ----

data <- clean_duindepe
for (i in which(sapply(data, is.numeric))) {
  for (j in which(is.na(data[, i]))) {
    data[j, i] <- mean(data[data[, "Site"] == data[j, "Site"], i],  na.rm = TRUE)
  }
}
library(dplyr)
data <- rename(data,
               CN = C_N)

data2 <- data[,c(-1)]

#By means:
data_mit <- data %>% 
  group_by(Site) %>%
  summarise_all("mean")
data_mit2 <- data_mit[,c(-1)]

site <- unique(my_data[,c(2,93)])
site <- site[,c(2)]
# site <- site$Site
# site <- factor(site, levels = site)

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

# ggsave(path = "Figures/1 GRADIENT","PCA_fsca_contributions.png", width = 10, height = 6, dpi = 300)


# mycolors2<-c("#427681","#3891A6","#9BBC79","#CCD263","#E5DD58",
#              "#FDE74C", "#EC9E57", "#E3655B", "#DB5461","#D84652",
#              "#7D1809", "#290500")


mycolors2<-c("#538D88","#9BBD7A",
             "#FDE74C", "#EC9E57","#D84652")


library(ggbiplot)
all_plot <- ggbiplot(all, obs.scale = 1, var.scale = 1, 
                     ellipse = TRUE, fill=site,
                     varname.adjust = 2.5,
                     circle = TRUE, alpha=0,
                     loadings.label.repel=TRUE) +
  theme_classic()+
  scale_fill_manual(values = mycolors2,
                    breaks=c('Temperate (No dry season)', 
                             'Temperate (Dry Warm Summer)', 
                             'Temperate (Dry Hot Summer)',
                             'Dry Semi-arid',
                             'Dry Arid'))+
  geom_point(aes(fill=site), colour= "black", pch=21, size = 6)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))+
  ggtitle("Physicochemical variables")+
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=12))
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

ggsave(path = "Figures/1 GRADIENT","PCA_Phsca_selected2.png", width = 7, height = 6, dpi = 300)



## >> Obtaining scores PC1 and PC2 ----
scores <- as.data.frame(all$x[,1:2])
scores$Site <- unique(data[,1])
scores <- scores[,c(3,1,2)] #Reorder to have Site as the first column

