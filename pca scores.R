# ANNA DOMÉNECH PASCUAL - 11/05/2023
# HOW TO OBTAIN SCORES FROM PCA


#1r. has de fer la teva PCA amb les rèpliques.

#Obtain scores:
scores <- as.data.frame(pc$x[,1:2])
scores$Site <- my_data[,2]

#Això mirat si ho necessites o no:
scores <- scores[,c(3,1,2)] #Reorder to have Site as the first column

scores_msd <- aggregate(cbind(PC1,PC2) ~ Site, data = scores, 
                        FUN = function(scores) c(mean = mean(scores), sd = sd(scores)))

#A partir d'aquí és on hi havia una mica de lio per culpa de la comanda anterior:

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

#Ara ja hauries de tenir un dataframe amb les mitjanes i SD per cada grup
#Només hauries de graficar amb un scatter plot normal i corrent i ja.