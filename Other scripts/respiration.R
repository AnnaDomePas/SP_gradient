#RESPIRATION ALL GRADIENTS -----

library(readxl)
data <- read_excel("Resultats SCRIPT Raz_FLUOR_summary_20230413.xlsx")

meta <- read_excel("Coordinates_alt.xlsx")

data <- merge(data, meta[,-c(2:4)], by.x=c("Site"), by.y=c("Site"), sort=T, all.x = TRUE)


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


# > EU ----
subdata <- data[data$Gradient == 'EU',]
subdata <- subdata[order(subdata$Latitude, decreasing = T),]

site_order <- meta[13:25,1] 
site_order <- merge(site_order, meta[13:25,-c(2:4)], by.x=c("Site"), by.y=c("Site"), sort=T, all.x = TRUE)
site_order <- site_order[order(site_order$Latitude, decreasing = T),]

files <- colnames(subdata[1:4])

for(i in 4){
  var = names(subdata)[i]
  print(
    ggplot(subdata, aes(x=factor(Site, level=c(site_order[,1]), ordered = TRUE), y=subdata[,i])) +
      geom_boxplot(outlier.shape = NA) +
      theme_classic() +
      labs(title=var) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.x = element_blank())
  )
}



# > GL ----
subdata <- data[data$Gradient == 'GL',]
subdata <- subdata[order(subdata$Latitude, decreasing = T),]

site_order <- meta[26:40,1] 
site_order <- merge(site_order, meta[26:40,-c(2:4)], by.x=c("Site"), by.y=c("Site"), sort=T, all.x = TRUE)
site_order <- site_order[order(site_order$Latitude, decreasing = T),]


files <- colnames(subdata[1:4])

for(i in 4){
  var = names(subdata)[i]
  print(
    ggplot(subdata, aes(x=factor(Site, level=c(site_order[,1]), ordered = TRUE), y=subdata[,i])) +
      geom_boxplot(outlier.shape = NA) +
      theme_classic() +
      labs(title=var) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.x = element_blank())
  )
}


# > SA ----
subdata <- data[data$Gradient == 'SA',]
subdata <- subdata[order(subdata$Latitude, decreasing = T),]

site_order <- meta[54:81,1] 
site_order <- merge(site_order, meta[54:81,-c(2:4)], by.x=c("Site"), by.y=c("Site"), sort=T, all.x = TRUE)
site_order <- site_order[order(site_order$Latitude, decreasing = T),]


files <- colnames(subdata[1:4])

for(i in 4){
  var = names(subdata)[i]
  print(
    ggplot(subdata, aes(x=factor(Site, level=c(site_order[,1]), ordered = TRUE), y=subdata[,i])) +
      geom_boxplot(outlier.shape = NA) +
      theme_classic() +
      labs(title=var) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylab(ifelse(files[i] =='phe', expression(paste("μmol DIQC·" ~  g ~ AFDW^-1,"·" ~ h^-1)),
                  (expression(paste("μmol MUF·" ~  g ~ AFDW^-1,"·" ~ h^-1))))) +
      theme(axis.title.x = element_blank())
  )
}


# > AL ----
subdata <- data[data$Gradient == 'AL',]
subdata <- subdata[order(subdata$Altitude_masl, decreasing = T),]

site_order <- meta[41:53,1] 
site_order <- merge(site_order, meta[41:53,-c(2:4)], by.x=c("Site"), by.y=c("Site"), sort=T, all.x = TRUE)
site_order <- site_order[order(site_order$Altitude_masl, decreasing = T),]


files <- colnames(subdata[1:4])

for(i in 4){
  var = names(subdata)[i]
  print(
    ggplot(subdata, aes(x=factor(Site, level=c(site_order[,1]), ordered = TRUE), y=subdata[,i])) +
      geom_boxplot(outlier.shape = NA) +
      theme_classic() +
      labs(title=var) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.x = element_blank())
  )
}
