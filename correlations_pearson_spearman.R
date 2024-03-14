# LOAD LIBRARY ----
library("readxl")
library("ggplot2")
library("dplyr")
library("ggpubr")
library("gridExtra")


# LOAD DATA ----
my_data <- read_excel("SP_metadata_2021.xlsx")
mit_16S <- read_excel("C:/Users/Anna/OneDrive - Universitat de Girona/GRADCATCH/ANALISIS/amplicon_sequencing_2021_JD/16S/Alphadiv_SP_16S_mit.xlsx")
mit_18S <- read_excel("C:/Users/Anna/OneDrive - Universitat de Girona/GRADCATCH/ANALISIS/amplicon_sequencing_2021_JD/ITS2/Alphadiv_SP_ITS2_mit.xlsx")

mit_16S <- rename(mit_16S,
                  Site = GxS,
                  Shannon = Shannon_mean,
                  Chao = Chao1_mean,
                  Sobs = Observed_mean)
mit_18S <- rename(mit_18S,
                  Site = GxS,
                  Shannon = Shannon_mean,
                  Chao = Chao1_mean,
                  Sobs = Observed_mean)

aridity <- my_data[,c(2,7)]
aridity <- aridity %>%
  group_by(Site) %>%
  summarise_all(mean)

mit_16S <- merge(aridity, mit_16S, by = "Site")
mit_18S <- merge(aridity, mit_18S, by = "Site")

mit_16S <- mit_16S[,-c(6:11)]
mit_18S <- mit_18S[,-c(6:11)]

rm(aridity)


#PHYSICOCHEMICAL ----

#Clean non-interesting variables:
cor_data <- my_data[,c(2,7,13,20:23,33,34,28,76,30,86,92)]

replace_na_with_mean <- function(data) {
  data %>%
    group_by(Site) %>%
    mutate_all(~ ifelse(is.na(.), mean(., na.rm = TRUE), .)) %>%
    ungroup()
}
cor_data <- replace_na_with_mean(cor_data)

cor_data <- cor_data %>%
  group_by(Site) %>%
  summarise_all(mean)


#PEARSON:
numerical_variables <- sapply(cor_data, is.numeric)
numerical_variables_names <- names(numerical_variables)[numerical_variables]

all_plots <- list()

for (var in numerical_variables_names) {
  if (var != "AI") {
    plot_title <- paste("AI vs", var)
    all_plots[[length(all_plots) + 1]] <- ggscatter(cor_data, x = "AI", y = var,
                                                    add = "reg.line", conf.int = TRUE,
                                                    cor.coef = TRUE, cor.method = "pearson",
                                                    xlab = "AI", ylab = var,
                                                    main = plot_title)
  }
}

grid.arrange(grobs = all_plots, nrow = 2)




#SPEARMAN:
numerical_variables <- sapply(cor_data, is.numeric)
numerical_variables_names <- names(numerical_variables)[numerical_variables]

all_plots <- list()

for (var in numerical_variables_names) {
  if (var != "AI") {
    plot_title <- paste("AI vs", var)
    all_plots[[length(all_plots) + 1]] <- ggscatter(cor_data, x = "AI", y = var,
                                                    add = "reg.line", conf.int = TRUE,
                                                    cor.coef = TRUE, cor.method = "spearman",
                                                    xlab = "AI", ylab = var,
                                                    main = plot_title)
  }
}

grid.arrange(grobs = all_plots, nrow = 2)



# FUNCTIONS ----
#Clean non-interesting variables:
cor_data <- my_data[,c(2,7,40,41,46:53,64)]

replace_na_with_mean <- function(data) {
  data %>%
    group_by(Site) %>%
    mutate_all(~ ifelse(is.na(.), mean(., na.rm = TRUE), .)) %>%
    ungroup()
}
cor_data <- replace_na_with_mean(cor_data)

cor_data <- cor_data %>%
  group_by(Site) %>%
  summarise_all(mean)

cor_data$aridity <- (1-cor_data$AI)
cor_data <- cor_data[,-c(2)]


#PEARSON:
numerical_variables <- sapply(cor_data, is.numeric)
numerical_variables_names <- names(numerical_variables)[numerical_variables]

all_plots <- list()

for (var in numerical_variables_names) {
  if (var != "aridity") {
    plot_title <- paste("Aridity vs", var)
    all_plots[[length(all_plots) + 1]] <- ggscatter(cor_data, x = "aridity", y = var,
                                                    add = "reg.line", conf.int = TRUE,
                                                    cor.coef = TRUE, cor.method = "pearson",
                                                    xlab = "Aridity", ylab = var,
                                                    main = plot_title)
  }
}

grid.arrange(grobs = all_plots, nrow = 2)


#SPEARMAN:
numerical_variables <- sapply(cor_data, is.numeric)
numerical_variables_names <- names(numerical_variables)[numerical_variables]

all_plots <- list()

for (var in numerical_variables_names) {
  if (var != "AI") {
    plot_title <- paste("AI vs", var)
    all_plots[[length(all_plots) + 1]] <- ggscatter(cor_data, x = "AI", y = var,
                                                    add = "reg.line", conf.int = TRUE,
                                                    cor.coef = TRUE, cor.method = "spearman",
                                                    xlab = "AI", ylab = var,
                                                    main = plot_title)
  }
}

grid.arrange(grobs = all_plots, nrow = 2)



# WITHOUT SP06:
cor_data <- subset(cor_data, Site != "SP06")

# PEARSON:
numerical_variables <- sapply(cor_data, is.numeric)
numerical_variables_names <- names(numerical_variables)[numerical_variables]

all_plots <- list()

for (var in numerical_variables_names) {
  if (var != "AI") {
    plot_title <- paste("AI vs", var)
    all_plots[[length(all_plots) + 1]] <- ggscatter(cor_data, x = "AI", y = var,
                                                    add = "reg.line", conf.int = TRUE,
                                                    cor.coef = TRUE, cor.method = "pearson",
                                                    xlab = "AI", ylab = var,
                                                    main = plot_title)
  }
}

grid.arrange(grobs = all_plots, nrow = 2)



#SPEARMAN:
numerical_variables <- sapply(cor_data, is.numeric)
numerical_variables_names <- names(numerical_variables)[numerical_variables]

all_plots <- list()

for (var in numerical_variables_names) {
  if (var != "AI") {
    plot_title <- paste("AI vs", var)
    all_plots[[length(all_plots) + 1]] <- ggscatter(cor_data, x = "AI", y = var,
                                                    add = "reg.line", conf.int = TRUE,
                                                    cor.coef = TRUE, cor.method = "spearman",
                                                    xlab = "AI", ylab = var,
                                                    main = plot_title)
  }
}

grid.arrange(grobs = all_plots, nrow = 2)







# ALPHA DIVERISTY ----

## 16S ----

#PEARSON:
numerical_variables <- sapply(mit_16S, is.numeric)
numerical_variables_names <- names(numerical_variables)[numerical_variables]

all_plots <- list()

for (var in numerical_variables_names) {
  if (var != "AI") {
    plot_title <- paste("AI vs", var)
    all_plots[[length(all_plots) + 1]] <- ggscatter(mit_16S, x = "AI", y = var,
                                                    add = "reg.line", conf.int = TRUE,
                                                    cor.coef = TRUE, cor.method = "pearson",
                                                    xlab = "AI", ylab = var,
                                                    main = plot_title)
  }
}

grid.arrange(grobs = all_plots, nrow = 1)


#SPEARMAN:
numerical_variables <- sapply(mit_16S, is.numeric)
numerical_variables_names <- names(numerical_variables)[numerical_variables]

all_plots <- list()

for (var in numerical_variables_names) {
  if (var != "AI") {
    plot_title <- paste("AI vs", var)
    all_plots[[length(all_plots) + 1]] <- ggscatter(mit_16S, x = "AI", y = var,
                                                    add = "reg.line", conf.int = TRUE,
                                                    cor.coef = TRUE, cor.method = "spearman",
                                                    xlab = "AI", ylab = var,
                                                    main = plot_title)
  }
}

grid.arrange(grobs = all_plots, nrow = 1)



## 18S ----
numerical_variables <- sapply(mit_18S, is.numeric)
numerical_variables_names <- names(numerical_variables)[numerical_variables]

all_plots <- list()

for (var in numerical_variables_names) {
  if (var != "AI") {
    plot_title <- paste("AI vs", var)
    all_plots[[length(all_plots) + 1]] <- ggscatter(mit_18S, x = "AI", y = var,
                                                    add = "reg.line", conf.int = TRUE,
                                                    cor.coef = TRUE, cor.method = "pearson",
                                                    xlab = "AI", ylab = var,
                                                    main = plot_title)
  }
}

grid.arrange(grobs = all_plots, nrow = 1)


#SPEARMAN:
numerical_variables <- sapply(mit_18S, is.numeric)
numerical_variables_names <- names(numerical_variables)[numerical_variables]

all_plots <- list()

for (var in numerical_variables_names) {
  if (var != "AI") {
    plot_title <- paste("AI vs", var)
    all_plots[[length(all_plots) + 1]] <- ggscatter(mit_18S, x = "AI", y = var,
                                                    add = "reg.line", conf.int = TRUE,
                                                    cor.coef = TRUE, cor.method = "spearman",
                                                    xlab = "AI", ylab = var,
                                                    main = plot_title)
  }
}

grid.arrange(grobs = all_plots, nrow = 1)

