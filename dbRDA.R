#dbRDA ----
rm(list = ls())
library("ggplot2")
library(ggrepel)
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
library(mlbench)
library(caret)
library(ggfortify)
library(readxl)

plot.theme1 <- theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "white",
                                                     colour = "black",
                                                     size = 0.5, linetype = "solid"),
                     panel.border= element_rect(fill=NA,size = 0.5, linetype = 'solid',colour = "black"),
                     axis.text.x = element_text(size=13),axis.text.y = element_text(size=13),legend.text = element_text(size=13),
                     axis.title = element_text(size=14),
                     legend.title = element_text(color = "black", size = 14),
                     strip.text.x = element_text(size=14),
                     strip.background = element_rect(colour="black", fill="white"))

plot.theme2 <- theme_classic() +
  theme(text=element_text(size=15),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        axis.title.x = element_text(size = rel(1.2), angle = 00, margin = margin(t=8)),
        axis.title.y = element_text(size = rel(1.2), angle = 90, margin = margin(t=8)),
        plot.title = element_text(size=22),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))


# IMPORT DATA ----
my_data <- read.csv("SP_metadata_2021.csv", sep=";")
my_data$aridity <- (1- my_data$AI)
#To replace NA values with a mean of the other values of the Site:
for (i in which(sapply(my_data, is.numeric))) {
  for (j in which(is.na(my_data[, i]))) {
    my_data[j, i] <- mean(my_data[my_data[, "Site"] == my_data[j, "Site"], i],  na.rm = TRUE)
  }
}

colnames(my_data)[colnames(my_data) == "alpha"] <- "ALPHA"
colnames(my_data)[colnames(my_data) == "beta"] <- "BETA"
colnames(my_data)[colnames(my_data) == "xyl"] <- "XYL"
colnames(my_data)[colnames(my_data) == "cbh"] <- "CBH"
colnames(my_data)[colnames(my_data) == "gla"] <- "NAG"
colnames(my_data)[colnames(my_data) == "fos"] <- "PHOS"
colnames(my_data)[colnames(my_data) == "leu"] <- "LEU"
colnames(my_data)[colnames(my_data) == "phe"] <- "PHE"
colnames(my_data)[colnames(my_data) == "aridity"] <- "Aridity"
colnames(my_data)[colnames(my_data) == "Name"] <- "Sample.name"



my_data <- my_data %>%
  mutate(Site = case_when(
    Site == "SP01" ~ "ATK",
    Site == "SP02" ~ "LHE",
    Site == "SP03" ~ "GAV",
    Site == "SP04" ~ "COY",
    Site == "SP05" ~ "TAB",
    Site == "SP06" ~ "VAL",
    Site == "SP07" ~ "ARZ",
    Site == "SP08" ~ "FDE",
    Site == "SP09" ~ "SAN",
    Site == "SP10" ~ "MAL",
    Site == "SP11" ~ "MON",
    Site == "SP12" ~ "ALB",
    TRUE ~ NA_character_  # In case there are other values not listed above
  ))



func <- my_data[,c(2,3,46:53,64,93)]


clima <- my_data[,c(2,3,93,5,6,76)]
soil_char <- my_data[,c(2,8:10,12:13,17:23,27:30,33:34,79:92)]
veg <- my_data[,c(2,77)]

env <- cbind(clima, veg)
env <- cbind(env, soil_char)
env <- env[,-c(7,9)]

rm(clima, soil_char, veg)


ASV_16S <- read_excel("ASV_table_16S_phylum.xlsx")
ASV_16S <- merge(ASV_16S, env[, c("Sample.name", "Aridity")], by = "Sample.name", all.x = TRUE)
ASV_16S <- ASV_16S[, c("Sample.name", "Aridity", setdiff(names(ASV_16S), c("Sample.name", "Aridity")))]


ASV_ITS2 <- read_excel("ASV_table_ITS2_class.xlsx")
ASV_ITS2 <- merge(ASV_ITS2, env[, c("Sample.name", "Aridity")], by = "Sample.name", all.x = TRUE)
ASV_ITS2 <- ASV_ITS2[, c("Sample.name", "Aridity", setdiff(names(ASV_ITS2), c("Sample.name", "Aridity")))]


# EDIT DATA ----

## Functions ----
func <- func[order(func$Aridity), ]

func$rep <- sub(".*(-\\w)", "\\1\\2", func$Sample.name)
func <- func[, c("rep", setdiff(names(func), "rep"))]

func$Site.rep <- paste(func$Site, func$rep, sep = "")
func <- func[, c("Site.rep", setdiff(names(func), "Site.rep"))]

func <- func[,-c(2:4,16)]

rownames(func) <- func$Site.rep
func <- func[,-c(1,11)]



## Environmental ----

env$rep <- sub(".*(-\\w)", "\\1\\2", env$Sample.name)
env <- env[, c("rep", setdiff(names(env), "rep"))]

env$Site.rep <- paste(env$Site, env$rep, sep = "")
env <- env[, c("Site.rep", setdiff(names(env), "Site.rep"))]

env <- env[,-c(2:4)]

rownames(env) <- env$Site.rep
env <- env[,-1]

env <- env[order(env$Aridity), ]

colnames(env)[colnames(env) == "Soil_Temp"] <- "Soil Temperature"
colnames(env)[colnames(env) == "Water_activity"] <- "Water activity"
colnames(env)[colnames(env) == "Water_content"] <- "Water content"
colnames(env)[colnames(env) == "C_N"] <- "C/N"

env <- env[,-5]



## 16S ----
ASV_16S$Site <- sub(".*-(SP)(\\d+).*", "\\1\\2", ASV_16S$Sample.name)  # First extract the part
ASV_16S$Site <- ifelse(as.numeric(gsub("SP", "", ASV_16S$Site)) < 10,
                       sub("SP(\\d+)", "SP0\\1", ASV_16S$Site),
                       ASV_16S$Site)
ASV_16S <- ASV_16S[, c("Site", setdiff(names(ASV_16S), "Site"))]
# ASV_16S <- ASV_16S[,-2]


ASV_16S <- ASV_16S %>%
  mutate(Site = case_when(
    Site == "SP01" ~ "ATK",
    Site == "SP02" ~ "LHE",
    Site == "SP03" ~ "GAV",
    Site == "SP04" ~ "COY",
    Site == "SP05" ~ "TAB",
    Site == "SP06" ~ "VAL",
    Site == "SP07" ~ "ARZ",
    Site == "SP08" ~ "FDE",
    Site == "SP09" ~ "SAN",
    Site == "SP10" ~ "MAL",
    Site == "SP11" ~ "MON",
    Site == "SP12" ~ "ALB",
    TRUE ~ NA_character_  # In case there are other values not listed above
  ))

ASV_16S$rep <- sub(".*(-\\w)", "\\1\\2", ASV_16S$Sample.name)
ASV_16S <- ASV_16S[, c("rep", setdiff(names(ASV_16S), "rep"))]

ASV_16S$Site.rep <- paste(ASV_16S$Site, ASV_16S$rep, sep = "")
ASV_16S <- ASV_16S[, c("Site.rep", setdiff(names(ASV_16S), "Site.rep"))]

ASV_16S <- ASV_16S[order(ASV_16S$Aridity), ]


ASV_16S <- ASV_16S[,-c(2:5)]


ASV_16S <- as.data.frame(ASV_16S)
# rownames(ASV_16S) <- ASV_16S$Sample.name
rownames(ASV_16S) <- ASV_16S$Site.rep
ASV_16S <- ASV_16S[,-1]




## ITS2 ----

ASV_ITS2$Site <- sub(".*-(SP)(\\d+).*", "\\1\\2", ASV_ITS2$Sample.name)  # First extract the part
ASV_ITS2$Site <- ifelse(as.numeric(gsub("SP", "", ASV_ITS2$Site)) < 10,
                        sub("SP(\\d+)", "SP0\\1", ASV_ITS2$Site),
                        ASV_ITS2$Site)
ASV_ITS2 <- ASV_ITS2[, c("Site", setdiff(names(ASV_ITS2), "Site"))]


ASV_ITS2 <- ASV_ITS2 %>%
  mutate(Site = case_when(
    Site == "SP01" ~ "ATK",
    Site == "SP02" ~ "LHE",
    Site == "SP03" ~ "GAV",
    Site == "SP04" ~ "COY",
    Site == "SP05" ~ "TAB",
    Site == "SP06" ~ "VAL",
    Site == "SP07" ~ "ARZ",
    Site == "SP08" ~ "FDE",
    Site == "SP09" ~ "SAN",
    Site == "SP10" ~ "MAL",
    Site == "SP11" ~ "MON",
    Site == "SP12" ~ "ALB",
    TRUE ~ NA_character_  # In case there are other values not listed above
  ))

ASV_ITS2$rep <- sub(".*(-\\w)", "\\1\\2", ASV_ITS2$Sample.name)
ASV_ITS2 <- ASV_ITS2[, c("rep", setdiff(names(ASV_ITS2), "rep"))]

ASV_ITS2$Site.rep <- paste(ASV_ITS2$Site, ASV_ITS2$rep, sep = "")
ASV_ITS2 <- ASV_ITS2[, c("Site.rep", setdiff(names(ASV_ITS2), "Site.rep"))]

ASV_ITS2 <- ASV_ITS2[order(ASV_ITS2$Aridity), ]

ASV_ITS2 <- ASV_ITS2[,-c(2:5)]



ASV_ITS2 <- as.data.frame(ASV_ITS2)
rownames(ASV_ITS2) <- ASV_ITS2$Site.rep
ASV_ITS2 <- ASV_ITS2[,-1]




## Scaling ----


#Standardazing (scale) for env and func
env_s <- scale(env)
func_s <- scale(func)
env_s <- as.data.frame(env_s)
func_s <- as.data.frame(func_s)

rm(env, func)




# 16s vs func ----
rankindex(func_s, ASV_16S, indices = c("euc", "man", "gow","bra", "kul"), stepacross= FALSE, method = "spearman")

# prok_h <- vegdist(ASV_16S, method = "bray")
# dbRDA = capscale(prok_h ~ ., func_s, dist="bray")
# dbRDA <- dbrda(ASV_16S ~ ., func_s, dist="bray")
# plot(dbRDA)


dbRDA = capscale(ASV_16S ~ ., func_s, dist="bray")
plot(dbRDA)

anova(dbRDA)  # overall test of the significant of the analysis
# anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=1000) # test for sign. FUNC variables





func_s_scores <- vegan::scores(dbRDA, display = "bp") %>% as.data.frame()
site_scores <- vegan::scores(dbRDA, display = "sites") %>% as.data.frame()



site_scores$SiteRowNames <- rownames(site_scores)
site_scores$Site <- rownames(site_scores)

site_scores <- site_scores %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores <- dplyr::left_join(site_scores, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores) <- site_scores$SiteRowNames

site_scores$Aridity <- (1 - site_scores$AI)
site_scores <- site_scores[,-c(3,5:6)]



species_scores <- vegan::scores(dbRDA, display = "species") %>% as.data.frame()

species_scores$var <- rownames(species_scores)



perc <- round(100*(summary(dbRDA)$cont$importance[2, 1:2]), 2)
perc



ggplot() +
  geom_point(data = site_scores, 
             aes(x = CAP1, y = CAP2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores, 
            aes(x = CAP1, y = CAP2, label = Site), 
            size = 3.5, color = "#404040", 
            nudge_y = -0.09,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores, 
                  aes(x = CAP1, y = CAP2, label = rownames(func_s_scores)), 
                  size = 5, color = "blue",
                  fontface = "bold") +
  # geom_text_repel(data = species_scores,  # Use rda_mcc_df for text labels
  #                 aes(x = CAP1, y = CAP2, label = var),  # Ensure 'var' is the correct column for labels
  #                 size = 4, color = "darkred", vjust = -1,
  #                 fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  labs(x = paste0("CAP1 (", perc[1], "%)"), y = paste0("CAP2 (", perc[2], "%)"), tag = "A")+
  # xlab(paste0("CAP1 (", perc[1], "%)")) + 
  # ylab(paste0("CAP2 (", perc[2], "%)")) +
  plot.theme2+
  annotate(geom="text", x=0.8, y=1.2, label="PERMANOVA ~ Site \n F = 6.21 \n p-value < 0.01",
           color="black", size = 5.5)


# ggsave(path = "Figures",
#        "dbRDA_16s_func.png", width = 10, height = 8, dpi = 300)



# summary(ord1.df)
# summary(site_scores)
# 
# combined_scores <- rbind(
#   data.frame(Axis1 = ord1.df$Axis.1, Axis2 = ord1.df$Axis.2, Method = "PCoA"),
#   data.frame(Axis1 = site_scores$CAP1, Axis2 = site_scores$CAP2, Method = "dbRDA")
# )
# 
# ggplot(combined_scores, aes(x = Axis1, y = Axis2, color = Method)) +
#   geom_point() +
#   labs(x = "Axis 1", y = "Axis 2", title = "Comparison of PCoA and dbRDA") +
#   theme_minimal()








# ITS2 vs func ----
rankindex(func_s, ASV_ITS2, indices = c("euc", "man", "gow","bra", "kul"), stepacross= FALSE, method = "spearman")

dbRDA = capscale(ASV_ITS2 ~ ., func_s, dist="bray")
plot(dbRDA)

anova(dbRDA)  # overall test of the significant of the analysis
# anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=1000) # test for sign. FUNC variables





func_s_scores <- vegan::scores(dbRDA, display = "bp") %>% as.data.frame()
site_scores <- vegan::scores(dbRDA, display = "sites") %>% as.data.frame()



site_scores$SiteRowNames <- rownames(site_scores)
site_scores$Site <- rownames(site_scores)

site_scores <- site_scores %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores <- dplyr::left_join(site_scores, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores) <- site_scores$SiteRowNames

site_scores$Aridity <- (1 - site_scores$AI)
site_scores <- site_scores[,-c(3,5:6)]



species_scores <- vegan::scores(dbRDA, display = "species") %>% as.data.frame()

species_scores$var <- rownames(species_scores)



perc <- round(100*(summary(dbRDA)$cont$importance[2, 1:2]), 2)
perc



ggplot() +
  geom_point(data = site_scores, 
             aes(x = CAP1, y = CAP2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores, 
            aes(x = CAP1, y = CAP2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.09,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores, 
                  aes(x = CAP1, y = CAP2, label = rownames(func_s_scores)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  # geom_text_repel(data = species_scores,  # Use rda_mcc_df for text labels
  #                 aes(x = CAP1, y = CAP2, label = var),  # Ensure 'var' is the correct column for labels
  #                 size = 4, color = "darkred", vjust = -1,
  #                 fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab(paste0("CAP1 (", perc[1], "%)")) + 
  ylab(paste0("CAP2 (", perc[2], "%)")) +
  plot.theme1+
  annotate(geom="text", x=-1.3, y=1.5, label= "PERMANOVA ~ Site \n F = 3.97 \n p-value < 0.01",
           color="black", size =  5.5)


# ggsave(path = "Figures", "dbRDA_ITS2_func.png", width = 10, height = 8, dpi = 300)




# Without NAG

func_s_scores2 <- func_s_scores[-5,]

ggplot() +
  geom_point(data = site_scores, 
             aes(x = CAP1, y = CAP2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores, 
            aes(x = CAP1, y = CAP2, label = Site), 
            size = 3.5, color = "#404040", 
            nudge_y = -0.09,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores2, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores2, 
                  aes(x = CAP1, y = CAP2, label = rownames(func_s_scores2)), 
                  size = 5, color = "blue",
                  fontface = "bold") +
  # geom_text_repel(data = species_scores,  # Use rda_mcc_df for text labels
  #                 aes(x = CAP1, y = CAP2, label = var),  # Ensure 'var' is the correct column for labels
  #                 size = 4, color = "darkred", vjust = -1,
  #                 fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  labs(x = paste0("CAP1 (", perc[1], "%)"), y = paste0("CAP2 (", perc[2], "%)"), tag = "B")+
  # xlab(paste0("CAP1 (", perc[1], "%)")) + 
  # ylab(paste0("CAP2 (", perc[2], "%)")) +
  plot.theme2+
  annotate(geom="text", x=-1.3, y=1.5, label= "PERMANOVA ~ Site \n F = 3.97 \n p-value < 0.01",
           color="black", size =  5.5)


# ggsave(path = "Figures", "dbRDA_ITS2_func_significative.png", width = 10, height = 8, dpi = 300)




# 16S vs env ----
rankindex(env_s, ASV_16S, indices = c("euc", "man", "gow","bra", "kul"), stepacross= FALSE, method = "spearman")

# prok_h <- vegdist(ASV_16S, method = "bray")
# dbRDA = capscale(prok_h ~ ., func_s, dist="bray")
# dbRDA <- dbrda(ASV_16S ~ ., func_s, dist="bray")
# plot(dbRDA)


dbRDA = capscale(ASV_16S ~ ., env_s, dist="bray")
plot(dbRDA)

anova(dbRDA)  # overall test of the significant of the analysis
# anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=1000) # test for sign. FUNC variables





env_s_scores <- vegan::scores(dbRDA, display = "bp") %>% as.data.frame()
site_scores <- vegan::scores(dbRDA, display = "sites") %>% as.data.frame()



site_scores$SiteRowNames <- rownames(site_scores)
site_scores$Site <- rownames(site_scores)

site_scores <- site_scores %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores <- dplyr::left_join(site_scores, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores) <- site_scores$SiteRowNames

site_scores$Aridity <- (1 - site_scores$AI)
site_scores <- site_scores[,-c(3,5:6)]



species_scores <- vegan::scores(dbRDA, display = "species") %>% as.data.frame()

species_scores$var <- rownames(species_scores)



perc <- round(100*(summary(dbRDA)$cont$importance[2, 1:2]), 2)
perc



ggplot() +
  geom_point(data = site_scores, 
             aes(x = CAP1, y = CAP2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores, 
            aes(x = CAP1, y = CAP2, label = Site), 
            size = 3.5, color = "#404040", 
            nudge_y = -0.09,  
            nudge_x = 0.01) + 
  geom_segment(data = env_s_scores, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = env_s_scores, 
                  aes(x = CAP1, y = CAP2, label = rownames(env_s_scores)), 
                  size = 5, color = "blue",
                  fontface = "bold") +
  # geom_text_repel(data = species_scores,  # Use rda_mcc_df for text labels
  #                 aes(x = CAP1, y = CAP2, label = var),  # Ensure 'var' is the correct column for labels
  #                 size = 4, color = "darkred", vjust = -1,
  #                 fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  labs(x = paste0("CAP1 (", perc[1], "%)"), y = paste0("CAP2 (", perc[2], "%)"), tag = "A")+
  # xlab(paste0("CAP1 (", perc[1], "%)")) + 
  # ylab(paste0("CAP2 (", perc[2], "%)")) +
  plot.theme2




# Without non-significative:

env_s_scores2 <- env_s_scores[c(1:9,11:12,15:17,19,21,24),]

rownames(env_s_scores2)[rownames(env_s_scores2) == "altitude"] <- "Altitude"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Soil Temperature`"] <- "STemp"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Water activity`"] <- "WA"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Water content`"] <- "WC"
rownames(env_s_scores2)[rownames(env_s_scores2) == "L_TC"] <- "LTC"
rownames(env_s_scores2)[rownames(env_s_scores2) == "s350.400"] <- "S350-400"



ggplot() +
  geom_point(data = site_scores, 
             aes(x = CAP1, y = CAP2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores, 
            aes(x = CAP1, y = CAP2, label = Site), 
            size = 3.5, color = "#404040", 
            nudge_y = -0.09,  
            nudge_x = 0.01) + 
  geom_segment(data = env_s_scores2, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = env_s_scores2, 
                  aes(x = CAP1, y = CAP2, label = rownames(env_s_scores2)), 
                  size = 5, color = "blue",
                  fontface = "bold") +
  # geom_text_repel(data = species_scores,  # Use rda_mcc_df for text labels
  #                 aes(x = CAP1, y = CAP2, label = var),  # Ensure 'var' is the correct column for labels
  #                 size = 4, color = "darkred", vjust = -1,
  #                 fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  labs(x = paste0("CAP1 (", perc[1], "%)"), y = paste0("CAP2 (", perc[2], "%)"), tag = "A")+
  # xlab(paste0("CAP1 (", perc[1], "%)")) + 
  # ylab(paste0("CAP2 (", perc[2], "%)")) +
  plot.theme2


# ggsave(path = "Figures", "dbRDA_16S_env_significative.png", width = 10, height = 8, dpi = 300)




# With significative by 0.01:

env_s_scores2 <- env_s_scores[c(1:9,11:12,15,17,21,24),]

rownames(env_s_scores2)[rownames(env_s_scores2) == "altitude"] <- "Altitude"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Soil Temperature`"] <- "STemp"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Water activity`"] <- "WA"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Water content`"] <- "WC"
rownames(env_s_scores2)[rownames(env_s_scores2) == "L_TC"] <- "LTC"
rownames(env_s_scores2)[rownames(env_s_scores2) == "s350.400"] <- "S350-400"



ggplot() +
  geom_point(data = site_scores, 
             aes(x = CAP1, y = CAP2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores, 
            aes(x = CAP1, y = CAP2, label = Site), 
            size = 3.5, color = "#404040", 
            nudge_y = -0.09,  
            nudge_x = 0.01) + 
  geom_segment(data = env_s_scores2, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = env_s_scores2, 
                  aes(x = CAP1, y = CAP2, label = rownames(env_s_scores2)), 
                  size = 5, color = "blue",
                  fontface = "bold") +
  # geom_text_repel(data = species_scores,  # Use rda_mcc_df for text labels
  #                 aes(x = CAP1, y = CAP2, label = var),  # Ensure 'var' is the correct column for labels
  #                 size = 4, color = "darkred", vjust = -1,
  #                 fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  labs(x = paste0("CAP1 (", perc[1], "%)"), y = paste0("CAP2 (", perc[2], "%)"), tag = "A")+
  # xlab(paste0("CAP1 (", perc[1], "%)")) + 
  # ylab(paste0("CAP2 (", perc[2], "%)")) +
  plot.theme2


# ggsave(path = "Figures", "dbRDA_16S_env_sig_0.01.png", width = 10, height = 8, dpi = 300)





# ITS2 vs env ----
rankindex(env_s, ASV_ITS2, indices = c("euc", "man", "gow","bra", "kul"), stepacross= FALSE, method = "spearman")

# prok_h <- vegdist(ASV_16S, method = "bray")
# dbRDA = capscale(prok_h ~ ., func_s, dist="bray")
# dbRDA <- dbrda(ASV_16S ~ ., func_s, dist="bray")
# plot(dbRDA)


dbRDA = capscale(ASV_ITS2 ~ ., env_s, dist="bray")
plot(dbRDA)

anova(dbRDA)  # overall test of the significant of the analysis
# anova(dbRDA, by="axis", perm.max=500) # test axes for significance
anova(dbRDA, by="terms", permu=1000) # test for sign. FUNC variables





env_s_scores <- vegan::scores(dbRDA, display = "bp") %>% as.data.frame()
site_scores <- vegan::scores(dbRDA, display = "sites") %>% as.data.frame()



site_scores$SiteRowNames <- rownames(site_scores)
site_scores$Site <- rownames(site_scores)

site_scores <- site_scores %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores <- dplyr::left_join(site_scores, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores) <- site_scores$SiteRowNames

site_scores$Aridity <- (1 - site_scores$AI)
site_scores <- site_scores[,-c(3,5:6)]



species_scores <- vegan::scores(dbRDA, display = "species") %>% as.data.frame()

species_scores$var <- rownames(species_scores)



perc <- round(100*(summary(dbRDA)$cont$importance[2, 1:2]), 2)
perc



ggplot() +
  geom_point(data = site_scores, 
             aes(x = CAP1, y = CAP2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores, 
            aes(x = CAP1, y = CAP2, label = Site), 
            size = 3.5, color = "#404040", 
            nudge_y = -0.09,  
            nudge_x = 0.01) + 
  geom_segment(data = env_s_scores, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = env_s_scores, 
                  aes(x = CAP1, y = CAP2, label = rownames(env_s_scores)), 
                  size = 5, color = "blue",
                  fontface = "bold") +
  # geom_text_repel(data = species_scores,  # Use rda_mcc_df for text labels
  #                 aes(x = CAP1, y = CAP2, label = var),  # Ensure 'var' is the correct column for labels
  #                 size = 4, color = "darkred", vjust = -1,
  #                 fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  labs(x = paste0("CAP1 (", perc[1], "%)"), y = paste0("CAP2 (", perc[2], "%)"), tag = "A")+
  # xlab(paste0("CAP1 (", perc[1], "%)")) + 
  # ylab(paste0("CAP2 (", perc[2], "%)")) +
  plot.theme2




# Without non-significative:

env_s_scores2 <- env_s_scores[-c(10,13,18,22:23,25:32,34:36),]

rownames(env_s_scores2)[rownames(env_s_scores2) == "altitude"] <- "Altitude"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Soil Temperature`"] <- "STemp"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Water activity`"] <- "WA"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Water content`"] <- "WC"
rownames(env_s_scores2)[rownames(env_s_scores2) == "L_TC"] <- "LTC"
rownames(env_s_scores2)[rownames(env_s_scores2) == "s350.400"] <- "S350-400"
rownames(env_s_scores2)[rownames(env_s_scores2) == "Peak_M"] <- "Peak M"



ggplot() +
  geom_point(data = site_scores, 
             aes(x = CAP1, y = CAP2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores, 
            aes(x = CAP1, y = CAP2, label = Site), 
            size = 3.5, color = "#404040", 
            nudge_y = -0.09,  
            nudge_x = 0.01) + 
  geom_segment(data = env_s_scores2, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = env_s_scores2, 
                  aes(x = CAP1, y = CAP2, label = rownames(env_s_scores2)), 
                  size = 5, color = "blue",
                  fontface = "bold") +
  # geom_text_repel(data = species_scores,  # Use rda_mcc_df for text labels
  #                 aes(x = CAP1, y = CAP2, label = var),  # Ensure 'var' is the correct column for labels
  #                 size = 4, color = "darkred", vjust = -1,
  #                 fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  labs(x = paste0("CAP1 (", perc[1], "%)"), y = paste0("CAP2 (", perc[2], "%)"), tag = "A")+
  # xlab(paste0("CAP1 (", perc[1], "%)")) + 
  # ylab(paste0("CAP2 (", perc[2], "%)")) +
  plot.theme2


# ggsave(path = "Figures", "dbRDA_ITS2_env_significative.png", width = 10, height = 8, dpi = 300)




# With significative by 0.01:

env_s_scores2 <- env_s_scores[-c(10,13,18,20,22:23,25:32,34:36),]

rownames(env_s_scores2)[rownames(env_s_scores2) == "altitude"] <- "Altitude"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Soil Temperature`"] <- "STemp"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Water activity`"] <- "WA"
rownames(env_s_scores2)[rownames(env_s_scores2) == "`Water content`"] <- "WC"
rownames(env_s_scores2)[rownames(env_s_scores2) == "L_TC"] <- "LTC"
rownames(env_s_scores2)[rownames(env_s_scores2) == "s350.400"] <- "S350-400"
rownames(env_s_scores2)[rownames(env_s_scores2) == "Peak_M"] <- "Peak M"



ggplot() +
  geom_point(data = site_scores, 
             aes(x = CAP1, y = CAP2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores, 
            aes(x = CAP1, y = CAP2, label = Site), 
            size = 3.5, color = "#404040", 
            nudge_y = -0.09,  
            nudge_x = 0.01) + 
  geom_segment(data = env_s_scores2, 
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = env_s_scores2, 
                  aes(x = CAP1, y = CAP2, label = rownames(env_s_scores2)), 
                  size = 5, color = "blue",
                  fontface = "bold") +
  # geom_text_repel(data = species_scores,  # Use rda_mcc_df for text labels
  #                 aes(x = CAP1, y = CAP2, label = var),  # Ensure 'var' is the correct column for labels
  #                 size = 4, color = "darkred", vjust = -1,
  #                 fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  labs(x = paste0("CAP1 (", perc[1], "%)"), y = paste0("CAP2 (", perc[2], "%)"), tag = "B")+
  # xlab(paste0("CAP1 (", perc[1], "%)")) + 
  # ylab(paste0("CAP2 (", perc[2], "%)")) +
  plot.theme2


# ggsave(path = "Figures", "dbRDA_ITS2_env_sig_0.01.png", width = 10, height = 8, dpi = 300)







#__________________________________________----

phy2 <- phyRar

phy2@sam_data

rownames(phy2@sam_data) <- phy2@sam_data$Sample.name

# Access the sample data and convert it to a data frame for manipulation
sample_data_phy2 <- data.frame(sample_data(phy2))

# Step 1: Extract the Site code from Sample.name and format single-digit numbers
sample_data_phy2$Site <- sub(".*-(SP)(\\d+).*", "\\1\\2", sample_data_phy2$Sample.name)
sample_data_phy2$Site <- ifelse(as.numeric(gsub("SP", "", sample_data_phy2$Site)) < 10,
                                sub("SP(\\d+)", "SP0\\1", sample_data_phy2$Site),
                                sample_data_phy2$Site)

# Step 2: Rename the Site codes
sample_data_phy2 <- sample_data_phy2 %>%
  mutate(Site = case_when(
    Site == "SP01" ~ "ATK",
    Site == "SP02" ~ "LHE",
    Site == "SP03" ~ "GAV",
    Site == "SP04" ~ "COY",
    Site == "SP05" ~ "TAB",
    Site == "SP06" ~ "VAL",
    Site == "SP07" ~ "ARZ",
    Site == "SP08" ~ "FDE",
    Site == "SP09" ~ "SAN",
    Site == "SP10" ~ "MAL",
    Site == "SP11" ~ "MON",
    Site == "SP12" ~ "ALB",
    TRUE ~ NA_character_
  ))

# Step 3: Extract the replicate information and add Site.rep column
sample_data_phy2$rep <- sub(".*(-\\w)", "\\1", sample_data_phy2$Sample.name)  # Corrected the regex to avoid capture group issues
sample_data_phy2$Site.rep <- paste(sample_data_phy2$Site, sample_data_phy2$rep, sep = "")

# # Step 4: Sort by Aridity and select desired columns
# sample_data_phy2$Aridity <- (1 - sample_data_phy2$AI)
# sample_data_phy2 <- sample_data_phy2[order(sample_data_phy2$Aridity), ]
# sample_data_phy2 <- sample_data_phy2[, c("Site.rep", setdiff(names(sample_data_phy2), c("Sample.name", "Site", "rep", "Site.rep")))]

# Step 5: Set rownames to Site.rep and remove unnecessary columns
rownames(sample_data_phy2) <- sample_data_phy2$Site.rep
sample_data_phy2 <- sample_data_phy2[, -1]

# Convert the modified data frame back to a sample_data object
new_sample_data_phy2 <- sample_data(sample_data_phy2)

# Update the sample_data in the phyloseq object with the transformed sample data
new_sample_data_phy2 <- sample_data(new_sample_data_phy2)

sample_data(phy2) <- new_sample_data_phy2



# Get sample names from the phyloseq object
existing_sample_names <- sample_names(phy2)

# Get sample names from the new sample data
new_sample_names <- rownames(new_sample_data_phy2)

# Print sample names to compare
print(existing_sample_names)
print(new_sample_names)

sample_names(phy2) <- new_sample_names

phy2@sam_data
phy2@otu_table



# 
# 
# 
# sample_data_phy2 <- sample_data(phy2)
# sample_data_phy2$Aridity <- (1 - sample_data_phy2$AI)
# sample_data_phy2 <- sample_data_phy2[order(sample_data_phy2$Aridity), ]
# 
# sample_data(phy2) <- sample_data_phy2
# 
# 
# phy2 <- phyloseq::prune_samples(sample_names(phy2)[order(sample_data_phy2$Aridity)], phy2)
# 
# phy2@sam_data




# Get sample names from the phyloseq object
phylo_sample_names <- sample_names(phy2)

# Get sample names from the OTU table and tax table
otu_sample_names <- colnames(otu_table(phy2))
tax_sample_names <- colnames(tax_table(phy2))

# Print the sample names to compare
print(phylo_sample_names)
print(otu_sample_names)
print(tax_sample_names)



# Extract the sample data
sample_data_phy2 <- sample_data(phy2)

# Create the Aridity column
sample_data_phy2$Aridity <- 1 - sample_data_phy2$AI

# Sort the sample data by Aridity in decreasing order
sorted_sample_data <- sample_data_phy2[order(sample_data_phy2$Aridity), ]





# Get the current OTU table and tax table
otu_table_phy2 <- otu_table(phy2)
tax_table_phy2 <- tax_table(phy2)

# Get the names of the samples from the sorted sample data
sorted_sample_names <- rownames(sorted_sample_data)

# Check if the sorted sample names are in the otu_table and tax_table
valid_otu_names <- intersect(sorted_sample_names, colnames(otu_table_phy2))
# valid_tax_names <- intersect(sorted_sample_names, colnames(tax_table_phy2))

# Print valid names to debug
print(valid_otu_names)
# print(valid_tax_names)

# Reorder the OTU table and tax table using valid sample names
otu_table_ordered <- otu_table_phy2[, valid_otu_names]
# tax_table_ordered <- tax_table_phy2[, valid_tax_names]

# Create a new phyloseq object with ordered sample data
phy2_ord <- phyloseq(
  otu_table(otu_table_ordered, taxa_are_rows = TRUE),
  tax_table(tax_table_phy2),
  sample_data(sorted_sample_data[match(colnames(otu_table_ordered), rownames(sorted_sample_data)), ])
)

# Validate the new phyloseq object
validObject(phy2_ord)


phy2_ord@sam_data
phy2_ord@otu_table
phy2_ord@tax_table





# dist_matrix <- phyloseq::distance(phy2, method = "bray")
# 
# 
# 
# 
# dbRDA_result <- capscale(dist_matrix ~ ., data = func_s)
# summary(dbRDA_result)
# 
# plot(dbRDA_result)
# 
# func_s_scores <- vegan::scores(dbRDA_result, display = "bp") %>% as.data.frame()
# site_scores <- vegan::scores(dbRDA_result, display = "sites") %>% as.data.frame()
# 
# 
# site_scores$SiteRowNames <- rownames(site_scores)
# site_scores$Site <- rownames(site_scores)
# 
# site_scores <- site_scores %>%
#   tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
# AI_df <- data.frame(
#   Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
#   AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
# )
# site_scores <- dplyr::left_join(site_scores, AI_df, by = "Site")
# rm(AI_df)
# rownames(site_scores) <- site_scores$SiteRowNames
# 
# site_scores$Aridity <- (1 - site_scores$AI)
# site_scores <- site_scores[,-c(3,5:6)]
# 
# 
# 
# species_scores <- vegan::scores(dbRDA_result, display = "species") %>% as.data.frame()
# 
# species_scores$var <- rownames(species_scores)
# 
# 
# perc <- round(100*(summary(dbRDA_result)$cont$importance[2, 1:2]), 2)
# perc
# 
# 
# 
# ggplot() +
#   geom_point(data = site_scores, 
#              aes(x = CAP1, y = CAP2, color = Aridity),
#              size = 3, show.legend = TRUE) + 
#   geom_text(data = site_scores, 
#             aes(x = CAP1, y = CAP2, label = Site), 
#             size = 3, color = "#404040", 
#             nudge_y = -0.12,  
#             nudge_x = 0.01) + 
#   geom_segment(data = func_s_scores, 
#                aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
#                arrow = arrow(length = unit(0.2, "inches")), 
#                color = "blue") +
#   geom_text_repel(data = func_s_scores, 
#                   aes(x = CAP1, y = CAP2, label = rownames(func_s_scores)), 
#                   size = 4, color = "blue",
#                   fontface = "bold") +
#   # geom_text_repel(data = species_scores,  # Use rda_mcc_df for text labels
#   #                 aes(x = CAP1, y = CAP2, label = var),  # Ensure 'var' is the correct column for labels
#   #                 size = 4, color = "darkred", vjust = -1,
#   #                 fontface = "bold") +
#   scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
#   xlab(paste0("CAP1 (", perc[1], "%)")) + 
#   ylab(paste0("CAP2 (", perc[2], "%)")) +
#   plot.theme1


bray_dist <- distance(phy2_ord, method = "bray")
dbRDA_result <- capscale(bray_dist ~ ., func_s)
summary(dbRDA_result)
plot(dbRDA_result)




