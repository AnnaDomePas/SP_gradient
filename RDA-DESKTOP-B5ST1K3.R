rm(list = ls())

library(ade4) 
library(vegan)  
library(gclus) 
library(cluster)
library(RColorBrewer)  
library(labdsv)
library(ape)
library(dplyr)
library(readxl)
# remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(ggplot2)
library(ggrepel)
library(tibble)
# install.packages("adespatial")
library(adespatial)


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



func <- my_data[,c(2,3,40,41,46:53,64,93)]

ASV_16S <- read_excel("ASV_table_16S_phylum.xlsx")

ASV_ITS2 <- read_excel("ASV_table_ITS2_class.xlsx")

clima <- my_data[,c(2,3,93,5,6,76)]
soil_char <- my_data[,c(2,8:10,12:13,17:23,27:30,33:34,79:92)]
veg <- my_data[,c(2,77)]

env <- cbind(clima, veg)
env <- cbind(env, soil_char)
env <- env[,-c(7,9)]

rm(clima, soil_char, veg,my_data)


ASV_16S <- merge(ASV_16S, env[, c("Sample.name", "Aridity")], by = "Sample.name", all.x = TRUE)
ASV_16S <- ASV_16S[, c("Sample.name", "Aridity", setdiff(names(ASV_16S), c("Sample.name", "Aridity")))]

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
func <- func[,-1]



# Functions WITHOUT BIOMASS:
func2 <- func[,-c(1:2)]





## Environmental ----

env$rep <- sub(".*(-\\w)", "\\1\\2", env$Sample.name)
env <- env[, c("rep", setdiff(names(env), "rep"))]

env$Site.rep <- paste(env$Site, env$rep, sep = "")
env <- env[, c("Site.rep", setdiff(names(env), "Site.rep"))]

env <- env[,-c(2:4)]

rownames(env) <- env$Site.rep
env <- env[,-1]

env <- env[order(env$Aridity), ]




env <- env %>%
  rename("Soil Temperature" = Soil_Temp,
         "Water activity" = Water_activity,
         "Water content" = Water_content,
         "C/N" = C_N)







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





# Step 1: Extract the phylum names from the column names
phylum_names <- sub("\\..*", "", colnames(ASV_16S))  # Remove everything after the first dot

# Step 2: Transpose the ASV table and add the phylum names as column names
ASV_16S_t <- t(ASV_16S)

# Step 3: Convert to a data frame for easier manipulation
ASV_16S_df <- as.data.frame(ASV_16S_t)

# Step 4: Add the phylum names as a new column for aggregation
ASV_16S_df$Phylum <- phylum_names

# Step 5: Summarize the ASV data by phylum (i.e., sum ASVs with the same phylum)
ASV_16S_grouped <- ASV_16S_df %>%
  group_by(Phylum) %>%
  summarise(across(everything(), sum))%>%
  as.data.frame()

# Step 6: Transpose the result back to the original shape
ASV_16S_grouped <- t(ASV_16S_grouped)
ASV_16S_grouped <- as.data.frame(ASV_16S_grouped)
colnames(ASV_16S_grouped) <- ASV_16S_grouped[1,]
ASV_16S_grouped <- ASV_16S_grouped[-1,]

# Convert all character columns to numeric
ASV_16S_grouped <- ASV_16S_grouped %>% 
  mutate(across(everything(), ~as.numeric(trimws(.))))


# Based on the barplot figure with the phylum abundances of > 1%, I remove
# manually the phylumns (columns) that I am not interested with.
ASV_16S_grouped_clean <- ASV_16S_grouped %>% 
  select(Acidobacteriota, Actinobacteriota, Chloroflexi,
         Crenarchaeota, Cyanobacteria, Firmicutes,
         Gemmatimonadota, Myxococcota, Patescibacteria,
         Planctomycetota, Proteobacteria, Verrucomicrobiota,
         `WPS-2`)



ASV_16S_grouped_clean <- ASV_16S_grouped_clean %>%
  rename(
    Actinomycetota = Actinobacteriota,
    Chloroflexota = Chloroflexi,
    Pseudomonadota = Proteobacteria,
    Bacillota = Firmicutes,
    Thermoproteota = Crenarchaeota
  )



rm(ASV_16S_grouped, ASV_16S, ASV_16S_df, ASV_16S_t)







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




# Step 1: Extract the phylum names from the column names
class_names <- sub("\\..*", "", colnames(ASV_ITS2))  # Remove everything after the first dot

# Step 2: Transpose the ASV table and add the phylum names as column names
ASV_ITS2_t <- t(ASV_ITS2)

# Step 3: Convert to a data frame for easier manipulation
ASV_ITS2_df <- as.data.frame(ASV_ITS2_t)

# Step 4: Add the phylum names as a new column for aggregation
ASV_ITS2_df$Class <- class_names

# Step 5: Summarize the ASV data by phylum (i.e., sum ASVs with the same phylum)
ASV_ITS2_grouped <- ASV_ITS2_df %>%
  group_by(Class) %>%
  summarise(across(everything(), sum))%>%
  as.data.frame()

# Step 6: Transpose the result back to the original shape
ASV_ITS2_grouped <- t(ASV_ITS2_grouped)
ASV_ITS2_grouped <- as.data.frame(ASV_ITS2_grouped)
colnames(ASV_ITS2_grouped) <- ASV_ITS2_grouped[1,]
ASV_ITS2_grouped <- ASV_ITS2_grouped[-1,]

# Convert all character columns to numeric
ASV_ITS2_grouped <- ASV_ITS2_grouped %>% 
  mutate(across(everything(), ~as.numeric(trimws(.))))


# Based on the barplot figure with the class abundances of > 1%, I remove
# manually the class (columns) that I am not interested with.
ASV_ITS2_grouped_clean <- ASV_ITS2_grouped %>% 
  select(Agaricomycetes, Dothideomycetes, Eurotiomycetes,
         Fungi_cls_Incertae_sedis, Geminibasidiomycetes,
         Geoglossomycetes, Lecanoromycetes, Leotiomycetes,
         Microbotryomycetes, Mortierellomycetes,
         Mucoromycetes, Mucoromycotina_cls_Incertae_sedis,
         Pezizomycetes, Pezizomycotina_cls_Incertae_sedis,
         Rozellomycota_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis,
         Saccharomycetes, Sordariomycetes, Tremellomycetes,
         Umbelopsidomycetes, `NA`)

rm(ASV_ITS2_grouped, ASV_ITS2, ASV_ITS2_df, ASV_ITS2_t)



ASV_ITS2_grouped_clean <- ASV_ITS2_grouped_clean %>%
  rename("Fungi (Incertae sedis)" = Fungi_cls_Incertae_sedis,
         "Mucoromycotina (Incertae sedis)" = Mucoromycotina_cls_Incertae_sedis,
         "Pezizomycotina (Incertae sedis)" = Pezizomycotina_cls_Incertae_sedis,
         "Rozellomycota (Incertae sedis)" = Rozellomycota_cls_Incertae_sedis,
         "Rozellomycotina (Incertae sedis)" = Rozellomycotina_cls_Incertae_sedis,
         "Others" = `NA`)








## Scaling and distances ----


#Standardazing (scale) for env and func, and hellinger distance for ASVs
env_s <- scale(env)
func_s <- scale(func)
func2_s <- scale(func2)


# prok_h <- decostand(ASV_16S, method = "hellinger")
prok_h <- decostand(ASV_16S_grouped_clean, method = "hellinger")
ITS2_h <- decostand(ASV_ITS2_grouped_clean, method = "hellinger")




rm(env, func, func2, ASV_16S_grouped_clean, ASV_ITS2_grouped_clean)













# Func vs 16S: Without biomass ----

rda13 <- rda(func2_s ~ ., data=data.frame(prok_h))
summary(rda13)

R2 <- RsquareAdj(rda13)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda13)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda13, permutations = how(nperm = 10000))
# significative
anova.cca(rda13, by = "axis", permutations = how(nperm = 1000))
# 2 first RDA significatives



plot(rda13, scaling=2, main='Triplot - scaling 2')



summary(rda13)


func2_s_scores <- vegan::scores(rda13, display = "bp")
func2_s_scores_df <- as.data.frame(func2_s_scores)
rm(func2_s_scores)

site_scores <- vegan::scores(rda13, display = "sites") 
site_scores_df <- as.data.frame(site_scores)
rm(site_scores)
site_scores_df$SiteRowNames <- rownames(site_scores_df)
site_scores_df$Site <- rownames(site_scores_df)
site_scores_df <- site_scores_df %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores_df <- dplyr::left_join(site_scores_df, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores_df) <- site_scores_df$SiteRowNames
site_scores_df$Aridity <- (1 - site_scores_df$AI)
site_scores_df <- site_scores_df[,-c(3,5:6)]

species_scores <- vegan::scores(rda13, display = "species")
species_scores_df <- as.data.frame(species_scores)
species_scores_df$var <- rownames(species_scores_df)
rm(species_scores)

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.12,  
            nudge_x = 0.01) + 
  geom_segment(data = func2_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func2_s_scores_df, 
                  aes(x = RDA1, y = RDA2, label = rownames(func2_s_scores_df)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()


# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_FUNC_16S_NOBIOMASS.png", width = 10, height = 8, dpi = 300)





# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(func2_s, prok_h,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, # This means that the forward selection will only add variables that explain at least 1% more variance to the model.
            alpha = 0.05) # The significance level used for testing each variable's inclusion into the model.
ordistep(rda13,direction="forward",perm.max = 200)


func2_s_scores_df_sig <- func2_s_scores_df[c(3,7,9,12:13),]
prok_h1 <- prok_h[,c(3,7,9,12:13)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.08,  
            nudge_x = 0.01) + 
  geom_segment(data = func2_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func2_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func2_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()

# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_FUNC_16S_sel_NOBIOMASS.png", width = 10, height = 8, dpi = 300)

rda131 <- rda(func2_s ~ ., data=data.frame(prok_h1))
summary(rda131)

R2 <- RsquareAdj(rda131)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda131)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda131, permutations = how(nperm = 10000))





# !!!!!! LESS RESTRICTED FORWARD SELECTION !!!!!!

# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(func2_s, prok_h,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, 
            alpha = 0.1) #INCREASED ALPHA


func2_s_scores_df_sig <- func2_s_scores_df[c(1:3,12:13),]
prok_h1 <- prok_h[,c(1:3,12:13)]


ggplot() +
  geom_point(data = site_scores_df,
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) +
  geom_text(data = site_scores_df,
            aes(x = RDA1, y = RDA2, label = Site),
            size = 3, color = "#404040",
            nudge_y = 0.12,
            nudge_x = 0.01) +
  geom_segment(data = func2_s_scores_df_sig,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "inches")),
               color = "blue") +
  geom_text_repel(data = func2_s_scores_df_sig,
                  aes(x = RDA1, y = RDA2, label = rownames(func2_s_scores_df_sig)),
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()

# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_FUNC_16S_sel2_NOBIOMASS.png", width = 10, height = 8, dpi = 300)

rda1312 <- rda(func2_s ~ ., data=data.frame(prok_h1))
summary(rda1312)

R2 <- RsquareAdj(rda1312)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda1312)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda1312, permutations = how(nperm = 10000))




rm(func2_s_scores_df, func2_s_scores_df_sig, site_scores_df, species_scores_df,
   species_scores_df_sig, prok_h1)

rm(rda13, rda131, rda1312)















# Func vs ITS2: Without biomass ----

rda222 <- rda(func2_s ~ ., data=data.frame(ITS2_h))
summary(rda222)

R2 <- RsquareAdj(rda222)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda222)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda222, permutations = how(nperm = 10000))
# significative
anova.cca(rda222, by = "axis", permutations = how(nperm = 1000))
# 2 significative RDA axis


plot(rda222, scaling=2, main='Triplot - scaling 2')

func2_s_scores <- vegan::scores(rda222, display = "bp")
func2_s_scores_df <- as.data.frame(func2_s_scores)
rm(func2_s_scores)

site_scores <- vegan::scores(rda222, display = "sites") 
site_scores_df <- as.data.frame(site_scores)
rm(site_scores)
site_scores_df$SiteRowNames <- rownames(site_scores_df)
site_scores_df$Site <- rownames(site_scores_df)
site_scores_df <- site_scores_df %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores_df <- dplyr::left_join(site_scores_df, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores_df) <- site_scores_df$SiteRowNames
site_scores_df$Aridity <- (1 - site_scores_df$AI)
site_scores_df <- site_scores_df[,-c(3,5:6)]

species_scores <- vegan::scores(rda222, display = "species")
species_scores_df <- as.data.frame(species_scores)
species_scores_df$var <- rownames(species_scores_df)
rm(species_scores)

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.12,  
            nudge_x = 0.01) + 
  geom_segment(data = func2_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text(data = func2_s_scores_df, 
            aes(x = RDA1, y = RDA2, label = rownames(func2_s_scores_df)), 
            size = 4, color = "blue",
            fontface = "bold") +
  geom_text(data = species_scores_df,  # Use rda_mcc_df for text labels
            aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
            size = 4, color = "darkred", vjust = -1,
            fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()


ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
       "RDA_FUNC_ITS2_NOBIOMASS.png", width = 10, height = 8, dpi = 300)




# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(func2_s, ITS2_h,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, # This means that the forward selection will only add variables that explain at least 1% more variance to the model.
            alpha = 0.05) # The significance level used for testing each variable's inclusion into the model.
ordistep(rda222,direction="forward",perm.max = 200)


func2_s_scores_df_sig <- func2_s_scores_df[c(1:2,7,11,17,20),]
ITS2_hs <- ITS2_h[,c(1:2,7,11,17,20)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.12,  
            nudge_x = 0.01) + 
  geom_segment(data = func2_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func2_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func2_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()
 
# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_FUNC_ITS2_sel_NOBIOMASS.png", width = 10, height = 8, dpi = 300)

rda2221 <- rda(func2_s ~ ., data=data.frame(ITS2_hs))
summary(rda2221)

R2 <- RsquareAdj(rda2221)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda2221)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda2221, permutations = how(nperm = 10000))






# !!!!!! LESS RESTRICTED FORWARD SELECTION !!!!!!

# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(func2_s, ITS2_h,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, 
            alpha = 0.1) #INCREASED ALPHA


func2_s_scores_df_sig <- func2_s_scores_df[c(13, 7, 19, 21),]
ITS2_hs <- ITS2_h[,c(13, 7, 19, 21)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.12,  
            nudge_x = 0.01) + 
  geom_segment(data = func2_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func2_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func2_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()

# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_FUNC_ITS2_sel2_NOBIOMASS.png", width = 10, height = 8, dpi = 300)

rda22221 <- rda(func2_s ~ ., data=data.frame(ITS2_hs))
summary(rda22221)

R2 <- RsquareAdj(rda22221)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda22221)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda22221, permutations = how(nperm = 10000))








rm(func2_s_scores_df, func2_s_scores_df_sig, site_scores_df, species_scores_df,
   species_scores_df_sig, ITS2_hs)
rm(rda2221, rda222, rda22221)










# 16S vs Env ----
rda3 <- rda(prok_h ~ ., data=data.frame(env_s))
summary(rda3)

R2 <- RsquareAdj(rda3)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda3)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda3, permutations = how(nperm = 10000))
# significative
anova.cca(rda3, by = "axis", permutations = how(nperm = 1000))
# First RDA significative



plot(rda3, scaling=2, main='Triplot - scaling 2')



summary(rda3)
# BIPLOT = FUNCTIONS (arrows)
# SITES = sample sites (points)
# SPECIES = MCC (text)


func_s_scores <- vegan::scores(rda3, display = "bp")
func_s_scores_df <- as.data.frame(func_s_scores)
rm(func_s_scores)

site_scores <- vegan::scores(rda3, display = "sites") 
site_scores_df <- as.data.frame(site_scores)
rm(site_scores)
site_scores_df$SiteRowNames <- rownames(site_scores_df)
site_scores_df$Site <- rownames(site_scores_df)
site_scores_df <- site_scores_df %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores_df <- dplyr::left_join(site_scores_df, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores_df) <- site_scores_df$SiteRowNames
site_scores_df$Aridity <- (1 - site_scores_df$AI)
site_scores_df <- site_scores_df[,-c(3,5:6)]

species_scores <- vegan::scores(rda3, display = "species")
species_scores_df <- as.data.frame(species_scores)
species_scores_df$var <- rownames(species_scores_df)
rm(species_scores)

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()


# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_16S_ENV.png", width = 10, height = 8, dpi = 300)




# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(prok_h, env_s,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, # This means that the forward selection will only add variables that explain at least 1% more variance to the model.
            alpha = 0.05) # The significance level used for testing each variable's inclusion into the model.

func_s_scores_df_sig <- func_s_scores_df[c(1,4:6,8,10,12:13,16,19,25,30:31),]
env_s1 <- env_s[,c(1,4:6,8,10,12:13,16,19,25,30:31)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()

# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_16S_ENV_sel.png", width = 10, height = 8, dpi = 300)

rda31 <- rda(prok_h ~ ., data=data.frame(env_s1))
summary(rda31)

R2 <- RsquareAdj(rda31)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda31)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda31, permutations = how(nperm = 10000))






# !!!!!! LESS RESTRICTED FORWARD SELECTION !!!!!!

# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(prok_h, env_s,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, 
            alpha = 0.1) #INCREASED ALPHA





rm(func_s_scores_df, func_s_scores_df_sig, site_scores_df, species_scores_df,
   species_scores_df_sig)




# ITS2 vs Env ----
rda4 <- rda(ITS2_h ~ ., data=data.frame(env_s))
summary(rda4)

R2 <- RsquareAdj(rda4)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda4)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda4, permutations = how(nperm = 10000))
# significative
anova.cca(rda4, by = "axis", permutations = how(nperm = 1000))
# 4 RDA significative



plot(rda4, scaling=2, main='Triplot - scaling 2')



summary(rda4)
# BIPLOT = FUNCTIONS (arrows)
# SITES = sample sites (points)
# SPECIES = MCC (text)


func_s_scores <- vegan::scores(rda4, display = "bp")
func_s_scores_df <- as.data.frame(func_s_scores)
rm(func_s_scores)

site_scores <- vegan::scores(rda4, display = "sites") 
site_scores_df <- as.data.frame(site_scores)
rm(site_scores)
site_scores_df$SiteRowNames <- rownames(site_scores_df)
site_scores_df$Site <- rownames(site_scores_df)
site_scores_df <- site_scores_df %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores_df <- dplyr::left_join(site_scores_df, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores_df) <- site_scores_df$SiteRowNames
site_scores_df$Aridity <- (1 - site_scores_df$AI)
site_scores_df <- site_scores_df[,-c(3,5:6)]

species_scores <- vegan::scores(rda4, display = "species")
species_scores_df <- as.data.frame(species_scores)
species_scores_df$var <- rownames(species_scores_df)
rm(species_scores)

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()


# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_ITS2_ENV.png", width = 10, height = 8, dpi = 300)




# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(ITS2_h, env_s,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, # This means that the forward selection will only add variables that explain at least 1% more variance to the model.
            alpha = 0.05) # The significance level used for testing each variable's inclusion into the model.

func_s_scores_df_sig <- func_s_scores_df[c(1:6,10,14,16:18,26,31,36:37),]
env_s1 <- env_s[,c(1:6,10,14,16:18,26,31,36:37)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()

# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_ITS2_ENV_sel.png", width = 10, height = 8, dpi = 300)

rda41 <- rda(ITS2_h ~ ., data=data.frame(env_s1))
summary(rda41)

R2 <- RsquareAdj(rda41)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda41)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda41, permutations = how(nperm = 10000))




# !!!!!! LESS RESTRICTED FORWARD SELECTION !!!!!!

# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(ITS2_h, env_s,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, 
            alpha = 0.1) #INCREASED ALPHA


func_s_scores_df_sig <- func_s_scores_df[c(1:2, 4:6, 14:17, 20:22, 36:37),]
env_s1 <- env_s[,c(1:2, 4:6, 14:17, 20:22, 36:37)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.12,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()

ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
       "RDA_ENV_ITS2_sel2.png", width = 10, height = 8, dpi = 300)

rda42 <- rda(func2_s ~ ., data=data.frame(env_s1))
summary(rda42)

R2 <- RsquareAdj(rda42)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda42)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda42, permutations = how(nperm = 10000))




rm(func_s_scores_df, func_s_scores_df_sig, site_scores_df, species_scores_df,
   species_scores_df_sig, env_s1)

rm(rda4, rda41)






# Func vs Env: Without biomass ----
rda6 <- rda(func2_s ~ ., data=data.frame(env_s))
summary(rda6)

R2 <- RsquareAdj(rda6)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda6)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda6, permutations = how(nperm = 10000))
# significative
anova.cca(rda6, by = "axis", permutations = how(nperm = 1000))
# 3 first RDA significatives



plot(rda6, scaling=2, main='Triplot - scaling 2')



summary(rda6)
# BIPLOT = FUNCTIONS (arrows)
# SITES = sample sites (points)
# SPECIES = MCC (text)


func2_s_scores <- vegan::scores(rda6, display = "bp")
func2_s_scores_df <- as.data.frame(func2_s_scores)
rm(func2_s_scores)

site_scores <- vegan::scores(rda6, display = "sites") 
site_scores_df <- as.data.frame(site_scores)
rm(site_scores)
site_scores_df$SiteRowNames <- rownames(site_scores_df)
site_scores_df$Site <- rownames(site_scores_df)
site_scores_df <- site_scores_df %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores_df <- dplyr::left_join(site_scores_df, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores_df) <- site_scores_df$SiteRowNames
site_scores_df$Aridity <- (1 - site_scores_df$AI)
site_scores_df <- site_scores_df[,-c(3,5:6)]

species_scores <- vegan::scores(rda6, display = "species")
species_scores_df <- as.data.frame(species_scores)
species_scores_df$var <- rownames(species_scores_df)
rm(species_scores)

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.12,  
            nudge_x = 0.01) + 
  geom_segment(data = func2_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func2_s_scores_df, 
                  aes(x = RDA1, y = RDA2, label = rownames(func2_s_scores_df)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()


# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_FUNC_ENV_NOBIOMASS.png", width = 10, height = 8, dpi = 300)



# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(func2_s, env_s,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, # This means that the forward selection will only add variables that explain at least 1% more variance to the model.
            alpha = 0.05) # The significance level used for testing each variable's inclusion into the model.
ordistep(rda6,direction="forward",perm.max = 200)


func2_s_scores_df_sig <- func2_s_scores_df[c(1,5:6,10:12,19,32),]
env_hs <- env_s[,c(1,5:6,10:12,19,32)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.12,  
            nudge_x = 0.01) + 
  geom_segment(data = func2_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func2_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func2_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()
 
# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_FUNC_ENV_sel_NOBIOMASS.png", width = 10, height = 8, dpi = 300)

rda61 <- rda(func2_s ~ ., data=data.frame(env_hs))
summary(rda61)

R2 <- RsquareAdj(rda61)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda61)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda61, permutations = how(nperm = 10000))









# !!!!!! LESS RESTRICTED FORWARD SELECTION !!!!!!

# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(func2_s, env_s,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, 
            alpha = 0.1) #INCREASED ALPHA


func2_s_scores_df_sig <- func2_s_scores_df[c(1, 3, 5:8, 10:12, 18:22, 26, 32, 37),]
env_s1 <- env_s[,c(1, 3, 5:8, 10:12, 18:22, 26, 32, 37)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = -0.12,  
            nudge_x = 0.01) + 
  geom_segment(data = func2_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func2_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func2_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()

# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_ENV_FUNC_sel2_NOBIOMASS.png", width = 10, height = 8, dpi = 300)

rda42 <- rda(func2_s ~ ., data=data.frame(env_s1))
summary(rda42)

R2 <- RsquareAdj(rda42)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda42)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda42, permutations = how(nperm = 10000))




rm(func2_s_scores_df, func2_s_scores_df_sig, site_scores_df, species_scores_df,
   species_scores_df_sig, env_s1)
rm(rda61, rda6, rda42)


# (Func vs 16S) ----

rda12 <- rda(func_s ~ ., data=data.frame(prok_h))
summary(rda12)

R2 <- RsquareAdj(rda12)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda12)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda12, permutations = how(nperm = 10000))
# significative
anova.cca(rda12, by = "axis", permutations = how(nperm = 1000))
# 2 first RDA significatives



plot(rda12, scaling=2, main='Triplot - scaling 2')



summary(rda12)


func_s_scores <- vegan::scores(rda12, display = "bp")
func_s_scores_df <- as.data.frame(func_s_scores)
rm(func_s_scores)

site_scores <- vegan::scores(rda12, display = "sites") 
site_scores_df <- as.data.frame(site_scores)
rm(site_scores)
site_scores_df$SiteRowNames <- rownames(site_scores_df)
site_scores_df$Site <- rownames(site_scores_df)
site_scores_df <- site_scores_df %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores_df <- dplyr::left_join(site_scores_df, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores_df) <- site_scores_df$SiteRowNames
site_scores_df$Aridity <- (1 - site_scores_df$AI)
site_scores_df <- site_scores_df[,-c(3,5:6)]

species_scores <- vegan::scores(rda12, display = "species")
species_scores_df <- as.data.frame(species_scores)
species_scores_df$var <- rownames(species_scores_df)
rm(species_scores)

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()


# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_FUNC_16S.png", width = 10, height = 8, dpi = 300)







# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(func_s, prok_h,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, # This means that the forward selection will only add variables that explain at least 1% more variance to the model.
            alpha = 0.05) # The significance level used for testing each variable's inclusion into the model.
ordistep(rda12,direction="forward",perm.max = 200)


func_s_scores_df_sig <- func_s_scores_df[c(1,2,12),]
prok_h1 <- prok_h[,c(1,2,12)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()

# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_FUNC_16S_sel.png", width = 10, height = 8, dpi = 300)

rda121 <- rda(func_s ~ ., data=data.frame(prok_h1))
summary(rda121)

R2 <- RsquareAdj(rda121)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda121)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda121, permutations = how(nperm = 10000))

rm(func_s_scores_df, func_s_scores_df_sig, site_scores_df, species_scores_df,
   species_scores_df_sig)




# (Func vs ITS2) ----

rda22 <- rda(func_s ~ ., data=data.frame(ITS2_h))
summary(rda22)

R2 <- RsquareAdj(rda22)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda22)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda22, permutations = how(nperm = 10000))
# significative
anova.cca(rda22, by = "axis", permutations = how(nperm = 1000))
# First RDA significatives


plot(rda22, scaling=2, main='Triplot - scaling 2')

func_s_scores <- vegan::scores(rda22, display = "bp")
func_s_scores_df <- as.data.frame(func_s_scores)
rm(func_s_scores)

site_scores <- vegan::scores(rda22, display = "sites") 
site_scores_df <- as.data.frame(site_scores)
rm(site_scores)
site_scores_df$SiteRowNames <- rownames(site_scores_df)
site_scores_df$Site <- rownames(site_scores_df)
site_scores_df <- site_scores_df %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores_df <- dplyr::left_join(site_scores_df, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores_df) <- site_scores_df$SiteRowNames
site_scores_df$Aridity <- (1 - site_scores_df$AI)
site_scores_df <- site_scores_df[,-c(3,5:6)]

species_scores <- vegan::scores(rda22, display = "species")
species_scores_df <- as.data.frame(species_scores)
species_scores_df$var <- rownames(species_scores_df)
rm(species_scores)

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text(data = func_s_scores_df, 
            aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df)), 
            size = 4, color = "blue",
            fontface = "bold") +
  geom_text(data = species_scores_df,  # Use rda_mcc_df for text labels
            aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
            size = 4, color = "darkred", vjust = -1,
            fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()


# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_FUNC_ITS2.png", width = 10, height = 8, dpi = 300)




# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(func_s, ITS2_h,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, # This means that the forward selection will only add variables that explain at least 1% more variance to the model.
            alpha = 0.05) # The significance level used for testing each variable's inclusion into the model.
ordistep(rda22,direction="forward",perm.max = 200)


func_s_scores_df_sig <- func_s_scores_df[c(19,21),]
ITS2_hs <- ITS2_h[,c(19,21)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.2) + 
  geom_segment(data = func_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()

ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
       "RDA_FUNC_ITS2_sel.png", width = 10, height = 8, dpi = 300)

rda221 <- rda(func_s ~ ., data=data.frame(ITS2_hs))
summary(rda221)

R2 <- RsquareAdj(rda221)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda221)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda221, permutations = how(nperm = 10000))

rm(func_s_scores_df, func_s_scores_df_sig, site_scores_df, species_scores_df,
   species_scores_df_sig)






# (16S vs Func) ----

rda1 <- rda(prok_h ~ ., data=data.frame(func_s))
summary(rda1)

R2 <- RsquareAdj(rda1)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda1)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda1, permutations = how(nperm = 10000))
# significative
anova.cca(rda1, by = "axis", permutations = how(nperm = 1000))
# 2 first RDA significatives



plot(rda1, scaling=2, main='Triplot - scaling 2')



summary(rda1)
# BIPLOT = FUNCTIONS (arrows)
# SITES = sample sites (points)
# SPECIES = MCC (text)


func_s_scores <- vegan::scores(rda1, display = "bp")
func_s_scores_df <- as.data.frame(func_s_scores)
rm(func_s_scores)

site_scores <- vegan::scores(rda1, display = "sites") 
site_scores_df <- as.data.frame(site_scores)
rm(site_scores)
site_scores_df$SiteRowNames <- rownames(site_scores_df)
site_scores_df$Site <- rownames(site_scores_df)
site_scores_df <- site_scores_df %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores_df <- dplyr::left_join(site_scores_df, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores_df) <- site_scores_df$SiteRowNames
site_scores_df$Aridity <- (1 - site_scores_df$AI)
site_scores_df <- site_scores_df[,-c(3,5:6)]

species_scores <- vegan::scores(rda1, display = "species")
species_scores_df <- as.data.frame(species_scores)
species_scores_df$var <- rownames(species_scores_df)
rm(species_scores)

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()


# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_16S_FUNC.png", width = 10, height = 8, dpi = 300)




# Selection of the significant variables from the response matrix:

# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(prok_h, func_s,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, # This means that the forward selection will only add variables that explain at least 1% more variance to the model.
            alpha = 0.05) # The significance level used for testing each variable's inclusion into the model.
ordistep(rda1,direction="forward",perm.max = 200)


func_s_scores_df_sig <- func_s_scores_df[c(1,6,10),]
func_s1 <- func_s[,c(1,6,10)]


ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df_sig, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df_sig, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df_sig)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()

ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
       "RDA_16S_FUNC_sel.png", width = 10, height = 8, dpi = 300)

rda11 <- rda(prok_h ~ ., data=data.frame(func_s1))
summary(rda11)

R2 <- RsquareAdj(rda11)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda11)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda11, permutations = how(nperm = 10000))

rm(func_s_scores_df, func_s_scores_df_sig, site_scores_df, species_scores_df,
   species_scores_df_sig)









# (ITS2 vs Func) ----

rda2 <- rda(ITS2_h ~ ., data=data.frame(func_s))
summary(rda2)

R2 <- RsquareAdj(rda2)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda2)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda2, permutations = how(nperm = 10000))
# significative
anova.cca(rda2, by = "axis", permutations = how(nperm = 1000))
# 1 RDA significative



plot(rda2, scaling=2, main='Triplot - scaling 2')



summary(rda2)
# BIPLOT = FUNCTIONS (arrows)
# SITES = sample sites (points)
# SPECIES = MCC (text)


func_s_scores <- vegan::scores(rda2, display = "bp")
func_s_scores_df <- as.data.frame(func_s_scores)
rm(func_s_scores)

site_scores <- vegan::scores(rda2, display = "sites") 
site_scores_df <- as.data.frame(site_scores)
rm(site_scores)
site_scores_df$SiteRowNames <- rownames(site_scores_df)
site_scores_df$Site <- rownames(site_scores_df)
site_scores_df <- site_scores_df %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores_df <- dplyr::left_join(site_scores_df, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores_df) <- site_scores_df$SiteRowNames
site_scores_df$Aridity <- (1 - site_scores_df$AI)
site_scores_df <- site_scores_df[,-c(3,5:6)]

species_scores <- vegan::scores(rda2, display = "species")
species_scores_df <- as.data.frame(species_scores)
species_scores_df$var <- rownames(species_scores_df)
rm(species_scores)

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text(data = func_s_scores_df, 
            aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df)), 
            size = 4, color = "blue",
            fontface = "bold") +
  geom_text(data = species_scores_df,  # Use rda_mcc_df for text labels
            aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
            size = 4, color = "darkred", vjust = -1,
            fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()


# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_ITS2_FUNC.png", width = 10, height = 8, dpi = 300)




# Selection of the significant variables from the response matrix:

install.packages("adespatial")
library(adespatial)
# order matrixes data: RESPONSE, EXPLANATORI MATRIXES
forward.sel(func_s, ITS2_h,
            R2thresh = 0.99,
            adjR2thresh = 0.99,
            nperm = 999,
            R2more = 0.01, # This means that the forward selection will only add variables that explain at least 1% more variance to the model.
            alpha = 0.05) # The significance level used for testing each variable's inclusion into the model.
ordistep(rda2, direction="forward",perm.max = 200)


species_scores_df_sig <- species_scores_df[c(19,21),]

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df_sig,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()




rm(func_s_scores_df, func_s_scores_df_sig, site_scores_df, species_scores_df,
   species_scores_df_sig)







# (Env vs. Func) ----
rda5 <- rda(env_s ~ ., data=data.frame(func_s))
summary(rda5)

R2 <- RsquareAdj(rda5)$r.squared # Aquest és el valor que ja havíem vist abans al summary
R2
R2adj <- RsquareAdj(rda5)$adj.r.squared # Aquests él el valor ajustat que ens interessa
R2adj


anova.cca(rda5, permutations = how(nperm = 10000))
# significative
anova.cca(rda5, by = "axis", permutations = how(nperm = 1000))
# 3 first RDA significatives



plot(rda5, scaling=2, main='Triplot - scaling 2')



summary(rda5)
# BIPLOT = FUNCTIONS (arrows)
# SITES = sample sites (points)
# SPECIES = MCC (text)


func_s_scores <- vegan::scores(rda5, display = "bp")
func_s_scores_df <- as.data.frame(func_s_scores)
rm(func_s_scores)

site_scores <- vegan::scores(rda5, display = "sites") 
site_scores_df <- as.data.frame(site_scores)
rm(site_scores)
site_scores_df$SiteRowNames <- rownames(site_scores_df)
site_scores_df$Site <- rownames(site_scores_df)
site_scores_df <- site_scores_df %>%
  tidyr::separate(Site, into = c("Site", "Replicate"), sep = "-")
AI_df <- data.frame(
  Site = c("FDE", "ATK", "LHE", "ARZ", "VAL", "GAV", "ALB", "MON", "COY", "MAL", "SAN", "TAB"),
  AI = c(1.33, 1.25, 1.09, 1.00, 0.70, 0.56, 0.43, 0.40, 0.25, 0.18, 0.18, 0.16)
)
site_scores_df <- dplyr::left_join(site_scores_df, AI_df, by = "Site")
rm(AI_df)
rownames(site_scores_df) <- site_scores_df$SiteRowNames
site_scores_df$Aridity <- (1 - site_scores_df$AI)
site_scores_df <- site_scores_df[,-c(3,5:6)]

species_scores <- vegan::scores(rda5, display = "species")
species_scores_df <- as.data.frame(species_scores)
species_scores_df$var <- rownames(species_scores_df)
rm(species_scores)

ggplot() +
  geom_point(data = site_scores_df, 
             aes(x = RDA1, y = RDA2, color = Aridity),
             size = 3, show.legend = TRUE) + 
  geom_text(data = site_scores_df, 
            aes(x = RDA1, y = RDA2, label = Site), 
            size = 3, color = "#404040", 
            nudge_y = 0.03,  
            nudge_x = 0.01) + 
  geom_segment(data = func_s_scores_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "inches")), 
               color = "blue") +
  geom_text_repel(data = func_s_scores_df, 
                  aes(x = RDA1, y = RDA2, label = rownames(func_s_scores_df)), 
                  size = 4, color = "blue",
                  fontface = "bold") +
  geom_text_repel(data = species_scores_df,  # Use rda_mcc_df for text labels
                  aes(x = RDA1, y = RDA2, label = var),  # Ensure 'var' is the correct column for labels
                  size = 4, color = "darkred", vjust = -1,
                  fontface = "bold") +
  scale_color_AI(discrete = FALSE, palette = "Sites", reverse = FALSE, name = "Aridity")+
  xlab("RDA1") + ylab("RDA2") +
  theme_minimal()


# ggsave(path = "C:/Users/ecologia.PCECO002/OneDrive - Universitat de Girona/GRADCATCH/Manuscripts/1 GRADIENT/Figures",
#        "RDA_ENV_FUNC.png", width = 10, height = 8, dpi = 300)




