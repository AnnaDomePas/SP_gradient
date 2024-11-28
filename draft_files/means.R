
library(writexl)
library(dplyr)
# means_by_site <- my_data %>%
#   group_by(Site) %>%
#   summarise(across(where(is.numeric), mean, na.rm = TRUE))

means_by_site <- my_data %>%
  group_by(Site) %>%
  summarise_if(is.numeric, 
               list(mean = ~mean(., na.rm = TRUE),
                    sd = ~sd(., na.rm = TRUE),
                    sem = ~sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))))
means_by_site_clean <- means_by_site %>%
  select(-contains("_sd"))


final <- means_by_site_clean[,-c(2:4,6,9,14:16,21:26, 30:75)]

final <- final[,-c(30:32,34,37,42:44,49:54,58:102,117)]



final <- final %>%
  arrange(aridity_mean)

final <- final %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))
print(final)

all_columns <- colnames(final)

# Identify mean and sem columns
mean_columns <- all_columns[grepl("_mean$", all_columns)]
sem_columns <- all_columns[grepl("_sem$", all_columns)]

# Create a vector of ordered column names
ordered_columns <- c("Site", unlist(mapply(function(mean, sem) c(mean, sem), mean_columns, sem_columns)))

# Reorder the dataset
final_reordered <- final %>%
  select(all_of(ordered_columns))

final <- final_reordered %>%
  mutate(site = case_when(
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


write_xlsx(final, "means_sem.xlsx")
