#CREATED BY: ANNA DOMÉNECH-PASCUAL (April 2024)

#Load packages ----
library(dplyr)
library(lme4)

# Load and edit data ----
my_data         = read.csv("SP_metadata_2021.csv", sep=";")

my_data$LTC_LTN <- my_data$L_TC/my_data$L_TN
my_data$xylcbh <- my_data$xyl + my_data$cbh
my_data$alphabeta <- my_data$alpha + my_data$beta

#NEW VARIABLE: ARIDITY
my_data$Aridity <- (1 - my_data$AI)

#To replace NA values with a mean of the other values of the Site:
for (i in which(sapply(my_data, is.numeric))) {
  for (j in which(is.na(my_data[, i]))) {
    my_data[j, i] <- mean(my_data[my_data[, "Site"] == my_data[j, "Site"], i],  na.rm = TRUE)
  }
}

max_values      = apply(my_data[,c(8:34,40,41,76,81:93)],2,max)
my_data.1       = as.data.frame(cbind(my_data[,c(2,96,46:53,64,94,95)],
                                      my_data[,c(8:34,40,41,76,81:93)] / as.list(max_values)))
a.1             = as.data.frame(colnames(my_data.1))

corr            = as.data.frame(cor(my_data.1[,14:56]))

my_data.1 <-my_data.1 %>%
  mutate(Cenz = select(., 3:6) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(MB = select(.,41:42) %>% rowSums(na.rm = TRUE))


#RESPIRATION ----
# > Clay:Aridity ----
# interaction.plot(
#   x.factor = my_data$Clay,
#   trace.factor = my_data$Respiration,
#   response = my_data$Aridity,
#   fun = mean,
#   ylab = "Aridity",
#   xlab = "Clay",
#   trace.label = "Respiration"
# )

my_data.resp = my_data.1[-c(9,21,35,39,48),]

respiration.full.2r = lmer(Respiration ~ 1 + Water_content + Clay + BIX + Aridity + BIX:Aridity + 
                             Clay:Aridity + (1 | Site),
                           data = my_data.resp)  




# Create a grid of values for Clay, Aridity, Water_content, and BIX
clay_values <- seq(min(my_data.resp$Clay), max(my_data.resp$Clay), length.out = 50)
aridity_values <- seq(min(my_data.resp$Aridity), max(my_data.resp$Aridity), length.out = 50)
water_content_values <- seq(min(my_data.resp$Water_content), max(my_data.resp$Water_content), length.out = 50)
bix_values <- seq(min(my_data.resp$BIX), max(my_data.resp$BIX), length.out = 50)

grid <- expand.grid(Clay = clay_values, Aridity = aridity_values, Water_content = water_content_values, BIX = bix_values)

# Predict Respiration values for the grid
pred <- predict(respiration.full.2r, newdata = grid, re.form = NA)

# Combine predicted values with the grid
grid$respiration <- pred

# Plot the interaction
interaction_plot <- plot(grid$Clay, grid$Aridity, type = "n", xlab = "Clay", ylab = "Aridity", main = "Interaction Plot: Respiration vs. Clay and Aridity")
points(grid$Clay, grid$Aridity, col = "blue", pch = 20, cex = 0.5)

# Reshape the data for contour plotting
contour_data <- with(grid, contour(x = sort(unique(Clay)), y = sort(unique(Aridity)), z = matrix(respiration, nrow = length(unique(Aridity)), byrow = TRUE)))

# Plot the contour
contour(x = contour_data$x, y = contour_data$y, z = contour_data$z, add = TRUE, labels = contour_data$level)# Add legend
legend("topright", legend = "Respiration", col = "blue", pch = 20)

# Display the plot
interaction_plot
