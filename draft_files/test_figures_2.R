library(buildmer)
library(lmerTest)

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
my_data.1       = as.data.frame(cbind(my_data[,c(2,7,46:53,64,94,95)],
                                      my_data[,c(8:34,40,41,76,81:93)] / as.list(max_values)))
a.1             = as.data.frame(colnames(my_data.1))

corr            = as.data.frame(cor(my_data.1[,14:56]))

my_data.1 <-my_data.1 %>%
  mutate(Cenz = select(., 3:6) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(MB = select(.,41:42) %>% rowSums(na.rm = TRUE))

# VIF ####

# A VIF greater than 10 is a signal that the model has a collinearity problem
# calculates VIF for all variables, excludes the one with the highest VIF 
# (if it is greater than the threshold), repeat the procedure untill no variables
# with a VIF greater than th remains

# install.packages('usdm')
library(usdm)
my_data_vif = as.data.frame(cbind(my_data.1[,c(2,14,16,18,19,23:29,33:36,39:56)]))
vif(my_data_vif)
vifcor(my_data_vif, th=0.9)
vifstep(my_data_vif, th=10) # I am using the results from this function 
                            # (https://www.rdocumentation.org/packages/usdm/versions/2.1-7/topics/vif)
                            # but I am adding some variables that we know are important
                            # like water content and aridity index and clay, FB

# Respiration ####

# Model 1 ####

respiration.full.I = lmer(Respiration ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                  TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                                  Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                                  BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI+
                            (1|Site), data = my_data.1)                                         

summary(respiration.full.I)
Anova(respiration.full.I)
cAIC(respiration.full.I)
AIC(respiration.full.I)
r.squaredGLMM(respiration.full.I)

qqnorm(residuals(respiration.full.I))
scatter.smooth(residuals(respiration.full.I) ~ fitted(respiration.full.I))

# Model 1 - reduction ####

max_mod = buildmer(Respiration ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                             TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                             Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                             BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(max_mod)
print(max_mod, correlation=TRUE)
vcov(max_mod)

# Model 2 - ####

respiration.full.2 = lmer(Respiration ~ 1 + Water_content + Clay + BIX + AI + BIX:AI + 
                            Clay:AI + (1 | Site),
                          data = my_data.1)                                         

summary(respiration.full.2)
Anova(respiration.full.2)
cAIC(respiration.full.2)
AIC(respiration.full.2)
r.squaredGLMM(respiration.full.2)

qqnorm(residuals(respiration.full.2))
scatter.smooth(residuals(respiration.full.2) ~ fitted(respiration.full.2))
# It seems that the assumptions are failing, so I need to reduce the statistical 
# outliers

# Model 2 - without outliers - WINNER ####

cooksd <- cooks.distance(respiration.full.2)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
my_data.resp = my_data.1[-c(9,21,35,39,48),]

respiration.full.2r = lmer(Respiration ~ 1 + Water_content + Clay + BIX + AI + BIX:AI + 
                            Clay:AI + (1 | Site),
                          data = my_data.resp)                                         

summary(respiration.full.2r)
Anova(respiration.full.2r)
cAIC(respiration.full.2r)
AIC(respiration.full.2r)
r.squaredGLMM(respiration.full.2r)

qqnorm(residuals(respiration.full.2r))
scatter.smooth(residuals(respiration.full.2r) ~ fitted(respiration.full.2r))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(Respiration ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.resp, 
                          sets = list("Water_content","Clay","BIX","AI",
                                      "BIX:AI","Clay:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Biomass - BB ####

# Model 1 ####

BB.full.I = lmer(log(BB) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                  TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                  Litter+L_TN+SR+E2.E3+E4.E6+
                                  BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + 
                   (1|Site), data = my_data.1)                                         

summary(BB.full.I)
Anova(BB.full.I)
cAIC(BB.full.I)
AIC(BB.full.I)
r.squaredGLMM(BB.full.I)

qqnorm(residuals(BB.full.I))
scatter.smooth(residuals(BB.full.I) ~ fitted(BB.full.I))

# Model 1 - reduction ####

max_BB = buildmer(log(BB) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                   TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                   Litter+L_TN+SR+E2.E3+E4.E6+
                                   BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + 
                    (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(max_BB)
print(max_BB, correlation=TRUE)
vcov(max_BB)

# Model 2 - WINNER MODEL ####

BB.full.2 = lmer(log(BB) ~ 1 + Water_content + Peak_B + L_TN + PO43 + E2.E3 +  
                   Peak_C + AI + E2.E3:AI + Peak_B:AI + Soil_Temp + Peak_C:AI + 
                   (1|Site), data = my_data.1)                                         

summary(BB.full.2)
Anova(BB.full.2)
cAIC(BB.full.2)
AIC(BB.full.2)
r.squaredGLMM(BB.full.2)

qqnorm(residuals(BB.full.2))
scatter.smooth(residuals(BB.full.2) ~ fitted(BB.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(log(BB) ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("Water_content","Peak_B","L_TN","PO43",
                                      "E2.E3","Peak_C","AI","E2.E3:AI","Peak_B:AI",
                                      "Soil_Temp","Peak_C:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Biomass - FB ####

# Model 1 ####

FB.full.I = lmer(log(FB) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                  TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                  Litter+L_TN+SR+E2.E3+E4.E6+
                                  BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + 
                   (1|Site), data = my_data.1)                                          

summary(FB.full.I)
Anova(FB.full.I)
cAIC(FB.full.I)
AIC(FB.full.I)
r.squaredGLMM(FB.full.I)

qqnorm(residuals(FB.full.I))
scatter.smooth(residuals(FB.full.I) ~ fitted(FB.full.I))

# Model 1 - reduction ####

max_FB = buildmer(log(FB) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                   TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                   Litter+L_TN+SR+E2.E3+E4.E6+
                                   BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + 
                    (1|Site), data = my_data.1,
                  buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                    ddf = "Satterthwaite"))
summary(max_FB)
print(max_FB, correlation=TRUE)
vcov(max_FB)

# Model 2 -  WINNER ####

FB.full.2 = lmer(log(FB) ~ 1 + TC + C_N + Peak_C + Litter + (1 | Site), data = my_data.1)                                         

summary(FB.full.2)
Anova(FB.full.2)
cAIC(FB.full.2)
AIC(FB.full.2)
r.squaredGLMM(FB.full.2)

qqnorm(residuals(FB.full.2))
scatter.smooth(residuals(FB.full.2) ~ fitted(FB.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(log(FB) ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("TC","C_N","Peak_C","Litter"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Biomass - MB ####

# Model 1 ####

MB.full.I = lmer(log(MB) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                  TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                  Litter+L_TN+SR+E2.E3+E4.E6+
                                  BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + 
                   (1|Site), data = my_data.1)                                          

summary(MB.full.I)
Anova(MB.full.I)
cAIC(MB.full.I)
AIC(MB.full.I)
r.squaredGLMM(MB.full.I)

qqnorm(residuals(MB.full.I))
scatter.smooth(residuals(MB.full.I) ~ fitted(MB.full.I))

# Model 1 - reduction ####

max_MB = buildmer(log(MB) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                   TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                   Litter+L_TN+SR+E2.E3+E4.E6+
                                   BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + 
                    (1|Site), data = my_data.1,
                  buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                    ddf = "Satterthwaite"))
summary(max_MB)
print(max_MB, correlation=TRUE)
vcov(max_MB)

# Model 2 - WINNER MODEL ####

MB.full.2 = lmer(log(MB) ~ 1 + TC + Peak_B + Water_content + C_N + Silt + Clay +  
                   PO43 + Peak_C + (1|Site), data = my_data.1)                                          

summary(MB.full.2)
Anova(MB.full.2)
cAIC(MB.full.2)
AIC(MB.full.2)
r.squaredGLMM(MB.full.2)

qqnorm(residuals(MB.full.2))
scatter.smooth(residuals(MB.full.2) ~ fitted(MB.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(log(MB) ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("TC","Peak_B","Water_content","C_N","Silt",
                                      "Clay","PO43","Peak_C"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - alpha ####

# Model 1 ####

alpha.full.I = lmer(alpha ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                        TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                        Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                        BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), 
                    data = my_data.1)                                         

summary(alpha.full.I)
Anova(alpha.full.I)
cAIC(alpha.full.I)
AIC(alpha.full.I)
r.squaredGLMM(alpha.full.I)

qqnorm(residuals(alpha.full.I))
scatter.smooth(residuals(alpha.full.I) ~ fitted(alpha.full.I))

# Model 1 - reduction ####

alpha_mod = buildmer(alpha ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                         TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                         Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                         BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), 
                     data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(alpha_mod)
print(alpha_mod, correlation=TRUE)
vcov(alpha_mod)

# Model 2 ####

alpha.full.2 = lmer(alpha ~ 1 + Soil_Temp + (1|Site), data = my_data.1)                                         

summary(alpha.full.2)
Anova(alpha.full.2)
cAIC(alpha.full.2)
AIC(alpha.full.2)
r.squaredGLMM(alpha.full.2)

qqnorm(residuals(alpha.full.2))
scatter.smooth(residuals(alpha.full.2) ~ fitted(alpha.full.2))

# Model 1 - without outliers ####

cooksd <- cooks.distance(alpha.full.2)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
my_data.alpha = my_data.1[-c(49,59),]

alpha.full.2r = lmer(alpha ~ 1 + Soil_Temp + (1|Site), data = my_data.alpha)                                         

summary(alpha.full.2r)
Anova(alpha.full.2r)
cAIC(alpha.full.2r)
AIC(alpha.full.2r)
r.squaredGLMM(alpha.full.2r)

qqnorm(residuals(alpha.full.2r))
scatter.smooth(residuals(alpha.full.2r) ~ fitted(alpha.full.2r))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(alpha ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.2, 
                          sets = list("Soil_Temp"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - beta ####

# Model 1 ####

beta.full.I = lmer(beta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                      TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                      Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                      BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site),
                   data = my_data.1)                                         

summary(beta.full.I)
Anova(beta.full.I)
cAIC(beta.full.I)
AIC(beta.full.I)
r.squaredGLMM(beta.full.I)

qqnorm(residuals(beta.full.I))
scatter.smooth(residuals(beta.full.I) ~ fitted(beta.full.I))

# Model 1 - reduction ####

beta_mod = buildmer(beta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                       Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                       BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                     buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                       ddf = "Satterthwaite"))
summary(beta_mod)
print(beta_mod, correlation=TRUE)
vcov(beta_mod)

# Model 2 ####

beta.full.2 = lmer(beta ~ 1 + Silt + HIX + Peak_C + (1|Site), data = my_data.1)                                         

summary(beta.full.2)
Anova(beta.full.2)
cAIC(beta.full.2)
AIC(beta.full.2)
r.squaredGLMM(beta.full.2)

qqnorm(residuals(beta.full.2))
scatter.smooth(residuals(beta.full.2) ~ fitted(beta.full.2))

# Model 2 - without outliers - WINNER ####

cooksd <- cooks.distance(beta.full.2)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
my_data.beta = my_data.1[-c(59),]

beta.full.2r = lmer(beta ~ 1 + Silt + HIX + Peak_C + (1|Site), data = my_data.beta)                                         

summary(beta.full.2r)
Anova(beta.full.2r)
cAIC(beta.full.2r)
AIC(beta.full.2r)
r.squaredGLMM(beta.full.2r)

qqnorm(residuals(beta.full.2r))
scatter.smooth(residuals(beta.full.2r) ~ fitted(beta.full.2r))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(beta ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.beta, 
                          sets = list("Silt","HIX","Peak_C"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - xyl ####

# Model 1 ####

xyl.full.I = lmer(xyl ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                    Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                    BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1)                                         

summary(xyl.full.I)
Anova(xyl.full.I)
cAIC(xyl.full.I)
AIC(xyl.full.I)
r.squaredGLMM(xyl.full.I)

qqnorm(residuals(xyl.full.I))
scatter.smooth(residuals(xyl.full.I) ~ fitted(xyl.full.I))

# Model 1 - reduction ####

xyl_mod = buildmer(xyl ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                     Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                     BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                    buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                      ddf = "Satterthwaite"))
summary(xyl_mod)
print(xyl_mod, correlation=TRUE)
vcov(xyl_mod)

# Model 2 ####

xyl.full.2 = lmer(xyl ~ 1 + L_TN + HIX + AI + Litter + Soil_Temp + Soil_Temp:AI + (1|Site),
                  data = my_data.1)                                         

summary(xyl.full.2)
Anova(xyl.full.2)
cAIC(xyl.full.2)
AIC(xyl.full.2)
r.squaredGLMM(xyl.full.2)

qqnorm(residuals(xyl.full.2))
scatter.smooth(residuals(xyl.full.2) ~ fitted(xyl.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(xyl ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("L_TN","HIX","AI","Litter","Soil_Temp",
                                      "Soil_Temp:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - cbh ####

# Model 1 ####

cbh.full.I = lmer(cbh ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                    Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                    BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site),
                  data = my_data.1)                                         

summary(cbh.full.I)
Anova(cbh.full.I)
cAIC(cbh.full.I)
AIC(cbh.full.I)
r.squaredGLMM(cbh.full.I)

qqnorm(residuals(cbh.full.I))
scatter.smooth(residuals(cbh.full.I) ~ fitted(cbh.full.I))

# Model 1 - reduction ####

cbh_mod = buildmer(cbh ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                     Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                     BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(cbh_mod)
print(cbh_mod, correlation=TRUE)
vcov(cbh_mod)

# Model 2 - WINNER MODEL ####

cbh.full.2 = lmer(cbh ~ 1 + (1|Site), data = my_data.1)                                         

summary(cbh.full.2)
Anova(cbh.full.2)
cAIC(cbh.full.2)
AIC(cbh.full.2)
r.squaredGLMM(cbh.full.2)

qqnorm(residuals(cbh.full.2))
scatter.smooth(residuals(cbh.full.2) ~ fitted(cbh.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
#dominance_output <- domin(cbh ~ 1, 
#                          lmer, 
#                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
#                          data = my_data.1, 
#                          sets = list("Peak_A","L_TN","FI","Clay","HIX","Peak_T",
#                                      "SR","SOM","pH","AI","Peak_T:AI","L_TN:AI",
#                                      "HIX:AI","TC"), 
#                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - gla ####

# Model 1 ####

gla.full.I = lmer(log(gla) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                    Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                    BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site),
                  data = my_data.1)                                         

summary(gla.full.I)
Anova(gla.full.I)
cAIC(gla.full.I)
AIC(gla.full.I)
r.squaredGLMM(gla.full.I)

qqnorm(residuals(gla.full.I))
scatter.smooth(residuals(gla.full.I) ~ fitted(gla.full.I))

# Model 1 - reduction ####

gla_mod = buildmer(log(gla) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                          TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                          Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                          BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(gla_mod)
print(gla_mod, correlation=TRUE)
vcov(gla_mod)

# Model 2 - WINNER ####

gla.full.2 = lmer(log(gla) ~ 1 + Water_content + SR + Soil_Temp + C_N + altitude + pH + TC + (1|Site), 
                  data = my_data.1)                                         

summary(gla.full.2)
Anova(gla.full.2)
cAIC(gla.full.2)
AIC(gla.full.2)
r.squaredGLMM(gla.full.2)

qqnorm(residuals(gla.full.2))
scatter.smooth(residuals(gla.full.2) ~ fitted(gla.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(log(gla) ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("Water_content","SR","Soil_Temp","C_N",
                                      "altitude","pH","TC"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - fos ####

# Model 1 ####

fos.full.I = lmer(log(fos) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                    Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                    BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1)                                         

summary(fos.full.I)
Anova(fos.full.I)
cAIC(fos.full.I)
AIC(fos.full.I)
r.squaredGLMM(fos.full.I)

qqnorm(residuals(fos.full.I))
scatter.smooth(residuals(fos.full.I) ~ fitted(fos.full.I))

# Model 1 - reduction ####

fos_mod = buildmer(log(fos) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                     Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                     BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(fos_mod)
print(fos_mod, correlation=TRUE)
vcov(fos_mod)

# Model 2 - WINNER MODEL ####

fos.full.2 = lmer(log(fos) ~ 1 + SOM + Litter + PO43 + BB + Clay + Silt + (1|Site), 
                  data = my_data.1)                                         

summary(fos.full.2)
Anova(fos.full.2)
cAIC(fos.full.2)
AIC(fos.full.2)
r.squaredGLMM(fos.full.2)

qqnorm(residuals(fos.full.2))
scatter.smooth(residuals(fos.full.2) ~ fitted(fos.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(log(fos) ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("SOM","Litter","PO43","BB","Clay","Silt"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - leu ####

# Model 1 ####

leu.full.I = lmer(log(leu) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                    Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                    BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1)                                         

summary(leu.full.I)
Anova(leu.full.I)
cAIC(leu.full.I)
AIC(leu.full.I)
r.squaredGLMM(leu.full.I)

qqnorm(residuals(leu.full.I))
scatter.smooth(residuals(leu.full.I) ~ fitted(leu.full.I))

# Model 1 - reduction ####

leu_mod = buildmer(log(leu) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                          TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                          Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                          BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(leu_mod)
print(leu_mod, correlation=TRUE)
vcov(leu_mod)

# Model 2 - ####

leu.full.2 = lmer(log(leu) ~ 1 + TC + pH + C_N + (1|Site), data = my_data.1)                                         

summary(leu.full.2)
Anova(leu.full.2)
cAIC(leu.full.2)
AIC(leu.full.2)
r.squaredGLMM(leu.full.2)

qqnorm(residuals(leu.full.2))
scatter.smooth(residuals(leu.full.2) ~ fitted(leu.full.2))

# Model 2 - without outliers - WINNER ####

cooksd <- cooks.distance(leu.full.2)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
my_data.leu = my_data.1[-c(40,43,45),]

leu.full.2r = lmer(log(leu) ~ 1 + TC + pH + C_N + (1|Site), data = my_data.leu)                                         

summary(leu.full.2r)
Anova(leu.full.2r)
cAIC(leu.full.2r)
AIC(leu.full.2r)
r.squaredGLMM(leu.full.2r)

qqnorm(residuals(leu.full.2r))
scatter.smooth(residuals(leu.full.2r) ~ fitted(leu.full.2r))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(log(leu) ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.leu, 
                          sets = list("TC","pH","C_N"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - phe ####

# Model 1 ####

phe.full.I = lmer(phe ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                    Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                    BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1)                                         

summary(phe.full.I)
Anova(phe.full.I)
cAIC(phe.full.I)
AIC(phe.full.I)
r.squaredGLMM(phe.full.I)

qqnorm(residuals(phe.full.I))
scatter.smooth(residuals(phe.full.I) ~ fitted(phe.full.I))

# Model 1 - reduction ####

phe_mod = buildmer(phe ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                     Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                     BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(phe_mod)
print(phe_mod, correlation=TRUE)
vcov(phe_mod)

# Model 2 - WINNER MODEL ####

phe.full.2 = lmer(phe ~ 1 + LTC_LTN + TC + SOM + Water_content + Peak_B + E4.E6 +  
                    HIX + Soil_Temp + (1|Site), data = my_data.1)                                         

summary(phe.full.2)
Anova(phe.full.2)
cAIC(phe.full.2)
AIC(phe.full.2)
r.squaredGLMM(phe.full.2)

qqnorm(residuals(phe.full.2))
scatter.smooth(residuals(phe.full.2) ~ fitted(phe.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(phe ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("LTC_LTN","TC","SOM","Water_content","Peak_B",
                                      "E4.E6","HIX","Soil_Temp"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - xylcbh ####

# Model 1 ####

xylcbh.full.I = lmer(xylcbh ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                          TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                          Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                          BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1)                                         

summary(xylcbh.full.I)
Anova(xylcbh.full.I)
cAIC(xylcbh.full.I)
AIC(xylcbh.full.I)
r.squaredGLMM(xylcbh.full.I)

qqnorm(residuals(xylcbh.full.I))
scatter.smooth(residuals(xylcbh.full.I) ~ fitted(xylcbh.full.I))

# Model 1 - reduction ####

xylcbh_mod = buildmer(xylcbh ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                           TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                           Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                           BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(xylcbh_mod)
print(xylcbh_mod, correlation=TRUE)
vcov(xylcbh_mod)

# Model 2 - WINNER MODEL ####

xylcbh.full.2 = lmer(xylcbh ~ 1 + L_TN + HIX + AI + L_TN:AI + C_N + Soil_Temp + E4.E6 +  
                       Silt + FB + TC + Silt:AI + (1|Site), data = my_data.1)                                         

summary(xylcbh.full.2)
Anova(xylcbh.full.2)
cAIC(xylcbh.full.2)
AIC(xylcbh.full.2)
r.squaredGLMM(xylcbh.full.2)

qqnorm(residuals(xylcbh.full.2))
scatter.smooth(residuals(xylcbh.full.2) ~ fitted(xylcbh.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(xylcbh ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("L_TN","HIX","AI","L_TN:AI","C_N","Soil_Temp",
                                      "E4.E6","Silt","FB","TC","Silt:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - alphabeta ####

# Model 1 ####

alphabeta.full.I = lmer(log(alphabeta) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                                Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                                BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1)                                         

summary(alphabeta.full.I)
Anova(alphabeta.full.I)
cAIC(alphabeta.full.I)
AIC(alphabeta.full.I)
r.squaredGLMM(alphabeta.full.I)

qqnorm(residuals(alphabeta.full.I))
scatter.smooth(residuals(alphabeta.full.I) ~ fitted(alphabeta.full.I))

# Model 1 - reduction ####

alphabeta_mod = buildmer(log(alphabeta) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                      TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                                      Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                                      BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                      buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                        ddf = "Satterthwaite"))
summary(alphabeta_mod)
print(alphabeta_mod, correlation=TRUE)
vcov(alphabeta_mod)

# Model 2 ####

alphabeta.full.2 = lmer(log(alphabeta) ~ 1 + Silt + L_TN + HIX + Soil_Temp + C_N + SR +  
                          FB + PO43 + Water_content + SO42 + AI + C_N:AI + PO43:AI +  
                          E4.E6 + L_TN:AI + BIX + LTC_LTN + HIX:AI + SO42:AI + Litter +  
                          AI:Litter + BB + AI:BIX + SR:AI + SOM + NH4 + (1|Site), data = my_data.1)                                         

summary(alphabeta.full.2)
Anova(alphabeta.full.2)
cAIC(alphabeta.full.2)
AIC(alphabeta.full.2)
r.squaredGLMM(alphabeta.full.2)

qqnorm(residuals(alphabeta.full.2))
scatter.smooth(residuals(alphabeta.full.2) ~ fitted(alphabeta.full.2))

# Model 1 - without outliers - WINNER ####

cooksd <- cooks.distance(alphabeta.full.2)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
my_data.alphabeta = my_data.1[-c(6,9,37,45),]

alphabeta.full.2r = lmer(log(alphabeta) ~ 1 + Silt + L_TN + HIX + Soil_Temp + C_N + SR +  
                          FB + PO43 + Water_content + SO42 + AI + C_N:AI + PO43:AI +  
                          E4.E6 + L_TN:AI + BIX + LTC_LTN + HIX:AI + SO42:AI + Litter +  
                          AI:Litter + BB + AI:BIX + SR:AI + SOM + NH4 + (1|Site), data = my_data.alphabeta)                                         

summary(alphabeta.full.2r)
Anova(alphabeta.full.2r)
cAIC(alphabeta.full.2r)
AIC(alphabeta.full.2r)
r.squaredGLMM(alphabeta.full.2r)

qqnorm(residuals(alphabeta.full.2r))
scatter.smooth(residuals(alphabeta.full.2r) ~ fitted(alphabeta.full.2r))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(log(alphabeta) ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.alphabeta, 
                          sets = list("Silt","L_TN","HIX","Soil_Temp","C_N","SR",
                                      "FB","PO43","Water_content","SO42","AI",
                                      "C_N:AI","PO43:AI","E4.E6","L_TN:AI","BIX",
                                      "LTC_LTN","HIX:AI","SO42:AI","Litter",
                                      "AI:Litter","BB","AI:BIX","SR:AI","SOM","NH4"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - Cenz ####

# Model 1 ####

Cenz.full.I = lmer(log(Cenz) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                      TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                      Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                      BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1)                                         

summary(Cenz.full.I)
Anova(Cenz.full.I)
cAIC(Cenz.full.I)
AIC(Cenz.full.I)
r.squaredGLMM(Cenz.full.I)

qqnorm(residuals(Cenz.full.I))
scatter.smooth(residuals(Cenz.full.I) ~ fitted(Cenz.full.I))

# Model 1 - reduction ####

Cenz_mod = buildmer(log(Cenz) ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TC+C_N+NH4+PO43+SO42+Clay+Silt+
                                       Litter+L_TN+BB+FB+SR+E2.E3+E4.E6+
                                       BIX+Peak_C+Peak_B+HIX+LTC_LTN)*AI + (1|Site), data = my_data.1,
                         buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                           ddf = "Satterthwaite"))
summary(Cenz_mod)
print(Cenz_mod, correlation=TRUE)
vcov(Cenz_mod)

# Model 2 - WINNER ####

Cenz.full.2 = lmer(log(Cenz) ~ 1 + Soil_Temp + C_N + TC + HIX + (1|Site), data = my_data.1)                                         

summary(Cenz.full.2)
Anova(Cenz.full.2)
cAIC(Cenz.full.2)
AIC(Cenz.full.2)
r.squaredGLMM(Cenz.full.2)

qqnorm(residuals(Cenz.full.2))
scatter.smooth(residuals(Cenz.full.2) ~ fitted(Cenz.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(Cenz ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("Soil_Temp","C_N","TC","HIX"), 
                          consmodel = "(1|Site)") # Replace with your actual function

