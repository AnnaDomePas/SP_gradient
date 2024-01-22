library(buildmer)
library(lmerTest)

max_values      = apply(my_data[,c(8:34,40,41,76,81:93)],2,max)
my_data.1       = as.data.frame(cbind(my_data[,c(2,7,46:53,64,94,95)],
                                      my_data[,c(8:34,40,41,76,81:93)] / as.list(max_values)))
a.1             = as.data.frame(colnames(my_data.1))

corr            = as.data.frame(cor(my_data.1[,14:56]))

my_data.1 <-my_data.1 %>%
  mutate(Cenz = select(., 3:6) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(MB = select(.,41:42) %>% rowSums(na.rm = TRUE))

# Respiration ####

# Model 1 ####

respiration.full.I = lmer(Respiration ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                    TOC+TC+TN+C_N+NH4+PO43+SO42+
                                                    Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                                    BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(respiration.full.I)
Anova(respiration.full.I)
cAIC(respiration.full.I)
AIC(respiration.full.I)
r.squaredGLMM(respiration.full.I)

qqnorm(residuals(respiration.full.I))
scatter.smooth(residuals(respiration.full.I) ~ fitted(respiration.full.I))

# Model 1 - reduction ####

max_mod = buildmer(Respiration ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                             TOC+TC+TN+C_N+NH4+PO43+SO42+
                                             Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                             BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                             Peak_T)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(max_mod)
print(max_mod, correlation=TRUE)
vcov(max_mod)

max_mod.1 = buildmer(Respiration ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                             TOC+TC+TN+C_N+NH4+PO43+SO42+
                                             Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                             BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                             Peak_T)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),
                                                     calc.anova = TRUE,direction='forward',
                                                     ddf = "Satterthwaite"))
summary(max_mod.1)
print(max_mod.1, correlation=TRUE)
vcov(max_mod.1)

# Model 2 - WINNER MODEL ####

respiration.full.2 = lmer(Respiration ~ 1 + Water_content + Clay + Silt + TC + SOM + TOC +  
                            AI + Clay:AI + Silt:AI + altitude + SR + AI:SR + Soil_Temp +  
                            L_TN + FI + AI:FI + C_N + AI:L_TN + Water_content:AI + SO42 +  
                            Litter + NH4 + AI:NH4 + AI:SO42 + SOM:AI + (1 | Site),
                          data = my_data.1)                                         

summary(respiration.full.2)
Anova(respiration.full.2)
cAIC(respiration.full.2)
AIC(respiration.full.2)
r.squaredGLMM(respiration.full.2)

qqnorm(residuals(respiration.full.2))
scatter.smooth(residuals(respiration.full.2) ~ fitted(respiration.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(Respiration ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("Water_content","Clay","Silt","TC","SOM",
                                      "TOC","AI","Clay:AI","Silt:AI","altitude",
                                      "SR","AI:SR","Soil_Temp","L_TN","FI",
                                      "AI:FI","C_N","AI:L_TN","Water_content:AI",
                                      "SO42","Litter","NH4","AI:NH4","AI:SO42","SOM:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 ####

respiration.full.3 = lmer(Respiration ~ 1 + Water_content + Clay + Silt + TC + SOM + TOC +  
                            AI + Clay:AI + SR + AI:SR + Soil_Temp + L_TN + E2.E3 + FI +  
                            AI:FI + SO42 + Litter + NH4 + SOM:AI + (1 | Site),
                          data = my_data.1)                                         

summary(respiration.full.3)
Anova(respiration.full.3)
cAIC(respiration.full.3)
AIC(respiration.full.3)
r.squaredGLMM(respiration.full.3)

qqnorm(residuals(respiration.full.3))
scatter.smooth(residuals(respiration.full.3) ~ fitted(respiration.full.3))

# Biomass - BB ####

# Model 1 ####

BB.full.I = lmer(BB ~ altitude+(Soil_Temp+Water_content+SOM+pH+TOC+TC+TN+C_N+NH4+PO43+SO42+
                                                    Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                                    FB+SR+E2.E3+FI+HIX+Peak_A+Peak_T)*AI + 
                   (1|Site), data = my_data.1)                                         

summary(BB.full.I)
Anova(BB.full.I)
cAIC(BB.full.I)
AIC(BB.full.I)
r.squaredGLMM(BB.full.I)

qqnorm(residuals(BB.full.I))
scatter.smooth(residuals(BB.full.I) ~ fitted(BB.full.I))

# Model 1 - reduction ####

max_BB = buildmer(BB ~ altitude+(Soil_Temp+Water_content+SOM+pH+TOC+TC+TN+C_N+NH4+PO43+SO42+
                                   Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                   FB+SR+E2.E3+FI+HIX+Peak_A+Peak_T)*AI + 
                    (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(max_BB)
print(max_BB, correlation=TRUE)
vcov(max_BB)

max_BB.1 = buildmer(BB ~ altitude+(Soil_Temp+Water_content+SOM+pH+TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     FB+SR+E2.E3+FI+HIX+Peak_A+Peak_T)*AI + 
                      (1|Site), data = my_data.1,
                     buildmerControl = buildmerControl(include = ~ (1|Site),
                                                       calc.anova = TRUE,direction='forward',
                                                       ddf = "Satterthwaite"))
summary(max_BB.1)
print(max_BB.1, correlation=TRUE)
vcov(max_BB.1)

# Model 2 - WINNER MODEL ####

BB.full.2 = lmer(BB ~ 1 + NH4 + SOM + pH + FB + E2.E3 + Water_content + AI + NH4:AI +  
                   Clay + (1|Site), data = my_data.1)                                         

summary(BB.full.2)
Anova(BB.full.2)
cAIC(BB.full.2)
AIC(BB.full.2)
r.squaredGLMM(BB.full.2)

qqnorm(residuals(BB.full.2))
scatter.smooth(residuals(BB.full.2) ~ fitted(BB.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(BB ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("NH4","SOM","pH","FB","E2.E3","Water_content",
                                      "AI","NH4:AI","Clay"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 ####

BB.full.3 = lmer(BB ~ 1 + TN + NH4 + Peak_A + Water_content + AI + NH4:AI + 
                   (1|Site), data = my_data.1)                                         

summary(BB.full.3)
Anova(BB.full.3)
cAIC(BB.full.3)
AIC(BB.full.3)
r.squaredGLMM(BB.full.3)

qqnorm(residuals(BB.full.3))
scatter.smooth(residuals(BB.full.3) ~ fitted(BB.full.3))

# Biomass - FB ####

# Model 1 ####

FB.full.I = lmer(FB ~ altitude+(Soil_Temp+Water_content+SOM+pH+TOC+TC+TN+C_N+NH4+PO43+SO42+
                                  Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                  BB+SR+E2.E3+FI+HIX+Peak_A+Peak_T)*AI + 
                   (1|Site), data = my_data.1)                                          

summary(FB.full.I)
Anova(FB.full.I)
cAIC(FB.full.I)
AIC(FB.full.I)
r.squaredGLMM(FB.full.I)

qqnorm(residuals(FB.full.I))
scatter.smooth(residuals(FB.full.I) ~ fitted(FB.full.I))

# Model 1 - reduction ####

max_FB = buildmer(FB ~ altitude+(Soil_Temp+Water_content+SOM+pH+TOC+TC+TN+C_N+NH4+PO43+SO42+
                                   Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                   BB+SR+E2.E3+FI+HIX+Peak_A+Peak_T)*AI + 
                    (1|Site), data = my_data.1,
                  buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                    ddf = "Satterthwaite"))
summary(max_FB)
print(max_FB, correlation=TRUE)
vcov(max_FB)

max_FB.1 = buildmer(FB ~ altitude+(Soil_Temp+Water_content+SOM+pH+TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     BB+SR+E2.E3+FI+HIX+Peak_A+Peak_T)*AI + 
                      (1|Site), data = my_data.1,
                    buildmerControl = buildmerControl(include = ~ (1|Site),
                                                      calc.anova = TRUE,direction='forward',
                                                      ddf = "Satterthwaite"))
summary(max_FB.1)
print(max_FB.1, correlation=TRUE)
vcov(max_FB.1)

# Model 2 -  ####

FB.full.2 = lmer(FB ~ 1 + TOC + Water_content + BB + SOM + NH4 + E2.E3 + Peak_A +  
                   SO42 + AI + pH + BB:AI + Soil_Temp + (1 | Site), data = my_data.1)                                         

summary(FB.full.2)
Anova(FB.full.2)
cAIC(FB.full.2)
AIC(FB.full.2)
r.squaredGLMM(FB.full.2)

qqnorm(residuals(FB.full.2))
scatter.smooth(residuals(FB.full.2) ~ fitted(FB.full.2))

# Model 3 - WINNER MODEL ####

FB.full.3 = lmer(FB ~ 1 + TOC + Water_content + BB + SOM + NH4 + E2.E3 + SR +  
                   SO42 + AI + pH + Soil_Temp + (1 | Site), data = my_data.1)                                         

summary(FB.full.3)
Anova(FB.full.3)
cAIC(FB.full.3)
AIC(FB.full.3)
r.squaredGLMM(FB.full.3)

qqnorm(residuals(FB.full.3))
scatter.smooth(residuals(FB.full.3) ~ fitted(FB.full.3))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(FB ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("TOC","Water_content","BB","SOM","NH4",
                                      "E2.E3","SR","SO42","AI","pH","Soil_Temp"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Biomass - MB ####

# Model 1 ####

MB.full.I = lmer(MB ~ altitude+(Soil_Temp+Water_content+SOM+pH+TOC+TC+TN+C_N+NH4+PO43+SO42+
                                  Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                  SR+E2.E3+FI+HIX+Peak_A+Peak_T)*AI + 
                   (1|Site), data = my_data.1)                                          

summary(MB.full.I)
Anova(MB.full.I)
cAIC(MB.full.I)
AIC(MB.full.I)
r.squaredGLMM(MB.full.I)

qqnorm(residuals(MB.full.I))
scatter.smooth(residuals(MB.full.I) ~ fitted(MB.full.I))

# Model 1 - reduction ####

max_MB = buildmer(MB ~ altitude+(Soil_Temp+Water_content+SOM+pH+TOC+TC+TN+C_N+NH4+PO43+SO42+
                                   Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                   SR+E2.E3+FI+HIX+Peak_A+Peak_T)*AI + 
                    (1|Site), data = my_data.1,
                  buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                    ddf = "Satterthwaite"))
summary(max_MB)
print(max_MB, correlation=TRUE)
vcov(max_MB)

max_MB.1 = buildmer(MB ~ altitude+(Soil_Temp+Water_content+SOM+pH+TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     SR+E2.E3+FI+HIX+Peak_A+Peak_T)*AI + 
                      (1|Site), data = my_data.1,
                    buildmerControl = buildmerControl(include = ~ (1|Site),
                                                      calc.anova = TRUE,direction='forward',
                                                      ddf = "Satterthwaite"))
summary(max_MB.1)
print(max_MB.1, correlation=TRUE)
vcov(max_MB.1)

# Model 2 - WINNER MODEL ####

MB.full.2 = lmer(MB ~ 1 + TOC + NH4 + (1|Site), data = my_data.1)                                          

summary(MB.full.2)
Anova(MB.full.2)
cAIC(MB.full.2)
AIC(MB.full.2)
r.squaredGLMM(MB.full.2)

qqnorm(residuals(MB.full.2))
scatter.smooth(residuals(MB.full.2) ~ fitted(MB.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(MB ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("TOC","NH4"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - alpha ####

# Model 1 ####

alpha.full.I = lmer(alpha ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                    TOC+TC+TN+C_N+NH4+PO43+SO42+
                                                    Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                                    BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(alpha.full.I)
Anova(alpha.full.I)
cAIC(alpha.full.I)
AIC(alpha.full.I)
r.squaredGLMM(alpha.full.I)

qqnorm(residuals(alpha.full.I))
scatter.smooth(residuals(alpha.full.I) ~ fitted(alpha.full.I))

# Model 1 - reduction ####

alpha_mod = buildmer(alpha ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                             TOC+TC+TN+C_N+NH4+PO43+SO42+
                                             Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                             BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                             Peak_T)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(alpha_mod)
print(alpha_mod, correlation=TRUE)
vcov(alpha_mod)

alpha_mod.1 = buildmer(alpha ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                               TOC+TC+TN+C_N+NH4+PO43+SO42+
                                               Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                               BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                               Peak_T)*AI + (1|Site), data = my_data.1,
                     buildmerControl = buildmerControl(include = ~ (1|Site),
                                                       calc.anova = TRUE,direction='forward',
                                                       ddf = "Satterthwaite"))
summary(alpha_mod.1)
print(alpha_mod.1, correlation=TRUE)
vcov(alpha_mod.1)

# Model 2 ####

alpha.full.2 = lmer(alpha ~ 1 + Soil_Temp + (1|Site), data = my_data.1)                                         

summary(alpha.full.2)
Anova(alpha.full.2)
cAIC(alpha.full.2)
AIC(alpha.full.2)
r.squaredGLMM(alpha.full.2)

qqnorm(residuals(alpha.full.2))
scatter.smooth(residuals(alpha.full.2) ~ fitted(alpha.full.2))
plot(cooks.distance(alpha.full.2))

# Model 3 ####

alpha.full.3 = lmer(alpha ~ 1 + Soil_Temp + FI + SR + SOM + PO43 + Litter + 
                      NH4 + E2.E3 + AI + (1|Site), data = my_data.1)                                         

summary(alpha.full.3)
Anova(alpha.full.3)
cAIC(alpha.full.3)
AIC(alpha.full.3)
r.squaredGLMM(alpha.full.3)

qqnorm(residuals(alpha.full.3))
scatter.smooth(residuals(alpha.full.3) ~ fitted(alpha.full.3))

# Model 1 - without outliers ####

plot(cooks.distance(alpha.full.2))
my_data.2 = my_data.1[-c(49,59),]

alpha.full.I.a = lmer(alpha ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                        TOC+TC+TN+C_N+NH4+PO43+SO42+
                                        Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                        BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                        Peak_T)*AI + (1|Site), data = my_data.2)                                         

summary(alpha.full.I.a)
Anova(alpha.full.I.a)
cAIC(alpha.full.I.a)
AIC(alpha.full.I.a)
r.squaredGLMM(alpha.full.I.a)

qqnorm(residuals(alpha.full.I.a))
scatter.smooth(residuals(alpha.full.I.a) ~ fitted(alpha.full.I.a))

# Model 1 - without outliers - reduction ####

alpha_mod.a = buildmer(alpha ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                         TOC+TC+TN+C_N+NH4+PO43+SO42+
                                         Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                         BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                         Peak_T)*AI + (1|Site), data = my_data.2,
                     buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                       ddf = "Satterthwaite"))
summary(alpha_mod.a)
print(alpha_mod.a, correlation=TRUE)
vcov(alpha_mod.a)

alpha_mod.a.1 = buildmer(alpha ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                           TOC+TC+TN+C_N+NH4+PO43+SO42+
                                           Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                           BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                           Peak_T)*AI + (1|Site), data = my_data.2,
                       buildmerControl = buildmerControl(include = ~ (1|Site),
                                                         calc.anova = TRUE,direction='forward',
                                                         ddf = "Satterthwaite"))
summary(alpha_mod.a.1)
print(alpha_mod.a.1, correlation=TRUE)
vcov(alpha_mod.a.1)

# Model 2 - without outliers - WINNER MODEL ####

alpha.full.2.a = lmer(alpha ~ 1 + Soil_Temp + TN + C_N + SOM + E2.E3 + BB + Peak_A +  
                        Clay + Peak_T + (1|Site), data = my_data.2)                                         

summary(alpha.full.2.a)
Anova(alpha.full.2.a)
cAIC(alpha.full.2.a)
AIC(alpha.full.2.a)
r.squaredGLMM(alpha.full.2.a)

qqnorm(residuals(alpha.full.2.a))
scatter.smooth(residuals(alpha.full.2.a) ~ fitted(alpha.full.2.a))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(alpha ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.2, 
                          sets = list("Soil_Temp","TN","C_N","SOM","E2.E3","BB",
                                      "Peak_A","Clay","Peak_T"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 - without outliers ####

alpha.full.3.a = lmer(alpha ~ 1 + Soil_Temp + TN + HIX + C_N + SOM + E2.E3 + 
                        (1|Site), data = my_data.2)                                         

summary(alpha.full.3.a)
Anova(alpha.full.3.a)
cAIC(alpha.full.3.a)
AIC(alpha.full.3.a)
r.squaredGLMM(alpha.full.3.a)

qqnorm(residuals(alpha.full.3.a))
scatter.smooth(residuals(alpha.full.3.a) ~ fitted(alpha.full.3.a))

# Enzyme - beta ####

# Model 1 ####

beta.full.I = lmer(beta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                        TOC+TC+TN+C_N+NH4+PO43+SO42+
                                        Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                        BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                        Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(beta.full.I)
Anova(beta.full.I)
cAIC(beta.full.I)
AIC(beta.full.I)
r.squaredGLMM(beta.full.I)

qqnorm(residuals(beta.full.I))
scatter.smooth(residuals(beta.full.I) ~ fitted(beta.full.I))

# Model 1 - reduction ####

beta_mod = buildmer(beta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                         TOC+TC+TN+C_N+NH4+PO43+SO42+
                                         Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                         BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                         Peak_T)*AI + (1|Site), data = my_data.1,
                     buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                       ddf = "Satterthwaite"))
summary(beta_mod)
print(beta_mod, correlation=TRUE)
vcov(beta_mod)

beta_mod.1 = buildmer(beta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                           TOC+TC+TN+C_N+NH4+PO43+SO42+
                                           Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                           BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                           Peak_T)*AI + (1|Site), data = my_data.1,
                       buildmerControl = buildmerControl(include = ~ (1|Site),
                                                         calc.anova = TRUE,direction='forward',
                                                         ddf = "Satterthwaite"))
summary(beta_mod.1)
print(beta_mod.1, correlation=TRUE)
vcov(beta_mod.1)

# Model 2 ####

beta.full.2 = lmer(beta ~ 1 + HIX + Soil_Temp + FI + AI + (1|Site), data = my_data.1)                                         

summary(beta.full.2)
Anova(beta.full.2)
cAIC(beta.full.2)
AIC(beta.full.2)
r.squaredGLMM(beta.full.2)

qqnorm(residuals(beta.full.2))
scatter.smooth(residuals(beta.full.2) ~ fitted(beta.full.2))

# Model 3 ####

beta.full.3 = lmer(beta ~ 1 + Silt + HIX + Peak_A + FI + AI + Clay + SOM + 
                     Litter + (1|Site), data = my_data.1)                                         

summary(beta.full.3)
Anova(beta.full.3)
cAIC(beta.full.3)
AIC(beta.full.3)
r.squaredGLMM(beta.full.3)

qqnorm(residuals(beta.full.3))
scatter.smooth(residuals(beta.full.3) ~ fitted(beta.full.3))

# Model 1 - without outliers ####

plot(cooks.distance(beta.full.3))
my_data.3 = my_data.1[-c(59),]

beta.full.I.a = lmer(beta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                          TOC+TC+TN+C_N+NH4+PO43+SO42+
                                          Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                          BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                          Peak_T)*AI + (1|Site), data = my_data.3)                                         

summary(beta.full.I.a)
Anova(beta.full.I.a)
cAIC(beta.full.I.a)
AIC(beta.full.I.a)
r.squaredGLMM(beta.full.I.a)

qqnorm(residuals(beta.full.I.a))
scatter.smooth(residuals(beta.full.I.a) ~ fitted(beta.full.I.a))

# Model 1 - without outliers - reduction ####

beta_mod.a = buildmer(beta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.3,
                    buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                      ddf = "Satterthwaite"))
summary(beta_mod.a)
print(beta_mod.a, correlation=TRUE)
vcov(beta_mod.a)

beta_mod.1.a = buildmer(beta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                         TOC+TC+TN+C_N+NH4+PO43+SO42+
                                         Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                         BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                         Peak_T)*AI + (1|Site), data = my_data.3,
                      buildmerControl = buildmerControl(include = ~ (1|Site),
                                                        calc.anova = TRUE,direction='forward',
                                                        ddf = "Satterthwaite"))
summary(beta_mod.1.a)
print(beta_mod.1.a, correlation=TRUE)
vcov(beta_mod.1.a)

# Model 2 - without outliers - WINNER MODEL ####

beta.full.2.a = lmer(beta ~ 1 + L_TN + Soil_Temp + TOC + FB + TC + Peak_A + 
                       Water_content + (1|Site), data = my_data.3)                                         

summary(beta.full.2.a)
Anova(beta.full.2.a)
cAIC(beta.full.2.a)
AIC(beta.full.2.a)
r.squaredGLMM(beta.full.2.a)

qqnorm(residuals(beta.full.2.a))
scatter.smooth(residuals(beta.full.2.a) ~ fitted(beta.full.2.a))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(beta ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.3, 
                          sets = list("L_TN","Soil_Temp","TOC","FB","TC","Peak_A",
                                      "Water_content"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 - without outliers ####

beta.full.3.a = lmer(beta ~ 1 + Silt + FB + (1|Site), data = my_data.3)                                         

summary(beta.full.3.a)
Anova(beta.full.3.a)
cAIC(beta.full.3.a)
AIC(beta.full.3.a)
r.squaredGLMM(beta.full.3.a)

qqnorm(residuals(beta.full.3.a))
scatter.smooth(residuals(beta.full.3.a) ~ fitted(beta.full.3.a))

# Enzyme - xyl ####

# Model 1 ####

xyl.full.I = lmer(xyl ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                      TOC+TC+TN+C_N+NH4+PO43+SO42+
                                      Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                      BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                      Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(xyl.full.I)
Anova(xyl.full.I)
cAIC(xyl.full.I)
AIC(xyl.full.I)
r.squaredGLMM(xyl.full.I)

qqnorm(residuals(xyl.full.I))
scatter.smooth(residuals(xyl.full.I) ~ fitted(xyl.full.I))

# Model 1 - reduction ####

xyl_mod = buildmer(xyl ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.1,
                    buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                      ddf = "Satterthwaite"))
summary(xyl_mod)
print(xyl_mod, correlation=TRUE)
vcov(xyl_mod)

xyl_mod.1 = buildmer(xyl ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                         TOC+TC+TN+C_N+NH4+PO43+SO42+
                                         Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                         BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                         Peak_T)*AI + (1|Site), data = my_data.1,
                      buildmerControl = buildmerControl(include = ~ (1|Site),
                                                        calc.anova = TRUE,direction='forward',
                                                        ddf = "Satterthwaite"))
summary(xyl_mod.1)
print(xyl_mod.1, correlation=TRUE)
vcov(xyl_mod.1)

# Model 2 ####

xyl.full.2 = lmer(xyl ~ 1 + L_TN + (1|Site), data = my_data.1)                                         

summary(xyl.full.2)
Anova(xyl.full.2)
cAIC(xyl.full.2)
AIC(xyl.full.2)
r.squaredGLMM(xyl.full.2)

qqnorm(residuals(xyl.full.2))
scatter.smooth(residuals(xyl.full.2) ~ fitted(xyl.full.2))

# Model 3 - WINNER MODEL ####

xyl.full.3 = lmer(xyl ~ 1 + L_TN + HIX + AI + Litter + Soil_Temp + Silt + 
                    SO42 + (1|Site), data = my_data.1)                                         

summary(xyl.full.3)
Anova(xyl.full.3)
cAIC(xyl.full.3)
AIC(xyl.full.3)
r.squaredGLMM(xyl.full.3)

qqnorm(residuals(xyl.full.3))
scatter.smooth(residuals(xyl.full.3) ~ fitted(xyl.full.3))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(xyl ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("L_TN","HIX","AI","Litter","Soil_Temp","Silt","SO42"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - cbh ####

# Model 1 ####

cbh.full.I = lmer(cbh ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TOC+TC+TN+C_N+NH4+PO43+SO42+
                                    Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                    BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(cbh.full.I)
Anova(cbh.full.I)
cAIC(cbh.full.I)
AIC(cbh.full.I)
r.squaredGLMM(cbh.full.I)

qqnorm(residuals(cbh.full.I))
scatter.smooth(residuals(cbh.full.I) ~ fitted(cbh.full.I))

# Model 1 - reduction ####

cbh_mod = buildmer(cbh ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                     Peak_T)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(cbh_mod)
print(cbh_mod, correlation=TRUE)
vcov(cbh_mod)

cbh_mod.1 = buildmer(cbh ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.1,
                     buildmerControl = buildmerControl(include = ~ (1|Site),
                                                       calc.anova = TRUE,direction='forward',
                                                       ddf = "Satterthwaite"))
summary(cbh_mod.1)
print(cbh_mod.1, correlation=TRUE)
vcov(cbh_mod.1)

# Model 2 - WINNER MODEL ####

cbh.full.2 = lmer(cbh ~ 1 + Peak_A + L_TN + FI + Clay + HIX + Peak_T + SR + SOM +  
                    pH + AI + Peak_T:AI + L_TN:AI + HIX:AI + TC + (1|Site), data = my_data.1)                                         

summary(cbh.full.2)
Anova(cbh.full.2)
cAIC(cbh.full.2)
AIC(cbh.full.2)
r.squaredGLMM(cbh.full.2)

qqnorm(residuals(cbh.full.2))
scatter.smooth(residuals(cbh.full.2) ~ fitted(cbh.full.2))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(cbh ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("Peak_A","L_TN","FI","Clay","HIX","Peak_T",
                                      "SR","SOM","pH","AI","Peak_T:AI","L_TN:AI",
                                      "HIX:AI","TC"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 ####

cbh.full.3 = lmer(cbh ~ 1 + Peak_A + L_TN + HIX + Peak_T + SO42 + AI + 
                    Peak_T:AI + L_TN:AI + (1|Site), data = my_data.1)                                         

summary(cbh.full.3)
Anova(cbh.full.3)
cAIC(cbh.full.3)
AIC(cbh.full.3)
r.squaredGLMM(cbh.full.3)

qqnorm(residuals(cbh.full.3))
scatter.smooth(residuals(cbh.full.3) ~ fitted(cbh.full.3))

# Enzyme - gla ####

# Model 1 ####

gla.full.I = lmer(gla ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TOC+TC+TN+C_N+NH4+PO43+SO42+
                                    Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                    BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(gla.full.I)
Anova(gla.full.I)
cAIC(gla.full.I)
AIC(gla.full.I)
r.squaredGLMM(gla.full.I)

qqnorm(residuals(gla.full.I))
scatter.smooth(residuals(gla.full.I) ~ fitted(gla.full.I))

# Model 1 - reduction ####

gla_mod = buildmer(gla ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                     Peak_T)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(gla_mod)
print(gla_mod, correlation=TRUE)
vcov(gla_mod)

gla_mod.1 = buildmer(gla ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.1,
                     buildmerControl = buildmerControl(include = ~ (1|Site),
                                                       calc.anova = TRUE,direction='forward',
                                                       ddf = "Satterthwaite"))
summary(gla_mod.1)
print(gla_mod.1, correlation=TRUE)
vcov(gla_mod.1)

# Model 2 ####

gla.full.2 = lmer(gla ~ 1 + AI + FB + SO42 + AI:SO42 + BB + FI + PO43 + TC + AI:FI +  
                    AI:FB + AI:PO43 + Soil_Temp + (1|Site), data = my_data.1)                                         

summary(gla.full.2)
Anova(gla.full.2)
cAIC(gla.full.2)
AIC(gla.full.2)
r.squaredGLMM(gla.full.2)

qqnorm(residuals(gla.full.2))
scatter.smooth(residuals(gla.full.2) ~ fitted(gla.full.2))

# Model 3 ####

gla.full.3 = lmer(gla ~ 1 + Water_content + SR + altitude + AI + FB + SO42 + NH4 +  
                    BB + FI + PO43 + AI:PO43 + Soil_Temp + HIX + (1|Site), data = my_data.1)                                         

summary(gla.full.3)
Anova(gla.full.3)
cAIC(gla.full.3)
AIC(gla.full.3)
r.squaredGLMM(gla.full.3)

qqnorm(residuals(gla.full.3))
scatter.smooth(residuals(gla.full.3) ~ fitted(gla.full.3))

# Model 1 - without outliers ####

plot(cooks.distance(gla.full.3))
my_data.4 = my_data.1[-c(39),]

gla.full.Ia = lmer(gla ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TOC+TC+TN+C_N+NH4+PO43+SO42+
                                    Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                    BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.4)                                         

summary(gla.full.Ia)
Anova(gla.full.Ia)
cAIC(gla.full.Ia)
AIC(gla.full.Ia)
r.squaredGLMM(gla.full.Ia)

qqnorm(residuals(gla.full.Ia))
scatter.smooth(residuals(gla.full.Ia) ~ fitted(gla.full.Ia))

# Model 1 - without outliers - reduction ####

gla_mod.a = buildmer(gla ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                     Peak_T)*AI + (1|Site), data = my_data.4,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(gla_mod.a)
print(gla_mod.a, correlation=TRUE)
vcov(gla_mod.a)

gla_mod.1.a = buildmer(gla ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.4,
                     buildmerControl = buildmerControl(include = ~ (1|Site),
                                                       calc.anova = TRUE,direction='forward',
                                                       ddf = "Satterthwaite"))
summary(gla_mod.1.a)
print(gla_mod.1.a, correlation=TRUE)
vcov(gla_mod.1.a)

# Model 2 - without outliers - WINNER MODEL ####

gla.full.2a = lmer(gla ~ 1 + Water_content + SR + altitude + TOC + Clay + 
                     AI + pH + (1|Site), data = my_data.4)                                         

summary(gla.full.2a)
Anova(gla.full.2a)
cAIC(gla.full.2a)
AIC(gla.full.2a)
r.squaredGLMM(gla.full.2a)

qqnorm(residuals(gla.full.2a))
scatter.smooth(residuals(c.full.2a) ~ fitted(gla.full.2a))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(gla ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.4, 
                          sets = list("Water_content","SR","altitude","TOC",
                                      "Clay","AI","pH"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 - without outliers ####

gla.full.3a = lmer(gla ~ 1 + Water_content + SR + Soil_Temp + altitude + TOC + AI +  
                     pH + E2.E3 + NH4 + (1|Site), data = my_data.4)                                         

summary(gla.full.3a)
Anova(gla.full.3a)
cAIC(gla.full.3a)
AIC(gla.full.3a)
r.squaredGLMM(gla.full.3a)

qqnorm(residuals(gla.full.3a))
scatter.smooth(residuals(gla.full.3a) ~ fitted(gla.full.3a))

# Enzyme - fos ####

# Model 1 ####

fos.full.I = lmer(fos ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TOC+TC+TN+C_N+NH4+PO43+SO42+
                                    Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                    BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(fos.full.I)
Anova(fos.full.I)
cAIC(fos.full.I)
AIC(fos.full.I)
r.squaredGLMM(fos.full.I)

qqnorm(residuals(fos.full.I))
scatter.smooth(residuals(fos.full.I) ~ fitted(fos.full.I))

# Model 1 - reduction ####

fos_mod = buildmer(fos ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                     Peak_T)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(fos_mod)
print(fos_mod, correlation=TRUE)
vcov(fos_mod)

fos_mod.1 = buildmer(fos ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.1,
                     buildmerControl = buildmerControl(include = ~ (1|Site),
                                                       calc.anova = TRUE,direction='forward',
                                                       ddf = "Satterthwaite"))
summary(fos_mod.1)
print(fos_mod.1, correlation=TRUE)
vcov(fos_mod.1)

# Model 2 - WINNER MODEL ####

fos.full.2 = lmer(fos ~ 1 + Water_content + SR + PO43 + Clay + BB + (1|Site), 
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
dominance_output <- domin(fos ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
                          sets = list("Water_content","SR","PO43","Clay","BB"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 ####

fos.full.3 = lmer(fos ~ 1 + Water_content + SR + altitude + Clay + BB + AI + C_N +
                    (1|Site), data = my_data.1)                                         

summary(fos.full.3)
Anova(fos.full.3)
cAIC(fos.full.3)
AIC(fos.full.3)
r.squaredGLMM(fos.full.3)

qqnorm(residuals(fos.full.3))
scatter.smooth(residuals(fos.full.3) ~ fitted(fos.full.3))

# Enzyme - leu ####

# Model 1 ####

leu.full.I = lmer(leu ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TOC+TC+TN+C_N+NH4+PO43+SO42+
                                    Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                    BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(leu.full.I)
Anova(leu.full.I)
cAIC(leu.full.I)
AIC(leu.full.I)
r.squaredGLMM(leu.full.I)

qqnorm(residuals(leu.full.I))
scatter.smooth(residuals(leu.full.I) ~ fitted(leu.full.I))

# Model 1 - reduction ####

leu_mod = buildmer(leu ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                     Peak_T)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(leu_mod)
print(leu_mod, correlation=TRUE)
vcov(leu_mod)

leu_mod.1 = buildmer(leu ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.1,
                     buildmerControl = buildmerControl(include = ~ (1|Site),
                                                       calc.anova = TRUE,direction='forward',
                                                       ddf = "Satterthwaite"))
summary(leu_mod.1)
print(leu_mod.1, correlation=TRUE)
vcov(leu_mod.1)

# Model 2 ####

leu.full.2 = lmer(leu ~ 1 + FB + SOM + Peak_A + (1|Site), data = my_data.1)                                         

summary(leu.full.2)
Anova(leu.full.2)
cAIC(leu.full.2)
AIC(leu.full.2)
r.squaredGLMM(leu.full.2)

qqnorm(residuals(leu.full.2))
scatter.smooth(residuals(leu.full.2) ~ fitted(leu.full.2))

# Model 3 ####

leu.full.3 = lmer(leu ~ 1 + FB + Peak_T + SOM + L_TC + (1|Site), data = my_data.1)                                         

summary(leu.full.3)
Anova(leu.full.3)
cAIC(leu.full.3)
AIC(leu.full.3)
r.squaredGLMM(leu.full.3)

qqnorm(residuals(leu.full.3))
scatter.smooth(residuals(leu.full.3) ~ fitted(leu.full.3))

# Model 1 - without outliers ####

plot(cooks.distance(leu.full.3))
my_data.5 = my_data.1[-c(30,43),]

leu.full.Ia = lmer(leu ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                     Peak_T)*AI + (1|Site), data = my_data.5)                                         

summary(leu.full.Ia)
Anova(leu.full.Ia)
cAIC(leu.full.Ia)
AIC(leu.full.Ia)
r.squaredGLMM(leu.full.Ia)

qqnorm(residuals(leu.full.Ia))
scatter.smooth(residuals(leu.full.Ia) ~ fitted(leu.full.Ia))

# Model 1 - without outliers - reduction ####

leu_mod.a = buildmer(leu ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.5,
                     buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                       ddf = "Satterthwaite"))
summary(leu_mod.a)
print(leu_mod.a, correlation=TRUE)
vcov(leu_mod.a)

leu_mod.1.a = buildmer(leu ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                         TOC+TC+TN+C_N+NH4+PO43+SO42+
                                         Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                         BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                         Peak_T)*AI + (1|Site), data = my_data.5,
                       buildmerControl = buildmerControl(include = ~ (1|Site),
                                                         calc.anova = TRUE,direction='forward',
                                                         ddf = "Satterthwaite"))
summary(leu_mod.1.a)
print(leu_mod.1.a, correlation=TRUE)
vcov(leu_mod.1.a)

# Model 2 - without outliers ####

leu.full.2a = lmer(leu ~ 1 + Clay + E2.E3 + SOM + TN + LTC_LTN + TC + L_TC + L_TN +  
                     TOC + AI + TC:AI + E2.E3:AI + SO42 + AI:SO42 + L_TC:AI + 
                     LTC_LTN:AI + Clay:AI + (1|Site), data = my_data.5)                                         

summary(leu.full.2a)
Anova(leu.full.2a)
cAIC(leu.full.2a)
AIC(leu.full.2a)
r.squaredGLMM(leu.full.2a)

qqnorm(residuals(leu.full.2a))
scatter.smooth(residuals(leu.full.2a) ~ fitted(leu.full.2a))

# Model 3 - without outliers - WINNER MODEL ####

leu.full.3a = lmer(leu ~ 1 + Clay + E2.E3 + Silt + SOM + TN + LTC_LTN + NH4 + pH +  
                     TC + L_TC + L_TN + TOC + AI + SO42 + SO42:AI + (1|Site), data = my_data.5)                                         

summary(leu.full.3a)
Anova(leu.full.3a)
cAIC(leu.full.3a)
AIC(leu.full.3a)
r.squaredGLMM(leu.full.3a)

qqnorm(residuals(leu.full.3a))
scatter.smooth(residuals(leu.full.3a) ~ fitted(leu.full.3a))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(leu ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.5, 
                          sets = list("Clay","E2.E3","Silt","SOM","TN","LTC_LTN",
                                      "NH4","pH","TC","L_TC","L_TN","TOC","AI",
                                      "SO42","SO42:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Enzyme - phe ####

# Model 1 ####

phe.full.I = lmer(phe ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TOC+TC+TN+C_N+NH4+PO43+SO42+
                                    Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                    BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(phe.full.I)
Anova(phe.full.I)
cAIC(phe.full.I)
AIC(phe.full.I)
r.squaredGLMM(phe.full.I)

qqnorm(residuals(phe.full.I))
scatter.smooth(residuals(phe.full.I) ~ fitted(phe.full.I))

# Model 1 - reduction ####

phe_mod = buildmer(phe ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                     Peak_T)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(phe_mod)
print(phe_mod, correlation=TRUE)
vcov(phe_mod)

phe_mod.1 = buildmer(phe ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.1,
                     buildmerControl = buildmerControl(include = ~ (1|Site),
                                                       calc.anova = TRUE,direction='forward',
                                                       ddf = "Satterthwaite"))
summary(phe_mod.1)
print(phe_mod.1, correlation=TRUE)
vcov(phe_mod.1)

# Model 2 - WINNER MODEL ####

phe.full.2 = lmer(phe ~ 1 + L_TC + TN + SOM + BB + Peak_A + Water_content + C_N + 
                    (1|Site), data = my_data.1)                                         

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
                          sets = list("L_TC","TN","SOM","BB","Peak_A","Water_content","C_N"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 ####

phe.full.3 = lmer(phe ~ 1 + L_TC + TN + SOM + C_N + (1|Site), data = my_data.1)                                         

summary(phe.full.3)
Anova(phe.full.3)
cAIC(phe.full.3)
AIC(phe.full.3)
r.squaredGLMM(phe.full.3)

qqnorm(residuals(phe.full.3))
scatter.smooth(residuals(phe.full.3) ~ fitted(phe.full.3))

# Enzyme - xylcbh ####

# Model 1 ####

xylcbh.full.I = lmer(xylcbh ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                    TOC+TC+TN+C_N+NH4+PO43+SO42+
                                    Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                    BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                    Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(xylcbh.full.I)
Anova(xylcbh.full.I)
cAIC(xylcbh.full.I)
AIC(xylcbh.full.I)
r.squaredGLMM(xylcbh.full.I)

qqnorm(residuals(xylcbh.full.I))
scatter.smooth(residuals(xylcbh.full.I) ~ fitted(xylcbh.full.I))

# Model 1 - reduction ####

xylcbh_mod = buildmer(xylcbh ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                     Peak_T)*AI + (1|Site), data = my_data.1,
                   buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                     ddf = "Satterthwaite"))
summary(xylcbh_mod)
print(xylcbh_mod, correlation=TRUE)
vcov(xylcbh_mod)

xylcbh_mod.1 = buildmer(xylcbh ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.1,
                     buildmerControl = buildmerControl(include = ~ (1|Site),
                                                       calc.anova = TRUE,direction='forward',
                                                       ddf = "Satterthwaite"))
summary(xylcbh_mod.1)
print(xylcbh_mod.1, correlation=TRUE)
vcov(xylcbh_mod.1)

# Model 2 - WINNER MODEL ####

xylcbh.full.2 = lmer(xylcbh ~ 1 + L_TN + Peak_A + AI + L_TN:AI + Water_content + 
                       C_N + (1|Site), data = my_data.1)                                         

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
                          sets = list("L_TN","Peak_A","AI","L_TN:AI","Water_content",
                                      "C_N"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 ####

xylcbh.full.3 = lmer(xylcbh ~ 1 + L_TN + Peak_A + SR + AI + Water_content + Soil_Temp + (1|Site), data = my_data.1)                                         

summary(xylcbh.full.3)
Anova(xylcbh.full.3)
cAIC(xylcbh.full.3)
AIC(xylcbh.full.3)
r.squaredGLMM(xylcbh.full.3)

qqnorm(residuals(xylcbh.full.3))
scatter.smooth(residuals(xylcbh.full.3) ~ fitted(xylcbh.full.3))

# Enzyme - alphabeta ####

# Model 1 ####

alphabeta.full.I = lmer(alphabeta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                          TOC+TC+TN+C_N+NH4+PO43+SO42+
                                          Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                          BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                          Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(alphabeta.full.I)
Anova(alphabeta.full.I)
cAIC(alphabeta.full.I)
AIC(alphabeta.full.I)
r.squaredGLMM(alphabeta.full.I)

qqnorm(residuals(alphabeta.full.I))
scatter.smooth(residuals(alphabeta.full.I) ~ fitted(alphabeta.full.I))

# Model 1 - reduction ####

alphabeta_mod = buildmer(alphabeta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                           TOC+TC+TN+C_N+NH4+PO43+SO42+
                                           Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                           BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                           Peak_T)*AI + (1|Site), data = my_data.1,
                      buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                        ddf = "Satterthwaite"))
summary(alphabeta_mod)
print(alphabeta_mod, correlation=TRUE)
vcov(alphabeta_mod)

alphabeta_mod.1 = buildmer(alphabeta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                             TOC+TC+TN+C_N+NH4+PO43+SO42+
                                             Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                             BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                             Peak_T)*AI + (1|Site), data = my_data.1,
                        buildmerControl = buildmerControl(include = ~ (1|Site),
                                                          calc.anova = TRUE,direction='forward',
                                                          ddf = "Satterthwaite"))
summary(alphabeta_mod.1)
print(alphabeta_mod.1, correlation=TRUE)
vcov(alphabeta_mod.1)

# Model 2 ####

alphabeta.full.2 = lmer(alphabeta ~ 1 + HIX + Soil_Temp + (1|Site), data = my_data.1)                                         

summary(alphabeta.full.2)
Anova(alphabeta.full.2)
cAIC(alphabeta.full.2)
AIC(alphabeta.full.2)
r.squaredGLMM(alphabeta.full.2)

qqnorm(residuals(alphabeta.full.2))
scatter.smooth(residuals(alphabeta.full.2) ~ fitted(alphabeta.full.2))

# Model 3 ####

alphabeta.full.3 = lmer(alphabeta ~ 1 + Silt + HIX + Soil_Temp + FI + AI + 
                          Clay + Litter + (1|Site), data = my_data.1)                                         

summary(alphabeta.full.3)
Anova(alphabeta.full.3)
cAIC(alphabeta.full.3)
AIC(alphabeta.full.3)
r.squaredGLMM(alphabeta.full.3)

qqnorm(residuals(alphabeta.full.3))
scatter.smooth(residuals(alphabeta.full.3) ~ fitted(alphabeta.full.3))

# Model 1 - without outliers ####

plot(cooks.distance(alphabeta.full.I))
my_data.6 = my_data.1[-c(59),]

alphabeta.full.Ia = lmer(alphabeta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                     TOC+TC+TN+C_N+NH4+PO43+SO42+
                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                     BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                     Peak_T)*AI + (1|Site), data = my_data.6)                                         

summary(alphabeta.full.Ia)
Anova(alphabeta.full.Ia)
cAIC(alphabeta.full.Ia)
AIC(alphabeta.full.Ia)
r.squaredGLMM(alphabeta.full.Ia)

qqnorm(residuals(alphabeta.full.Ia))
scatter.smooth(residuals(alphabeta.full.Ia) ~ fitted(alphabeta.full.Ia))

# Model 1 - without outliers - reduction ####

alphabeta_mod.a = buildmer(alphabeta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                       TOC+TC+TN+C_N+NH4+PO43+SO42+
                                       Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                       BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                       Peak_T)*AI + (1|Site), data = my_data.6,
                     buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                       ddf = "Satterthwaite"))
summary(alphabeta_mod.a)
print(alphabeta_mod.a, correlation=TRUE)
vcov(alphabeta_mod.a)

alphabeta_mod.1.a = buildmer(alphabeta ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                         TOC+TC+TN+C_N+NH4+PO43+SO42+
                                         Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                         BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                         Peak_T)*AI + (1|Site), data = my_data.6,
                       buildmerControl = buildmerControl(include = ~ (1|Site),
                                                         calc.anova = TRUE,direction='forward',
                                                         ddf = "Satterthwaite"))
summary(alphabeta_mod.1.a)
print(alphabeta_mod.1.a, correlation=TRUE)
vcov(alphabeta_mod.1.a)

# Model 2 - without outliers - WINNER MODEL ####

alphabeta.full.2a = lmer(alphabeta ~ 1 + FI + Soil_Temp + TOC + TC + Peak_A + AI + FB +  
                           Clay + L_TN + AI:L_TN + Litter + Litter:AI + (1|Site), data = my_data.6)                                         

summary(alphabeta.full.2a)
Anova(alphabeta.full.2a)
cAIC(alphabeta.full.2a)
AIC(alphabeta.full.2a)
r.squaredGLMM(alphabeta.full.2a)

qqnorm(residuals(alphabeta.full.2a))
scatter.smooth(residuals(alphabeta.full.2a) ~ fitted(alphabeta.full.2a))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(alphabeta ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.6, 
                          sets = list("FI","Soil_Temp","TOC","TC","Peak_A","AI",
                                      "FB","Clay","L_TN","AI:L_TN","Litter",
                                      "Litter:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 - without outliers ####

alphabeta.full.3a = lmer(alphabeta ~ 1 + FI + Soil_Temp + C_N + TOC + TC + Peak_A + TN +  
                           NH4 + AI + SOM + Water_content + Clay + L_TN + SO42 + 
                           Litter + L_TC + LTC_LTN + (1|Site), data = my_data.6)                                         

summary(alphabeta.full.3a)
Anova(alphabeta.full.3a)
cAIC(alphabeta.full.3a)
AIC(alphabeta.full.3a)
r.squaredGLMM(alphabeta.full.3a)

qqnorm(residuals(alphabeta.full.3a))
scatter.smooth(residuals(alphabeta.full.3a) ~ fitted(alphabeta.full.3a))

# Enzyme - Cenz ####

# Model 1 ####

Cenz.full.I = lmer(Cenz ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                TOC+TC+TN+C_N+NH4+PO43+SO42+
                                                Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                                BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                                Peak_T)*AI + (1|Site), data = my_data.1)                                         

summary(Cenz.full.I)
Anova(Cenz.full.I)
cAIC(Cenz.full.I)
AIC(Cenz.full.I)
r.squaredGLMM(Cenz.full.I)

qqnorm(residuals(Cenz.full.I))
scatter.smooth(residuals(Cenz.full.I) ~ fitted(Cenz.full.I))

# Model 1 - reduction ####

Cenz_mod = buildmer(Cenz ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                 TOC+TC+TN+C_N+NH4+PO43+SO42+
                                                 Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                                 BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                                 Peak_T)*AI + (1|Site), data = my_data.1,
                         buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                           ddf = "Satterthwaite"))
summary(Cenz_mod)
print(Cenz_mod, correlation=TRUE)
vcov(Cenz_mod)

Cenz_mod.1 = buildmer(Cenz ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                   TOC+TC+TN+C_N+NH4+PO43+SO42+
                                                   Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                                   BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                                   Peak_T)*AI + (1|Site), data = my_data.1,
                           buildmerControl = buildmerControl(include = ~ (1|Site),
                                                             calc.anova = TRUE,direction='forward',
                                                             ddf = "Satterthwaite"))
summary(Cenz_mod.1)
print(Cenz_mod.1, correlation=TRUE)
vcov(Cenz_mod.1)

# Model 2 ####

Cenz.full.2 = lmer(Cenz ~ 1 + HIX + Soil_Temp + TN + C_N + (1|Site), data = my_data.1)                                         

summary(Cenz.full.2)
Anova(Cenz.full.2)
cAIC(Cenz.full.2)
AIC(Cenz.full.2)
r.squaredGLMM(Cenz.full.2)

qqnorm(residuals(Cenz.full.2))
scatter.smooth(residuals(Cenz.full.2) ~ fitted(Cenz.full.2))

# Model 3 ####

Cenz.full.3 = lmer(Cenz ~ 1 + HIX + Soil_Temp + TN + C_N + BB + FI + AI + Clay + Litter + (1|Site), data = my_data.1)                                         

summary(Cenz.full.3)
Anova(Cenz.full.3)
cAIC(Cenz.full.3)
AIC(Cenz.full.3)
r.squaredGLMM(Cenz.full.3)

qqnorm(residuals(Cenz.full.3))
scatter.smooth(residuals(Cenz.full.3) ~ fitted(Cenz.full.3))

# Model 1 - without outliers ####

plot(cooks.distance(Cenz.full.I))
my_data.8 = my_data.1[-c(59),]

Cenz.full.Ia = lmer(Cenz ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                 TOC+TC+TN+C_N+NH4+PO43+SO42+
                                                 Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                                 BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                                 Peak_T)*AI + (1|Site), data = my_data.8)                                         

summary(Cenz.full.Ia)
Anova(Cenz.full.Ia)
cAIC(Cenz.full.Ia)
AIC(Cenz.full.Ia)
r.squaredGLMM(Cenz.full.Ia)

qqnorm(residuals(Cenz.full.Ia))
scatter.smooth(residuals(Cenz.full.Ia) ~ fitted(Cenz.full.Ia))

# Model 1 - without outliers - reduction ####

Cenz_mod.a = buildmer(Cenz ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                   TOC+TC+TN+C_N+NH4+PO43+SO42+
                                                   Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                                   BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                                   Peak_T)*AI + (1|Site), data = my_data.8,
                           buildmerControl = buildmerControl(include = ~ (1|Site),calc.anova = TRUE,
                                                             ddf = "Satterthwaite"))
summary(Cenz_mod.a)
print(Cenz_mod.a, correlation=TRUE)
vcov(Cenz_mod.a)

Cenz_mod.1.a = buildmer(Cenz ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                                     TOC+TC+TN+C_N+NH4+PO43+SO42+
                                                     Silt+Clay+Litter+LTC_LTN+L_TN+L_TC+
                                                     BB+FB+SR+E2.E3+FI+HIX+Peak_A+
                                                     Peak_T)*AI + (1|Site), data = my_data.8,
                             buildmerControl = buildmerControl(include = ~ (1|Site),
                                                               calc.anova = TRUE,direction='forward',
                                                               ddf = "Satterthwaite"))
summary(Cenz_mod.1.a)
print(Cenz_mod.1.a, correlation=TRUE)
vcov(Cenz_mod.1.a)

# Model 2 - without outliers - WINNER MODEL ####

Cenz.full.2a = lmer(Cenz ~ 1 + FI + L_TN + Peak_A + C_N + Soil_Temp + TOC + FB +  
                      BB + SOM + Clay + AI + Peak_A:AI + BB:AI + L_TN:AI + 
                      FI:AI + (1|Site), data = my_data.8)                                         

summary(Cenz.full.2a)
Anova(Cenz.full.2a)
cAIC(Cenz.full.2a)
AIC(Cenz.full.2a)
r.squaredGLMM(Cenz.full.2a)

qqnorm(residuals(Cenz.full.Ia))
scatter.smooth(residuals(Cenz.full.Ia) ~ fitted(Cenz.full.Ia))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(Cenz ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.8, 
                          sets = list("FI","L_TN","Peak_A","C_N","Soil_Temp",
                                      "TOC","FB","BB","SOM","Clay","AI","Peak_A:AI",
                                      "BB:AI","L_TN:AI","FI:AI"), 
                          consmodel = "(1|Site)") # Replace with your actual function

# Model 3 - without outliers ####

Cenz.full.3a = lmer(Cenz ~ 1 + FI + L_TN + Peak_A + C_N + Soil_Temp + TN + TC + SO42 +  
                      BB + LTC_LTN + Clay + E2.E3 + L_TC + AI + LTC_LTN:AI + 
                      Peak_A:AI + Litter + (1|Site), data = my_data.8)                                         

summary(Cenz.full.3a)
Anova(Cenz.full.3a)
cAIC(Cenz.full.3a)
AIC(Cenz.full.3a)
r.squaredGLMM(Cenz.full.3a)

qqnorm(residuals(Cenz.full.3a))
scatter.smooth(residuals(Cenz.full.3a) ~ fitted(Cenz.full.3a))
