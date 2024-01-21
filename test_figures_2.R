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
                          data = my_data.1, 
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
dominance_output <- domin(alpha ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
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

# Model 3 - MODEL WINNER ####

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
dominance_output <- domin(alpha ~ 1, 
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
dominance_output <- domin(alpha ~ 1, 
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
scatter.smooth(residuals(gla.full.2a) ~ fitted(gla.full.2a))

# Parameter ranking plot ----
#Change the signs of the variables
dominance_output <- domin(alpha ~ 1, 
                          lmer, 
                          list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
                          data = my_data.1, 
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

