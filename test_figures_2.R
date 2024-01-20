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

# Model 2 ####

respiration.full.II = lmer(Respiration ~ altitude+(Soil_Temp+Water_content+SOM+TOC+TC+TN+
                                            C_N+PO43+SR+Clay+L_TN+FB+FI+Silt+L_TC+LTC_LTN)*AI + (1|Site),
                           data = my_data.1)                                         

summary(respiration.full.II)
Anova(respiration.full.II)
cAIC(respiration.full.II)
AIC(respiration.full.II)
r.squaredGLMM(respiration.full.II)

qqnorm(residuals(respiration.full.II))
scatter.smooth(residuals(respiration.full.II) ~ fitted(respiration.full.II))

# Model 3 ####

respiration.full.III = lmer(Respiration ~ altitude + (Soil_Temp+Water_content+SR+Clay+
                                                        Silt+FI)*AI + (1|Site),data = my_data.1)                                         

summary(respiration.full.III)
Anova(respiration.full.III)
cAIC(respiration.full.III)
AIC(respiration.full.III)
r.squaredGLMM(respiration.full.III)

qqnorm(residuals(respiration.full.III))
scatter.smooth(residuals(respiration.full.III) ~ fitted(respiration.full.III))

# Model 4 ####

respiration.full.4 = lmer(Respiration ~ altitude + (Water_content+SR+Clay+Silt)*AI + 
                            (1|Site),data = my_data.1)                                         

summary(respiration.full.4)
Anova(respiration.full.4)
cAIC(respiration.full.4)
AIC(respiration.full.4)
r.squaredGLMM(respiration.full.4)

qqnorm(residuals(respiration.full.4))
scatter.smooth(residuals(respiration.full.4) ~ fitted(respiration.full.4))

# Model 5 ####

respiration.full.5 = lmer(Respiration ~ altitude + Water_content + Clay + Silt + 
                            SR:AI + Clay:AI + Silt:AI + 
                            (1|Site),data = my_data.1)                                         

summary(respiration.full.5)
Anova(respiration.full.5)
cAIC(respiration.full.5)
AIC(respiration.full.5)
r.squaredGLMM(respiration.full.5)

qqnorm(residuals(respiration.full.5))
scatter.smooth(residuals(respiration.full.5) ~ fitted(respiration.full.5))

# Importance assessment ####

domin(Respiration ~ 1, 
      lmer, 
      list(\(x) list(R2m = MuMIn::r.squaredGLMM(x)[[1]]), "R2m"), 
      data = my_data.1, 
      sets = list("altitude","Water_content","Silt","Silt:AI","Clay","Clay:AI","SR:AI"), 
      consmodel = "(1|Site)")

# Enzymes ####

# alpha ####

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

# Model 2 ####

alpha.full.2 = lmer(alpha ~ altitude+(Soil_Temp+Water_content+SOM+pH+
                                        TOC+TC+TN+C_N+PO43+
                                        Silt+Clay+L_TN+L_TC+LTC_LTN+
                                        FB+SR+E2.E3+FI+HIX+Peak_A)*AI + 
                      (1|Site), data = my_data.1)                                         

summary(alpha.full.2)
Anova(alpha.full.2)
cAIC(alpha.full.2)
AIC(alpha.full.2)
r.squaredGLMM(alpha.full.2)

qqnorm(residuals(alpha.full.2))
scatter.smooth(residuals(alpha.full.2) ~ fitted(alpha.full.2))

# Model 3 ####

alpha.full.3 = lmer(alpha ~ (Soil_Temp+SOM+TN+C_N+L_TN+L_TC+LTC_LTN+
                                        FB+SR+HIX)*AI + 
                      (1|Site), data = my_data.1)                                         

summary(alpha.full.3)
Anova(alpha.full.3)
cAIC(alpha.full.3)
AIC(alpha.full.3)
r.squaredGLMM(alpha.full.3)

qqnorm(residuals(alpha.full.3))
scatter.smooth(residuals(alpha.full.3) ~ fitted(alpha.full.3))

# Model 4 ####

alpha.full.4 = lmer(alpha ~ Soil_Temp + TN + HIX + SOM:AI + TN:AI + C_N:AI + FB + 
                      (1|Site), data = my_data.1)                                         

summary(alpha.full.4)
Anova(alpha.full.4)
cAIC(alpha.full.4)
AIC(alpha.full.4)
r.squaredGLMM(alpha.full.4)

qqnorm(residuals(alpha.full.4))
scatter.smooth(residuals(alpha.full.4) ~ fitted(alpha.full.4))

