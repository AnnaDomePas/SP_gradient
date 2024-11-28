gam_mod = gam(alpha ~ s(FB,AI) + s(as.factor(my_data.1$Site), bs = "re"),data = my_data.1, method = "REML")
gam.check(gam_mod)
summary(gam_mod)
plot(resid(gam_mod)~fitted(gam_mod))

gam_mod.1 = gam(alpha ~ s(FB) + s(AI) + s(as.factor(my_data.1$Site), bs = "re"),data = my_data.1, method = "REML")
gam.check(gam_mod.1)
summary(gam_mod.1)
plot(resid(gam_mod.1)~fitted(gam_mod.1))
k.check(gam_mod.1)

gam_mod.2 = gam(alpha ~ s(BB) + s(AI) + s(L_TC) + s(as.factor(my_data.1$Site), bs = "re"),data = my_data.1, method = "REML")
gam.check(gam_mod.2)
summary(gam_mod.2)
plot(resid(gam_mod.2)~fitted(gam_mod.2))
k.check(gam_mod.2)

gam_mod.3 = gam(alpha ~ s(AI) + s(as.factor(my_data.1$Site), bs = "re"),data = my_data.1, method = "REML")
gam.check(gam_mod.3)
summary(gam_mod.3)
plot(resid(gam_mod.3)~fitted(gam_mod.3))
k.check(gam_mod.3)

