# 14 March 2023

rm(list=ls(all=TRUE))

# Packages----
library(here)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(piecewiseSEM)
library(multcomp)


sessionInfo()

# R version 4.2.2 (2022-10-31 ucrt)

# [1] emmeans_1.8.4-1    piecewiseSEM_2.3.0 MASS_7.3-58.1      r2glmm_0.1.2       lmerTest_3.1-3     lme4_1.1-31       
# [7] Matrix_1.5-1       lubridate_1.9.2    forcats_1.0.0      stringr_1.5.0      dplyr_1.1.0        purrr_1.0.1       
# [13] readr_2.1.4        tidyr_1.3.0        tibble_3.1.8       ggplot2_3.4.1      tidyverse_2.0.0   


# local working directory
path <- here::here()
path

# Data----
library(tidyverse)
k.dat<-read_csv ("Data/Aphid_Data.csv")
str(k.dat) 
names (k.dat)

# check missing data
which(is.na(k.dat))

k.dat$Pot <- factor(k.dat$Pot)


# Check the relationships
library(ggplot2)
ggplot(k.dat, aes(Plant_weight,Aphid_Number)) + geom_point()
ggplot(k.dat, aes(Predtr, Aphid_Number)) + geom_boxplot()
ggplot(k.dat, aes(Hst_div, Aphid_Number)) + geom_boxplot()
ggplot(k.dat, aes(Garlic_weight, Aphid_Number)) + geom_point()
ggplot(k.dat, aes(Garlic, Aphid_Number)) + geom_boxplot()


# Mixed models----

library(lme4)
library(lmerTest)

## Plant weight------

mm1 <- lmer(Plant_weight  ~ Hst_div + Garlic_weight +
              (1|Pot), data = k.dat)
plot(mm1)
qqnorm(resid(mm1))
qqline(resid(mm1))

# check multicolinearity
car::vif(mm1)

# estimates
summary(mm1)
car::Anova(mm1)

# Marginal means and pairwise differences of Hst_div levels
emmeans::emmeans(mm1, list(pairwise ~ Hst_div))
# to add letters for the post-hoc test:
library(multcomp)
library(emmeans)
cld(emmeans(mm1, list(pairwise ~ Hst_div)),  
    #  type="response",
    Letters = letters, adjust = "none")

# check in psem
summary(piecewiseSEM::psem(mm1))

#  model R2
MuMIn::r.squaredGLMM(mm1)
# R2m and R2c are marginal (for fixed predictors) and conditional (for fixed and random predictors) coefficients of determination

# Partial R2 for fixed effects
library(r2glmm)
R2part_mm1 <- r2beta(mm1, method = 'nsj', partial = T, data = k.dat)
R2part_mm1

## Aphid density-----
mm2 <- glmer(Aphid_Number ~ Hst_div + Plant_weight + Predtr  + 
                    Garlic_weight + (1|Pot), data = k.dat, 
             family = "poisson") 

# check model
summary(mm2)

# check multicolinearity
car::vif(mm2)

plot(mm2)
qqnorm(resid(mm2))
qqline(resid(mm2))

# check overdispersion

resid_pearson <- residuals(mm2, type = "pearson")
SSQ <- sum(resid_pearson^2)
SSQ/df.residual(mm2) 
# [1] 3.404141 -> clear overdispersion

# change family
library(MASS)

mm2b <- glmmPQL(Aphid_Number ~ Hst_div + Plant_weight + Predtr  + 
                    Garlic_weight, 
                  random = ~ 1 | Pot,  data = k.dat,
                family = "quasipoisson") 
  
plot(mm2b)
qqnorm(resid(mm2b))
qqline(resid(mm2b))

# check multicolinearity
car::vif(mm2b)

# estimates
summary(mm2b)
car::Anova(mm2b)

# Marginal means and pairwise differences of Hst_div levels and of Predtr levels

emmeans::emmeans(mm2b, list(pairwise ~ Hst_div, 
                            pairwise ~ Predtr))
# to add letters for post-hoc test:
multcomp::cld(emmeans::emmeans(mm2b, list(pairwise ~ Hst_div, 
                       pairwise ~ Predtr)),  
    #  type="response",
    Letters = letters, adjust = "none")

# check in psem
summary(piecewiseSEM::psem(mm2b), standardize = "scale", standardize.type = "Menard.OE")


# R2 for the entire model
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(mm2b)

# Partial R2 for fixed effects
# Only the SGV method is compatible with glmmPQL object
R2part <- r2glmm::r2beta(mm2b, method = 'SGV', partial = T, data = k.dat)

# plot the partial R2

R2<- R2part %>% 
  filter(Effect!="Model") %>% 
  mutate(Effect=fct_recode(Effect, "Predator"="Predtrpresent", 
                      "Host mass" = "Plant_weight", 
                      "Garlic mass"  ="Garlic_weight",
                      "Host intercropping" ="Hst_divbarley intercropped")) %>% 
  mutate(Effect=fct_relevel(Effect, "Garlic mass", "Host mass",
                            "Host intercropping","Predator"))
                           
R2

library(ggplot2)

col <- c("#57B1C9","#81C5AD","#DB9B99" , "#FFFF4B" )
fill <- c("#B7DEE8","#AEDACA","#E6B9B8","#FFFFB9" )

ggplot(R2, aes(y=Rsq, x=Effect)) + 
geom_bar(position="stack", stat="identity", colour = col, fill=fill)+
  coord_flip() +
  #  ylab(bquote('Explained variance (partial R'^'2'*')'))+
  ylab(bquote('partial R'^'2'))+
  labs(x =element_blank()) +
   scale_x_discrete(labels= c( "Garlic mass",
                               "Host mass",
                               "Host 
                               intercropping",
                               "Predator
                               presence"))+
  theme(axis.text.y=element_text(colour = "black", size=15),
        axis.text.x=element_text(colour = "black", size=13),
        axis.title=element_text(size=15),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12) ,
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks =  element_line(colour = "black"),
        axis.ticks.y = element_blank())

###############################-

# Structural equation model:----
library(piecewiseSEM)
psem_model <- psem (mm1, mm2b)

# Testing missing links by D-separation test

dSep(psem_model, .progressBar = FALSE)


# Model fit:

# Fisher’s C statistic
fisherC(psem_model)
# χ2 statistic
# LLchisq(psem_model) # χ2 statistic can not be calculated, we rely on Fisher’s C

# sample size
nrow(k.dat)

# Estimates:
unstdCoefs(psem_model)
# model R2
rsquared(psem_model)

# summary( psem_model,  standardize = "none", .progressBar = FALSE)

# standardized coefficients are not returned
# Family (or link) function not supported, e.g. for quasipoisson 

coefs(psem_model, standardize = "scale", standardize.type = "latent.linear") # default
coefs(psem_model, standardize = "scale", standardize.type = "Menard.OE")
coefs(psem_model, standardize = "range")

#relying on unstandardized coefficients:
coefs(psem_model,  standardize = "none")

path
write_csv(unstdCoefs(psem_model), file="coefs.csv")



#################################################-

# Supplementary ----
## Fig S1.A----
# Effects of plant species on plant mass

mS1 <-lmer(Plant_weight~Plant +(1|Pot), data = k.dat) 

k.dat$Plant_weight
plot(mS1)
qqnorm(resid(mS1))
qqline(resid(mS1))

summary(mS1)
car::Anova(mS1)

# Marginal means and pairwise differences of Plant species
emmeans::emmeans(mS1, list(pairwise ~ Plant))

ggplot(k.dat, aes(y=Plant_weight, x=Plant)) + 
  geom_boxplot() +
  #  ylab(bquote('Explained variance (partial R'^'2'*')'))+
  #ylab(bquote('Plant mass, kg'^'2'))+
  labs(x ="Plant species", y="Plant mass, g") +
  theme_bw()+
  theme(axis.text.y=element_text(colour = "black", size=13),
        axis.text.x=element_text(colour = "black", size=13),
        axis.title=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.ticks =  element_line(colour = "black"))



## Fig S1.B----
# Effects of plant species on aphid density

mS2 <- glmer(Aphid_Number ~ Plant +   Predtr  + 
               Garlic_weight + (1|Pot), data = k.dat, 
             family = "poisson") 
# check model
summary(mS2)


# check multicolinearity
car::vif(mS2)

plot(mS2)
qqnorm(resid(mS2))
qqline(resid(mS2))

# check overdispersion

resid_pearson <- residuals(mS2, type = "pearson")
SSQ <- sum(resid_pearson^2)
SSQ/df.residual(mS2) 
# clear overdispersion

# change family
library(MASS)

mS2b <- glmmPQL(Aphid_Number ~ Plant +   Predtr  + 
                  Garlic_weight  , 
                random = ~ 1 | Pot,  data = k.dat,
                family = "quasipoisson") 

plot(mS2b)
qqnorm(resid(mS2b))
qqline(resid(mS2b))

# check multicolinearity
car::vif(mS2b)

# estimates
summary(mS2b)
car::Anova(mS2b)

# Marginal means and pairwise differences of Plant species
emmeans::emmeans(mS2b, list(pairwise ~ Plant))


ggplot(k.dat, aes(y=Aphid_Number, x=Plant)) + 
  geom_boxplot() +
   labs(x ="Plant species", y="Aphid density") +
  theme_bw()+
  theme(axis.text.y=element_text(colour = "black", size=13),
        axis.text.x=element_text(colour = "black", size=13),
        axis.title=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.ticks =  element_line(colour = "black"))


#End