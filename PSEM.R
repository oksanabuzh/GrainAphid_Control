# Analysis

rm(list=ls(all=TRUE))

# Packages----
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MASS)
library(car)
library(emmeans)
library(multcomp)
library(r2glmm)
library(MuMIn)
library(piecewiseSEM)
library(sjPlot)
library(ggeffects)


sessionInfo()

# R version 4.2.2 (2022-10-31 ucrt)

# [1] emmeans_1.8.4-1    piecewiseSEM_2.3.0 MASS_7.3-58.1      r2glmm_0.1.2       lmerTest_3.1-3     lme4_1.1-31       
# [7] Matrix_1.5-1       lubridate_1.9.2    forcats_1.0.0      stringr_1.5.0      dplyr_1.1.0        purrr_1.0.1       
# [13] readr_2.1.4        tidyr_1.3.0        tibble_3.1.8       ggplot2_3.4.1      tidyverse_2.0.0   


set_theme(base = theme_bw(), axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4)


# Data----
#  library(tidyverse)
k.dat<-read_csv ("Data/Aphid_Data.csv")
str(k.dat) 
names (k.dat)
# check missing data
which(is.na(k.dat))

k.dat$Pot <- factor(k.dat$Pot)



# Check the relationships
# library(ggplot2)
ggplot(k.dat, aes(Plant_weight,Aphid_Number)) + geom_point()
ggplot(k.dat, aes(Predtr, Aphid_Number)) + geom_boxplot()
ggplot(k.dat, aes(Hst_div, Aphid_Number)) + geom_boxplot()
ggplot(k.dat, aes(Garlic_weight, Aphid_Number)) + geom_point()
ggplot(k.dat, aes(Garlic, Aphid_Number)) + geom_boxplot()


# Mixed models----

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
# library(emmeans)
# library(multcomp)
emmeans::emmeans(mm1, list(pairwise ~ Hst_div))
# to add letters for the post-hoc test:
cld(emmeans(mm1, list(pairwise ~ Hst_div)),  
    #  type="response",
    Letters = letters, adjust = "none")

# check in psem
summary(piecewiseSEM::psem(mm1))

#  model R2
MuMIn::r.squaredGLMM(mm1)
# R2m and R2c are marginal (for fixed predictors) and conditional (for fixed and random predictors) coefficients of determination

# Partial R2 for fixed effects
# library(r2glmm)

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

### Distribution of the Aphid_Number ------
# Aphid_Number is a count variable with the right-skewed distribution
summary(k.dat$Aphid_Number)
ks.test(k.dat$Aphid_Number,"pnorm")

hist(k.dat$Aphid_Number, breaks = 50, freq=FALSE)
plot(density(k.dat$Aphid_Number))

mean(k.dat$Aphid_Number)
var(k.dat$Aphid_Number)

ks.test(k.dat$Aphid_Number,"ppois", lambda=mean(k.dat$Aphid_Number))
data1 <- rpois(n=120, lambda=17.6)
ks.test(data1, k.dat$Aphid_Number)

# change family
# library(MASS)
mm2b <- glmmPQL(Aphid_Number ~ Hst_div + Plant_weight + Predtr  + 
                    Garlic_weight , 
                  random = ~ 1 | Pot,  data = k.dat,
                family = "quasipoisson") 

PearsonResiduals <- resid(mm2b, type = "pearson")
n_cases <- nrow(k.dat) # extract number of samples
n_predictors <- length(fixef(mm2b)) + 1 # extract number of estimates (plus intercept)
Overdispersion <- sum(PearsonResiduals^2) / (n_cases-n_predictors) # calculate overdispersion
Overdispersion # The overdispersion has decreased and is no longer an issue.

plot(mm2b) #  the plot exhibits a slight funnel shape, but not drastically, and thus indicates heteroscedasticity.
qqnorm(resid(mm2b))
qqline(resid(mm2b))

# check multicolinearity
car::vif(mm2b)

# estimates
summary(mm2b)
car::Anova(mm2b)



# try negative binomial for comparison of the results to the quasipoisson
mm2_c <- glmer.nb(Aphid_Number ~ Hst_div + Plant_weight + Predtr  + 
                    Garlic_weight + (1|Pot), data = k.dat) 

getME(mm2_c, "glmer.nb.theta")


plot(mm2_c)
car::Anova(mm2_c)



# Quasi-Poisson vs Negative binomial----

k.dat$pred_qpoiss <- predict(mm2b, k.dat, type = "response",)  # ,level=0
k.dat$resid_qpoiss <- residuals(mm2b)


k.dat$pred_nb <- predict(mm2_c, k.dat, type = "response")
k.dat$resid_nb <- residuals(mm2_c)



new_dat <- k.dat %>% 
  dplyr::select(Plant_N, pred_qpoiss,  resid_qpoiss,
        pred_nb, resid_nb) %>% 
   pivot_longer(cols=c(pred_qpoiss,  resid_qpoiss,
                       pred_nb, resid_nb),
                names_to = c("property", "familly"),
                names_pattern = "(.*)_(.*)") %>% 
  pivot_wider(names_from = property, values_from = value) %>% 
  mutate(pred_groups   = cut_number(pred, n = 10)) %>% 
  group_by(familly, pred_groups) %>% 
  summarise(
    count = n(),
    mean=mean(pred),
    sq_resid=mean(resid^2),
    sum=sum(resid^2)) %>% 
  ungroup()
 



new_dat %>% 
  mutate(familly= fct_recode(familly, "negative binomial"="nb", "quasi-Poisson"="qpoiss")) %>% 
  ggplot(aes(mean, sum))+
 # geom_smooth(aes(color=familly, fill=familly), method='lm', alpha=0.1, se = FALSE) +
    geom_point(aes(color=familly, size=count, fill=familly), pch=21, alpha=0.6) +
  labs(x="Predicted mean", y="Sum of Squared Residuals", colour="GLMM", fill="GLMM")+
  guides(size = "none")+ theme(legend.key=element_blank())


new_dat <- k.dat %>% 
  dplyr::select(Plant_N, pred_qpoiss,  resid_qpoiss,
                pred_nb, resid_nb) %>% 
  pivot_longer(cols=c(pred_qpoiss,  resid_qpoiss,
                      pred_nb, resid_nb),
               names_to = c("property", "familly"),
               names_pattern = "(.*)_(.*)") %>% 
  pivot_wider(names_from = property, values_from = value) %>% 
  mutate(sq_resid=resid^2)


new_dat %>% 
  mutate(familly= fct_recode(familly, "negative binomial"="nb", "quasi-Poisson"="qpoiss")) %>% 
  ggplot(aes(pred, sq_resid))+
  geom_smooth(aes(color=familly, fill=familly), method='lm', alpha=0.1, se = FALSE) +
  geom_point(aes(color=familly, fill=familly), pch=21, alpha=0.6, size=2) +
labs(x="Predicted mean", y="Squared residuals", colour="GLMM", fill="GLMM")+ theme(legend.key=element_blank())


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

col <- c("#57B1C9","#81C5AD", "#FFFF4B", "#DB9B99" )
fill <- c("#57B1C9","#81C5AD", "#FFFF4B", "#DB9B99" )


ggplot(R2, aes(y=Rsq, x=Effect)) + 
geom_bar(position="stack", stat="identity", colour = "black", fill="#DCE6F2")+
  coord_flip() +
  #  ylab(bquote('Explained variance (partial R'^'2'*')'))+
  ylab(bquote('partial R'^'2'))+
  labs(x =element_blank()) +
   scale_x_discrete(labels= c( "Garlic biomass",
                               "Host-plant
                               biomass",
                               "Host-plant 
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
# Figure 3 -----

## Figure 3A----
### Interactive effects  ----

# glmmTMB for fitting quasipoisson
# https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf
# https://stackoverflow.com/questions/75799875/glmm-with-quasi-poisson-distribution

mS4 <- glmmTMB(Aphid_Number ~  Hst_div + 
                 Plant_weight + Predtr  + Garlic + 
                 Predtr:Plant_weight +
                 Predtr:Hst_div +
                 Predtr:Garlic + (1|Pot),
        data=k.dat, family=nbinom1) # nbinom1 fits quasipoisson familly




library(performance)
check_convergence(mS4)
check_collinearity(mS4)

# estimates
summary(mS4)
car::Anova(mS4)
write.csv(Anova(mS4),  file = "results/Table_S4.csv")


# Plot "Plant_weight" interaction with predator

mS4a <- glmmTMB(Aphid_Number ~   Hst_div + 
                 Plant_weight + Predtr  + Garlic + 
                 Predtr:Plant_weight +
                 #  Predtr:Hst_div +
                 #  Predtr:Garlic + 
                 (1|Pot),
               data=k.dat , # to include zero inflation use ziformula=~1
               family=nbinom1) # nbinom1 fits quasipoisson


plot_model(mS4a,type = "pred", terms=c("Plant_weight[0.004:0.063, by=.001]", "Predtr"),   show.data=T,
           title = "", line.size=1)

pred_host_Pred <-get_model_data(mS4a, type = "pred", 
                                terms=c("Plant_weight[0.004:0.063, by=.001]", "Predtr"))
colr=c("#FC4E07", "#00AFBB")

ggplot(pred_host_Pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, 
              data=as.data.frame(pred_host_Pred) %>%
                filter(group_col=="absent"), fill="#FC4E07")+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, 
               data=as.data.frame(pred_host_Pred) %>%
                filter(group_col=="present"), fill="#00AFBB")+
  geom_point(data=k.dat, aes(Plant_weight, Aphid_Number, fill = Predtr),
             size=2, alpha=0.6, pch=21)+
  scale_fill_manual(values = colr)+  scale_color_manual(values="gray" , guide = 'none')  +  
  labs(y='Aphid density', x="Host-plant biomass, g", 
       fill="Predator") +
    geom_line(aes(x, predicted), data=as.data.frame(pred_host_Pred) %>%
               filter(group_col=="present"), col="#00AFBC", size=1) +
   geom_line(aes(x, predicted), data=as.data.frame(pred_host_Pred) %>%
              filter(group_col=="absent"), col="#FC4E07", size=1) 
 ## geom_line(aes(x, predicted), data=as.data.frame(pred_host_Pred),
  #          col="black", size=1)




## Figure 3B----
# test the interaction of predator presence with garlic biomass

mS5 <- glmer.nb(Aphid_Number ~ Hst_div + 
                  Plant_weight + Predtr  + Garlic_weight + 
                  Predtr:Plant_weight+
                  Predtr:Hst_div +
                  Predtr:Garlic_weight +  (1|Pot), data = k.dat) 

check_convergence(mS5)

# estimates
summary(mS5)
car::Anova(mS5)
write.csv(Anova(mS5),  file = "results/Table_S5.csv")

# Plot "Predtr:Garlic_weight"

plot_model(mS5,type = "pred", terms=c("Garlic_weight[0:0.872, by=.001]", "Predtr"),   show.data=T,
           title = "", line.size=1)

pred_GarlWght_Pred <-get_model_data(mS5, type = "pred", 
                                    terms=c("Garlic_weight[0:0.872, by=.001]", "Predtr"))
colr=c("#FC4E07", "#00AFBB")

ggplot(pred_GarlWght_Pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, 
              data=as.data.frame(pred_GarlWght_Pred) %>%
                filter(group_col=="present"), fill="#00AFBB")+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, 
              data=as.data.frame(pred_GarlWght_Pred) %>%
                filter(group_col=="absent"), fill="#FC4E07")+
  geom_point(data=k.dat, aes(Garlic_weight, Aphid_Number, fill = Predtr),
             size=2, alpha=0.7, pch=21)+
  scale_fill_manual(values = colr)+  scale_color_manual(values="gray" , guide = 'none')  +  
  labs(y='Aphid density', x="Garlic biomass, g", 
       fill="Predator") +
  geom_line(aes(x, predicted), data=as.data.frame(pred_GarlWght_Pred) %>%
              filter(group_col=="present"), col="#00AFBC", size=1) +
  geom_line(aes(x, predicted), data=as.data.frame(pred_GarlWght_Pred) %>%
              filter(group_col=="absent"), col="#FC4E07", size=1) 



#################################################-

## Figure 3C ----
## effects of treatment
k.dat$Treatment_specific

mean.data <-k.dat %>% 
  group_by(Treatment_general,Pot) %>% 
  summarise(Aphid_Number=mean(Aphid_Number))


m <- lm(log(Aphid_Number) ~ Treatment_general, data = mean.data) 
summary(m)
anova(m)
plot(m)

emmeans::emmeans(m, list(pairwise ~ Treatment_general))
# to add letters for post-hoc test:
multcomp::cld(emmeans::emmeans(m, list(pairwise ~ Treatment_general)),  
              #  type="response",
              Letters = letters, adjust = "none")


mm3 <- glmer(Aphid_Number ~ Treatment_general + Plant_weight + (1|Pot), data = k.dat, 
            family = "poisson") 

plot(mm3)
qqnorm(resid(mm3))
qqline(resid(mm3))

# check overdispersion
resid_pearson <- residuals(mm3, type = "pearson")
SSQ <- sum(resid_pearson^2)
SSQ/df.residual(mm2) 


# change family
library(MASS)

mm3b <- glmmPQL(Aphid_Number ~ Treatment_general  + Plant_weight, 
                random = ~ 1 | Pot,  data = k.dat,
                family = "quasipoisson") 



plot(mm3b)
qqnorm(resid(mm3b))
qqline(resid(mm3b))


# estimates
summary(mm3b)
car::Anova(mm3b)

# Marginal means and pairwise differences of Hst_div levels and of Predtr levels

emmeans::emmeans(mm3b, list(pairwise ~ Treatment_general  ))
# to add letters for post-hoc test:
fig1_emmeans <-multcomp::cld(emmeans::emmeans(mm3b, list(pairwise ~ Treatment_general)),  
              #  type="response",
              Letters = letters, adjust = "none")


fig1_emmeans <- fig1_emmeans[order(fig1_emmeans$Treatment_general), ]
fig1_emmeans



fig1_summarized = k.dat %>% 
  group_by(Treatment_general, Group) %>% 
  summarize(Max.Aphid_Number=max(Aphid_Number),
            mean=mean(Aphid_Number),
            sd=sd(Aphid_Number),
            se = sd(Aphid_Number, na.rm = T)/sqrt(length(Aphid_Number))) 

fig1_summarized


ggplot(fig1_summarized, aes(y=mean, x=Treatment_general, fill=Group), col="black") + 
  geom_bar(stat="identity", position=position_dodge(), col="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9))+
  geom_text(data=fig1_summarized,aes(x=Treatment_general,y=0.5+mean+se,
                                       label=fig1_emmeans$.group),vjust=0, size=5)+
 
 # geom_point(data=fig1_emmeans, aes(y=emmean, x=Treatment_general),
 #            shape = 23, color = "black", fill="blue",  size=3) +
  ylim(0,33)+
  labs(x =" ", y="Aphid density") +
  scale_fill_brewer(palette="Pastel2") +
  theme_bw()+
  theme(axis.text.y=element_text(colour = "black", size=11),
        axis.text.x=element_blank(),
        axis.title=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.ticks =  element_line(colour = "black"),
        legend.position="none")

#################################################-
## Figure 3D ----
### Combination of control agents----

mS6 <- glmer.nb(Aphid_Number ~ Group +  (1|Pot), data = k.dat)

check_convergence(mS6)
# check_collinearity(mS6)
plot(mS6)

# estimates
summary(mS6)
car::Anova(mS6)


# Marginal means and pairwise differences of Hst_div levels and of Predtr levels

emmeans::emmeans(mS6, list(pairwise ~ Group))

# Number_agents

# to add letters for post-hoc test:
fig_emmeans_S6 <- multcomp::cld(emmeans::emmeans(mS6, list(pairwise ~ Group)),  
                                #  type="response",
                                Letters = letters, adjust = "none")

fig_emmeans_S6 <- fig_emmeans_S6[order(fig_emmeans_S6$Group),]
fig_emmeans_S6


figS6_summarized = k.dat %>% 
  group_by(Group) %>% 
  summarize(Max.Aphid_Number=max(Aphid_Number),
            mean=mean(Aphid_Number),
            sd=sd(Aphid_Number),
            se = sd(Aphid_Number, na.rm = T)/sqrt(length(Aphid_Number))) 

figS6_summarized


ggplot(figS6_summarized, aes(y=mean, x=Group, fill=Group), col="black") + 
  geom_bar(stat="identity", position=position_dodge(), col="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) +
  geom_text(data=figS6_summarized,aes(x=Group,y=1+mean+se,
                                      label=c("a", "a", "b", "c")),vjust=0, size=5,
            inherit.aes = TRUE , position = position_dodge(width = 0.9))+
  ylim(0,31) + 
  labs(x =" ", y="Aphid density") +
  scale_fill_brewer(palette="Pastel2") +
  theme_bw()+
  theme(axis.text.y=element_text(colour = "black", size=11),
        axis.text.x=element_blank(),
        axis.title=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.ticks =  element_line(colour = "black"),
        legend.title=element_text(face="bold",size=11),
        legend.text=element_text(size=11),
        legend.position = "none")



#################################################-
## Figure 3E----
### Host biomass as response variable  ----

mod_host_mass <- lmer(Plant_weight  ~ Group +
                        # Treatment_general +
                        Hst_div + #Garlic +
                        (1|Pot), data = k.dat)
plot(mod_host_mass)
qqnorm(resid(mod_host_mass))
qqline(resid(mod_host_mass))

# check multicolinearity
car::vif(mod_host_mass)

# estimates
summary(mod_host_mass)
car::Anova(mod_host_mass)

# Marginal means and pairwise differences of Hst_div levels
emmeans::emmeans(mod_host_mass, list(pairwise ~ Group),
                 type="response", adjust = "Tukey")
# to add letters for the post-hoc test:
library(multcomp)
library(emmeans)
figS2_emmeans <- cld(emmeans(mod_host_mass, list(pairwise ~ Group)),  
                     type="response",
                     Letters = letters, adjust = "Tukey")

figS2_emmeans <- figS2_emmeans[order(figS2_emmeans$Group), ]
figS2_emmeans

#  model R2
MuMIn::r.squaredGLMM(mod_host_mass)
# R2m and R2c are marginal (for fixed predictors) and conditional (for fixed and random predictors) coefficients of determination

# Partial R2 for fixed effects
library(r2glmm)
R2part_mm1 <- r2beta(mod_host_mass, method = 'nsj', partial = T, data = k.dat)
R2part_mm1

figS2_summarized = k.dat %>% 
  group_by(Group) %>% 
  summarize(Max.host.nass=max(Plant_weight),
            mean=mean(Plant_weight),
            sd=sd(Plant_weight),
            se = sd(Plant_weight, na.rm = T)/sqrt(length(Plant_weight))) 

figS2_summarized


ggplot(figS2_summarized, aes(y=mean, x=Group, fill=Group), col="black") + 
  geom_bar(stat="identity", position=position_dodge(), col="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9))+
  geom_text(data=figS2_summarized,aes(x=Group,y=0.001+mean+se,
                                      label=figS2_emmeans$.group),vjust=0, size=5)+
  
  # geom_point(data=fig1_emmeans, aes(y=emmean, x=Treatment_general),
  #            shape = 23, color = "black", fill="blue",  size=3) +
  
  labs(x =" ", y="Host-plant biomass, g") +
  scale_fill_brewer(palette="Pastel2") +
  theme_bw()+ ylim(0,0.04)+
  theme(axis.text.y=element_text(colour = "black", size=11),
        axis.text.x=element_blank(),
        axis.title=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.ticks =  element_line(colour = "black"),
        legend.position="none")






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

figS1A_emmeans <- multcomp::cld(emmeans::emmeans(mS1, list(pairwise ~ Plant)),  
                                #  type="response",
                                Letters = letters, adjust = "none")

figS1A_emmeans <- figS1A_emmeans[order(figS1A_emmeans$Plant), ]
figS1A_emmeans


figS1A_summarized = k.dat %>% 
  group_by(Plant) %>% 
  summarize(Max.Plant_weight=max(Plant_weight),
            Mean.Plant_weight=mean(Plant_weight))


ggplot(k.dat, aes(y=Plant_weight, x=Plant)) + 
  geom_boxplot(outlier.shape=NA)+
  geom_text(data=figS1A_summarized,aes(x=Plant,y=0.002+Max.Plant_weight,
                                      label=figS1A_emmeans$.group),vjust=0, , size=5)+
  geom_jitter(data = k.dat, aes(y = Plant_weight, x = Plant),
    alpha = 0.2,width = 0.2) +
  geom_point(data=figS1A_summarized, aes(y=Mean.Plant_weight, x=Plant),  
             position = position_nudge(x=0.1),
             shape = 23, color = "black", fill="red", size=3) +
  geom_point(data=figS1A_emmeans, aes(y=emmean, x=Plant), 
             position = position_nudge(x=-0.1),
             shape = 23, color = "black", fill="blue",  size=3) +
  
  labs(x ="Plant species", y="Plant biomass (g)") +
  ylim(0.01, 0.07)+
  theme_bw()+
  theme(axis.text.y=element_text(colour = "black", size=13),
        axis.text.x=element_text(colour = "black", size=13),
        axis.title=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.ticks =  element_line(colour = "black"))



## Fig S1.B----
# Effects of plant species on aphid density

mS2 <- glmer(Aphid_Number ~  Plant +   Predtr  + 
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

mS2b <- glmmPQL(Aphid_Number ~ Plant  +   Predtr  +  Garlic_weight, 
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

figS1B_emmeans <- multcomp::cld(emmeans::emmeans(mS2b, list(pairwise ~ Plant)),  
              #  type="response",
              Letters = letters, adjust = "none")

figS1B_emmeans <- figS1B_emmeans[order(figS1B_emmeans$Plant), ]
figS1B_emmeans


figS1B.summarized = k.dat %>% 
  group_by(Plant) %>% 
  summarize(Max.Aphid_Number=max(Aphid_Number),
            Mean.Aphid_Number=mean(Aphid_Number))


 
  
ggplot(k.dat, aes(x=Plant,y=Aphid_Number)) + 
  geom_boxplot(outlier.shape=NA)+
  geom_text(data=figS1B.summarized,aes(x=Plant,y=3+Max.Aphid_Number,
                                     label=figS1B_emmeans$.group),vjust=0, , size=5)+
  geom_jitter(data = k.dat, aes(y = Aphid_Number, x = Plant),
              alpha = 0.2,width = 0.2) +
  geom_point(data=figS1B.summarized, aes(y=Mean.Aphid_Number, x=Plant), 
             position = position_nudge(x=0.1),
             shape = 23, color = "black", fill="red", , size=3) +
  geom_point(data=figS1B_emmeans, aes(y=emmean, x=Plant), 
             position = position_nudge(x=-0.1),
             shape = 23, color = "black", fill="blue",  size=3) +
  
  labs(x ="Plant species", y="Aphid density (number)") +
  ylim(0, 65)+
    theme_bw()+
  theme(axis.text.y=element_text(colour = "black", size=13),
        axis.text.x=element_text(colour = "black", size=13),
        axis.title=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.ticks =  element_line(colour = "black"))


##Fig S1.C ----------

# Effects of plant species on aphid density/host biomass

k.dat <- k.dat %>% 
  mutate(Aphid_load=Aphid_Number/Plant_weight)

mS3 <- lmer(sqrt(Aphid_load) ~ Plant +   Predtr  + 
               Garlic_weight + (1|Pot), data = k.dat) 
# check model
summary(mS3)


# check multicolinearity
car::vif(mS3)

plot(mS3)
qqnorm(resid(mS3))

qqline(resid(mS3))


car::Anova(mS3)

# Marginal means and pairwise differences of Plant species
emmeans::emmeans(mS3, list(pairwise ~ Plant))

figS1C_emmeans <- multcomp::cld(emmeans::emmeans(mS3, list(pairwise ~ Plant)),  
                                #  type="response",
                                Letters = letters, adjust = "none")

figS1C_emmeans <- figS1C_emmeans[order(figS1C_emmeans$Plant), ]
figS1C_emmeans


figS1C.summarized = k.dat %>% 
  group_by(Plant) %>% 
  summarize(Max.Aphid_load=max(sqrt(Aphid_load)),
            Mean.Aphid_load=mean(sqrt(Aphid_load)))




ggplot(k.dat, aes(x=Plant,y=sqrt(Aphid_load))) + 
  geom_boxplot(outlier.shape=NA)+
  geom_text(data=figS1C.summarized,aes(x=Plant,y=c(1550, 1260, 1100),
                                       label=figS1C_emmeans$.group),vjust=0, , size=5)+
  geom_jitter(data = k.dat, aes(y = Aphid_load, x = Plant),
              alpha = 0.2,width = 0.2) +
  geom_point(data=figS1C.summarized, aes(y=Mean.Aphid_load, x=Plant), 
             position = position_nudge(x=0.1),
             shape = 23, color = "black", fill="red", , size=3) +
  geom_point(data=figS1C_emmeans, aes(y=emmean, x=Plant), 
             position = position_nudge(x=-0.1),
             shape = 23, color = "black", fill="blue",  size=3) +
  
  labs(x ="Plant species", y="Aphid load (number/g)") +
#  ylim(0, 65)+
  theme_bw()+
  theme(axis.text.y=element_text(colour = "black", size=13),
        axis.text.x=element_text(colour = "black", size=13),
        axis.title=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.ticks =  element_line(colour = "black"))



#################################################-
