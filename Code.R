##### PACKAGES ##### 
library(lcmm)
library(NormPsy)
library(lattice)
library(ggpubr)
library(Rmisc)
library(splines)
library(splines2)
library(MASS)
library(sqldf)
library(dplyr)
library(ggplot2)
library(viridis)
library(lmtest)
library(moments)
library(dplyr, warn = FALSE)
library(ggplot2)
library(patchwork)
library(readxl)
library(xlsx)
library(survminer)
library(foreach) 
library(doParallel)


tab_tot3  = read.table(file="C:/Users/mwagner4/Downloads/LIBRA score/3C/tab_tot3.txt",header=T,sep="\t",dec=".")
tab_tot3$DIPNIV_3cat = as.factor(tab_tot3$DIPNIV_3cat)
tab_tot3$CENTRE = as.factor(tab_tot3$CENTRE)
tab_tot3$SEXE = as.factor(tab_tot3$SEXE)
dim(tab_tot3)
tab_apoe_genotype   = read.table(file="C:/Users/mwagner4/Downloads/LIBRA score/3C/Apoe_genotype.txt",header=T,sep="\t",dec=".")
str(tab_apoe_genotype)
dim(tab_apoe_genotype) # N = 9294
tab_tot3 = merge(tab_tot3, tab_apoe_genotype, by = "NUM")
dim(tab_tot3)
tab_tot3 = tab_tot3[order(tab_tot3$NUM, tab_tot3$SUIVI),]
## merge of different CR measure ##
tab_cov = subset(tab_tot3, select = c("NUM",
                                      "genotype",
                                      "scoreLIBRAF",
                                      "AlcMod.libra",
                                      "CoronaryDisease.libra",
                                      "PhysicalInactivity.libra",
                                      "RenalDisease.libra",
                                      "Diabetes.libra",
                                      "Hyperchol.libra",
                                      "Tabac.libra",
                                      "Obesity.libra",
                                      "HTA.libra",
                                      "MediDiet.libra",
                                      "Depression.libra",
                                      "CogAct.libra",
                                      "follow",
                                      "nb_tot_mmse",
                                      "SUIVI0","time_y",
                                      "SEXE","CENTRE","AGE0","AGE0_75",
                                      "DIPNIV_3cat", "APOE4","MMSTT","norm_MMSE",
                                      "BENTON","ISA_15","TMT_A","TMT_B"))
tab_cov = tab_cov[order(tab_cov$NUM, tab_cov$SUIVI),]
tab_APOE = as.data.frame(tab_cov %>% group_by(NUM) %>% filter(row_number(NUM) == 1))


load(file='C:/Users/mwagner4/Downloads/LIBRA score/3C/m1_global_apoe4_splinesMMSE_bis.R')
summary(m1_global_apoe4_splinesMMSE_bis) # without posfix




#############################
####       BOOTSTRAPS   #####
#############################
m = m1_global_apoe4_splinesMMSE_bis
nboot = 1000
mu    = as.matrix(estimates(m)) 
Sigma = as.matrix(VarCov(m)) 
boot  = as.matrix(mvrnorm(nboot, mu, Sigma))
head(boot)
mtemp = m

########################################
##  INDIVIDUAL SLOPES - Fixed effects ##
########################################
summary(m1_global_apoe4_splinesMMSE_bis)
#
beta_apoe4       = -0.08446
beta_time        = -0.16728
beta_apoe4_time  = -0.02662

# mtemp = multlcmm(fixed = norm_MMSE + BENTON + ISA_15 + TMT_A + TMT_B ~ 
#                    SUIVI0 + SEXE + CENTRE + AGE0_75 + DIPNIV_3cat + APOE4 + 
#                    time_y + APOE4 * time_y + SEXE * time_y + DIPNIV_3cat * 
#                    time_y + AGE0_75 * time_y + CENTRE * time_y, random = ~time_y, 
#                  subject = "NUM", randomY = T, link = c("linear", "3-quant-splines", 
#                                                         "3-quant-splines", "3-quant-splines", "3-quant-splines"), 
#                  data = tab_tot3, maxiter = 0)
# 
# mtemp$best = boot[2,]
# RE_boot = predictRE(mtemp, newdata = tab_tot3) # Vector B should be of length 41
# head(RE_boot)
# 
# mtemp$best = boot[3,]
# RE_boot = predictRE(mtemp, newdata = tab_tot3) # Vector B should be of length 41
# head(RE_boot)




for (i in 1:nboot){
  
  
  mtemp$best = boot[i,]
  RE_boot = predictRE(mtemp, newdata = tab_tot3) # Vector B should be of length 41
  # head(RE_boot)
  
  ##########################################
  ### CLASSIFICATION (modele avec APOE4) ###
  ##########################################
  tab_apoe = subset(tab_tot3, select = c("NUM", "APOE4"))
  tempo    = as.data.frame(tab_apoe %>% group_by(NUM) %>% filter(row_number(NUM) == 1))
  #dim(tempo) # n = 6774
  # model with linear function of time #
  tempo1 = as.data.frame(RE_boot)
  tempo1$RE_int   = tempo1$intercept
  tempo1$RE_slope = tempo1$time_y
  tempo1 = subset(tempo1, select = c("NUM","RE_int","RE_slope"))
  #
  tab_RE0 = merge(tempo1, tempo, by = c("NUM"), all = T)
  #dim(tab_RE0) # n = 6774
  

  
  ##########################
  ### INTRA-DISTRIBUTION ###
  ##########################
  
  # linear function of time #
  tab_RE0$CR_level = ifelse(tab_RE0$APOE4 == 1, tab_RE0$RE_int + beta_apoe4, tab_RE0$RE_int)
  tab_RE0$CR_slope = ifelse(tab_RE0$APOE4 == 1, tab_RE0$RE_slope + beta_time + beta_apoe4_time, tab_RE0$RE_slope + beta_time)
  
  ###########################
  ## Cut-off for CR status ##
  ###########################
  # Linear time #
  P50_level = quantile(tab_RE0$CR_level, probs = c(0.50))
  P75_level = quantile(tab_RE0$CR_level, probs = c(0.75))
  P50_slope = quantile(tab_RE0$CR_slope, probs = c(0.50))
  P75_slope = quantile(tab_RE0$CR_slope, probs = c(0.75))
  # c(P50_level, P75_level, P50_slope, P75_slope)
  
  
  ##############
  ## CR_SLOPE ##
  ##############
  ##############
  ### LINEAR ###
  ##############
  # <P75 vs. P>=75 (25% slowest)
  tab_RE0$CR_slope_P75_all = as.factor(ifelse(tab_RE0$CR_slope  >= P75_slope, 1, 0))
  
  # <P50 vs. P>=50 (50% slowest)
  tab_RE0$CR_slope_P50_all = as.factor(ifelse(tab_RE0$CR_slope  >= P50_slope, 1, 0))
  
  # <P50 vs. P>=75 (25% slowest)
  tab_RE0$CR_slope_P50_P75 = NA
  tab_RE0$CR_slope_P50_P75 = ifelse(tab_RE0$CR_slope >= P75_slope, 1, tab_RE0$CR_slope_P50_P75)
  tab_RE0$CR_slope_P50_P75 = as.factor(ifelse(tab_RE0$CR_slope < P50_slope, 0, tab_RE0$CR_slope_P50_P75))
  
  # summary(tab_RE0$CR_slope_P75_all)
  # summary(tab_RE0$CR_slope_P50_all)
  # summary(tab_RE0$CR_slope_P50_P75) # verif OK
  
  
  ###########################
  ####   MERGING DATA    ####
  ###########################
  #dim(tab_RE0)
  tab_CR = subset(tab_RE0, select= c("NUM",
                                     "CR_level",
                                     "RE_slope",
                                     "CR_slope",
                                     "CR_slope_P75_all", # <P75 vs. P>=75 (25% slowest)
                                     "CR_slope_P50_all", # <P50 vs. P>=50 (50% slowest)
                                     "CR_slope_P50_P75"))  # <P50 vs. P>=75 (25% slowest)
  
  # dim(tab_APOE) # N = 6774 #
  # dim(tab_CR)    # N = 6774 #
  #
  tab_CR_slope = merge(tab_APOE, tab_CR, by = c("NUM"), all = T)
  tab_CR_slope$CogAct.libra     = as.factor(tab_CR_slope$CogAct.libra)
  tab_CR_slope$MediDiet.libra   = as.factor(tab_CR_slope$MediDiet.libra)
  tab_CR_slope$AlcMod.libra     = as.factor(tab_CR_slope$AlcMod.libra)
  tab_CR_slope$CoronaryDisease.libra = as.factor(tab_CR_slope$CoronaryDisease.libra)
  tab_CR_slope$PhysicalInactivity.libra = as.factor(tab_CR_slope$PhysicalInactivity.libra)
  tab_CR_slope$RenalDisease.libra = as.factor(tab_CR_slope$RenalDisease.libra)
  tab_CR_slope$Diabetes.libra   = as.factor(tab_CR_slope$Diabetes.libra)
  tab_CR_slope$Hyperchol.libra  = as.factor(tab_CR_slope$Hyperchol.libra)
  tab_CR_slope$Tabac.libra      = as.factor(tab_CR_slope$Tabac.libra)
  tab_CR_slope$HTA.libra        = as.factor(tab_CR_slope$HTA.libra)
  tab_CR_slope$Obesity.libra    = as.factor(tab_CR_slope$Obesity.libra)
  tab_CR_slope$Depression.libra = as.factor(tab_CR_slope$Depression.libra)
  # dim(tab_CR_slope) # n = 6774
  
  t = subset(tab_CR_slope, APOE4 == 1)
  t$LIBRA2 = -t$scoreLIBRAF
  
  ## LIBRA ## 
  a = glm(CR_slope_P50_P75 ~ LIBRA2 + AGE0_75 + SEXE + DIPNIV_3cat + CENTRE, data=t, family = "binomial")
  y = c(i,coef(a))
  t1 = as.data.frame(y)
  p1 = as.data.frame(coef(summary(a))[,'Pr(>|z|)'])
  s1 = as.data.frame(coef(summary(a))[,'Std. Error'])
  
  a = glm(CR_slope_P75_all ~ LIBRA2 + AGE0_75 + SEXE + DIPNIV_3cat + CENTRE, data=t, family = "binomial")
  y = c(i,coef(a))
  t2 = as.data.frame(y)
  p2 = as.data.frame(coef(summary(a))[,'Pr(>|z|)'])
  s2 = as.data.frame(coef(summary(a))[,'Std. Error'])
  
  a = glm(CR_slope_P50_all ~ LIBRA2 + AGE0_75 + SEXE + DIPNIV_3cat + CENTRE, data=t, family = "binomial")
  y = c(i,coef(a))
  t3 = as.data.frame(y)
  p3 = as.data.frame(coef(summary(a))[,'Pr(>|z|)'])
  s3 = as.data.frame(coef(summary(a))[,'Std. Error'])
  
  # excluding e2/e4 
  tempo = subset(t, genotype!="24")
  a = glm(CR_slope_P50_P75 ~ LIBRA2 + AGE0_75 + SEXE + DIPNIV_3cat + CENTRE, data=tempo, family = "binomial")
  y = c(i,coef(a))
  t4 = as.data.frame(y)
  p4 = as.data.frame(coef(summary(a))[,'Pr(>|z|)'])
  s4 = as.data.frame(coef(summary(a))[,'Std. Error'])
  
  # excluding e4/e4 
  tempo = subset(t, genotype!="44")
  a = glm(CR_slope_P50_P75 ~ LIBRA2 + AGE0_75 + SEXE + DIPNIV_3cat + CENTRE, data=tempo, family = "binomial")
  y = c(i,coef(a))
  t5 = as.data.frame(y)
  p5 = as.data.frame(coef(summary(a))[,'Pr(>|z|)'])
  s5 = as.data.frame(coef(summary(a))[,'Std. Error'])
  
  if (i==1){ 
    tab_boot_P50_P75_libra = t1
    tab_boot_P75_all_libra = t2
    tab_boot_P50_all_libra = t3 
    tab_boot_P50_P75_libra_e2e4 = t4
    tab_boot_P50_P75_libra_e4e4 = t5
    #
    sd_boot_P50_P75_libra = s1
    sd_boot_P75_all_libra = s2
    sd_boot_P50_all_libra = s3 
    sd_boot_P50_P75_libra_e2e4 = s4 
    sd_boot_P50_P75_libra_e4e4 = s5
    #
    p_boot_P50_P75_libra = p1
    p_boot_P75_all_libra = p2
    p_boot_P50_all_libra = p3 
    p_boot_P50_P75_libra_e2e4 = p4 
    p_boot_P50_P75_libra_e4e4 = p5
    
  } else if (i>1){ 
    tab_boot_P50_P75_libra = cbind(tab_boot_P50_P75_libra, t1) 
    tab_boot_P75_all_libra = cbind(tab_boot_P75_all_libra, t2) 
    tab_boot_P50_all_libra = cbind(tab_boot_P50_all_libra, t3)
    tab_boot_P50_P75_libra_e2e4 = cbind(tab_boot_P50_P75_libra_e2e4, t4) 
    tab_boot_P50_P75_libra_e4e4 = cbind(tab_boot_P50_P75_libra_e4e4, t5)  
    #
    sd_boot_P50_P75_libra = cbind(sd_boot_P50_P75_libra, s1) 
    sd_boot_P75_all_libra = cbind(sd_boot_P75_all_libra, s2) 
    sd_boot_P50_all_libra = cbind(sd_boot_P50_all_libra, s3)  
    sd_boot_P50_P75_libra_e2e4 = cbind(sd_boot_P50_P75_libra_e2e4, s4) 
    sd_boot_P50_P75_libra_e4e4 = cbind(sd_boot_P50_P75_libra_e4e4, s5) 
    #
    p_boot_P50_P75_libra = cbind(p_boot_P50_P75_libra, p1) 
    p_boot_P75_all_libra = cbind(p_boot_P75_all_libra, p2) 
    p_boot_P50_all_libra = cbind(p_boot_P50_all_libra, p3)
    p_boot_P50_P75_libra_e2e4 = cbind(p_boot_P50_P75_libra_e2e4, p4) 
    p_boot_P50_P75_libra_e4e4 = cbind(p_boot_P50_P75_libra_e4e4, p5)
  } 
  
  
  ## Components ##
  a = glm(CR_slope_P50_P75 ~ CogAct.libra +
            AlcMod.libra +
            MediDiet.libra +
            #
            CoronaryDisease.libra +
            Tabac.libra+
            Obesity.libra+
            Depression.libra +
            HTA.libra+
            Diabetes.libra+
            PhysicalInactivity.libra +
            RenalDisease.libra+
            Hyperchol.libra+
            AGE0_75 + SEXE + DIPNIV_3cat + CENTRE, data=t, family = "binomial")
  y = c(i,coef(a))
  t1 = as.data.frame(y)
  p1 = as.data.frame(coef(summary(a))[,'Pr(>|z|)'])
  s1 = as.data.frame(coef(summary(a))[,'Std. Error'])
  
  a = glm(CR_slope_P75_all ~ CogAct.libra +
            AlcMod.libra +
            MediDiet.libra +
            #
            CoronaryDisease.libra +
            Tabac.libra+
            Obesity.libra+
            Depression.libra +
            HTA.libra+
            Diabetes.libra+
            PhysicalInactivity.libra +
            RenalDisease.libra+
            Hyperchol.libra+
            AGE0_75 + SEXE + DIPNIV_3cat + CENTRE, data=t, family = "binomial")
  y = c(i,coef(a))
  t2 = as.data.frame(y)
  p2 = as.data.frame(coef(summary(a))[,'Pr(>|z|)'])
  s2 = as.data.frame(coef(summary(a))[,'Std. Error'])
  
  a = glm(CR_slope_P50_all ~ CogAct.libra +
            AlcMod.libra +
            MediDiet.libra +
            #
            CoronaryDisease.libra +
            Tabac.libra+
            Obesity.libra+
            Depression.libra +
            HTA.libra+
            Diabetes.libra+
            PhysicalInactivity.libra +
            RenalDisease.libra+
            Hyperchol.libra+
            AGE0_75 + SEXE + DIPNIV_3cat + CENTRE, data=t, family = "binomial")
  y = c(i,coef(a))
  t3 = as.data.frame(y)
  p3 = as.data.frame(coef(summary(a))[,'Pr(>|z|)'])
  s3 = as.data.frame(coef(summary(a))[,'Std. Error'])
  
  # excluding e2/e4 
  tempo = subset(t, genotype!="24")
  a = glm(CR_slope_P50_all ~ CogAct.libra +
            AlcMod.libra +
            MediDiet.libra +
            #
            CoronaryDisease.libra +
            Tabac.libra+
            Obesity.libra+
            Depression.libra +
            HTA.libra+
            Diabetes.libra+
            PhysicalInactivity.libra +
            RenalDisease.libra+
            Hyperchol.libra+
            AGE0_75 + SEXE + DIPNIV_3cat + CENTRE, data=tempo, family = "binomial")
  y = c(i,coef(a))
  t4 = as.data.frame(y)
  p4 = as.data.frame(coef(summary(a))[,'Pr(>|z|)'])
  s4 = as.data.frame(coef(summary(a))[,'Std. Error'])
  
  # excluding e4/e4 
  tempo = subset(t, genotype!="44")
  a = glm(CR_slope_P50_all ~ CogAct.libra +
            AlcMod.libra +
            MediDiet.libra +
            #
            CoronaryDisease.libra +
            Tabac.libra+
            Obesity.libra+
            Depression.libra +
            HTA.libra+
            Diabetes.libra+
            PhysicalInactivity.libra +
            RenalDisease.libra+
            Hyperchol.libra+
            AGE0_75 + SEXE + DIPNIV_3cat + CENTRE, data=tempo, family = "binomial")
  y = c(i,coef(a))
  t5 = as.data.frame(y)
  p5 = as.data.frame(coef(summary(a))[,'Pr(>|z|)'])
  s5 = as.data.frame(coef(summary(a))[,'Std. Error'])
  
  if (i==1){ 
    
    tab_boot_P50_P75_compo = t1     
    tab_boot_P75_all_compo = t2 
    tab_boot_P50_all_compo = t3 
    tab_boot_P50_P75_compo_e2e4 = t4
    tab_boot_P50_P75_compo_e4e4 = t5
    #
    sd_boot_P50_P75_compo = s1
    sd_boot_P75_all_compo = s2
    sd_boot_P50_all_compo = s3 
    sd_boot_P50_P75_compo_e2e4 = s4 
    sd_boot_P50_P75_compo_e4e4 = s5
    #
    p_boot_P50_P75_compo = p1
    p_boot_P75_all_compo = p2
    p_boot_P50_all_compo = p3 
    p_boot_P50_P75_compo_e2e4 = p4 
    p_boot_P50_P75_compo_e4e4 = p5
    
  } else if (i>1){ 
    
    tab_boot_P50_P75_compo = cbind(tab_boot_P50_P75_compo, t1) 
    tab_boot_P75_all_compo = cbind(tab_boot_P75_all_compo, t2) 
    tab_boot_P50_all_compo = cbind(tab_boot_P50_all_compo, t3) 
    tab_boot_P50_P75_compo_e2e4 = cbind(tab_boot_P50_P75_compo_e2e4, t4) 
    tab_boot_P50_P75_compo_e4e4 = cbind(tab_boot_P50_P75_compo_e4e4, t5) 
    #
    sd_boot_P50_P75_compo = cbind(sd_boot_P50_P75_compo, s1) 
    sd_boot_P75_all_compo = cbind(sd_boot_P75_all_compo, s2) 
    sd_boot_P50_all_compo = cbind(sd_boot_P50_all_compo, s3)  
    sd_boot_P50_P75_compo_e2e4 = cbind(sd_boot_P50_P75_compo_e2e4, s4) 
    sd_boot_P50_P75_compo_e4e4 = cbind(sd_boot_P50_P75_compo_e4e4, s5) 
    #
    p_boot_P50_P75_compo = cbind(p_boot_P50_P75_compo, p1) 
    p_boot_P75_all_compo = cbind(p_boot_P75_all_compo, p2) 
    p_boot_P50_all_compo = cbind(p_boot_P50_all_compo, p3)
    p_boot_P50_P75_compo_e2e4 = cbind(p_boot_P50_P75_compo_e2e4, p4) 
    p_boot_P50_P75_compo_e4e4 = cbind(p_boot_P50_P75_compo_e4e4, p5)
    
  } 
  
  i = i + 1
  print(i)
  
  
}

# LIBRA

write.table(tab_boot_P50_P75_libra, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/tab_boot_P50_P75_libra.txt",sep="\t",  row.names=F)
write.table(tab_boot_P75_all_libra, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/tab_boot_P75_all_libra.txt",sep="\t",  row.names=F)
write.table(tab_boot_P50_all_libra, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/tab_boot_P50_all_libra.txt",sep="\t",  row.names=F)
write.table(tab_boot_P50_P75_libra_e2e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/tab_boot_P50_P75_libra_e2e4.txt",sep="\t",  row.names=F)
write.table(tab_boot_P50_P75_libra_e4e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/tab_boot_P50_P75_libra_e4e4.txt",sep="\t",  row.names=F)

write.table(sd_boot_P50_P75_libra, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/sd_boot_P50_P75_libra.txt",sep="\t",  row.names=F)
write.table(sd_boot_P75_all_libra, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/sd_boot_P75_all_libra.txt",sep="\t",  row.names=F)
write.table(sd_boot_P50_all_libra, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/sd_boot_P50_all_libra.txt",sep="\t",  row.names=F)
write.table(sd_boot_P50_P75_libra_e2e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/sd_boot_P50_P75_libra_e2e4.txt",sep="\t",  row.names=F)
write.table(sd_boot_P50_P75_libra_e4e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/sd_boot_P50_P75_libra_e4e4.txt",sep="\t",  row.names=F)

write.table(p_boot_P50_P75_libra, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/p_boot_P50_P75_libra.txt",sep="\t",  row.names=F)
write.table(p_boot_P75_all_libra, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/p_boot_P75_all_libra.txt",sep="\t",  row.names=F)
write.table(p_boot_P50_all_libra, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/p_boot_P50_all_libra.txt",sep="\t",  row.names=F)
write.table(p_boot_P50_P75_libra_e2e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/p_boot_P50_P75_libra_e2e4.txt",sep="\t",  row.names=F)
write.table(p_boot_P50_P75_libra_e4e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/p_boot_P50_P75_libra_e4e4.txt",sep="\t",  row.names=F)

# Compo

write.table(tab_boot_P50_P75_compo, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/tab_boot_P50_P75_compo.txt",sep="\t",  row.names=F)
write.table(tab_boot_P75_all_compo, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/tab_boot_P75_all_compo.txt",sep="\t",  row.names=F)
write.table(tab_boot_P50_all_compo, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/tab_boot_P50_all_compo.txt",sep="\t",  row.names=F)
write.table(tab_boot_P50_P75_compo_e2e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/tab_boot_P50_P75_compo_e2e4.txt",sep="\t",  row.names=F)
write.table(tab_boot_P50_P75_compo_e4e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/tab_boot_P50_P75_compo_e4e4.txt",sep="\t",  row.names=F)

write.table(sd_boot_P50_P75_compo, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/sd_boot_P50_P75_compo.txt",sep="\t",  row.names=F)
write.table(sd_boot_P75_all_compo, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/sd_boot_P75_all_compo.txt",sep="\t",  row.names=F)
write.table(sd_boot_P50_all_compo, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/sd_boot_P50_all_compo.txt",sep="\t",  row.names=F)
write.table(sd_boot_P50_P75_compo_e2e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/sd_boot_P50_P75_compo_e2e4.txt",sep="\t",  row.names=F)
write.table(sd_boot_P50_P75_compo_e4e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/sd_boot_P50_P75_compo_e4e4.txt",sep="\t",  row.names=F)

write.table(p_boot_P50_P75_compo, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/p_boot_P50_P75_compo.txt",sep="\t",  row.names=F)
write.table(p_boot_P75_all_compo, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/p_boot_P75_all_compo.txt",sep="\t",  row.names=F)
write.table(p_boot_P50_all_compo, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/p_boot_P50_all_compo.txt",sep="\t",  row.names=F)
write.table(p_boot_P50_P75_compo_e2e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/p_boot_P50_P75_compo_e2e4.txt",sep="\t",  row.names=F)
write.table(p_boot_P50_P75_compo_e4e4, "C:/Users/mwagner4/Downloads/LIBRA score/3C/Boot_e4_MMSEsplines/p_boot_P50_P75_compo_e4e4.txt",sep="\t",  row.names=F)


