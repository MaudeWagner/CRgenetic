######################################################################################################
#####    R CODE FOR THE ASSOCIATION BETWEEN LIBRA & COGNITIVE RESILIENCE IN APOE-e4 CARRIERS     #####
######################################################################################################

##### PACKAGES ##### 
library(lcmm)
library(NormPsy)
library(lattice)
library(ggpubr)
library(Rmisc)
library(MASS)
library(sqldf)
library(dplyr)
library(ggplot2)
library(viridis)
library(lmtest)
library(moments)
library(dplyr, warn = FALSE)
library(patchwork)
library(survminer)
library(foreach) 
library(doParallel)
library(Cairo)

## This code applies to the dataset named tab_cog, which contains the following longitudinal 
## data (one row per visit):

tab_cog = read.table(file="~/tab_cog.txt", header=T, sep="\t", dec=".")
tab_V0  = subset(tab_cog, visit0 == 1)

# time: follow-up time (years)
# MMSE: Mini-Mental State Examination (points)
# BENTON: Benton Visual Retention Test (points)
# ISAAC_15: Isaacs’ Set Test (points) 
# TMT_A: Trail Making Test – Part A (points)
# TMT_B: Trail Making Test – Part B (points) 
# visit0: 0=follow-up visit, 1=baseline

## tab_cog also contains the time-invariant exposure and covariates:

# LIBRA: LIBRA risk score at baseline (points)
# sex: 0=men, 1=women
# apoe4: 0=ApoE4 non-carriers, 1=carriers 
# center: 0=Bordeaux, 1=Dijon, 2=Montpellier
# ageBL_74: age at baseline given in decades and centered around 74 (mean age at baseline)
# education: 0=no graduation, 1=high school graduate, 2=college graduate

library("NormPsy","lcmm","ggplot2","Cairo")

### Normalizing the MMSE using the NormPsy R package (doi:10.1159/000365637)
tab_cog$norm_MMSE = normMMSE(tab_cog$MMSE)

### Modeling trajectories of global cognition using a latent process mixed model for 
### multivariate longitudinal outcomes adjusted for demographics and ApoE4. We used the 
### multlcmm of the lcmm R package (doi:10.18637/jss.v078.i02).
model_cognition = multlcmm(norm_MMSE + BENTON + ISAAC_15 + TMT_A + TMT_B ~ visit0 + 
                     time + apoe4 + apoe4*time + sex + sex*time + ageBL_74 + ageBL_74*time + 
                     education + education*time + center + center*time, random =~ time, subject = 'projid', randomY = T, 
                     data = tab_cog, link = c('3-quant-splines','3-quant-splines', '3-quant-splines','3-quant-splines','3-quant-splines'))

### Calculating the individual adjusted cognitive slopes in ApoE4 carriers and non-carriers
beta_time       = model_cognition$best[9]    # fixed effect _time_   
beta_apoe4_time = model_cognition$best[10]   # fixed effect _apoe4*time_   
tempo1          = as.data.frame(model_cognition$predRE) 
tempo1$RE_slope = tempo1$time # random slope parameters = var _time_ in model_cognition$predRE
tempo1          = subset(tempo1, select = c("projid","RE_slope"))
tempo2          = subset(tab_cog, visit0 == 1, select = c("projid", "apoe4"))
tab_cr          = merge(tempo1, tempo2, by = c("projid"), all = T)
tab_cr$adj_slope = ifelse(tab_cr$apoe4 == 1, beta_time + beta_apoe4_time + tab_cr$RE_slope, beta_time + tab_cr$RE_slope)

### Superimposing distributions of the adjusted cognitive slopes in e4-carriers and non-carriers 
### with thresholds 50th vs. 75th percentiles depicted using dashed lines (cf. Figure 1)
P50 = quantile(tab_cr$cor_slope, probs = c(0.50))   
P75 = quantile(tab_cr$cor_slope, probs = c(0.75))     
ggplot(tab_cr, aes(x = adj_slope)) + 
  labs(y = "Frequency", x = "Individual adjusted slopes of global cognition") +
  geom_histogram(aes(color = apoe4, fill = apoe4), colour = "grey20", alpha = 0.7) +
  scale_fill_manual(values = c("#238A8DFF","#DCE319FF"), labels=c("ApoE-\u03B54 non-carriers","ApoE-\u03B54 carriers")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = c(P50, P75), linetype="dashed")


### Assigning CR_gen status to e4-carriers based on their adjusted cognitive slopes estimated
### in the overall population (<50th percentile [non-resilient] vs. >=75th perc. [resilient])
tab_cr$CR_gen = NA
tab_cr$CR_gen = ifelse(tab_cr$adj_slope >= P75, 1, tab_cr$CR_gen)
tab_cr$CR_gen = as.factor(ifelse(tab_cr$adj_slope < P50, 0, tab_cr$CR_gen))

### Fitting association between LIBRA (reverse scale) and CR_gen using a logistic regression model 
tab_libra = merge(tab_V0, tab_cr, by = c("projid"))
tab_libra$LIBRA_reverse = tab_libra$LIBRA*(-1)
model_libra  = glm(CR_gen ~ LIBRA_reverse + ageBL_74 + sex + education + center, family="binomial", data=tab_libra)


##################################################
###   PARAMETRIC BOOTSTRAPS, 1000 REPLICATES   ###
##################################################
m = mtemp = model_cognition
nboot = 1000
mu    = as.matrix(estimates(m)) 
Sigma = as.matrix(VarCov(m)) 
boot  = as.matrix(mvrnorm(nboot, mu, Sigma))
beta_time        = model_cognition$best[9]
beta_apoe4_time  = model_cognition$best[10]

for (i in 1:nboot){ 
  
  # Boot
  mtemp$best = boot[i,]
  RE_boot  = predictRE(mtemp, newdata = tab_cog) 
  tab_apoe = subset(tab_cog, select = c("projid", "apoe4"))
  tempo    = as.data.frame(tab_apoe %>% group_by(projid) %>% filter(row_number(projid) == 1))
  tempo1 = as.data.frame(RE_boot)
  tempo1$RE_int   = tempo1$intercept
  tempo1$RE_slope = tempo1$time_y
  tempo1 = subset(tempo1, select = c("projid","RE_int","RE_slope"))
  tab_RE0 = merge(tempo1, tempo, by = c("projid"), all = T)
  
  # CR slope <P50 vs. P>=75 (25% slowest)
  tab_RE0$CR_gen = ifelse(tab_RE0$APOE4 == 1, tab_RE0$RE_slope + beta_time + beta_apoe4_time, tab_RE0$RE_slope + beta_time)
  P50 = quantile(tab_RE0$CR_gen, probs = c(0.50))
  P75 = quantile(tab_RE0$CR_gen, probs = c(0.75))
  tab_RE0$CR_gen_P50_P75 = NA
  tab_RE0$CR_gen_P50_P75 = ifelse(tab_RE0$CR_gen >= P75, 1, tab_RE0$CR_gen_P50_P75)
  tab_RE0$CR_gen_P50_P75 = as.factor(ifelse(tab_RE0$CR_gen < P50, 0, tab_RE0$CR_gen_P50_P75))

  # merging  
  tab_CR = subset(tab_RE0, select= c("projid", "CR_gen", "CR_gen_P50_P75"))    
  tab_CR_slope = merge(tab_V0, tab_CR, by = c("projid"), all = T)  

  
  ## Association with reversed LIBRA ## 
  tab_libra = subset(tab_CR_slope, apoe4 == 1)
  tab_libra$LIBRA_reverse = tab_libra$LIBRA*(-1) 
  model_libra = glm(CR_slope_P50_P75 ~ LIBRA_reverse + ageBL_74 + sex + education + center, data=tab_libra, family = "binomial")
  b = as.data.frame(coef(summary(model_libra))[,'Estimate'])
  se = as.data.frame(coef(summary(model_libra))[,'Std. Error'])
  
  if (i==1){ 
    
    b_boot_P50_P75_libra  = b
    se_boot_P50_P75_libra = se
    
  } else if (i>1){ 
    
    b_boot_P50_P75_libra  = cbind(b_boot_P50_P75_libra, b) 
    se_boot_P50_P75_libra = cbind(se_boot_P50_P75_libra, se)
  
  } 
    
  i = i + 1
  
}

write.table(b_boot_P50_P75_libra, "~/b_boot_P50_P75_libra.txt", sep="\t",  row.names=F)
write.table(se_boot_P50_P75_libra, "~/se_boot_P50_P75_libra.txt", sep="\t",  row.names=F)

######################################################
###   Total Variance  and 95% confidence intervals ###
######################################################
mean_b    = rowSums(b_boot_P50_P75_libra[2,])/nboot
V_within  = rowSums(se_boot_P50_P75_libra[2,]**2)/nboot
V_between = (rowSums(b_boot_P50_P75_libra[2,]-mean_b)**2)/(nboot-1)
V_tot     = V_within + V_between + V_between/nboot

# calculation of 95%CI based on original beta and SE_boot
binf = coef(model_libra)[2] - 1.96*sqrt(V_tot)
bsup = coef(model_libra)[2] + 1.96*sqrt(V_tot)


##########################
###    Forest Plot     ###
##########################

c = round(exp(coef(model_libra)[2]),3) # boot
d = round(exp(binf),3)
e = round(exp(bsup),3)
f = "("
g = ","
h = ")"
space = " "
i = paste(c,space,f,d,g,space,e,h, sep="")
j = paste(d,g,space,e, sep="")
#
fig1a = data.frame(treatmentgroup = "LIBRA score (reversed scale)",
                   rr     = c,
                   low_ci = d,
                   up_ci  = e,
                   RR_ci = i,
                   ci = j,
                   X = "COPD",
                   no = 1)

##### LIBRA score alone #####

forest <-  ggplot(
  data = fig1a,
  aes(x = treatmentgroup, y = rr, ymin = low_ci, ymax = up_ci)) +
  geom_pointrange(aes(col = treatmentgroup), size=0.4) +
  geom_hline(yintercept = 1, colour = "grey", lty = 2, size=0.7) +
  xlab("") +
  ylab("Odd Ratio (95%CI)") +
  geom_errorbar(aes(ymin = low_ci, ymax = up_ci, col = treatmentgroup), width = 0, cex = 1.5) +
  facet_wrap(~X, strip.position = "top", nrow = 9, scales = "free_y") +
  theme_classic() +
  theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
    strip.text.y = element_text(face = "bold", size = 10),
    panel.grid.major.y = element_line(colour = col_grid, size = 0.5),
    strip.text = element_text(color = "white"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none",
    axis.text = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 10),
    plot.title = element_text("none")) +
  coord_flip()

dat_table <- fig1a %>%
  select(treatmentgroup, X, RR_ci) %>%
  tidyr::pivot_longer(c(RR_ci), names_to = "stat") %>%
  mutate(stat = factor(stat, levels = c("RR_ci")))

table_base <- ggplot(dat_table, aes(stat, treatmentgroup, label = value)) +
  geom_text(size = 3.5) +
  scale_x_discrete(position = "top", labels = c("Odd Ratio (95%CI)")) +
  facet_wrap(~X, strip.position = "top", ncol = 1, scales = "free_y",
             labeller = labeller(X = c(Cancer = "", COPD = ""))) +
  labs(y = NULL, x = NULL) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title = element_text(face = "bold"),
  )

forest + table_base + plot_layout(widths = c(15, 8))




