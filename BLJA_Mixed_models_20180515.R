# Statistical analysis of Blue Jay alarm call results
# Last updated March 21, 2018

# Use 'BLJA_encoding_results_201801018.csv'
# AND Use 'BLJA_encoding_results_201801018_FeedingRate.csv'

# Mixed effects models to test for significance of factors on the following dependent variables:
#**Results are significant.
# 1. **Ave number of calls (as the difference from Pre phase)
      # Significant effect Stimuli
      # CONCLUSION: SIGNIFICANCE between GHOW and AMRO and GHOW and EASO
# 2. **Ave number of elements/call (as the difference from Pre phase)
      # Significant effect of Stimuli
      # CONCLUSION: Significant difference between GHOW and EASO, and a trend of difference between GHOW and AMRO
# 3. **Ave duration of call (as the difference from Pre phase)
      # Significant effect of Stimuli
      # CONCLUSION: Trend of a difference between GHOW and AMRO
# 4. Duty cycle = ratio of total duration of elements to total duration of call (as the difference from Pre phase)
      # Data are sqrt transformed for normality
      # CONCLUSION: No significance
# 5. Ave feeding rate (not as difference, but for each phase; data are rank-transformed)
      # Significant effect of Stimuli and Phase (without interaction)
      # CONCLUSION: There is a significant difference between EASO and AMRO and between Pre and Present phases
# Each dependent variable (#1-5) is tested with the following steps:
    #     a. Normality
    #     b. Interaction
    #     c. Interaction Fixed Effects Bootstrap
    #     d. Fixed Effect - Stimuli
    #     e. Fixed Effect - Phase
    #     f. Stimuli Fixed Effects Bootstrap
    #     g. Phase Fixed Effects Bootstrap
    #     h. Tukey Bootstrap

# Chi-squared test for significance
# 7. Ave closest distance of approach
# Conclusion: No significant differences

# 8. Correlation Matrix

###############
# System Set-up
################
# Janelle used:
# R 3.4.3 GUI 1.70 El Capitan build (7463)

# Alexis used:
# R

# Any differences between us are noted.

##########
# Set-up #
##########

library(Matrix)
library(lattice)
library(lme4)
library(plotrix)
library(ggplot2)
library(multcomp)
library(nortest)
library(pbkrtest)
library(emmeans)
library(lmerTest)

# Janelle
setwd("/Users/jlm394/Desktop/R_working/NSFAlarm")
# Alexis
# setwd("/Users/alexisbillings/Desktop/PROJECTS/Blue_jay_project/From_Janelle_1-19-18")

encd_blja <- read.csv('BLJA_encoding_results_20180118.csv', header = T)


#Set reference level of AMRO
encd_blja$Stimuli<-relevel(encd_blja$Stimuli, ref = "AMRO")
#Make Site a Factor
encd_blja$Site<-as.factor(encd_blja$Site) 


str(encd_blja)
  # 'data.frame':	36 obs. of  45 variables:
  # $ X                   : int  1 2 3 4 5 6 7 8 9 10 ...
  # $ Date                : int  20141119 20141119 20141203 20141203 20141205 20141222 20150122 20150122 20150123 20150204 ...
  # $ Site                : Factor w/ 5 levels "1","2","3","4",..: 3 3 3 3 1 2 5 5 5 5 ...
  # $ Stimuli             : Factor w/ 3 levels "AMRO","EASO",..: 2 3 1 3 1 3 1 2 3 2 ...
  # $ Phase               : Factor w/ 2 levels "Post","Present": 1 1 1 1 1 1 1 1 1 1 ...
  # $ Phase_dur           : num  287 292 286 296 312 ...
  # $ Total_calls         : int  50 0 28 10 10 35 4 0 15 23 ...
  # $ Total_ditonal       : int  61 0 41 12 0 44 5 0 24 28 ...
  # $ Avedur_ditonal      : num  0.387 0 0.296 0.398 0 ...
  # $ Total_monotonal     : int  NA NA NA NA NA NA NA NA NA NA ...
  # $ Total_bell          : int  NA NA NA NA NA NA NA NA NA NA ...
  # $ Total_alert         : int  NA NA NA NA NA NA NA NA NA NA ...
  # $ Total_redtail       : int  NA NA NA 3 NA NA NA NA NA NA ...
  # $ Total_trains        : int  10 NA 10 2 2 8 1 NA 3 4 ...
  # $ Avenum_elements     : num  1.22 0 1.46 1.2 0 ...
  # $ Avedur_trains       : num  1.042 NA 0.916 1.173 NA ...
  # $ Avedur_silence      : num  0.243 NA 0.221 0.275 NA ...
  # $ Avedur_delta        : num  0.634 NA 0.489 0.739 NA ...
  # $ Pre.Phase_dur       : num  310 331 326 326 315 ...
  # $ Pre.Total_calls     : int  0 0 91 2 2 0 19 7 3 0 ...
  # $ Pre.Total_ditonal   : int  0 0 99 1 6 0 23 15 4 0 ...
  # $ Pre.Avedur_ditonal  : num  0 0 0.364 0.351 0.353 ...
  # $ Pre.Total_monotonal : int  NA NA NA 8 NA NA NA NA NA NA ...
  # $ Pre.Total_bell      : int  NA NA NA NA NA NA NA NA NA NA ...
  # $ Pre.Total_alert     : int  NA NA NA NA NA NA NA NA NA NA ...
  # $ Pre.Total_redtail   : int  NA NA NA 2 NA NA NA NA NA NA ...
  # $ Pre.Total_trains    : int  NA NA 7 NA NA NA 3 5 1 NA ...
  # $ Pre.Avenum_elements : num  NA NA 1.09 0.5 3 ...
  # $ Pre.Avedur_trains   : num  NA NA 1.02 NA NA ...
  # $ Pre.Avedur_silence  : num  NA NA 0.263 NA NA ...
  # $ Pre.Avedur_delta    : num  NA NA 0.654 NA NA ...
  # $ Approach            : int  NA NA NA NA NA NA NA NA NA NA ...
  # $ Video_dur           : int  288 292 285 294 300 285 309 307 322 327 ...
  # $ Feed_dur            : int  90 1 51 27 153 0 2 0 0 1 ...
  # $ Pre.Approach        : logi  NA NA NA NA NA NA ...
  # $ Pre.Video_dur       : int  293 313 316 306 300 309 302 317 309 300 ...
  # $ Pre.Feed_dur        : int  0 128 90 26 90 0 0 10 21 10 ...
  # $ Corr.Total_calls    : int  50 0 -63 8 8 35 -15 -7 12 23 ...
  # $ Corr.Avenum_elements: num  NA NA 0.376 0.7 -3 ...
  # $ Corr.Avedur_ditonal : num  0.3869 0 -0.068 0.0465 -0.3532 ...
  # $ Feed_rate           : num  18.75 0.205 10.737 5.51 30.6 ...
  # $ Pre.Feed_rate       : num  0 24.5 17.1 5.1 18 ...
  # $ Corr.Feed_rate      : num  18.75 -24.331 -6.352 0.412 12.6 ...
  # $ Corr.Avgnum_calls   : num  10.44 0 -13.22 1.62 1.54 ...
  # $ DutyCycle           : num  0.0821 0 0.0424 0.0161 0 ...

############
# N values #
############
library(dplyr)
# Number of sites
nsite <- count(encd_blja, Site)
nsite
# Site      n
# <fct> <int>
#   1 1         8
# 2 2         6
# 3 3         8
# 4 4         6
# 5 6         8

nstimuli <- count(encd_blja, Stimuli)
nstimuli
# Stimuli     n
# <fct>   <int>
# 1 AMRO       12
# 2 EASO       12
# 3 GHOW       12

aveelem <- select(encd_blja, Stimuli, Avenum_elements) %>%
  group_by(Stimuli) %>%
  summarise(avg = mean(Avenum_elements))
# Stimuli   avg
# <fct>   <dbl>
# 1 AMRO    1.02 
# 2 EASO    0.971
# 3 GHOW    1.31 

###########################################################################################
# 1. Corr.Avgnum_calls     # Ave number of calls (as the difference from Pre phase)
##########################################################################################

#Get rid of NAs 
encd_blja = encd_blja[!is.na(encd_blja$Corr.Avgnum_calls),]

#CHECK FOR ZERO-INFLATED (likely not relevant)
(length(encd_blja$Corr.Avgnum_calls[encd_blja$Corr.Avgnum_calls==0])/length(encd_blja$Corr.Avgnum_calls))*100 
#[1] 8.333333
# Not zero-inflated

######################### 
# 1.a. Normality        # 
######################### 

Corr.Avgnum_calls.full<-lmer(Corr.Avgnum_calls ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
Corr.Avgnum_calls.attributes<-attributes(Corr.Avgnum_calls.full)
qqnorm(resid(Corr.Avgnum_calls.full))
# Looks good

sf.test(resid(Corr.Avgnum_calls.full))

## If p-value >0.05, then data did not violate the normality test and can assume normality
# 	Shapiro-Francia normality test

# data:  resid(Corr.Avgnum_calls.full)
# W = 0.97737, p-value = 0.5653
####CONCLUSION: Normally distributed


#########################
# 1.b. INTERACTION      #
#########################
Corr.Avgnum_calls.int<-lmer(Corr.Avgnum_calls ~ Stimuli * Phase + (1|Site),  data=encd_blja, REML="F")
#REML off, restricted maximum likelihood. You just want maximum likelihood. 
# * interaction between Stimuli and Phase

Corr.Avgnum_calls.add<-lmer(Corr.Avgnum_calls ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
# Stimuli and Phase are main effects

anova(Corr.Avgnum_calls.int, Corr.Avgnum_calls.add)
# Data: encd_blja
# Models:
#   Corr.Avgnum_calls.add: Corr.Avgnum_calls ~ Stimuli + Phase + (1 | Site)
# Corr.Avgnum_calls.int: Corr.Avgnum_calls ~ Stimuli * Phase + (1 | Site)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Corr.Avgnum_calls.add  6 239.46 248.97 -113.73   227.46                         
# Corr.Avgnum_calls.int  8 242.88 255.55 -113.44   226.88 0.5796      2     0.7484
#####CONCLUSION: No significance.

############################################
# 1.c. Interaction Fixed Effects Bootstrap # 
############################################
# Bootstrapping fixed effects models, comparing log liklelihood to null model
# Bootstrapping is resampling data, to look at how likely there is a significance with "different data"
# more robust than just running the mixed effects model alone
# p-value should be similar to confirm the results of the mixed effects model

#### Data are normally distributed; therefore, bootstrapping is not necessary.

###############################
# 1.d. Fixed Effect - Stimuli #
###############################

Corr.Avgnum_calls.full.1<-lmer(Corr.Avgnum_calls ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")

Corr.Avgnum_calls.null.1<-lmer(Corr.Avgnum_calls ~ Phase + (1|Site),  data=encd_blja, REML="F")

Corr.Avgnum_calls.full.1
anova(Corr.Avgnum_calls.full.1, Corr.Avgnum_calls.null.1)	


# Data: encd_blja
# Models:
#   Corr.Avgnum_calls.null.1: Corr.Avgnum_calls ~ Phase + (1 | Site)
# Corr.Avgnum_calls.full.1: Corr.Avgnum_calls ~ Stimuli + Phase + (1 | Site)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# Corr.Avgnum_calls.null.1  4 243.22 249.55 -117.61   235.22                           
# Corr.Avgnum_calls.full.1  6 239.46 248.97 -113.73   227.46 7.7516      2    0.02074 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### CONCLUSION: SIGNIFICANT EFFECT OF STIMULI

################################
# 1.e. Fixed Effect - Phase #
################################

Corr.Avgnum_calls.full.2<-lmer(Corr.Avgnum_calls ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")

Corr.Avgnum_calls.null.2<-lmer(Corr.Avgnum_calls ~ Stimuli + (1|Site),  data=encd_blja, REML="F")

Corr.Avgnum_calls.full.2
anova(Corr.Avgnum_calls.full.2, Corr.Avgnum_calls.null.2)  

#Data: encd_blja
# Models:
#   Corr.Avgnum_calls.null.2: Corr.Avgnum_calls ~ Stimuli + (1 | Site)
# Corr.Avgnum_calls.full.2: Corr.Avgnum_calls ~ Stimuli + Phase + (1 | Site)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# Corr.Avgnum_calls.null.2  5 237.82 245.74 -113.91   227.82                         
# Corr.Avgnum_calls.full.2  6 239.46 248.97 -113.73   227.46 0.3559      1     0.5508

#####CONCLUSION: Not significant


#######################################
# 1.f. Stimuli Fixed Effects Bootstrap# 
#######################################

# not necessary, since data are normally distributed


########################################
# 1.g. Phase Fixed Effects Bootstrap# 
########################################

# not necessary, since data are normally distributed


############################
# 1.h. Tukey Bootstrapping # 
############################
# Run on pairwise comparisons for significant Stimuli only
#Model from above with Phase
Corr.Avgnum_calls.full.3<-lmer(Corr.Avgnum_calls ~ Stimuli + Phase + (1|Site), REML="F", data=encd_blja)

lsmeans(Corr.Avgnum_calls.full.3, list(pairwise ~ Stimuli), adjust = "tukey") 
# $`lsmeans of Stimuli`
# Stimuli    lsmean       SE df  lower.CL upper.CL
# AMRO    -2.306112 2.231769 15 -7.062960 2.450737
# EASO     3.259327 2.231769 15 -1.497521 8.016175
# GHOW     3.649129 2.231769 15 -1.107719 8.405977

# Results are averaged over the levels of: Phase 
# Confidence level used: 0.95 

# $`pairwise differences of contrast`
# contrast      estimate       SE    df t.ratio p.value
# AMRO - EASO -5.5654388 2.302724 35.12  -2.417  0.0534
# AMRO - GHOW -5.9552405 2.302724 35.12  -2.586  0.0364
# EASO - GHOW -0.3898018 2.302724 35.12  -0.169  0.9843

# Results are averaged over the levels of: Phase 
# P value adjustment: tukey method for comparing a family of 3 estimates

#### CONCLUSION: SIGNIFICANCE between GHOW and AMRO and GHOW and EASO

####################################################################################################
# 2. Corr.Avenum_elements     # Number of elements in jeer calls (as the difference from Pre phase)
####################################################################################################

#Get rid of NAs
encd_blja = encd_blja[!is.na(encd_blja$Corr.Avenum_elements),]

#CHECK FOR ZERO-INFLATED (likely not relevant)
(length(encd_blja$Corr.Avenum_elements[encd_blja$Corr.Avenum_elements==0])/length(encd_blja$Corr.Avenum_elements))*100 
#0% ARE 0'S

######################### 
# 2.a. Normality        # 
######################### 

Corr.Avenum_elements.full<-lmer(Corr.Avenum_elements ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
Corr.Avenum_elements.attributes<-attributes(Corr.Avenum_elements.full)
qqnorm(resid(Corr.Avenum_elements.full))
# looks good

sf.test(resid(Corr.Avenum_elements.full))

## If p-value >0.05, then data did not violate the normality test and can assume normality
# Shapiro-Francia normality test
# 
# data:  resid(Corr.Avenum_elements.full)
# W = 0.92689, p-value = 0.06308      
####CONCLUSION: NORMALLY DISTRIBUTED                  


#########################
# 2.b. INTERACTION      #
#########################
Corr.Avenum_elements.int<-lmer(Corr.Avenum_elements ~ Stimuli * Phase + (1|Site),  data=encd_blja, REML="F")
#REML off, restricted maximum likelihood. You just want maximum likelihood. 
# * interaction

Corr.Avenum_elements.add<-lmer(Corr.Avenum_elements ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
# Stimuli and Phase are main effects

anova(Corr.Avenum_elements.int, Corr.Avenum_elements.add)
# Data: encd_blja
# Models:
# ..1: Corr.Avenum_elements ~ Stimuli + Phase + (1 | Site)
# object: Corr.Avenum_elements ~ Stimuli * Phase + (1 | Site)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# ..1     6 95.876 103.42 -41.938   83.876                         
# object  8 98.053 108.12 -41.026   82.053 1.8233      2     0.4019
#####Conclusion: No significance

###########################################
# 2.c. Interaction Fixed Effects Bootstrap# 
###########################################
# Data are normally distributed, so not necessary


###############################
# 2.d. Fixed Effect - Stimuli #
###############################

Corr.Avenum_elements.full.1<-lmer(Corr.Avenum_elements ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")

Corr.Avenum_elements.null.1<-lmer(Corr.Avenum_elements ~ Phase + (1|Site),  data=encd_blja, REML="F")

Corr.Avenum_elements.full.1
anova(Corr.Avenum_elements.full.1, Corr.Avenum_elements.null.1)	

# Data: encd_blja
# Models:
# ..1: Corr.Avenum_elements ~ Phase + (1 | Site)
# object: Corr.Avenum_elements ~ Stimuli + Phase + (1 | Site)
# Df     AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# ..1     4 102.883 107.92 -47.442   94.883                            
# object  6  95.876 103.42 -41.938   83.876 11.007      2   0.004072 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#### CONCLUSION: SIGNIFICANT effect of Stimuli

################################
# 2.e. Fixed Effect - Phase #
################################

Corr.Avenum_elements.full.2<-lmer(Corr.Avenum_elements ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")

Corr.Avenum_elements.null.2<-lmer(Corr.Avenum_elements ~ Stimuli + (1|Site),  data=encd_blja, REML="F")

Corr.Avenum_elements.full.2
anova(Corr.Avenum_elements.full.2, Corr.Avenum_elements.null.2)  

# Data: encd_blja
# Models:
# ..1: Corr.Avenum_elements ~ Stimuli + (1 | Site)
# object: Corr.Avenum_elements ~ Stimuli + Phase + (1 | Site)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# ..1     5 94.612 100.90 -42.306   84.612                         
# object  6 95.876 103.42 -41.938   83.876 0.7359      1      0.391
#### CONCLUSION: not significant


########################################
# 2.f. Stimuli Fixed Effects Bootstrap # 
########################################
# Not necessary because data are normally distributed.


########################################
# 2.g. Phase Fixed Effects Bootstrap# 
########################################
# Not necessary because data are normally distributed.


############################
# 2.h. Tukey Bootstrapping # 
############################
### We should only look at Stimuli, not interaction?
#Model from above with Phase
Corr.Avenum_elements.full.4<-lmer(Corr.Avenum_elements ~ Stimuli + Phase + (1|Site), REML="F", data=encd_blja)

lsmeans(Corr.Avenum_elements.full.4, list(pairwise ~ Stimuli), adjust = "tukey")
# $`lsmeans of Stimuli`
# Stimuli     lsmean        SE    df   lower.CL   upper.CL
# AMRO    -0.5362012 0.4174031 23.58 -1.3984992  0.3260967
# EASO    -1.3585852 0.4831401 23.58 -2.3566868 -0.3604836
# GHOW     0.8603793 0.4831401 23.58 -0.1377224  1.8584809

# Results are averaged over the levels of: Phase 
# Confidence level used: 0.95 

# $`pairwise differences of contrast`
# contrast      estimate        SE    df t.ratio p.value
# AMRO - EASO  0.8223839 0.6384745 27.43   1.288  0.4137
# AMRO - GHOW -1.3965805 0.6384745 27.43  -2.187  0.0913
# EASO - GHOW -2.2189644 0.6862110 29.47  -3.234  0.0082

# Results are averaged over the levels of: Phase 
# P value adjustment: tukey method for comparing a family of 3 estimates 

#### CONCLUSION: Significant difference between GHOW and EASO, and a trend of difference between GHOW and AMRO


####################################################################################################
# 3. Corr.Avedur_ditonal     #  Ave duration of call (as the difference from Pre phase)
####################################################################################################
#Get rid of NAs
encd_blja = encd_blja[!is.na(encd_blja$Corr.Avedur_ditonal),]

#CHECK FOR ZERO-INFLATED (likely not relevant)
(length(encd_blja$Corr.Avedur_ditonal[encd_blja$Corr.Avedur_ditonal==0])/length(encd_blja$Corr.Avedur_ditonal))*100 
# [1] 5.555556

######################### 
# 3.a. Normality        # 
######################### 
Corr.Avedur_ditonal.full<-lmer(Corr.Avedur_ditonal ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
Corr.Avedur_ditonal.attributes<-attributes(Corr.Avedur_ditonal.full)
qqnorm(resid(Corr.Avedur_ditonal.full))
# Looks good

sf.test(resid(Corr.Avedur_ditonal.full))

## If p-value >0.05, then data did not violate the normality test and can assume normality
# Shapiro-Francia normality test
# 
# data:  resid(Corr.Avedur_ditonal.full)
# W = 0.95715, p-value = 0.2882
#####Conclusion: Normally distributed



##############################
# 3.b. INTERACTION           #
##############################
Corr.Avedur_ditonal.int<-lmer(Corr.Avedur_ditonal ~ Stimuli * Phase + (1|Site),  data=encd_blja, REML="F")
#REML off, restricted maximum likelihood. You just want maximum likelihood. 
# * interaction

Corr.Avedur_ditonal.add<-lmer(Corr.Avedur_ditonal ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
# Stimuli and Phase are main effects

anova(Corr.Avedur_ditonal.int, Corr.Avedur_ditonal.add)
# Data: encd_blja
# Models:
#   Corr.Avedur_ditonal.add: Corr.Avedur_ditonal ~ Stimuli + Phase + (1 | Site)
# Corr.Avedur_ditonal.int: Corr.Avedur_ditonal ~ Stimuli * Phase + (1 | Site)
# Df     AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# Corr.Avedur_ditonal.add  6 -3.9845 3.5641 7.9922  -15.985                         
# Corr.Avedur_ditonal.int  8 -1.2459 8.8189 8.6229  -17.246 1.2614      2     0.5322
#####Conclusion: No significance

###########################################
# 3.c. Interaction Fixed Effects Bootstrap# 
###########################################
# Not necessary because data are normally distributed.


###############################
# 3.d. Fixed Effect - Stimuli #
###############################

Corr.Avedur_ditonal.full.1<-lmer(Corr.Avedur_ditonal ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")

Corr.Avedur_ditonal.null.1<-lmer(Corr.Avedur_ditonal ~ Phase + (1|Site),  data=encd_blja, REML="F")

Corr.Avedur_ditonal.full.1
anova(Corr.Avedur_ditonal.full.1, Corr.Avedur_ditonal.null.1)	

#### Original results
# Data: encd_blja
# Models:
#   Corr.Avedur_ditonal.null.1: Corr.Avedur_ditonal ~ Phase + (1 | Site)
# Corr.Avedur_ditonal.full.1: Corr.Avedur_ditonal ~ Stimuli + Phase + (1 | Site)
# Df     AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# Corr.Avedur_ditonal.null.1  4 -1.3008 3.7316 4.6504  -9.3008                           
# Corr.Avedur_ditonal.full.1  6 -3.9845 3.5641 7.9922 -15.9845 6.6837      2    0.03537 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#####CONCLUSION: Significant effect of Stimuli

################################
# 3.e. Fixed Effect - Phase #
################################

Corr.Avedur_ditonal.full.2<-lmer(Corr.Avedur_ditonal ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")

Corr.Avedur_ditonal.null.2<-lmer(Corr.Avedur_ditonal ~ Stimuli + (1|Site),  data=encd_blja, REML="F")

Corr.Avedur_ditonal.full.2
anova(Corr.Avedur_ditonal.full.2, Corr.Avedur_ditonal.null.2)  

# Data: encd_blja
# Models:
#   Corr.Avedur_ditonal.null.2: Corr.Avedur_ditonal ~ Stimuli + (1 | Site)
# Corr.Avedur_ditonal.full.2: Corr.Avedur_ditonal ~ Stimuli + Phase + (1 | Site)
# Df     AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# Corr.Avedur_ditonal.null.2  5 -5.8510 0.4395 7.9255  -15.851                         
# Corr.Avedur_ditonal.full.2  6 -3.9845 3.5641 7.9922  -15.985 0.1335      1     0.7148

#####Conclusion: Not significant


#######################################
# 3.f. Stimuli Fixed Effects Bootstrap# 
#######################################
# Not necessary because data are normally distributed


########################################
# 3.g. Phase Fixed Effects Bootstrap# 
########################################
# Not necessary because data are normally distributed


############################
# 3.h. Tukey Bootstrapping # 
############################
#Model from above with Phase
Corr.Avedur_ditonal.full.4<-lmer(Corr.Avedur_ditonal ~ Stimuli + Phase + (1|Site), REML="F", data=encd_blja)

lsmeans(Corr.Avedur_ditonal.full.4, list(pairwise ~ Stimuli), adjust = "tukey")
# $`lsmeans of Stimuli`
# Stimuli      lsmean         SE    df   lower.CL   upper.CL
# AMRO    -0.08481802 0.06117005 23.58 -0.2111870 0.04155096
# EASO    -0.07742503 0.07080375 23.58 -0.2236959 0.06884586
# GHOW     0.12709279 0.07080375 23.58 -0.0191781 0.27336367

# Results are averaged over the levels of: Phase 
# Confidence level used: 0.95 

# $`pairwise differences of contrast`
# contrast        estimate         SE    df t.ratio p.value
# AMRO - EASO -0.007392991 0.09356787 27.43  -0.079  0.9966
# AMRO - GHOW -0.211910805 0.09356787 27.43  -2.265  0.0781
# EASO - GHOW -0.204517814 0.10056360 29.47  -2.034  0.1218

# Results are averaged over the levels of: Phase 
# P value adjustment: tukey method for comparing a family of 3 estimates

#### CONCLUSION: Trend of a difference between GHOW and AMRO

####################################################################################################
# 4. DutyCycle	 Duty cycle = (Avedur_ditonal * Total_ditonal)/ Phase_dur)
####################################################################################################  
  
#Get rid of NAs
encd_blja = encd_blja[!is.na(encd_blja$DutyCycle),]

#CHECK FOR ZERO-INFLATED 
(length(encd_blja$DutyCycle[encd_blja$DutyCycle==0])/length(encd_blja$DutyCycle))*100 
#[1] 15.38462
# Not zero-inflated
  
  ######################### 
  # 4.a. Normality        # 
  ######################### 
  
  DutyCycle.full<-lmer(DutyCycle ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
  DutyCycle.attributes<-attributes(DutyCycle.full)
  qqnorm(resid(DutyCycle.full))
  
  
  sf.test(resid(DutyCycle.full))
  
  ## If p-value >0.05, then data did not violate the normality test and can assume normality
  # Shapiro-Francia normality test
  
  # data:  resid(DutyCycle.full)
  # W = 0.7198, p-value = 3.595e-05
  ####Conclusion: THIS IS NOT NORMALLY DISTRIBUTED 
  hist(encd_blja$DutyCycle)
  
  # Need to corrected skewed data
  # Trying square root transformation
  encd_blja$Trans.DutyCycle <- sqrt(encd_blja$DutyCycle)
  
  Trans.DutyCycle.full<-lmer(Trans.DutyCycle ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
  Trans.DutyCycle.attributes<-attributes(Trans.DutyCycle.full)
  qqnorm(resid(Trans.DutyCycle.full))
  
  
  sf.test(resid(Trans.DutyCycle.full))
  # Shapiro-Francia normality test
  # 
  # data:  resid(Trans.DutyCycle.full)
  # W = 0.96359, p-value = 0.2367
  ####Conclusion: Sqrt transformed data is normally distributed!

  
  #########################
  # 4.b. INTERACTION      #
  #########################
  DutyCycle.int<-lmer(Trans.DutyCycle ~ Stimuli * Phase + (1|Site),  data=encd_blja, REML="F")
  #REML off, restricted maximum likelihood. You just want maximum likelihood. 
  # * interaction
  
  DutyCycle.add<-lmer(Trans.DutyCycle ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
  # Stimuli and Phase are main effects
  
  anova(DutyCycle.int, DutyCycle.add)
  # Data: encd_blja
  # Models:
  #   DutyCycle.add: Trans.DutyCycle ~ Stimuli + Phase + (1 | Site)
  # DutyCycle.int: Trans.DutyCycle ~ Stimuli * Phase + (1 | Site)
  # Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
  # DutyCycle.add  6 -48.102 -38.601 30.051  -60.102                         
  # DutyCycle.int  8 -47.782 -35.114 31.891  -63.782 3.6795      2     0.1589
  ##### CONCLUSION: Not significant
  
  ############################################
  # 4.c. Interaction Fixed Effects Bootstrap # 
  ############################################
  # Not necessary, data are normally distributed
  
  
  ###############################
  # 4.d. Fixed Effect - Stimuli #
  ###############################
  
  DutyCycle.full.1<-lmer(Trans.DutyCycle ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
  
  DutyCycle.null.1<-lmer(Trans.DutyCycle ~ Phase + (1|Site),  data=encd_blja, REML="F")
  
  DutyCycle.full.1
  anova(DutyCycle.full.1, DutyCycle.null.1)	
  
  # Data: encd_blja
  # Data: encd_blja
  # Models:
  #   DutyCycle.null.1: Trans.DutyCycle ~ Phase + (1 | Site)
  # DutyCycle.full.1: Trans.DutyCycle ~ Stimuli + Phase + (1 | Site)
  # Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
  # DutyCycle.null.1  4 -50.172 -43.838 29.086  -58.172                         
  # DutyCycle.full.1  6 -48.102 -38.601 30.051  -60.102 1.9302      2     0.3809
  #####Conclusion: Not significant
  
  ################################
  # 4.e. Fixed Effect - Phase #
  ################################
  #Because this is not significant, we can just look at Stimuli
  
  DutyCycle.full.2<-lmer(Trans.DutyCycle ~ Stimuli + Phase + (1|Site),  data=encd_blja, REML="F")
  
  DutyCycle.null.2<-lmer(Trans.DutyCycle ~ Stimuli + (1|Site),  data=encd_blja, REML="F")
  
  DutyCycle.full.2
  anova(DutyCycle.full.2, DutyCycle.null.2)  
  
  # Data: encd_blja
  # Data: encd_blja
  # Models:
  #   DutyCycle.null.2: Trans.DutyCycle ~ Stimuli + (1 | Site)
  # DutyCycle.full.2: Trans.DutyCycle ~ Stimuli + Phase + (1 | Site)
  # Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
  # DutyCycle.null.2  5 -50.067 -42.149 30.034  -60.067                         
  # DutyCycle.full.2  6 -48.102 -38.601 30.051  -60.102 0.0354      1     0.8507
  #####Conclusion: Not significant
  
  
  #######################################
  # 4.f. Stimuli Fixed Effects Bootstrap# 
  #######################################
  # Not necessary

  
  ########################################
  # 4.g. Phase Fixed Effects Bootstrap# 
  ########################################
  # Not necessary
  
  
  ############################
  # 4.h. Tukey Bootstrapping # 
  ############################
  # Not necessary, no significance.
 
  
####################################################################################################
# 5. Feed_rate				Average feeding rate (as the difference from Pre phase)
####################################################################################################
  
  # The results of Feed_rate are zero-inflated (50%), and when correcting for Pre phase, 
  # there are negative values. 
  # Therefore, it is inappropriate to apply a Poisson-distribution zero-inflated model.
  # To avoid the negative values, we are not correcting for Pre phase, but including Pre phase when testing for Phase as a factor.
  
  #Use 'BLJA_encoding_results_201801018_FeedingRate.csv'
  feed_blja <- read.csv('BLJA_encoding_results_20180118_FeedingRate.csv', header = T)
  
  #Set reference level of AMRO
  feed_blja$Stimuli<-relevel(feed_blja$Stimuli, ref = "AMRO")
  #Make Site a Factor
  feed_blja$Site<-as.factor(feed_blja$Site) 
  
  ggplot(feed_blja, aes(x=Phase, y=Feed_rate, fill=Stimuli)) + 
    geom_bar(stat="summary", fun.y="mean", position="dodge") +
    geom_errorbar(stat="summary", fun.data=mean_se, position="dodge")
  
  par(mfrow=c(2,3))
  
  # square root transformation
  Feed_rate.full <- lmer(sqrt(Feed_rate) ~ Stimuli * Phase + (1|Site),  data=feed_blja, REML="F")
  qqnorm(resid(Feed_rate.full))
  qqline(resid(Feed_rate.full))
  hist(resid(Feed_rate.full))
  plot(predict(Feed_rate.full), resid(Feed_rate.full))
  sf.test(resid(Feed_rate.full))
  # Shapiro-Francia normality test
  # 
  # data:  resid(Feed_rate.full)
  # W = 0.95849, p-value = 0.05703
  # Conclusion: normally distributed
  
  anova(Feed_rate.full)
  # Analysis of Variance Table of type III  with  Satterthwaite 
  # approximation for degrees of freedom
  # Sum Sq Mean Sq NumDF  DenDF F.value   Pr(>F)   
  # Stimuli       16.0186  8.0093     2 49.306  5.3399 0.007958 **
  #   Phase          5.5092  2.7546     2 49.003  1.8365 0.170174   
  # Stimuli:Phase  4.8870  1.2217     4 49.003  0.8146 0.522052   
 
  
  # rank transformation
  Feed_rate.full2 <- lmer(rank(Feed_rate) ~ Stimuli * Phase + (1|Site),  data=feed_blja, REML="F")
  qqnorm(resid(Feed_rate.full2))
  qqline(resid(Feed_rate.full2))
  hist(resid(Feed_rate.full2))
  plot(predict(Feed_rate.full2), resid(Feed_rate.full2))
  # sf.test(resid(Feed_rate.full2)) 
  # Shapiro-Francia normality test
  # 
  # data:  resid(Feed_rate.full2)
  # W = 0.97726, p-value = 0.3311
  # doesn't fail sf
  
  anova(Feed_rate.full2) 
  # Analysis of Variance Table of type III  with  Satterthwaite 
  # approximation for degrees of freedom
  # Sum Sq Mean Sq NumDF  DenDF F.value  Pr(>F)  
  # Stimuli       634.82  317.41     2 49.254  3.8321 0.02839 *
  #   Phase         607.00  303.50     2 48.988  3.6642 0.03289 *
  #   Stimuli:Phase 316.89   79.22     4 48.988  0.9565 0.43978  
  # CONCLUSION: no significant interaction, but significant effect of Phase and Stimuli
  
  Feed_rate.full2_reduced <- lmer(rank(Feed_rate) ~ Stimuli + Phase + (1|Site),  data=feed_blja, REML="F")
  anova(Feed_rate.full2_reduced)
  # Analysis of Variance Table of type III  with  Satterthwaite 
  # approximation for degrees of freedom
  # Sum Sq Mean Sq NumDF  DenDF F.value  Pr(>F)  
  # Stimuli  637.6   318.8     2 49.274  3.5701 0.03565 *
  #   Phase    607.0   303.5     2 48.987  3.3988 0.04147 *
  # CONCLUSION: Significant effect of Stimuli and Phase
  
  # pairwise comparisons
  emmeans(Feed_rate.full2_reduced, pairwise ~ Stimuli)
  # Stimuli   emmean       SE  df  lower.CL upper.CL
  # AMRO    31.15145 5.627345 8.3 18.256762 44.04613
  # EASO    22.64408 5.627345 8.3  9.749393 35.53876
  # GHOW    25.89776 5.627345 8.3 13.003078 38.79244
  # 
  # Results are averaged over the levels of: Phase 
  # Degrees-of-freedom method: kenward-roger 
  # Results are given on the rank (not the response) scale. 
  # Confidence level used: 0.95 
  # 
  # $contrasts
  # contrast     estimate       SE    df t.ratio p.value
  # AMRO - EASO  8.507369 3.357803 53.65   2.534  0.0373
  # AMRO - GHOW  5.253684 3.357803 53.65   1.565  0.2697
  # EASO - GHOW -3.253684 3.357803 53.65  -0.969  0.5995
  
  # Results are averaged over the levels of: Phase 
  # P value adjustment: tukey method for comparing a family of 3 estimates
  
  ####CONCLUSION: There is a significant difference between EASO and AMRO
  
  emmeans(Feed_rate.full2_reduced, pairwise ~ Phase)
  # $emmeans
  # Phase     emmean      SE   df  lower.CL upper.CL
  # Post    27.06443 5.61333 8.23 14.181726 39.94713
  # Pre     30.39776 5.61333 8.23 17.515059 43.28046
  # Present 22.23109 5.61333 8.23  9.348392 35.11380
  # 
  # Results are averaged over the levels of: Stimuli 
  # Degrees-of-freedom method: kenward-roger 
  # Results are given on the rank (not the response) scale. 
  # Confidence level used: 0.95 
  # 
  # $contrasts
  # contrast        estimate       SE    df t.ratio p.value
  # Post - Pre     -3.333333 3.286675 53.34  -1.014  0.5712
  # Post - Present  4.833333 3.286675 53.34   1.471  0.3130
  # Pre - Present   8.166667 3.286675 53.34   2.485  0.0420
  # 
  # Results are averaged over the levels of: Stimuli 
  # P value adjustment: tukey method for comparing a family of 3 estimates 
  
  ####CONCLUSION: There is a significant difference between Pre and Present
  
  
  # but there's no significant interaction, so this doesn't make sense to do.
  emmeans(Feed_rate.full2, pairwise ~ Phase + Stimuli)
  # Phase   Stimuli   emmean      SE    df  lower.CL upper.CL
  # Post    AMRO    30.08103 6.53697 14.64 16.118107 44.04395
  # Pre     AMRO    32.24769 6.53697 14.64 18.284773 46.21062
  # Present AMRO    31.08103 6.53697 14.64 17.118107 45.04395
  # Post    EASO    24.36943 6.53697 14.64 10.406512 38.33235
  # Pre     EASO    26.03610 6.53697 14.64 12.073178 39.99902
  # Present EASO    17.53610 6.53697 14.64  3.573178 31.49902
  # Post    GHOW    26.72523 6.53697 14.64 12.762309 40.68815
  # Pre     GHOW    32.89190 6.53697 14.64 18.928976 46.85482
  # Present GHOW    18.05856 6.53697 14.64  4.095643 32.02149
  # 
  # Degrees-of-freedom method: kenward-roger 
  # Results are given on the rank (not the response) scale. 
  # Confidence level used: 0.95 
  # 
  # $contrasts
  # contrast                       estimate       SE    df t.ratio p.value
  # Post,AMRO - Pre,AMRO         -2.1666667 5.743863 58.54  -0.377  1.0000
  # Post,AMRO - Present,AMRO     -1.0000000 5.743863 58.54  -0.174  1.0000
  # Post,AMRO - Post,EASO         5.7115950 5.785623 58.65   0.987  0.9857
  # Post,AMRO - Pre,EASO          4.0449283 5.785623 58.65   0.699  0.9986
  # Post,AMRO - Present,EASO     12.5449283 5.785623 58.65   2.168  0.4392
  # Post,AMRO - Post,GHOW         3.3557975 5.785623 58.65   0.580  0.9996
  # Post,AMRO - Pre,GHOW         -2.8108692 5.785623 58.65  -0.486  0.9999
  # Post,AMRO - Present,GHOW     12.0224641 5.785623 58.65   2.078  0.4977
  # Pre,AMRO - Present,AMRO       1.1666667 5.743863 58.54   0.203  1.0000
  # Pre,AMRO - Post,EASO          7.8782616 5.785623 58.65   1.362  0.9074
  # Pre,AMRO - Pre,EASO           6.2115950 5.785623 58.65   1.074  0.9758
  # Pre,AMRO - Present,EASO      14.7115950 5.785623 58.65   2.543  0.2336
  # Pre,AMRO - Post,GHOW          5.5224641 5.785623 58.65   0.955  0.9885
  # Pre,AMRO - Pre,GHOW          -0.6442025 5.785623 58.65  -0.111  1.0000
  # Pre,AMRO - Present,GHOW      14.1891308 5.785623 58.65   2.452  0.2764
  # Present,AMRO - Post,EASO      6.7115950 5.785623 58.65   1.160  0.9616
  # Present,AMRO - Pre,EASO       5.0449283 5.785623 58.65   0.872  0.9936
  # Present,AMRO - Present,EASO  13.5449283 5.785623 58.65   2.341  0.3355
  # Present,AMRO - Post,GHOW      4.3557975 5.785623 58.65   0.753  0.9977
  # Present,AMRO - Pre,GHOW      -1.8108692 5.785623 58.65  -0.313  1.0000
  # Present,AMRO - Present,GHOW  13.0224641 5.785623 58.65   2.251  0.3881
  # Post,EASO - Pre,EASO         -1.6666667 5.743863 58.54  -0.290  1.0000
  # Post,EASO - Present,EASO      6.8333333 5.743863 58.54   1.190  0.9556
  # Post,EASO - Post,GHOW        -2.3557975 5.785623 58.65  -0.407  1.0000
  # Post,EASO - Pre,GHOW         -8.5224641 5.785623 58.65  -1.473  0.8632
  # Post,EASO - Present,GHOW      6.3108692 5.785623 58.65   1.091  0.9734
  # Pre,EASO - Present,EASO       8.5000000 5.743863 58.54   1.480  0.8602
  # Pre,EASO - Post,GHOW         -0.6891308 5.785623 58.65  -0.119  1.0000
  # Pre,EASO - Pre,GHOW          -6.8557975 5.785623 58.65  -1.185  0.9566
  # Pre,EASO - Present,GHOW       7.9775359 5.785623 58.65   1.379  0.9012
  # Present,EASO - Post,GHOW     -9.1891308 5.785623 58.65  -1.588  0.8072
  # Present,EASO - Pre,GHOW     -15.3557975 5.785623 58.65  -2.654  0.1874
  # Present,EASO - Present,GHOW  -0.5224641 5.785623 58.65  -0.090  1.0000
  # Post,GHOW - Pre,GHOW         -6.1666667 5.743863 58.54  -1.074  0.9758
  # Post,GHOW - Present,GHOW      8.6666667 5.743863 58.54   1.509  0.8469
  # Pre,GHOW - Present,GHOW      14.8333333 5.743863 58.54   2.582  0.2163
  # 
  # P value adjustment: tukey method for comparing a family of 9 estimates

  ####CONCLUSION:

############################################################################
# 6.  Approach     # Ave closest distance of approach
############################################################################
# Get rid of NAs
encd_blja = encd_blja[!is.na(encd_blja$Approach),]

# Chi-Squared Test

# If there isn't much variablity in site (review above results; doesn't seem to be),
# then just do chi-squared
# The sample size only causes an issue for a test of association when you are using the chi-square approximation.
# You can use the exact version of the test of association, called Fisher’s exact test, in this case. R should give 
# you a warning message when you use the chisq.test function and your counts are too low, in which case you can just 
# use the fisher.test function.  

# Yes, error messages
# chisq.test(table(encd_blja$Approach,encd_blja$Stimuli))
# Pearson's Chi-squared test
# data:  table(encd_blja$Approach, encd_blja$Stimuli)
# X-squared = 12.169, df = 8, p-value = 0.1438
# 
# Warning message:
#   In chisq.test(table(encd_blja$Approach, encd_blja$Stimuli)) :
#   Chi-squared approximation may be incorrect
chisq.test(table(encd_blja$Approach,encd_blja$Stimuli), simulate.p.value = TRUE)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# 
# data:  table(encd_blja$Approach, encd_blja$Stimuli)
# X-squared = 12.169, df = NA, p-value = 0.1289
# Conclusion: No significant differences

fisher.test(table(encd_blja$Approach,encd_blja$Stimuli))
# Fisher's Exact Test for Count Data
# 
# data:  table(encd_blja$Approach, encd_blja$Stimuli)
# p-value = 0.2424
# alternative hypothesis: two.sided
# Conclusion: No significant differences


####################
#CORRELATION MATRIX#
####################
str(encd_blja)

PCA_dataframe <- data.frame(encd_blja$Corr.Avenum_elements, encd_blja$Corr.Avedur_ditonal, encd_blja$Corr.Feed_rate, encd_blja$Corr.Avgnum_calls, encd_blja$Approach, encd_blja$DutyCycle)

PCA_correlations<- cor(PCA_dataframe, use = "pairwise.complete.obs")
PCA_correlations

#Alexis RE-RAN 3/19/18
                                # encd_blja.Corr.Avenum_elements encd_blja.Corr.Avedur_ditonal encd_blja.Corr.Feed_rate
# encd_blja.Corr.Avenum_elements                    1.000000000                     0.2272727              -0.09090909
# encd_blja.Corr.Avedur_ditonal                     0.227272727                     1.0000000               0.17272727
# encd_blja.Corr.Feed_rate                         -0.090909091                     0.1727273               1.00000000
# encd_blja.Corr.Avgnum_calls                       0.018181818                     0.4636364               0.15454545
# encd_blja.Approach                               -0.009656091                     0.1882938               0.84490796
# encd_blja.DutyCycle                               0.636363636                     0.1000000               0.16363636
                               # encd_blja.Corr.Avgnum_calls encd_blja.Approach encd_blja.DutyCycle
# encd_blja.Corr.Avenum_elements                  0.01818182       -0.009656091           0.6363636
# encd_blja.Corr.Avedur_ditonal                   0.46363636        0.188293774           0.1000000
# encd_blja.Corr.Feed_rate                        0.15454545        0.844907962           0.1636364
# encd_blja.Corr.Avgnum_calls                     1.00000000        0.188293774           0.5727273
# encd_blja.Approach                              0.18829377        1.000000000           0.2510584
# encd_blja.DutyCycle                             0.57272727        0.251058366           1.0000000


# Janelle ran 3/20/18 -- I think the difference is when NAs are removed (see str below)
                                # encd_blja.Corr.Avenum_elements encd_blja.Corr.Avedur_ditonal encd_blja.Corr.Feed_rate
# encd_blja.Corr.Avenum_elements                     1.00000000                     0.5367053             -0.208499199
# encd_blja.Corr.Avedur_ditonal                      0.53670525                     1.0000000              0.123957682
# encd_blja.Corr.Feed_rate                          -0.20849920                     0.1239577              1.000000000
# encd_blja.Corr.Avgnum_calls                       -0.07152543                     0.3430239              0.003860715
# encd_blja.Approach                                -0.08651946                     0.3204781             -0.253231436
# encd_blja.DutyCycle                                0.63885561                     0.3194309              0.045690901
                                # encd_blja.Corr.Avgnum_calls encd_blja.Approach encd_blja.DutyCycle
# encd_blja.Corr.Avenum_elements                -0.071525428        -0.08651946          0.63885561
# encd_blja.Corr.Avedur_ditonal                  0.343023895         0.32047810          0.31943092
# encd_blja.Corr.Feed_rate                       0.003860715        -0.25323144          0.04569090
# encd_blja.Corr.Avgnum_calls                    1.000000000         0.30251767          0.38028516
# encd_blja.Approach                             0.302517669         1.00000000         -0.08454219
# encd_blja.DutyCycle                            0.380285159        -0.08454219          1.00000000


# 'data.frame':	16 obs. of  46 variables:
#   $ X                   : int  19 20 21 22 23 24 25 26 27 28 ...
# $ Date                : int  20141119 20141119 20141203 20141203 20141205 20141222 20150122 20150122 20150123 20150204 ...
# $ Site                : int  3 3 3 3 1 2 6 6 6 6 ...
# $ Stimuli             : Factor w/ 3 levels "AMRO","EASO",..: 2 3 1 3 1 3 1 2 3 2 ...
# $ Phase               : Factor w/ 2 levels "Post","Present": 2 2 2 2 2 2 2 2 2 2 ...
# $ Phase_dur           : num  279 274 289 281 284 ...
# $ Total_calls         : int  22 25 10 46 46 0 7 1 27 2 ...
# $ Total_ditonal       : int  25 35 11 59 12 0 11 2 35 3 ...
# $ Avedur_ditonal      : num  0.37 0.352 0.355 0.359 0.431 ...
# $ Total_monotonal     : int  NA NA NA NA NA NA NA NA NA NA ...
# $ Total_bell          : int  NA NA NA NA NA NA NA NA NA NA ...
# $ Total_alert         : int  NA NA NA NA NA NA NA NA NA NA ...
# $ Total_redtail       : int  NA NA NA NA NA NA NA NA NA NA ...
# $ Total_trains        : int  3 5 1 8 8 NA 4 1 7 1 ...
# $ Avenum_elements     : num  1.136 1.4 1.1 1.283 0.261 ...
# $ Avedur_trains       : num  1.075 1.448 0.831 1.397 NA ...
# $ Avedur_silence      : num  0.298 0.196 0.096 0.203 NA ...
# $ Avedur_delta        : num  0.697 0.537 0.473 0.598 NA ...
# $ Pre.Phase_dur       : num  310 331 326 326 315 ...
# $ Pre.Total_calls     : int  0 0 91 2 2 0 19 7 3 0 ...
# $ Pre.Total_ditonal   : int  0 0 99 1 6 0 23 15 4 0 ...
# $ Pre.Avedur_ditonal  : num  0 0 0.364 0.351 0.353 ...
# $ Pre.Total_monotonal : int  NA NA NA 8 NA NA NA NA NA NA ...
# $ Pre.Total_bell      : int  NA NA NA NA NA NA NA NA NA NA ...
# $ Pre.Total_alert     : int  NA NA NA NA NA NA NA NA NA NA ...
# $ Pre.Total_redtail   : int  NA NA NA 2 NA NA NA NA NA NA ...
# $ Pre.Total_trains    : int  NA NA 7 NA NA NA 3 5 1 NA ...
# $ Pre.Avenum_elements : num  NA NA 1.09 0.5 3 ...
# $ Pre.Avedur_trains   : num  NA NA 1.02 NA NA ...
# $ Pre.Avedur_silence  : num  NA NA 0.263 NA NA ...
# $ Pre.Avedur_delta    : num  NA NA 0.654 NA NA ...
# $ Approach            : int  3 5 0 2 2 5 3 3 2 5 ...
# $ Video_dur           : int  280 273 288 281 284 154 225 238 252 232 ...
# $ Feed_dur            : int  0 0 49 5 75 0 22 0 0 0 ...
# $ Pre.Approach        : logi  NA NA NA NA NA NA ...
# $ Pre.Video_dur       : int  293 313 316 306 300 309 302 317 309 300 ...
# $ Pre.Feed_dur        : int  0 128 90 26 90 0 0 10 21 10 ...
# $ Corr.Total_calls    : int  22 25 -81 44 44 0 -12 -6 24 2 ...
# $ Corr.Avenum_elements: num  NA NA 0.0121 0.7826 -2.7391 ...
# $ Corr.Avedur_ditonal : num  0.37016 0.35189 -0.00933 0.00783 0.07731 ...
# $ Feed_rate           : num  0 0 10.21 1.07 15.85 ...
# $ Pre.Feed_rate       : num  0 24.5 17.1 5.1 18 ...
# $ Corr.Feed_rate      : num  0 -24.54 -6.88 -4.03 -2.15 ...
# $ Corr.Avgnum_calls   : num  4.73 5.48 -16.83 9.4 9.3 ...
# $ DutyCycle           : num  0.0331 0.045 0.0135 0.0754 0.0182 ...
# $ Trans.DutyCycle     : num  0.182 0.212 0.116 0.275 0.135 ...
