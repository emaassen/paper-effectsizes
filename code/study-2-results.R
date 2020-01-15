# rm(list = ls()) # clear workspace
options(scipen=999)

datm <- read.table("../codebooks/codebook-meta-analyses-complete.csv", header=T, sep = '')

# How many meta-analyses cannot be reproduced based on the reported results only?
sum(is.na(datm$model)) # 7 meta-analyses do not explicitly report whether they use RE or FE models 
sum(is.na(datm$software)) # 19 meta-analyses do not explicitly report on software
sum(is.na(datm$estimator)) # 13 meta-analyses do not explicitly report on estimator

# sum of all:
sum(!is.na(datm$model) & !is.na(datm$software) & !is.na(datm$estimator))
# this number includes MA26 and MA27, which refer to Lipsey and Wilson SPSS Macros.
# From that information we inferred which estimator the authors used, but they do not
# explicitly report it. As such, we don't count these two when we talk about
# clear reporting of information on estimation.

sum(datm$model == "FE", na.rm=T) # 1
sum(datm$model == "FE+RE", na.rm=T) # 3
sum(datm$model == "RE", na.rm=T) # 22
 
sum(datm$software == "CMA", na.rm=T)

# For how many meta-analyses couldn't we reproduce the original results as reported?
# Calculated this by hand because not every meta-analysis reports their method (FE/RE/Mixed)
# We were conservative and always chose the smallest estimate

# Meta-analyses with correlations
# No 6 small discrepancy in both FE and RE
# No 12 small discrepancy in RE estimate
# No 15 small discrepancy in RE estimate
# No 27 large discrepancy in RE estimate

# Meta-analyses with hedges' g
# No 5 small discrepancy in RE estimate
# No 17 small discrepancy in RE estimate
# No 19 small discrepancy in RE estimate

# Meta-analyses with cohen's d
# No 16 moderate discrepancy in RE estimate
# No 23 small discrepancy in RE estimate

# In total 9 meta-analyses for which we couldn't reproduce the reported result based on 
# their reported individual effect sizes

# How many meta-analyses had a small / med / large  change in effect size estimate?
sum(datm$disccat.s == 0) # 24 no discrepancy
sum(datm$disccat.s == 1) # 9 a small discrepancy
sum(datm$disccat.s == 2) # 0
sum(datm$disccat.s == 3) # 0

# How many meta-analyses larger or smaller mean estimate?
# Based on the random-effects meta-analytic estimates reported
sum(datm$recalc.re > datm$effest.re, na.rm=T) # 19 larger
sum(datm$recalc.re < datm$effest.re, na.rm=T) # 13 smaller
sum(datm$recalc.re == datm$effest.re, na.rm=T)

datm$recalc.fe[31] > datm$effest.fe[31] # MA 31 is the only one where only FE was reported
datm$recalc.fe[31] < datm$effest.fe[31] # total 19 larger 14 smaller (13 RE, 1 FE)

sum(datm$recalc.fe > datm$effest.fe, na.rm=T) # 5 larger
sum(datm$recalc.fe < datm$effest.fe, na.rm=T) # 4 smaller
sum(datm$recalc.fe == datm$effest.fe, na.rm=T) # 0

# How many meta-analyses had a small / med / large  change in CI estimate?
sum(datm$disccat.ci.s == 0) # 21 no discrepancy
sum(datm$disccat.ci.s == 1) # 9 a small discrepancy
sum(datm$disccat.ci.s == 2) # 3 a moderate discrepancy
sum(datm$disccat.ci.s == 3) # 0

# How many confidence intervals were wider/smaller?
ci.diff.so <- datm$ciub.so - datm$cilb.so
ci.diff.sc <- datm$ciub.sc - datm$cilb.sc
sum(ci.diff.sc > ci.diff.so) # 18
sum(ci.diff.sc < ci.diff.so) # 15
sum(ci.diff.sc == ci.diff.so) # 0

# How many changes in statistical significance? - none
round(datm$pval.sc - datm$pval.so,2)
datm$pval.sc[1];datm$pval.so[1]
datm$pval.sc[2];datm$pval.so[2]
datm$pval.sc[4];datm$pval.so[4]
datm$pval.sc[8];datm$pval.so[8]
datm$pval.sc[11];datm$pval.so[11]
datm$pval.sc[14];datm$pval.so[14]
datm$pval.sc[18];datm$pval.so[18]
datm$pval.sc[24];datm$pval.so[24]
datm$pval.sc[28];datm$pval.so[28]
datm$pval.sc[32];datm$pval.so[32]

# Tau 2 discrepancies
round(datm$pval.het.so, 3)
round(datm$pval.het.sc, 3)
round(datm$pval.het.so - datm$pval.het.sc, 3)

round(datm$pval.het.so[5], 3);round(datm$pval.het.sc[5], 3)
round(datm$pval.het.so[29], 3);round(datm$pval.het.sc[29], 3)

# How many tau2 estimates were larger/smaller?
sum(datm$tau2.sc > datm$tau2.so) # 17
sum(datm$tau2.sc < datm$tau2.so) # 12
sum(datm$tau2.sc == datm$tau2.so) # 4

# how many heterogeneity estimates increased?
sum(datm$tau2.sc < datm$tau2.so) # 12

# how many heterogeneity estimates decreased?
sum(datm$tau2.so < datm$tau2.sc) # 16

# how many heterogeneity estimates stayed the same?
sum(datm$tau2.sc == datm$tau2.so) # 5


# How often articles cited in 2018?
# Checked in Web of Science on Oct. 16 2019
# Adesope 6
# Alfieri 44
# Babbage 15
# Balliet 16
# Benish 29
# Berry1 37
# Berry2 6
# Card 15
# Crook 50
# De Wit 66
# Else-Quest 18
# Farber 6
# Fischer 40
# Fox 21
# Freund 15
# Green 2
# Hallion 51
# Ihle 7
# Koenig 86
# Kolden 8
# Lucassen 13 
# Mol 74
# Morgan 12
# Munder 5
# Piet 25
# Smith 47
# Tillman 1
# Toosi 8
# Van Iddekinge 8
# Webb 100
# Woodin 4
# Woodley 3
# Yoon 8
sum(c(6,44,15,16,29,37,6,15,50,66,18,6,40,21,15,2,51,7,86,8,13,74,12,5,25,47,1,8,8,100,4,3,8))
length(c(6,44,15,16,29,37,6,15,50,66,18,6,40,21,15,2,51,7,86,8,13,74,12,5,25,47,1,8,8,100,4,3,8))
