rm(list = ls())
setwd("C:/Users/s421506/tiu/research/effectsizes/codebooks/")
packages <- c("readxl")
options(scipen=999)
#lapply(packages,install.packages(packages),character.only=T)     # if packages are not yet installed
lapply(packages,library,character.only=T)

datm <- read.table("codebook-meta-analyses-final-complete.csv", header=T, sep = '')

# How many meta-analyses cannot be reproduced based on the reported results only?

sum(is.na(datm$model)) # 7 meta-analyses do not explicitly report whether they use RE or FE models 
sum(is.na(datm$software)) # 19 meta-analyses do not explicitly report on software
sum(is.na(datm$estimator)) # 13 meta-analyses do not explicitly report on estimator

for (i in 1:nrow(datm)) {
  
  datm$diff.fe[i] <- abs(datm$effest.fe[i]) - abs(datm$recalc.fe[i])
  datm$diff.re[i] <- abs(datm$effest.re[i]) - abs(datm$recalc.re[i])
  
  if (is.na(datm$diff.fe[i])) {
    
    datm$diff.fe[i] <- 0
    
  }
  
  if (is.na(datm$diff.re[i])) {
    
    datm$diff.re[i] <- 0
    
  }
  
}

for (i in 1:nrow(datm)) {
  
  if (datm$efftype[i] == "g") {
    
    # Fill in discrepancy category for CI of effect sizes on subset checked MAs
    if (abs(datm$diff.fe[i]) < 0.050 & abs(datm$diff.re[i]) < 0.050) {
      datm$recalc.cat[i] = 0
    } else if (abs(datm$diff.fe[i]) > 0.249 | abs(datm$diff.re[i]) > 0.249) {
      datm$recalc.cat[i] = 3
    } else if (abs(datm$diff.fe[i]) > 0.149 & abs(datm$diff.fe[i]) <= 0.249 | abs(datm$diff.re[i]) > 0.149 & abs(datm$diff.re[i]) <= 0.249) {
      datm$recalc.cat[i] = 2
    } else if (abs(datm$diff.fe[i]) >= 0.050 & abs(datm$diff.fe[i]) <= 0.149 | abs(datm$diff.re[i]) >= 0.050 & abs(datm$diff.re[i]) <= 0.149) {
      datm$recalc.cat[i] = 1
    } else {
      datm$recalc.cat[i] == "check"
    }
  }
  
  if (datm$efftype[i] == "d") {
    
    # Fill in discrepancy category for CI of effect sizes on subset checked MAs
    if (abs(datm$diff.fe[i]) < 0.051 & abs(datm$diff.re[i]) < 0.051) {
      datm$recalc.cat[i] = 0
    } else if (abs(datm$diff.fe[i]) > 0.252 | abs(datm$diff.re[i]) > 0.252) {
      datm$recalc.cat[i] = 3
    } else if (abs(datm$diff.fe[i]) > 0.151 & abs(datm$diff.fe[i]) <= 0.252 | abs(datm$diff.re[i]) > 0.151 & abs(datm$diff.re[i]) <= 0.252) {
      datm$recalc.cat[i] = 2
    } else if (abs(datm$diff.fe[i]) >= 0.051 & abs(datm$diff.fe[i]) <= 0.151 | abs(datm$diff.re[i]) >= 0.051 & abs(datm$diff.re[i]) <= 0.151) {
      datm$recalc.cat[i] = 1
    } else {
      datm$recalc.cat[i] == "check"
    }
  }
  
  if (datm$efftype[i] == "r" | datm$efftype[i] == "z") {
    
    # Fill in discrepancy category for CI of effect sizes on subset checked MAs
    if (abs(datm$diff.fe[i]) < 0.025 & abs(datm$diff.re[i]) < 0.025) {
      datm$recalc.cat[i] = 0
    } else if (abs(datm$diff.fe[i]) > 0.125 | abs(datm$diff.re[i]) > 0.125) {
      datm$recalc.cat[i] = 3
    } else if (abs(datm$diff.fe[i]) > 0.075 & abs(datm$diff.fe[i]) <= 0.125 | abs(datm$diff.re[i]) > 0.075 & abs(datm$diff.re[i]) <= 0.125) {
      datm$recalc.cat[i] = 2
    } else if (abs(datm$diff.fe[i]) >= 0.025 & abs(datm$diff.fe[i]) <= 0.075 | abs(datm$diff.re[i]) >= 0.025 & abs(datm$diff.re[i]) <= 0.075) {
      datm$recalc.cat[i] = 1
    } else {
      datm$recalc.cat[i] == "check"
    }
  }
}

sum(datm$recalc.cat != 0) # for 11 meta-analyes we cannot reproduce the reported effect size based on the reported primary study effect sizes

# How many meta-analyses had a small / med / large  change in effect size estimate?
sum(datm$disccat.s == 0) # 24 no discrepancy
sum(datm$disccat.s == 1) # 9 a small discrepancy
sum(datm$disccat.s == 2) # 0
sum(datm$disccat.s == 3) # 0

# How many meta-analyses had a small / med / large  change in CI estimate?
sum(datm$disccat.ci.s == 0) # 21 no discrepancy
sum(datm$disccat.ci.s == 1) # 9 a small discrepancy
sum(datm$disccat.ci.s == 2) # 3 a moderate discrepancy
sum(datm$disccat.ci.s == 3) # 0

# How many meta-analyses had a small / med / large  change in tau2 estimate?
sum(datm$disccat.tau2.s == 0) # 25 no discrepancy
sum(datm$disccat.tau2.s == 1) # 6 a small discrepancy
sum(datm$disccat.tau2.s == 2) # 2
sum(datm$disccat.tau2.s == 3) # 0

sum(datm$disccat.tau2.s == 0 & datm$disccat.ci.s == 0 & datm$disccat.s == 0) # 19 MA no discrepancies at all 

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

