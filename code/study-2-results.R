rm(list = ls())
options(scipen=999)

datm <- read.table("../codebooks/codebook-meta-analyses-final-complete.csv", header=T, sep = '')

# How many meta-analyses cannot be reproduced based on the reported results only?

sum(is.na(datm$model)) # 7 meta-analyses do not explicitly report whether they use RE or FE models 
sum(is.na(datm$software)) # 19 meta-analyses do not explicitly report on software
sum(is.na(datm$estimator)) # 13 meta-analyses do not explicitly report on estimator

sum(datm$model == "FE", na.rm=T)
sum(datm$model == "FE+RE", na.rm=T)
sum(datm$model == "RE", na.rm=T)

sum(datm$software == "CMA", na.rm=T)

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


write.table(datm, "../codebooks/codebook-meta-analyses-final-complete.csv", row.names=F, col.names=T, sep = "")
