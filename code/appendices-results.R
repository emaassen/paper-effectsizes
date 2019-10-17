# rm(list = ls()) clear workspace
options(scipen=999)

datm <- read.table("../codebooks/codebook-meta-analyses-complete.csv", header=T, sep = '')
ret <- read.table("../codebooks/nonretrieved-primary-studies.csv", header=T, sep = ';')

# Supplement A: Additional Results Part 2 ---------------------------------

# How many meta-analyses had a small / med / large  change in effect size estimate?
sum(datm$disccat.c == 0) # 30 no discrepancy
sum(datm$disccat.c == 1) # 3 a small discrepancy
sum(datm$disccat.c == 2) # 0
sum(datm$disccat.c == 3) # 0


# How many meta-analyses had a small / med / large  change in CI estimate?
sum(datm$disccat.ci.c == 0) # 26 no discrepancy
sum(datm$disccat.ci.c == 1) # 6 a small discrepancy
sum(datm$disccat.ci.c == 2) # 1 a moderate discrepancy
sum(datm$disccat.ci.c == 3) # 0

# Tau2 discrepancy
(abs(datm$tau2.co) - abs(datm$tau2.cc)) != 0
abs(datm$disc.tau2.c) > 0.1

# Statistical significance changed? 
round(datm$pval.het.co, 3)
round(datm$pval.het.cc, 3)
round(datm$pval.het.co - datm$pval.het.cc, 3)

round(datm$pval.het.co[5], 3);round(datm$pval.het.cc[5], 3)
round(datm$pval.het.co[17], 3);round(datm$pval.het.cc[17], 3)
round(datm$pval.het.co[31], 3);round(datm$pval.het.cc[31], 3)

round(datm$pval.co, 3)
round(datm$pval.cc, 3)
round(datm$pval.co - datm$pval.cc, 3)

# Reporting standards: heterogeneity --------------------------------------
# Which models were used?
datm$model
sum(datm$model == "FE+RE", na.rm=T) 
sum(datm$model == "RE", na.rm=T)
sum(datm$model == "FE", na.rm=T) 
sum(is.na(datm$model))

# Which software was used?
datm$software
sum(datm$software == "CMA", na.rm=T)
sum(datm$software == "SL", na.rm=T)
sum(datm$software == "SPSS", na.rm=T)
sum(datm$software == "Stata", na.rm=T)
sum(is.na(datm$software))

# Which estimators were used?
datm$estimator
sum(datm$estimator == "DL", na.rm=T)
sum(datm$estimator == "HO", na.rm=T)
sum(datm$estimator == "HS", na.rm=T)
sum(datm$estimator == "LW", na.rm=T)
sum(datm$estimator == "ML", na.rm=T)
sum(datm$estimator == "MoM", na.rm=T)
sum(datm$estimator == "REML", na.rm=T)
sum(is.na(datm$estimator))

# Was subgroup analysis performed?
datm$subgroup
sum(datm$subgroup)


# Reporting standards: outliers -------------------------------------------
# How many report on outliers and how many do separate analyses
datm$out
datm$outsep
sum(datm$out)
sum(datm$outsep)
sum(datm$out & datm$outsep)

# Reporting standards: publication bias -----------------------------------
# How many mention publication bias and how many find publication bias?
datm$pub
datm$pubfound

sum(datm$pub)
sum(datm$pubfound == 1, na.rm=T)
sum(datm$pubfound == 0, na.rm=T)
sum(is.na(datm$pubfound))

# Which publication bias methods were used?
levels(datm$pubmet1)
levels(datm$pubmet2)
levels(datm$pubmet3)
levels(datm$pubmet4)

sum(sum(datm$pubmet1 == "Egger's test", na.rm=T),sum(datm$pubmet2 == "Egger's test", na.rm=T),sum(datm$pubmet3 == "Egger's test", na.rm=T),sum(datm$pubmet4 == "Egger's test", na.rm=T))
sum(sum(datm$pubmet1 == "Fail-safe N", na.rm=T),sum(datm$pubmet2 == "Fail-safe N", na.rm=T),sum(datm$pubmet3 == "Fail-safe N", na.rm=T),sum(datm$pubmet4 == "Fail-safe N", na.rm=T))
sum(sum(datm$pubmet1 == "Adjusted Fail-safe N", na.rm=T),sum(datm$pubmet2 == "Adjusted Fail-safe N", na.rm=T),sum(datm$pubmet3 == "Adjusted Fail-safe N", na.rm=T),sum(datm$pubmet4 == "Adjusted Fail-safe N", na.rm=T))
sum(sum(datm$pubmet1 == "Funnel plot", na.rm=T),sum(datm$pubmet2 == "Funnel plot", na.rm=T),sum(datm$pubmet3 == "Funnel plot", na.rm=T),sum(datm$pubmet4 == "Funnel plot", na.rm=T))
sum(sum(datm$pubmet1 == "Publication status used in meta-regression", na.rm=T),sum(datm$pubmet2 == "Publication status used in meta-regression", na.rm=T),sum(datm$pubmet3 == "Publication status used in meta-regression", na.rm=T),sum(datm$pubmet4 == "Publication status used in meta-regression", na.rm=T))
sum(sum(datm$pubmet1 == "Trim and fill", na.rm=T),sum(datm$pubmet2 == "Trim and fill", na.rm=T),sum(datm$pubmet3 == "Trim and fill", na.rm=T),sum(datm$pubmet4 == "Trim and fill", na.rm=T))

sum(datm$contact)

# Dependency of primary study effect sizes --------------------------------
datm$depreport
datm$depmeasure

sum(datm$depreport)
sum(datm$depmeasure == 0)

# I did not count #22 as one of the measures being used since they are talking about independent samples
datm$depmeasure[22]

# Effect of outliers on pooled effect size estimates ----------------------
# Number of primary studies that had a small (1), moderate (2) or large (3) effect on pooled effect size estimates
datm$outeff1.s;datm$outeff2.s;datm$outeff3.s;
sum(datm$outeff1.s, na.rm=T)
sum(datm$outeff2.s, na.rm=T)
sum(datm$outeff3.s, na.rm=T)


# Nonretrieval ------------------------------------------------------------
nrow(ret)
levels(ret$type)
# In total 154 primary papers not found in 26 out of 33 meta-analyses = 79%.

sum(ret$type == "article");sum(ret$type == "article")/nrow(ret)*100                                                  # of which 70 (45%) are articles
sum((ret$type == "dissertation")+(ret$type == "unpublished doctoral dissertation"))
sum((ret$type == "dissertation")+(ret$type == "unpublished doctoral dissertation"))/nrow(ret)*100                    # of which 48 (31%) are dissertations
sum(ret$type == "paper presented at conference");sum(ret$type == "paper presented at conference")/nrow(ret)*100      # of which 8 (5%) are papers presented at conferences
sum((ret$type == "book")+(ret$type == "book chapter"));sum((ret$type == "book")+(ret$type == "book chapter"))/nrow(ret)*100  # of which 7 (5%) are books or book chapters
sum((ret$type == "manuscript submitted for publication")+(ret$type == "unpublished manuscript"))
sum((ret$type == "manuscript submitted for publication")+(ret$type == "unpublished manuscript"))/nrow(ret)*100       # of which 7 (5%) are unpublished manuscripts
sum(ret$type == "poster");sum(ret$type == "poster")/nrow(ret)*100                                                    # of which 7 (5%) are posters
sum(ret$type == "unpublished master thesis");sum(ret$type == "unpublished master thesis")/nrow(ret)*100              # of which 3 (2%) are unpublished master theses
sum(ret$type == "unpublished raw data");sum(ret$type == "unpublished raw data")/nrow(ret)*100                        # of which 3 (2%) are unpublished raw data
sum(ret$type == "report");sum(ret$type == "report")/nrow(ret)*100                                                    # of which 1 (1%) is a report


# Replicating meta-analytic effect size estimate --------------------------
datm$recalc.cat
33 - sum(datm$recalc.cat == 0)
datm$effest.fe[27]
datm$effest.re[27]
datm$recalc.re[27]
datm$diff.re[27]
datm$efftype[27]

datm$effest.fe[16]
datm$effest.re[16]
datm$recalc.re[16]
datm$diff.re[16]
datm$efftype[16]


