rm(list = ls())
x <- c("metafor","readxl","dplyr")
#lapply(x,install.packages(x),character.only=T)
lapply(x,library,character.only=T)

est.fe <- est.re <- k.tot <- outlier.total <- regular.total <- c() # empty vector to store results 
eff.so <- cilb.so <- ciub.so <- tau2.so <- c() #empty vectors for subset original meta-analyses
eff.sc <- cilb.sc <- ciub.sc <- tau2.sc <- c() #empty vectors for subset checked meta-analyses
krep1 <- percent <- krep2 <- percent2 <- pval.so <- pval.sc <- c() #empty vectors for misc results

setwd("C:/Users/s421506/tiu/research/effectsizes/codebooks/")
dat <- read_excel("codebook-meta-analyses-final-empty.xlsx", 2)

# Functions ---------------------------------------------------------------
dependency <- function(y,z) {
  
  dep <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 
  colnames(dep) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") 
  dependent <- c()
  df$eff <- y
  df$effvar <- z
  df <- df[order(df$eff),]                                # order df by effect size
  
  for (i in 1:nrow(df)) {  
    
    subset <- df[i:(i+1),]                              # make subset of two adjacent studies (e.g., 1 with 2, 2 with 3 etc.)
    res <- rma(eff, effvar, data = subset, method="FE") # FE MA
    
    dep[i,] <- c(i, (i+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, etc.) 
                 round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
                 df[i,]$n,df[i+1,]$n,                   # sample size of two studies
                 df[i,]$eff,df[i+1,]$eff,               # effect size of two studies  
                 as.character(df[i,]$a1),               # authors study 1
                 as.character(df[i,]$a2),
                 as.character(df[i,]$a3),
                 as.character(df[i,]$a4),
                 as.character(df[i,]$a5),
                 as.character(df[i,]$a6),
                 as.character(df[i,]$a7),
                 as.character(df[i+1,]$a1),             # authors study 2
                 as.character(df[i+1,]$a2),
                 as.character(df[i+1,]$a3),
                 as.character(df[i+1,]$a4),
                 as.character(df[i+1,]$a5),
                 as.character(df[i+1,]$a6),
                 as.character(df[i+1,]$a7))
  } 
  
  dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
  dep[dep == ""] <- NA
  dep <- dep[dep[, "QEp"] >= 0.80,]                   
  
  for (i in 1:nrow(dep)) {                             # if any authors are duplicates in both studies, keep the row with two studies being compared
    if (any(duplicated(dep[i,9:22], incomparables=NA) == TRUE)) {
      dependent <- rbind(dependent,dep[i,])
    } 
  }
  
  return(dependent)
}

setwd("C:/Users/s421506/tiu/research/effectsizes/data-per-ma")

# Adesope -----------------------------------------------------------------
df <- read.table("adesope_complete.csv", header=T, sep=';') # recalculated effect size based on reported effect sizes
J <- 1 - (3 / (4 * df$n - 9))
df$d <- df$g / J                                        # cohen's d
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd                                    # variance hedges' g 

# select outliers from dataset
res <- rma(g, vg, data=df)                            # random effects meta-analysis
l1o <- leave1out(res)                                 # leave1out analysis
q <- res$QE - l1o$Q                                   # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                     # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                      # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg))
write.table(df, file = "adesope_complete.csv", row.names=FALSE, sep=";")

# complete analysis with reported values
res.fe <- rma(g, vg, data=df, method="FE")
res.re <- rma(g, vg, data=df, method="DL") 
est.fe <- c(est.fe,res.fe$b)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))

rm(df)

# subset
df <- read.table("adesope_subset.csv", header=T, sep=';') # reproduced results subset
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$d.o <- df$effest.exp / J.o                              # cohen's d
df$d.c <- df$effestnew.exp / J.c                           # cohen's d
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n))) # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.c^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 
write.table(df, file = "adesope_subset.csv", row.names=FALSE, sep=";")

# subset original - the authors are unclear but do estimate heterogeneity, so we compare with the random effects model
res.so <- rma(effest.exp, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked

krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[1]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Alfieri -----------------------------------------------------------------
df <- read.table("alfieri_complete.csv", header=T, sep=";")

df$n <- df$n1 + df$n2
dfs <- df$n1 + df$n2 - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J                                        # cohen's d
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd                                    # variance hedges' g 

# select outliers from dataset
res <- rma(g, vg, data=df)                            # random effects meta-analysis
l1o <- leave1out(res)                                 # leave1out analysis
q <- res$QE - l1o$Q                                   # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                     # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                      # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg))
write.table(df, file = "alfieri_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.fe <- rma(g, vg, data=df, method="FE")
res.re <- rma(g, vg, data=df, method="DL") 
est.fe <- c(est.fe,res.fe$b)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))

rm(df)

# subset
df <- read.table("alfieri_subset.csv", header=T, sep=";")

df$n <- df$n1 + df$n2
dfs <- df$n1 + df$n2 - 2
dfs.c <- df$ncnew + df$ntnew - 2
J.o <- 1 - (3 / (4 * dfs - 1))
J.c <- 1 - (3 / (4 * dfs.c - 1))
df$d.o <- df$effest.exp / J.o                          # cohen's d
df$d.c <- df$effestnew.exp / J.c 
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n)))           # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.o^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                             # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                             # variance hedges' g 

write.table(df, file = "alfieri_subset.csv", row.names=FALSE, sep=";")

# subset original - the authors are unclear but do estimate heterogeneity, so we compare with the random effects model
res.so <- rma(effest.exp, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked

krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[2]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Babbage -----------------------------------------------------------------
df <- read.table("babbage_complete.csv", header=T, sep=";")

J <- 1 - (3 / (4 * df$n - 9))
df$d <- df$g / J                                        # cohen's d
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd                                    # variance hedges' g 

# select outliers from dataset
res <- rma(g, vg, data=df)                            # random effects meta-analysis
l1o <- leave1out(res)                                 # leave1out analysis
q <- res$QE - l1o$Q                                   # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                     # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                      # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg))        # no possible dependent studies
write.table(df, file = "babbage_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(g, vg, data=df, method="DL") 
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))

rm(df)

df <- read.table("babbage_subset.csv", header=T, sep=";")

J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$d.o <- df$effest.exp / J.o                              # cohen's d
df$d.c <- df$effestnew.exp / J.c                           # cohen's d
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n))) # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.c^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 

write.table(df, file = "babbage_subset.csv", row.names=FALSE, sep=";")

# subset original - the authors are unclear but do estimate heterogeneity, so we compare with the random effects model
res.so <- rma(effest.exp, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked

krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[3]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Balliet -----------------------------------------------------------------
df <- read.table("balliet_complete.csv", header=T, sep=";")

J <- 1 - (3 / (4 * df$n - 9))
df$d <- df$g / J                                        # cohen's d
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd                                    # variance hedges' g 

# select outliers from dataset
res <- rma(g, vg, data=df)                            # random effects meta-analysis
l1o <- leave1out(res)                                 # leave1out analysis
q <- res$QE - l1o$Q                                   # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                     # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                      # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg))
write.table(df, file = "balliet_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.fe <- rma(g, vg, data=df, method="FE")
res.re <- rma(g, vg, data=df, method="DL") 
est.fe <- c(est.fe,res.fe$b)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("balliet_subset.csv", header=T, sep=";")

J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$d.o <- df$effest.exp / J.o                              # cohen's d
df$d.c <- df$effestnew.exp / J.c                           # cohen's d
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n))) # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.c^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 

write.table(df, file = "balliet_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[4]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Benish ------------------------------------------------------------------
df <- read.table("benish_complete.csv", header=T, sep=";")

J <- 1 - (3 / (4 * df$n - 9))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 

res <- rma(g, vg, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg))
write.table(df, file = "benish_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(g, vg, data=df, method="HE") 
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

df <- read.table("benish_subset.csv", header=T, sep=";")

J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$d.o <- df$effest.exp / J.o                              # cohen's d
df$d.c <- df$effestnew.exp / J.c                           # cohen's d
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n))) # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.c^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c  

write.table(df, file = "benish_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vg.o, data=df, method="HE")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="HE")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[5]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Berry1 ------------------------------------------------------------------
df <- read.table("berry1_complete.csv", header=T, sep=";")

df$vr <- ((1 - (df$r^2))^2) / (df$n-1)             # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))               
write.table(df, file = "berry1_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.fe <- rma(r, vr, data=df, method="FE")
res.re <- rma(r, vr, data=df, method="HS") 
est.fe <- c(est.fe,res.fe$b)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

#subset
df <- read.table("berry1_subset.csv", header=T, sep=";")

df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1)             # variance correlation r

write.table(df, file = "berry1_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(effest.exp, vr.o, data=df, method="HS")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="HS")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[6]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Berry2 ------------------------------------------------------------------
df <- read.table("berry2_complete.csv", header=T, sep=";")

df$vr <- ((1 - (df$r^2))^2) / (df$n-1)             # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))   
write.table(df, file = "berry2_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.fe <- rma(r, vr, data=df, method="FE")
res.re <- rma(r, vr, data=df, method="HS") 
est.fe <- c(est.fe,res.fe$b)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

#subset
df <- read.table("berry2_subset.csv", header=T, sep=";")

df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1)             # variance correlation r

write.table(df, file = "berry2_subset.csv", row.names=FALSE, sep=";")

# subset original - the authors are unclear but do check for moderators, so we compare with the random effects model
res.so <- rma(effest.exp, vr.o, data=df, method="HS")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="HS")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[7]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Card --------------------------------------------------------------------
df <- read.table("card_complete.csv", header=T, sep=";")

df$vr <- ((1 - (df$r^2))^2) / (df$n-1)             # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))   
write.table(df, file = "card_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
z <- 0.5 * log((1 + df$r) / (1 - df$r))
vz <- 1 / (df$n - 3)    
res.re <- rma(z, vz, data=df, method="DL") 
# transform estimate back to correlation
res.re <-tanh(res.re$b)
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("card_subset.csv", header=T, sep=";")

df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1)             # variance correlation r

write.table(df, file = "card_subset.csv", row.names=FALSE, sep=";")

# subset original 
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$vz.o <- 1 / (df$n - 3)    
res.so <- rma(z.o, vz.o, data=df, method="DL") 

# subset reproduced
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.c <- 1 / (df$nnew - 3) 
res.sc <- rma(z.c, vz.c, data=df, method="DL")               

# transform estimate back to correlation
eff.so <- c(eff.so,tanh(res.so$b)) # MA subset original effect size estimate
cilb.so <- c(cilb.so,tanh(res.so$ci.lb)) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,tanh(res.so$ci.ub)) # upperbound
tau2.so <- c(tau2.so,tanh(res.so$tau2)) # tau2 estimate

eff.sc <- c(eff.sc,tanh(res.sc$b)) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,tanh(res.sc$ci.lb))
ciub.sc <- c(ciub.sc,tanh(res.sc$ci.ub))
tau2.sc <- c(tau2.sc,tanh(res.sc$tau2))

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[8]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Crook -------------------------------------------------------------------
df <- read.table("crook_complete.csv", header=T, sep=";")

df$vr <- ((1 - (df$r^2))^2) / (df$n-1)             # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                  # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id,  study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))   
write.table(df, file = "crook_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.fe <- res.re <- sum(df$n * df$r)/sum(df$n)

# estimate confidence interval (see formula sheet):
r.bar <- res.re
N <- sum(df$n)
K <- nrow(df)
sd2.res <- var(df$r)
se <- ((1 - r.bar^2)^2 / (N - K)) + (sd2.res / K)^1/2
ci.lb <- r.bar - (1.96 * se)
ci.ub <- r.bar + (1.96 * se)
# mean estimate corrected for attenuation (see formula sheet):
res.fe <- res.re <- (res.fe / 0.82)
est.fe <- c(est.fe,res.fe)
est.re <- c(est.re,res.re)
k.tot <- c(k.tot,nrow(df))
rm(df)
# insert hier tau2 formula from Hunter & Schmidt

# subset
df <- read.table("crook_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1)       # variance correlation r
write.table(df, file = "crook_subset.csv", row.names=FALSE, sep=";")

# subset original
est.o <- sum(df$n * df$effest.exp)/sum(df$n)
# estimate confidence interval (see formula sheet):
r.bar <- est.o
N <- sum(df$n)
K <- nrow(df)
sd2.res <- var(df$effest.exp)
se <- ((1 - r.bar^2)^2 / (N - K)) + (sd2.res / K)^1/2
ci.lb.o <- r.bar - (1.96 * se)
ci.ub.o <- r.bar + (1.96 * se)
# mean estimate corrected for attenuation (see formula sheet):
est.o <- est.o / 0.82
# insert hier tau2 formula from Hunter & Schmidt
# k <- nrow(df)
# Y <- as.matrix(df$r)
# wi <- 1/df$vr.o
# W <- diag(wi, nrow = k, ncol = k)
# stXWX <- .invcalc(X = X, W = W, k = k)
# P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
# RSS <- crossprod(Y, P) %*% Y
# tau2 <- ifelse(tau2.fix, tau2.val, (RSS - k)/sum(wi))


# subset reproduced
est.c <- sum(df$nnew * df$effestnew.exp)/sum(df$nnew)
# estimate confidence interval (see formula sheet):
r.bar <- est.c
N <- sum(df$nnew)
K <- nrow(df)
sd2.res <- var(df$effestnew.exp)
se <- ((1 - r.bar^2)^2 / (N - K)) + (sd2.res / K)^1/2
ci.lb.c <- r.bar - (1.96 * se)
ci.ub.c <- r.bar + (1.96 * se)
# mean estimate corrected for attenuation (see formula sheet):
est.c <- est.c / 0.82
# insert hier tau2 formula from Hunter & Schmidt

eff.so <- c(eff.so,est.o) # MA subset original effect size estimate
cilb.so <- c(cilb.so,ci.lb.o) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,ci.ub.o) # upperbound
tau2.so <- c(tau2.so,0) # tau2 estimate

eff.sc <- c(eff.sc,est.c) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,ci.lb.c)
ciub.sc <- c(ciub.sc,ci.ub.c)
tau2.sc <- c(tau2.sc,0)

pval.so <- c(pval.so,NA) # p value original
pval.sc <- c(pval.sc,NA) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[9]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# DeWit -------------------------------------------------------------------
df <- read.table("dewit_complete.csv", header=T, sep=";")

df$vr <- ((1 - (df$r^2))^2) / (df$n-1)             # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))  
write.table(df, file = "dewit_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.fe <- rma(r, vr, data=df, method="FE")
res.re <- rma(r, vr, data=df, method="HS") 
est.fe <- c(est.fe,res.fe$b)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("dewit_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1)       
write.table(df, file = "dewit_subset.csv", row.names=FALSE, sep=";")

# subset original - the authors are unclear but do check for moderators, so we compare with the random effects model
res.so <- rma(effest.exp, vr.o, data=df, method="HS")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="HS")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[10]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Elsequest ---------------------------------------------------------------
df <- read.table("elsequest_complete.csv", header=T, sep=";")

df$n <- df$n1 + df$n2
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
res <- rma(d, vd, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                  # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, d))
regular.total <- bind_rows(regular.total,select(regular, id, study, d))

# dependency check with d
depdf <- as.data.frame(dependency(df$d,df$vd))
write.table(df, file = "elsequest_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(d, vd, data=df, method="DL") 
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

df <- read.table("elsequest_subset.csv", header=T, sep=";")

df$n <- df$n1 + df$n2
df$vd.o <- df$n / ((df$n / 2)^2 + (df$effest.exp^2 / (2 * df$n)))  # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$effestnew.exp^2 / (2 * df$nnew)))  # variance cohen's d

write.table(df, file = "elsequest_subset.csv", row.names=FALSE, sep=";")

# subset original - the authors are unclear but do check for moderators, so we compare with the random effects model
res.so <- rma(effest.exp, vd.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vd.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[11]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Farber ------------------------------------------------------------------
df <- read.table("farber_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))        # no possible dependent studies
write.table(df, file = "farber_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported estimates
res.re <- rma(r, vr, data=df, method="DL") 
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

#subset
df <- read.table("farber_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1)             # variance correlation r
write.table(df, file = "farber_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vr.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[12]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Fisher -----------------------------------------------------------------
df <- read.table("fischer_complete.csv", header=T, sep=";")
df$n <- df$n1 + df$n2
dfs <- df$n1 + df$n2 - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 

res <- rma(g, vg, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "fischer_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.fe <- rma(g, vg, data=df, method="FE")
res.re <- rma(g, vg, data=df, method="DL")
est.fe <- c(est.fe,res.fe$b)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

df <- read.table("fischer_subset.csv", header=T, sep=";")

df$n <- df$n1 + df$n2
dfs <- df$n1 + df$n2 - 2
dfs.c <- df$ncnew + df$ntnew - 2
J.o <- 1 - (3 / (4 * dfs - 1))
J.c <- 1 - (3 / (4 * dfs.c - 1))
df$d.o <- df$effest.exp / J.o                          # cohen's d
df$d.c <- df$effestnew.exp / J.c 
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n)))           # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.o^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                             # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                             # variance hedges' g 

write.table(df, file = "fischer_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[13]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Fox ---------------------------------------------------------------------
df <- read.table("fox_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "fox_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
z <- 0.5 * log((1 + df$r) / (1 - df$r))
vz <- 1 / (df$n - 3)    
res.re <- rma(z, vz, data=df, method="ML") 
# transform estimate back to correlation
res.re <-tanh(res.re$b)
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("fox_subset.csv", header=T, sep=";")

df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1)             # variance correlation r

write.table(df, file = "fox_subset.csv", row.names=FALSE, sep=";")

# subset original 
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$vz.o <- 1 / (df$n - 3)    
res.so <- rma(z.o, vz.o, data=df, method="ML") 

# subset reproduced
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.c <- 1 / (df$nnew - 3) 
res.sc <- rma(z.c, vz.c, data=df, method="ML")               

# transform estimate back to correlation
eff.so <- c(eff.so,tanh(res.so$b)) # MA subset original effect size estimate
cilb.so <- c(cilb.so,tanh(res.so$ci.lb)) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,tanh(res.so$ci.ub)) # upperbound
tau2.so <- c(tau2.so,tanh(res.so$tau2)) # tau2 estimate

eff.sc <- c(eff.sc,tanh(res.sc$b)) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,tanh(res.sc$ci.lb))
ciub.sc <- c(ciub.sc,tanh(res.sc$ci.ub))
tau2.sc <- c(tau2.sc,tanh(res.sc$tau2))

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[14]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Freund ------------------------------------------------------------------
df <- read.table("freund_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "freund_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma.mv(r, vr, data=df)
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("freund_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1)             # variance correlation r
write.table(df, file = "freund_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma.mv(effest.exp, vr.o, data=df)  

# subset reproduced
res.sc <- rma.mv(effestnew.exp, vr.c, data=df)               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[15]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Green -------------------------------------------------------------------
df <- read.table("green_complete.csv", header=T, sep=";")
df$n <- df$n1 + df$n2
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
res <- rma(d, vd, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, d))
regular.total <- bind_rows(regular.total,select(regular, id, study, d))

# dependency check with d
depdf <- as.data.frame(dependency(df$d,df$vd))
write.table(df, file = "green_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(d, vd, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("green_subset.csv", header=T, sep=";")

df$n <- df$n1 + df$n2
df$vd.o <- df$n / ((df$n / 2)^2 + (df$effest.exp^2 / (2 * df$n)))  # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$effestnew.exp^2 / (2 * df$nnew)))  # variance cohen's d

write.table(df, file = "green_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vd.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vd.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[16]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Hallion -----------------------------------------------------------------
df <- read.table("hallion_complete.csv", header=T, sep=";")
J <- 1 - (3 / (4 * (df$n - 1)))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
res <- rma(g, vg, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "hallion_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(g, vg, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("hallion_subset.csv", header=T, sep=";")
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$d.o <- df$effest.exp / J.o                              # cohen's d
df$d.c <- df$effestnew.exp / J.c                           # cohen's d
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n))) # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.c^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c       
write.table(df, file = "hallion_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[17]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Ihle --------------------------------------------------------------------
df <- read.table("ihle_complete.csv", header=T, sep=";")

df$n <- df$n1 + df$n2
dfs <- df$n1 + df$n2 - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
res <- rma(g, vg, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "ihle_complete.csv", row.names=FALSE, sep=";")

# meta-analysis on reported values
res.re <- rma(g, vg, data=df, method="HE")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("ihle_subset.csv", header=T, sep=";")
df$n <- df$n1 + df$n2
dfs <- df$n1 + df$n2 - 2
dfs.c <- df$ncnew + df$ntnew - 2
J.o <- 1 - (3 / (4 * dfs - 1))
J.c <- 1 - (3 / (4 * dfs.c - 1))
df$d.o <- df$effest.exp / J.o                          # cohen's d
df$d.c <- df$effestnew.exp / J.c 
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n)))           # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.o^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                             # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                             # variance hedges' g 

write.table(df, file = "ihle_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vg.o, data=df, method="HE")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="HE")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[18]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Koenig ------------------------------------------------------------------
df <- read.table("koenig_complete.csv", header=T, sep=";")
dfs <- df$n - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
res <- rma(g, vg, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "koenig_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(g, vg, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("koenig_subset.csv", header=T, sep=";")
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$d.o <- df$effest.exp / J.o                              # cohen's d
df$d.c <- df$effestnew.exp / J.c                           # cohen's d
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n))) # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.c^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c   
write.table(df, file = "koenig_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[19]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Kolden ------------------------------------------------------------------
df <- read.table("kolden_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "kolden_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(r, vr, data=df, method="HE")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("kolden_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1) 
write.table(df, file = "kolden_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vr.o, data=df, method="HE")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="HE")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[20]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Lucassen ----------------------------------------------------------------
df <- read.table("lucassen_complete.csv", header=T, sep = ";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "lucassen_complete.csv", row.names=FALSE, sep=";")

# meta-analysis on reported values
res.re <- rma(r, vr, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("lucassen_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1)
write.table(df, file = "lucassen_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vr.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[21]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Mol ---------------------------------------------------------------------
df <- read.table("mol_complete.csv", header=T, sep=";")
df$vz <- 1 / (df$n - 3)                              # variance fisher's z 
res <- rma(z, vz, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, z))
regular.total <- bind_rows(regular.total,select(regular, id, study, z))

# dependency check with r
depdf <- as.data.frame(dependency(df$z,df$vz))
write.table(df, file = "mol_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(z, vz, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("mol_subset.csv", header=T, sep=";")
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "mol_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vz.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vz.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[22]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Morgan ------------------------------------------------------------------
df <- read.table("morgan_complete.csv", header=T, sep=";")
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
res <- rma(d, vd, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, d))
regular.total <- bind_rows(regular.total,select(regular, id, study, d))

# dependency check with d
depdf <- as.data.frame(dependency(df$d,df$vd))
write.table(df, file = "morgan_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(d, vd, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("morgan_subset.csv", header=T, sep=";")
df$vd.o <- df$n / ((df$n / 2)^2 + (df$effest.exp^2 / (2 * df$n)))  # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$effestnew.exp^2 / (2 * df$nnew)))  # variance cohen's d
write.table(df, file = "morgan_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vd.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vd.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[23]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)



# Munder ------------------------------------------------------------------
df <- read.table("munder_complete.csv", header=T, sep=";")
dfs <- df$n - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
res <- rma(g, vg, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier <- outlier[order(outlier$study),]
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "munder_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(g, vg, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("munder_subset.csv", header=T, sep=";")
dfs <- df$n - 2
dfs.c <- df$nnew - 2
J.o <- 1 - (3 / (4 * dfs - 1))
J.c <- 1 - (3 / (4 * dfs.c - 1))
df$d.o <- df$effest.exp / J.o                          # cohen's d
df$d.c <- df$effestnew.exp / J.c 
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n)))           # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.o^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                             # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                             # variance hedges' g 
write.table(df, file = "munder_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[24]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Piet --------------------------------------------------------------------
df <- read.table("piet_complete.csv", header=T, sep=";")
dfs <- df$n - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
res <- rma(g, vg, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "piet_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(g, vg, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("piet_subset.csv", header=T, sep=";")
dfs <- df$n - 2
dfs.c <- df$nnew - 2
J.o <- 1 - (3 / (4 * dfs - 1))
J.c <- 1 - (3 / (4 * dfs.c - 1))
df$d.o <- df$effest.exp / J.o                          # cohen's d
df$d.c <- df$effestnew.exp / J.c 
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n)))           # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.o^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                             # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c 
write.table(df, file = "piet_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[25]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Smith -------------------------------------------------------------------
df <- read.table("smith_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "smith_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(r, vr, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("smith_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1) 
write.table(df, file = "smith_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vr.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[26]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Tillman -----------------------------------------------------------------
df <- read.table("tillman_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "tillman_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(r, vr, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("tillman_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1) 
write.table(df, file = "tillman_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vr.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[27]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Toosi -------------------------------------------------------------------
df <- read.table("toosi_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "toosi_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.re <- rma(r, vr, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("toosi_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1) 
write.table(df, file = "toosi_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vr.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[28]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# VanIddekinge ------------------------------------------------------------
df <- read.table("vaniddekinge_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "vaniddekinge_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.fe <- rma(r, vr, data=df, method="FE")
res.re <- rma(r, vr, data=df, method="HS")
est.fe <- c(est.fe,res.fe$b)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("vaniddekinge_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1) 
write.table(df, file = "vaniddekinge_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vr.o, data=df, method="HS")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="HS")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[29]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Webb --------------------------------------------------------------------
df <- read.table("webb_complete.csv", header=T, sep=";")
df$n <- df$n1 + df$n2
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
res <- rma(d, vd, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, d))
regular.total <- bind_rows(regular.total,select(regular, id, study, d))

# dependency check with d
depdf <- as.data.frame(dependency(df$d,df$vd))
write.table(df, file = "webb_complete.csv", row.names=FALSE, sep=";")

# meta-analysis on reported values
res.re <- rma(d, vd, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("webb_subset.csv", header=T, sep=";")
df$n <- df$n1 + df$n2
df$vd.o <- df$n / ((df$n / 2)^2 + (df$effest.exp^2 / (2 * df$n)))           # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$effestnew.exp^2 / (2 * df$nnew)))
write.table(df, file = "webb_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vd.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vd.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[30]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Woodin ------------------------------------------------------------------
df <- read.table("woodin_complete.csv", header=T, sep=";")
J <- 1 - (3 / (4 * df$n - 9))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
res <- rma(g, vg, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier <- outlier[order(outlier$study),]
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, g))
regular.total <- bind_rows(regular.total,select(regular, id, study, g))

# dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "woodin_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
res.fe <- rma(d, vd, data=df, method="FE")
est.fe <- c(est.fe,res.fe$b)
est.re <- c(est.re,NA)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("woodin_subset.csv", header=T, sep=";")
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$d.o <- df$effest.exp / J.o                              # cohen's d
df$d.c <- df$effestnew.exp / J.c                           # cohen's d
df$vd.o <- df$n / ((df$n / 2)^2 + (df$d.o^2 / (2 * df$n))) # variance cohen's d
df$vd.c <- df$nnew / ((df$nnew / 2)^2 + (df$d.c^2 / (2 * df$nnew)))  # variance cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c 
write.table(df, file = "woodin_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vg.o, data=df, method="FE")  

# subset reproduced
res.sc <- rma(effestnew.exp, vg.c, data=df, method="FE")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[31]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

# Woodley -----------------------------------------------------------------
df <- read.table("woodley_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "woodley_complete.csv", row.names=FALSE, sep=";")

# meta-analysis on reported values
res.fe <- rma(r, vr, data=df, method="FE")
res.re <- rma(r, vr, data=df, method="DL")
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re$b)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("woodley_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1) 
write.table(df, file = "woodley_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(effest.exp, vr.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(effestnew.exp, vr.c, data=df, method="DL")               

eff.so <- c(eff.so,res.so$b) # MA subset original effect size estimate
cilb.so <- c(cilb.so,res.so$ci.lb) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,res.so$ci.ub) # upperbound
tau2.so <- c(tau2.so,res.so$tau2) # tau2 estimate

eff.sc <- c(eff.sc,res.sc$b) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,res.sc$ci.lb)
ciub.sc <- c(ciub.sc,res.sc$ci.ub)
tau2.sc <- c(tau2.sc,res.sc$tau2)

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[32]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)

 # Yoon --------------------------------------------------------------------
df <- read.table("yoon_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                  # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "yoon_complete.csv", row.names=FALSE, sep=";")

# meta-analysis with reported values
df$z <- 0.5 * log((1 + df$r) / (1 - df$r))
df$vz <- 1 / (df$n - 3)    
res.re <- rma(z, vz, data=df, method="DL") 
# transform estimate back to correlation
res.re <- tanh(res.re$b)
est.fe <- c(est.fe,NA)
est.re <- c(est.re,res.re)
k.tot <- c(k.tot,nrow(df))
rm(df)

# subset
df <- read.table("yoon_subset.csv", header=T, sep=";")
df$vr.o <- ((1 - (df$effest.exp^2))^2) / (df$n-1)             # variance correlation r
df$vr.c <- ((1 - (df$effestnew.exp^2))^2) / (df$nnew-1) 
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$vz.o <- 1 / (df$n - 3) 
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.c <- 1 / (df$nnew - 3) 
write.table(df, file = "yoon_subset.csv", row.names=FALSE, sep=";")

# subset original 
res.so <- rma(z.o, vz.o, data=df, method="DL") 

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="DL")               

# transform estimate back to correlation
eff.so <- c(eff.so,tanh(res.so$b)) # MA subset original effect size estimate
cilb.so <- c(cilb.so,tanh(res.so$ci.lb)) # MA subset original effect size CI lowerbound
ciub.so <- c(ciub.so,tanh(res.so$ci.ub)) # upperbound
tau2.so <- c(tau2.so,tanh(res.so$tau2)) # tau2 estimate

eff.sc <- c(eff.sc,tanh(res.sc$b)) # MA subset checked effect size estimate
cilb.sc <- c(cilb.sc,tanh(res.sc$ci.lb))
ciub.sc <- c(ciub.sc,tanh(res.sc$ci.ub))
tau2.sc <- c(tau2.sc,tanh(res.sc$tau2))

pval.so <- c(pval.so,res.so$pval) # p value original
pval.sc <- c(pval.sc,res.sc$pval) # p value checked
krep1 <- c(krep1,nrow(df))
percent <- c(percent,(nrow(df)/k.tot[33]*100))
krep2 <- c(krep2,sum(df$disccat.eff == 0 & df$info == 0))
percent2 <- c(percent2,sum(df$disccat.eff == 0 & df$info == 0) / nrow(df) * 100)

rm(df)


# Calculate discrepancies and fill in table -------------------------------
dat$k <- k.tot
dat$recalc.fe <- est.fe
dat$recalc.re <- est.re
dat$eff.so <- eff.so
dat$cilb.so <- cilb.so
dat$ciub.so <- ciub.so
dat$tau2.so <- tau2.so
dat$eff.sc <- eff.sc
dat$cilb.sc <- cilb.sc
dat$ciub.sc <- ciub.sc
dat$tau2.sc <- tau2.sc
dat$krep1 <- krep1
dat$krep2 <- krep2
dat$percent <- percent
dat$percent2 <- percent2

# Discrepancies
dat$disc.s <- abs(eff.so) - abs(eff.sc)
dat$disc.cilb.s <- abs(cilb.so) - abs(cilb.sc)
dat$disc.ciub.s <- abs(ciub.so) - abs(ciub.sc)
dat$disc.tau2.s <- abs(tau2.so) - abs(tau2.sc)

# Are there any p-values which change from significant to non-sig
round(pval.so,2)
round(pval.sc,2)


# Effect size discrepancies -----------------------------------------------
for (i in 1:nrow(dat)) {
  
  if (dat$efftype[i] == "g") {
    
    if (abs(dat$disc.s[i]) < 0.05) {
      dat$disccat.s[i] = 0
    } else if (abs(dat$disc.s[i]) >= 0.050 & abs(dat$disc.s[i]) <= 0.149) {
      dat$disccat.s[i] = 1
    } else if (abs(dat$disc.s[i]) > 0.149 & abs(dat$disc.s[i]) <= 0.249) {
      dat$disccat.s[i] = 2
    } else if (abs(dat$disc.s[i]) > 0.249) {
      dat$disccat.s[i] = 3
    } else {
      dat$disccat.s[i] = "check"
    }
  } 
  
  if (dat$efftype[i] == "d") {
    
    if (abs(dat$disc.s[i]) < 0.051) {
      dat$disccat.s[i] = 0
    } else if (abs(dat$disc.s[i]) >= 0.051 & abs(dat$disc.s[i]) <= 0.151) {
      dat$disccat.s[i] = 1
    } else if (abs(dat$disc.s[i]) > 0.151 & abs(dat$disc.s[i]) <= 0.252) {
      dat$disccat.s[i] = 2
    } else if (abs(dat$disc.s[i]) > 0.252) {
      dat$disccat.s[i] = 3
    } else {
      dat$disccat.s[i] = "check"
    }
    
  }
  
  if (dat$efftype[i] == "r" | dat$efftype[i] == "z") {
    
    if (abs(dat$disc.s[i]) < 0.025) {
      dat$disccat.s[i] = 0
    } else if (abs(dat$disc.s[i]) >= 0.025 & abs(dat$disc.s[i]) <= 0.075) {
      dat$disccat.s[i] = 1
    } else if (abs(dat$disc.s[i]) > 0.075 & abs(dat$disc.s[i]) <= 0.125) {
      dat$disccat.s[i] = 2
    } else if (abs(dat$disc.s[i]) > 0.125) {
      dat$disccat.s[i] = 3
    } else {
      dat$disccat.s[i] = "check"
    }
  }
}


# Confidence Interval discrepancies ---------------------------------------
for (i in 1:nrow(dat)) {
  
  if (dat$efftype[i] == "g") {
    
    # Fill in discrepancy category for CI of effect sizes on subset checked MAs
    if (abs(dat$disc.cilb.s[i]) < 0.050 & abs(dat$disc.ciub.s[i]) < 0.050) {
      dat$disccat.ci.s[i] = 0
    } else if (abs(dat$disc.cilb.s[i]) > 0.249 | abs(dat$disc.ciub.s[i]) > 0.249) {
      dat$disccat.ci.s[i] = 3
    } else if (abs(dat$disc.cilb.s[i]) > 0.149 & abs(dat$disc.cilb.s[i]) <= 0.249 | abs(dat$disc.ciub.s[i]) > 0.149 & abs(dat$disc.ciub.s[i]) <= 0.249) {
      dat$disccat.ci.s[i] = 2
    } else if (abs(dat$disc.cilb.s[i]) >= 0.050 & abs(dat$disc.cilb.s[i]) <= 0.149 | abs(dat$disc.ciub.s[i]) >= 0.050 & abs(dat$disc.ciub.s[i]) <= 0.149) {
      dat$disccat.ci.s[i] = 1
    } else {
      dat$disccat.ci.s[i] == "check"
    }
  }
    
  if (dat$efftype[i] == "d") {
    
    # Fill in discrepancy category for CI of effect sizes on subset checked MAs
    if (abs(dat$disc.cilb.s[i]) < 0.051 & abs(dat$disc.ciub.s[i]) < 0.051) {
      dat$disccat.ci.s[i] = 0
    } else if (abs(dat$disc.cilb.s[i]) > 0.252 | abs(dat$disc.ciub.s[i]) > 0.252) {
      dat$disccat.ci.s[i] = 3
    } else if (abs(dat$disc.cilb.s[i]) > 0.151 & abs(dat$disc.cilb.s[i]) <= 0.252 | abs(dat$disc.ciub.s[i]) > 0.151 & abs(dat$disc.ciub.s[i]) <= 0.252) {
      dat$disccat.ci.s[i] = 2
    } else if (abs(dat$disc.cilb.s[i]) >= 0.051 & abs(dat$disc.cilb.s[i]) <= 0.151 | abs(dat$disc.ciub.s[i]) >= 0.051 & abs(dat$disc.ciub.s[i]) <= 0.151) {
      dat$disccat.ci.s[i] = 1
    } else {
      dat$disccat.ci.s[i] == "check"
    }
  }
  
  if (dat$efftype[i] == "r" | dat$efftype[i] == "z") {
    
    # Fill in discrepancy category for CI of effect sizes on subset checked MAs
    if (abs(dat$disc.cilb.s[i]) < 0.025 & abs(dat$disc.ciub.s[i]) < 0.025) {
      dat$disccat.ci.s[i] = 0
    } else if (abs(dat$disc.cilb.s[i]) > 0.125 | abs(dat$disc.ciub.s[i]) > 0.125) {
      dat$disccat.ci.s[i] = 3
    } else if (abs(dat$disc.cilb.s[i]) > 0.075 & abs(dat$disc.cilb.s[i]) <= 0.125 | abs(dat$disc.ciub.s[i]) > 0.075 & abs(dat$disc.ciub.s[i]) <= 0.125) {
      dat$disccat.ci.s[i] = 2
    } else if (abs(dat$disc.cilb.s[i]) >= 0.025 & abs(dat$disc.cilb.s[i]) <= 0.075 | abs(dat$disc.ciub.s[i]) >= 0.025 & abs(dat$disc.ciub.s[i]) <= 0.075) {
      dat$disccat.ci.s[i] = 1
    } else {
      dat$disccat.ci.s[i] == "check"
    }
  }
}


# Tau2 discrepancies ------------------------------------------------------
for (i in 1:nrow(dat)) {
  
  if (dat$efftype[i] == "g") {
    
    if (abs(dat$disc.tau2.s[i]) < 0.05) {
      dat$disccat.tau2.s[i] = 0
    } else if (abs(dat$disc.tau2.s[i]) >= 0.050 & abs(dat$disc.tau2.s[i]) <= 0.149) {
      dat$disccat.tau2.s[i] = 1
    } else if (abs(dat$disc.tau2.s[i]) > 0.149 & abs(dat$disc.tau2.s[i]) <= 0.249) {
      dat$disccat.tau2.s[i] = 2
    } else if (abs(dat$disc.tau2.s[i]) > 0.249) {
      dat$disccat.tau2.s[i] = 3
    } else {
      dat$disccat.tau2.s[i] = "check"
    }
  } 
  
  if (dat$efftype[i] == "d") {
    
    if (abs(dat$disc.tau2.s[i]) < 0.051) {
      dat$disccat.tau2.s[i] = 0
    } else if (abs(dat$disc.tau2.s[i]) >= 0.051 & abs(dat$disc.tau2.s[i]) <= 0.151) {
      dat$disccat.tau2.s[i] = 1
    } else if (abs(dat$disc.tau2.s[i]) > 0.151 & abs(dat$disc.tau2.s[i]) <= 0.252) {
      dat$disccat.tau2.s[i] = 2
    } else if (abs(dat$disc.tau2.s[i]) > 0.252) {
      dat$disccat.tau2.s[i] = 3
    } else {
      dat$disccat.tau2.s[i] = "check"
    }
    
  }
  
  if (dat$efftype[i] == "r" | dat$efftype[i] == "z") {
    
    if (abs(dat$disc.tau2.s[i]) < 0.025) {
      dat$disccat.tau2.s[i] = 0
    } else if (abs(dat$disc.tau2.s[i]) >= 0.025 & abs(dat$disc.tau2.s[i]) <= 0.075) {
      dat$disccat.tau2.s[i] = 1
    } else if (abs(dat$disc.tau2.s[i]) > 0.075 & abs(dat$disc.tau2.s[i]) <= 0.125) {
      dat$disccat.tau2.s[i] = 2
    } else if (abs(dat$disc.tau2.s[i]) > 0.125) {
      dat$disccat.tau2.s[i] = 3
    } else {
      dat$disccat.tau2.s[i] = "check"
    }
  }
}

# dataset with all studies
regular.total[,"type"] <- "r"
outlier.total[,"type"] <- "o"
dftot <- rbind(regular.total,outlier.total)
dftot <- dftot[order(dftot$id),]

# Further calculations ----------------------------------------------------
k.tot <- nrow(outlier.total)+nrow(regular.total);k.tot     # number of total studies 1951
nrow(outlier.total);nrow(regular.total)                    # 581 outliers 1370 regulars
k.out <- nrow(outlier.total)/k.tot*100; k.out              # perc of outlier studies of total k 30.1%
k.reg <- nrow(regular.total)/k.tot*100; k.reg              # perc of regular studies of total k 69.9%

# save results ------------------------------------------------------------
setwd("C:/Users/s421506/tiu/research/effectsizes/codebooks/")
write.table(dat, file = "codebook-meta-analyses-final-complete.csv", row.names=F, col.names=T, sep=' ')
write.table(outlier.total, file = "studies-outliers.csv", row.names=FALSE, sep=";")
write.table(regular.total, file = "studies-nonoutliers.csv", row.names=FALSE, sep=";")
write.table(dftot, file = "studies-all.csv", row.names=FALSE, sep=";")

# Save file in xlsx -------------------------------------------------------
write_xlsx(dat,"codebook-meta-analyses-final-complete.xlsx",col_names=T)



# modyop <- rma(znewexp, newvz, data=df, mods = yop) # year of publication moderator analysis
# dat$yop[i] <- modyop$b[2]                           # year of publication parameter estimate
# dat$yopse[i] <- modyop$se[2]                        # year of publication parameter standard error
# dat$yopsig[i] <- as.numeric(modyop$pval[2] < .05)   # year of publication significant
# 
# modfunnel <- rma(znewexp, newvz, data=df, mods = sqrt(newvz)) # standard error moderator analysis
# dat$funnel[i] <- modfunnel$b[2]                     # standard error parameter estimate
# dat$funnelse[i] <- modfunnel$se[2]                  # standard error parameter standard error
# dat$funnelsig[i] <- as.numeric(modfunnel$pval[2] < .05) # standard error significant
# 
# modboth <- rma(znewexp, newvz, data=df, mods = cbind(yop, sqrt(newvz))) # year of publication and standard error moderator analysis
# dat$bothyop[i] <- modboth$b[2]                      # year of publication parameter estimate
# dat$bothyopse[i] <- modboth$se[2]                   # year of publication parameter standard error
# dat$bothfunnel[i] <- modboth$b[3]                   # year of publication significant
# dat$bothfunnelse[i] <- modboth$se[3]                # standard errorparameter estimate
# dat$bothyopsig[i] <- as.numeric(modboth$pval[2] < .05) # standard error parameter standard error
# dat$bothfunnelsig[i] <- as.numeric(modboth$pval[3] < .05) # standard error significant
# 
# cortest <- cor.test(df$yop, sqrt(df$newvz))      # correlation year of publication and standard error
# dat$cor[i] <- cortest$estimate                      # correlation year of publication and standard error
#
#l1o <- leave1out(res.sc)              # leave-one-out analysis
#estdiff <- l1o$estimate - c(res.sc$b)  # MA pooled effect size estimates if ith study was omitted
#dat$outeff1.s[i] <- length(estdiff[abs(estdiff) >= 0.020 & abs(estdiff) <= 0.056])  # no of studies that have a small effect on MA effect size
#dat$outeff2.s[i] <- length(estdiff[abs(estdiff) > 0.056 & abs(estdiff) <= 0.093])  # no of studies that have a medium effect on MA effect size
#dat$outeff3.s[i] <- length(estdiff[abs(estdiff) > 0.093])                          # no of studies that have a large effect on MA effect size
#
#ci.lb.diff <- l1o$ci.lb - c(res.sc$ci.lb) # CI LB of effect size estimates if ith study was omitted   
#ci.ub.diff <- l1o$ci.ub - c(res.sc$ci.ub) # CI UB of effect size estimates if ith study was omitted 
#dat$outci1.s[i] <- max(length(ci.lb.diff[abs(ci.lb.diff) >= 0.020 & abs(ci.lb.diff) <= 0.056]),length(ci.ub.diff[abs(ci.ub.diff) >= 0.020 & abs(ci.ub.diff) <= 0.056])) # no of studies that have a small effect on MA effect size CI
#dat$outci2.s[i] <- max(length(ci.lb.diff[abs(ci.lb.diff) > 0.056 & abs(ci.lb.diff) <= 0.093]),length(ci.ub.diff[abs(ci.ub.diff) > 0.056 & abs(ci.ub.diff) <= 0.093])) # no of studies that have a medium effect on MA effect size CI
#dat$outci3.s[i] <- max(length(ci.lb.diff[abs(ci.lb.diff) > 0.093]),length(ci.ub.diff[abs(ci.ub.diff) > 0.093]))                                             # no of studies that have a large effect on MA effect size CI
#
#taudiff <- l1o$tau2 - c(res.sc$tau2) # tau2 estimate if ith study was omitted
#dat$outtau1.s[i] <- length(taudiff[abs(taudiff)>= 0.020 & abs(taudiff)<= 0.056])  # no of studies that have a small effect on tau2 estimate
#dat$outtau2.s[i] <- length(taudiff[abs(taudiff)> 0.056 & abs(taudiff)<= 0.093])  # no of studies that have a medium effect on tau2 estimate
#dat$outtau3.s[i] <- length(taudiff[abs(taudiff)> 0.093])                         # no of studies that have a large effect on tau2 estimate
      
      