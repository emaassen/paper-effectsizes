rm(list = ls())
setwd("C:/Users/s421506/tiu/research/effectsizes/codebooks/")
x <- c("metafor","readxl")
#lapply(x,install.packages(x),character.only=T)
lapply(x,library,character.only=T)

#dat <- read.table("codebook-meta-analyses-final.csv", header=T, sep = '')
dat <- read_excel("codebook-meta-analyses-final-empty.xlsx", 2)

# Recalculated effect size based on original estimates --------------------
setwd("C:/Users/s421506/tiu/research/effectsizes/data-per-ma")

est <- k <- c() # empty vector to store results

# Adesope
df <- read.table("adesope_subset.csv", header=T, sep=';')
res.fe <- rma(g, vg, data=df, method="FE")
res.re <- rma(g, vg, data=df) 
est <- c(est,paste("FE ",res.fe$b,"RE",res.re$b))
k <- c(k,nrow(df))
rm(df)

#Alfieri
df <- read.table("alfieri_subset.csv", header=T, sep=';')
res <- rma(g, vg, data=df)
est <- c(est,res$b)
k <- c(k,nrow(df))
rm(df)







#Babbage
df <- read.table("babbage_complete.csv", header=T, sep=';')
res <- rma(g, vg, data=df)
dat[3,"recalc"] <- res$b
dat[3, "k"] <- nrow(df)
rm(df)

#Balliet
df <- read.table("balliet_complete.csv", header=T, sep=';')
res <- rma(d, vd, data=df)
dat[4,"recalc"] <- res$b
dat[4, "k"] <- nrow(df)
rm(df)

#Benish
df <- read.table("benish_complete.csv", header=T, sep=';')
res <- rma(d, vd, data=df)
dat[5,"recalc"] <- res$b
dat[5, "k"] <- nrow(df)
rm(df)

#Berry1
df <- read.table("berry1_complete.csv", header=T, sep=';')
res.fe <- rma(r, vr, data=df, method="FE")
res.re <- rma(r, vr, data=df)
est <- paste("FE ",res.fe$b,"RE",res.re$b)
dat[6,"recalc"] <- est
dat[6, "k"] <- nrow(df)
rm(df)

#Berry2
df <- read.table("berry2_complete.csv", header=T, sep=';')
res.fe <- rma(r, vr, data=df, method="FE")
res.re <- rma(r, vr, data=df)
est <- paste("FE ",res.fe$b,"RE",res.re$b)
dat[7,"recalc"] <- est
dat[7, "k"] <- nrow(df)
rm(df)

#Card
df <- read.table("card_complete.csv", header=T, sep=';')
res <- rma(r, vr, data=df)
dat[8,"recalc"] <- res$b
dat[8, "k"] <- nrow(df)
rm(df)

#Crook
df <- read.table("crook_complete.csv", header=T, sep=';')
res.fe <- rma(r, vr, data=df,method="FE")
res.re <- rma(r, vr, data=df)
est <- paste("FE ",res.fe$b,"RE",res.re$b)
dat[9,"recalc"] <- est
dat[9, "k"] <- nrow(df)
rm(df)

#DeWit
df <- read.table("dewit_complete.csv", header=T, sep=';')
res.fe <- rma(r, vr, data=df,method="FE")
res.re <- rma(r, vr, data=df)
est <- paste("FE ",res.fe$b,"RE",res.re$b)
dat[10,"recalc"] <- est
dat[10, "k"] <- nrow(df)
rm(df)

#Elsequest
df <- read.table("elsequest_complete.csv", header=T, sep=';')
res <- rma(d, vd, data=df)
dat[11,"recalc"] <- res$b
dat[11, "k"] <- nrow(df)
rm(df)

#Farber
df <- read.table("farber_complete.csv", header=T, sep=';')
res <- rma(r, vr, data=df)
dat[12,"recalc"] <- res$b
dat[12, "k"] <- nrow(df)
rm(df)

#Fischer
df <- read.table("fischer_complete.csv", header=T, sep=';')
res.fe <- rma(g, vg, data=df, method="FE")
res.re <- rma(g, vg, data=df)
est <- paste("FE ",res.fe$b,"RE",res.re$b)
dat[13,"recalc"] <- est
dat[13, "k"] <- nrow(df)
rm(df)

#Fox
df <- read.table("fox_complete.csv", header=T, sep=';')
res <- rma(r, vr, data=df)
dat[14,"recalc"] <- res$b
dat[14, "k"] <- nrow(df)
rm(df)

#Freund
df <- read.table("freund_complete.csv", header=T, sep=';')
res <- rma(r, vr, data=df)
dat[15,"recalc"] <- res$b
dat[15, "k"] <- nrow(df)
rm(df)

#Green
df <- read.table("green_complete.csv", header=T, sep=';')
res <- rma(d, vd, data=df)
dat[16,"recalc"] <- res$b
dat[16, "k"] <- nrow(df)
rm(df)

#Hallion
df <- read.table("hallion_complete.csv", header=T, sep=';')
res <- rma(g, vg, data=df)
dat[17,"recalc"] <- res$b
dat[17, "k"] <- nrow(df)
rm(df)

#Ihle
df <- read.table("ihle_complete.csv", header=T, sep=';')
res <- rma(d, vd, data=df)
dat[18,"recalc"] <- res$b
dat[18, "k"] <- nrow(df)
rm(df)

#Koenig
df <- read.table("koenig_complete.csv", header=T, sep=';')
res <- rma(g, vg, data=df)
dat[19,"recalc"] <- res$b
dat[19, "k"] <- nrow(df)
rm(df)

#Kolden
df <- read.table("kolden_complete.csv", header=T, sep=';')
res <- rma (r, vr, data=df)
dat[20,"recalc"] <- res$b
dat[20, "k"] <- nrow(df)
rm(df)

#Lucassen
df <- read.table("lucassen_complete.csv", header=T, sep=';')
res <- rma (r, vr, data=df)
dat[21,"recalc"] <- res$b
dat[21, "k"] <- nrow(df)
rm(df)

#Mol
df <- read.table("mol_complete.csv", header=T, sep=';')
res <- rma (z, vz, data=df)
dat[22,"recalc"] <- res$b
dat[22, "k"] <- nrow(df)
rm(df)

#Morgan
df <- read.table("morgan_complete.csv", header=T, sep=';')
res <- rma(d, vd, data=df)
dat[23,"recalc"] <- res$b
dat[23, "k"] <- nrow(df)
rm(df)

#Munder
df <- read.table("munder_complete.csv", header=T, sep=';')
res <- rma(d, vd, data=df)
dat[24,"recalc"] <- res$b
dat[24, "k"] <- nrow(df)
rm(df)

#Piet
df <- read.table("piet_complete.csv", header=T, sep=';')
res <- rma(g, vg, data=df)
dat[25,"recalc"] <- res$b
dat[25, "k"] <- nrow(df)
rm(df)

#Smith
df <- read.table("smith_complete.csv", header=T, sep=';')
res <- rma (r, vr, data=df)
dat[26,"recalc"] <- res$b
dat[26, "k"] <- nrow(df)
rm(df)

#Tillman
df <- read.table("tillman_complete.csv", header=T, sep=';')
res <- rma(r, vr, data=df)
dat[27,"recalc"] <- res$b
dat[27, "k"] <- nrow(df)
rm(df)

#Toosi
df <- read.table("toosi_complete.csv", header=T, sep=';')
res <- rma(r, vr, data=df)
dat[28,"recalc"] <- res$b
dat[28, "k"] <- nrow(df)
rm(df)

#vanIddekinge
df <- read.table("vaniddekinge_complete.csv", header=T, sep=';')
res.fe <- rma(r, vr, data=df, method="FE")
res.re <- rma(r, vr, data=df)
est <- paste("FE ",res.fe$b,"RE",res.re$b)
dat[29,"recalc"] <- est
dat[29, "k"] <- nrow(df)
rm(df)

#Webb
df <- read.table("webb_complete.csv", header=T, sep=';')
res <- rma(d, vd, data=df)
dat[30,"recalc"] <- res$b
dat[30, "k"] <- nrow(df)
rm(df)

#Woodin
df <- read.table("woodin_complete.csv", header=T, sep=';')
res <- rma(d, vd, data=df, method="FE")
dat[31,"recalc"] <- res$b
dat[31, "k"] <- nrow(df)
rm(df)

#Woodley
df <- read.table("woodley_complete.csv", header=T, sep=';')
res.fe <- rma(r, vr, data=df, method="FE")
res.re <- rma(r, vr, data=df)
est <- paste("FE ",res.fe$b,"RE",res.re$b)
dat[32,"recalc"] <- est
dat[32, "k"] <- nrow(df)
rm(df)

#Yoon
df <- read.table("yoon_complete.csv", header=T, sep=';')
res <- rma(r, vr, data=df)
dat[33,"recalc"] <- res$b
dat[33, "k"] <- nrow(df)
rm(df)




dat[[i],"recalc"] <- est

# Save results k + recalc column
setwd("C:/Users/s421506/tiu/research/effectsizes/")
write.table(dat, file = "codebook-meta-analyses-final2.csv", row.names=F, col.names=T, sep=' ')

# Run all meta-analyses 
setwd("C:/Users/s421506/tiu/research/effectsizes/data-per-ma")
authors_complete <- list.files("C:/Users/s421506/tiu/research/effectsizes/data-per-ma", pattern="_complete.csv")
authors_subset <- list.files("C:/Users/s421506/tiu/research/effectsizes/data-per-ma", pattern="_subset.csv")

pval.co <- pval.cc <- pval.so <- pval.sc <- c()

for (i in 1:nrow(dat)) {
  
  # Complete original MAs ---------------------------------------------------
  df <- read.table(paste(authors_complete[i]), header=T, sep=';')
  
  #trycatch because Card doesn't converge with method REML, so ML is used. All other meta-analyses use method=REML
  res.co = tryCatch ({ res.co = rma(zexp, vz, data=df) }, error = function() { res.co = rma(zexp, vz, data=df, method="ML") })
  
  dat$z.co[i] <- res.co$b           # MA complete original effect size estimate
  dat$ci.li.co[i] <- res.co$ci.lb   # MA complete original effect size LB
  dat$ci.ul.co[i] <- res.co$ci.ub   # MA complete original effect size UB
  dat$tau.co[i] <- res.co$tau2      # MA complete original tau2 estimate
  
  # Complete checked MAs ----------------------------------------------------
  df$newvz <- 1 / (df$nnew-3) 
  write.table(df, file = paste(authors_complete[i]), row.names=FALSE, sep=";")
  
  res.cc = tryCatch ({ res.cc = rma(znewexp, newvz, data=df) }, error = function(err) { res.cc = rma(znewexp, newvz, data=df, method="ML") })
  
  dat$z.cc[i] <- res.cc$b           # MA complete checked effect size estimate
  dat$ci.li.cc[i] <- res.cc$ci.lb   # MA complete checked effect size LB
  dat$ci.ul.cc[i] <- res.cc$ci.ub   # MA complete checked effect size LB
  dat$tau.cc[i] <- res.cc$tau2      # MA complete checked tau2 estimate
  dat$disc.c[i] <- res.co$b - res.cc$b # Difference complete original and complete checked effect size estimate
  
  pval.co <- c(pval.co,res.co$pval)
  pval.cc <- c(pval.cc,res.cc$pval)
  
  # Fill in discrepancy category for effect sizes on complete checked MAs
  if (abs(dat$disc.c[i]) < 0.020) {
    dat$disccat.c[i] = 0
  } else if (abs(dat$disc.c[i]) >= 0.020 & abs(dat$disc.c[i]) <= 0.056) {
    dat$disccat.c[i] = 1
  } else if (abs(dat$disc.c[i]) > 0.056 & abs(dat$disc.c[i]) <= 0.093) {
    dat$disccat.c[i] = 2
  } else if (abs(dat$disc.c[i]) > 0.093) {
    dat$disccat.c[i] = 3
  } else {
    dat$disccat.c[i] = "check"
  }
  
  dat$disc.ci.li.c[i] <- res.co$ci.lb - res.cc$ci.lb  # difference complete original and complete checked CI LB of effect size estimate 
  dat$disc.ci.ul.c[i] <- res.co$ci.ub - res.cc$ci.ub  # difference complete original and complete checked CI UB of effect size estimate
  
  # Fill in discrepancy category for CI of effect sizes on complete checked MAs
  if (abs(dat$disc.ci.li.c[i]) < 0.020 & abs(dat$disc.ci.ul.c[i]) < 0.020) {
    dat$disccat.ci.c[i] = 0
  } else if (abs(dat$disc.ci.li.c[i]) > 0.093 | abs(dat$disc.ci.ul.c[i]) > 0.093) {
    dat$disccat.ci.c[i] = 3
  } else if (abs(dat$disc.ci.li.c[i]) > 0.056 & abs(dat$disc.ci.li.c[i]) <= 0.093 | abs(dat$disc.ci.ul.c[i]) > 0.056 & abs(dat$disc.ci.ul.c[i]) <= 0.093) {
    dat$disccat.ci.c[i] = 2
  } else if (abs(dat$disc.ci.li.c[i]) >= 0.020 & abs(dat$disc.ci.li.c[i]) <= 0.056 | abs(dat$disc.ci.ul.c[i]) >= 0.020 & abs(dat$disc.ci.ul.c[i]) <= 0.056) {
    dat$disccat.ci.c[i] = 1
  } else {
    dat$disccat.ci.c[i] == "check"
  }
  
  dat$disc.tau.c[i] <- res.co$tau2 - res.cc$tau2 # difference complete original and complete checked tausquared
  
  # Fill in discrepancy category for tau2 on complete checked MAs
  if (abs(dat$disc.tau.c[i]) < 0.020) {
    dat$disccat.tau.c[i] = 0
  } else if (abs(dat$disc.tau.c[i]) >= 0.020 & abs(dat$disc.tau.c[i]) <= 0.056) {
    dat$disccat.tau.c[i] = 1
  } else if (abs(dat$disc.tau.c[i]) > 0.056 & abs(dat$disc.tau.c[i]) <= 0.093) {
    dat$disccat.tau.c[i] = 2
  } else if (abs(dat$disc.tau.c[i]) > 0.093) {
    dat$disccat.tau.c[i] = 3
  } else {
    dat$disccat.tau.c[i] = "check"
  }
  
  l1o <- leave1out(res.cc)              # leave-one-out analysis
  estdiff <- l1o$estimate - c(res.cc$b)  # MA pooled effect size estimates if ith study was omitted
  dat$outeff1.c[i] <- length(estdiff[abs(estdiff) >= 0.020 & abs(estdiff) <= 0.056])  # no of studies that have a small effect on MA effect size
  dat$outeff2.c[i] <- length(estdiff[abs(estdiff) > 0.056 & abs(estdiff) <= 0.093])  # no of studies that have a medium effect on MA effect size
  dat$outeff3.c[i] <- length(estdiff[abs(estdiff) > 0.093])                     # no of studies that have a large effect on MA effect size
  
  ci.lb.diff <- l1o$ci.lb - c(res.cc$ci.lb) # CI LB of effect size estimates if ith study was omitted   
  ci.ub.diff <- l1o$ci.ub - c(res.cc$ci.ub) # CI UB of effect size estimates if ith study was omitted 
  dat$outci1.c[i] <- max(length(ci.lb.diff[abs(ci.lb.diff) >= 0.020 & abs(ci.lb.diff) <= 0.056]),length(ci.ub.diff[abs(ci.ub.diff) >= 0.020 & abs(ci.ub.diff) <= 0.056])) # no of studies that have a small effect on MA effect size CI
  dat$outci2.c[i] <- max(length(ci.lb.diff[abs(ci.lb.diff) > 0.056 & abs(ci.lb.diff) <= 0.093]),length(ci.ub.diff[abs(ci.ub.diff) > 0.056 & abs(ci.ub.diff) <= 0.093])) # no of studies that have a medium effect on MA effect size CI
  dat$outci3.c[i] <- max(length(ci.lb.diff[abs(ci.lb.diff) > 0.093]),length(ci.ub.diff[abs(ci.ub.diff) > 0.093]))                                             # no of studies that have a large effect on MA effect size CI
  
  taudiff <- l1o$tau2 - c(res.cc$tau2) # tau2 estimate if ith study was omitted
  dat$outtau1.c[i] <- length(taudiff[abs(taudiff)>= 0.020 & abs(taudiff)<= 0.056])  # no of studies that have a small effect on tau2 estimate
  dat$outtau2.c[i] <- length(taudiff[abs(taudiff)> 0.056 & abs(taudiff)<= 0.093])  # no of studies that have a medium effect on tau2 estimate
  dat$outtau3.c[i] <- length(taudiff[abs(taudiff)> 0.093])                     # no of studies that have a large effect on tau2 estimate
  
  modyop <- rma(znewexp, newvz, data=df, mods = yop) # year of publication moderator analysis
  dat$yop[i] <- modyop$b[2]                           # year of publication parameter estimate
  dat$yopse[i] <- modyop$se[2]                        # year of publication parameter standard error
  dat$yopsig[i] <- as.numeric(modyop$pval[2] < .05)   # year of publication significant
  
  modfunnel <- rma(znewexp, newvz, data=df, mods = sqrt(newvz)) # standard error moderator analysis
  dat$funnel[i] <- modfunnel$b[2]                     # standard error parameter estimate
  dat$funnelse[i] <- modfunnel$se[2]                  # standard error parameter standard error
  dat$funnelsig[i] <- as.numeric(modfunnel$pval[2] < .05) # standard error significant
  
  modboth <- rma(znewexp, newvz, data=df, mods = cbind(yop, sqrt(newvz))) # year of publication and standard error moderator analysis
  dat$bothyop[i] <- modboth$b[2]                      # year of publication parameter estimate
  dat$bothyopse[i] <- modboth$se[2]                   # year of publication parameter standard error
  dat$bothfunnel[i] <- modboth$b[3]                   # year of publication significant
  dat$bothfunnelse[i] <- modboth$se[3]                # standard errorparameter estimate
  dat$bothyopsig[i] <- as.numeric(modboth$pval[2] < .05) # standard error parameter standard error
  dat$bothfunnelsig[i] <- as.numeric(modboth$pval[3] < .05) # standard error significant
  
  cortest <- cor.test(df$yop, sqrt(df$newvz))      # correlation year of publication and standard error
  dat$cor[i] <- cortest$estimate                      # correlation year of publication and standard error

  # Subset original MAs -----------------------------------------------------
  df <- read.table(paste(authors_subset[i]), header=T, sep=';')
  
  res.so = tryCatch ({ res.so = rma(zexp, vz, data=df) }, error = function(err) { res.so = rma(zexp, vz, data=df, method="ML") })
  
  dat$z.so[i] <- res.so$b           # MA subset original effect size estimate
  dat$ci.li.so[i] <- res.so$ci.lb   # MA subset original effect size LB
  dat$ci.ul.so[i] <- res.so$ci.ub   # MA subset original effect size UB
  dat$tau.so[i] <- res.so$tau2      # MA subset original tau2 estimate
  
  dat$krep1[i] <- nrow(df)             # no of reproduced studies
  dat$percent[i] <- nrow(df)/dat$k[i]  # percentage of studies from original MA that were reproduced
  dat$krep2[i] <- sum(as.numeric(df$disccat == 0)) # no of reproduced studies that could be reproduced without a discrepancy
  dat$percent2[i] <- sum(as.numeric(df$disccat == 0)) / nrow(df) # percentage of studies that could be reproduced without a discrepancy, out of all reproduced studies
  
  # Subset checked MAs ------------------------------------------------------
  df$newvz <- 1 / (df$nnew-3) 
  write.table(df, file = paste(authors_subset[i]), row.names=FALSE, sep=";")
  
  res.sc = tryCatch ({ res.sc = rma(znewexp, newvz, data=df) }, error = function(err) { res.sc = rma(znewexp, newvz, data=df, method="ML") })
  
  dat$z.sc[i] <- res.sc$b           # MA subset checked effect size estimate
  dat$ci.li.sc[i] <- res.sc$ci.lb   # MA subset checked effect size LB
  dat$ci.ul.sc[i] <- res.sc$ci.ub   # MA subset checked effect size LB
  dat$tau.sc[i] <- res.sc$tau2      # MA subset checked tau2 estimate
  dat$disc.s[i] <- res.so$b - res.sc$b # Difference subset original and subset checked effect size estimate
  
  pval.so <- c(pval.so,res.so$pval)
  pval.sc <- c(pval.sc,res.sc$pval)
  
  # Fill in discrepancy category for effect sizes on subset checked MAs
  if (abs(dat$disc.s[i]) < 0.020) {
    dat$disccat.s[i] = 0
  } else if (abs(dat$disc.s[i]) >= 0.020 & abs(dat$disc.s[i]) <= 0.056) {
    dat$disccat.s[i] = 1
  } else if (abs(dat$disc.s[i]) > 0.056 & abs(dat$disc.s[i]) <= 0.093) {
    dat$disccat.s[i] = 2
  } else if (abs(dat$disc.s[i]) > 0.093) {
    dat$disccat.s[i] = 3
  } else {
    dat$disccat.s[i] = "check"
  }
  
  dat$disc.ci.li.s[i] <- res.so$ci.lb - res.sc$ci.lb  # difference subset original and complete checked CI LB of effect size estimate 
  dat$disc.ci.ul.s[i] <- res.so$ci.ub - res.sc$ci.ub  # difference subset original and complete checked CI UB of effect size estimate
  
  # Fill in discrepancy category for CI of effect sizes on subset checked MAs
  if (abs(dat$disc.ci.li.s[i]) < 0.020 & abs(dat$disc.ci.ul.s[i]) < 0.020) {
    dat$disccat.ci.s[i] = 0
  } else if (abs(dat$disc.ci.li.s[i]) > 0.093 | abs(dat$disc.ci.ul.s[i]) > 0.093) {
    dat$disccat.ci.s[i] = 3
  } else if (abs(dat$disc.ci.li.s[i]) > 0.056 & abs(dat$disc.ci.li.s[i]) <= 0.093 | abs(dat$disc.ci.ul.s[i]) > 0.056 & abs(dat$disc.ci.ul.s[i]) <= 0.093) {
    dat$disccat.ci.s[i] = 2
  } else if (abs(dat$disc.ci.li.s[i]) >= 0.020 & abs(dat$disc.ci.li.s[i]) <= 0.056 | abs(dat$disc.ci.ul.s[i]) >= 0.020 & abs(dat$disc.ci.ul.s[i]) <= 0.056) {
    dat$disccat.ci.s[i] = 1
  } else {
    dat$disccat.ci.s[i] == "check"
  }
  
  dat$disc.tau.s[i] <- res.so$tau2 - res.sc$tau2 # difference subset original and subset checked tausquared
  
  # Fill in discrepancy category for tau2 on subset checked MAs
  if (abs(dat$disc.tau.s[i]) < 0.020) {
    dat$disccat.tau.s[i] = 0
  } else if (abs(dat$disc.tau.s[i]) >= 0.020 & abs(dat$disc.tau.s[i]) <= 0.056) {
    dat$disccat.tau.s[i] = 1
  } else if (abs(dat$disc.tau.s[i]) > 0.056 & abs(dat$disc.tau.s[i]) <= 0.093) {
    dat$disccat.tau.s[i] = 2
  } else if (abs(dat$disc.tau.s[i]) > 0.093) {
    dat$disccat.tau.s[i] = 3
  } else {
    dat$disccat.tau.s[i] = "check"
  }
  
  
  l1o <- leave1out(res.sc)              # leave-one-out analysis
  estdiff <- l1o$estimate - c(res.sc$b)  # MA pooled effect size estimates if ith study was omitted
  dat$outeff1.s[i] <- length(estdiff[abs(estdiff) >= 0.020 & abs(estdiff) <= 0.056])  # no of studies that have a small effect on MA effect size
  dat$outeff2.s[i] <- length(estdiff[abs(estdiff) > 0.056 & abs(estdiff) <= 0.093])  # no of studies that have a medium effect on MA effect size
  dat$outeff3.s[i] <- length(estdiff[abs(estdiff) > 0.093])                          # no of studies that have a large effect on MA effect size
  
  ci.lb.diff <- l1o$ci.lb - c(res.sc$ci.lb) # CI LB of effect size estimates if ith study was omitted   
  ci.ub.diff <- l1o$ci.ub - c(res.sc$ci.ub) # CI UB of effect size estimates if ith study was omitted 
  dat$outci1.s[i] <- max(length(ci.lb.diff[abs(ci.lb.diff) >= 0.020 & abs(ci.lb.diff) <= 0.056]),length(ci.ub.diff[abs(ci.ub.diff) >= 0.020 & abs(ci.ub.diff) <= 0.056])) # no of studies that have a small effect on MA effect size CI
  dat$outci2.s[i] <- max(length(ci.lb.diff[abs(ci.lb.diff) > 0.056 & abs(ci.lb.diff) <= 0.093]),length(ci.ub.diff[abs(ci.ub.diff) > 0.056 & abs(ci.ub.diff) <= 0.093])) # no of studies that have a medium effect on MA effect size CI
  dat$outci3.s[i] <- max(length(ci.lb.diff[abs(ci.lb.diff) > 0.093]),length(ci.ub.diff[abs(ci.ub.diff) > 0.093]))                                             # no of studies that have a large effect on MA effect size CI
  
  taudiff <- l1o$tau2 - c(res.sc$tau2) # tau2 estimate if ith study was omitted
  dat$outtau1.s[i] <- length(taudiff[abs(taudiff)>= 0.020 & abs(taudiff)<= 0.056])  # no of studies that have a small effect on tau2 estimate
  dat$outtau2.s[i] <- length(taudiff[abs(taudiff)> 0.056 & abs(taudiff)<= 0.093])  # no of studies that have a medium effect on tau2 estimate
  dat$outtau3.s[i] <- length(taudiff[abs(taudiff)> 0.093])                         # no of studies that have a large effect on tau2 estimate
  
}

# Save all MA results
setwd("C:/Users/s421506/tiu/research/effectsizes/")
write.table(dat, file = "codebook-meta-analyses-final2.csv", row.names=F, col.names=T, sep=' ')


