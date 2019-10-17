#rm(list = ls()) # clear workspace
options(scipen=999) # no scientific notation
packages <- c("psych","irr","readxl") # library packages
# sapply(packages,install.packages(packages),character.only=T)     # if packages are not yet installed
sapply(packages,library,character.only=T)

datp <- read.table("../codebooks/codebook-primary-studies-complete.csv", header=T, sep = '')
datm <- read.table("../codebooks/codebook-meta-analyses-complete.csv", header=T, sep = '')
df.studies <- read.table("../codebooks/studies-all.csv", header=T, sep = ';')
df.out <- read.table("../codebooks/studies-outliers.csv", header=T, sep = ';')
df.reg <- read.table("../codebooks/studies-nonoutliers.csv", header=T, sep = ';')
datret <- read.table("../codebooks/nonretrieved-primary-studies.csv", header=T, sep = ';')


# PREREGISTRATION RESULTS -------------------------------------------------

# Discrepancies primary studies -------------------------------------------
# What is the proportion of primary studies for which the reproduced effect size differs from the effect size estimate of the authors of the meta-analysis 
# with a small, moderate, or large discrepancy? [It is expected that several primary and meta-analytic effect sizes will be incorrectly reproduced but there is no specific 
# prediction as to how large this proportion will be].
attach(datp)
table(info)   # based on original effect size estimate
# category 0 = 276
# category 1 = 74
# category 2 = 54
# category 3 = 96

table(disccat.eff)
# none: 386
# small: 62
# moderate: 21
# large: 31

# How many studies were larger than reported?
sum(datp$effestnew < datp$effest) # 165
sum(datp$effestnew > datp$effest) # 162
sum(datp$effest == datp$effestnew) # 173

# Crosstabs effect type / discrepancy -------------------------------------
# For which types of effect sizes are errors found most frequently?  [It is expected that more data extraction errors will be found in primary studies that report SMDs compared to correlations].
attach(datp)
typeof.info <- table(efftype,info)
sum(typeof.info["d",2:4]);sum(typeof.info["d",1:4]);sum(typeof.info["d",2:4])/sum(typeof.info["d",1:4])   # effect size d = 44 either different, omitted, or unclear out of 62 = 70.96774%
sum(typeof.info["g",2:4]);sum(typeof.info["g",1:4]);sum(typeof.info["g",2:4])/sum(typeof.info["g",1:4])   # effect size g = 96 either different, omitted, or unclear out of 185 = 51.89189%
sum(typeof.info["r",2:4]);sum(typeof.info["r",1:4]);sum(typeof.info["r",2:4])/sum(typeof.info["r",1:4])   # effect size r = 84 either different, omitted, or unclear out of 245 = 3428571%
sum(typeof.info["z",2:4]);sum(typeof.info["z",1:4]);sum(typeof.info["z",2:4])/sum(typeof.info["z",1:4])   # effect size z = 0 either different, omitted, or unclear out of 8 = 0%

typeof <- table(efftype,disccat.eff) # original effect sizes
sum(typeof["d",2:4]);sum(typeof["d",1:4]);sum(typeof["d",2:4])/sum(typeof["d",1:4])   # effect size d = 12 discrepancies out of 62 = 19.35484%
sum(typeof["g",2:4]);sum(typeof["g",1:4]);sum(typeof["g",2:4])/sum(typeof["g",1:4])   # effect size g = 49 discrepancies out of 185 = 26.48649%
sum(typeof["r",2:4]);sum(typeof["r",1:4]);sum(typeof["r",2:4])/sum(typeof["r",1:4])   # effect size r = 53 discrepancies out of 245 = 21.63265%
sum(typeof["z",2:4]);sum(typeof["z",1:4]);sum(typeof["z",2:4])/sum(typeof["z",1:4])   # effect size z = 0 discrepancies out of 8 = 0%

# Overall, more data extraction errors are expected in unpublished studies compared to published ones, due to more varying reporting standards in the former. [To determine in which types of articles more data extraction errors occur, 
# the frequency of errors will be compared for published and unpublished studies.]
attach(datp)

typecat.info <- table(typecat,info)
sum(typecat.info[1,2:4]);sum(typecat.info[1,1:4]);sum(table(typecat,info)[1,2:4])/sum(table(typecat,info)[1,])        # 216 out of 455 published papers either different, omitted, or unclear = 47.47253%
sum(typecat.info[2,2:4]);sum(typecat.info[2,1:4]);sum(table(typecat,info)[2,2:4])/sum(table(typecat,info)[2,])        # 8 out of 45 non-published papers either different, omitted, or unclear = 17.7777%

table(typecat,disccat.eff) # original effect sizes
sum(table(typecat,disccat.eff)[1,2:4]);sum(table(typecat,disccat.eff)[1,]);sum(table(typecat,disccat.eff)[1,2:4])/sum(table(typecat,disccat.eff)[1,])  # 109 out of 455 published papers contained a discrepancy = 23.95604%
sum(table(typecat,disccat.eff)[2,2:4]);sum(table(typecat,disccat.eff)[2,]);sum(table(typecat,disccat.eff)[2,2:4])/sum(table(typecat,disccat.eff)[2,])  # 5 out of 45 non-published papers contained a discrepancy = 11.11%


# NON PREREGISTRATION / EXPLORATORY RESULTS -------------------------------

# Outliers ----------------------------------------------------------------
# How many of the outlying studies contained discrepancies
attach(datp)

table(typestudy,info)
sum(table(typestudy,info)[1,2:4]);sum(table(typestudy,info)[1,]);sum(table(typestudy,info)[1,2:4])/sum(table(typestudy,info)[1,])        # 77 out of 179 outliers either different, omitted, or unclear = 0.3908629%
sum(table(typestudy,info)[2,2:4]);sum(table(typestudy,info)[2,]);sum(table(typestudy,info)[2,2:4])/sum(table(typestudy,info)[2,])        # 147 out of 303 non-outliers either different, omitted, or unclear = 0.4851485%

table(typestudy,disccat.eff)
sum(table(typestudy,disccat.eff)[1,2:4]);sum(table(typestudy,disccat.eff)[1,]);sum(table(typestudy,disccat.eff)[1,2:4])/sum(table(typestudy,disccat.eff)[1,])  # 44 out of 197 outliers contained a discrepancy = 22.33503%
sum(table(typestudy,disccat.eff)[2,2:4]);sum(table(typestudy,disccat.eff)[2,]);sum(table(typestudy,disccat.eff)[2,2:4])/sum(table(typestudy,disccat.eff)[2,])  # 70 out of 303 non-outliers contained a discrepancy = 23.10231%

# Discrepancies in sample size --------------------------------------------
attach(datp)
length(which(!abs(datp$disc.n) == 0))    # 67

# INTERRATER RELIABILITY --------------------------------------------------
dfe <- read_excel("../codebooks/individual/codebook-primary-studies-em.xlsx", 2)
dfa <- read_excel("../codebooks/individual/codebook-primary-studies-aoc.xlsx", 2)

# Whereas unweighted kappa does not distinguish among degrees of disagreement, 
# Weighted kappa incorporates the magnitude of each disagreement
# Since our categories are not ordered, we need unweighted
infoe <- dfe$info
infoa <- dfa$info
dfinfo <- cbind(infoe,infoa)

cohen.kappa(dfinfo) # unweighted = 0.71

# now by hand to check:
tab <- table(infoe,infoa);tab
# simple agreement among raters
pr.a <- (tab[1,1] + tab[2,2] + tab[3,3] + tab[4,4])  / sum(tab);pr.a
# expected agreement among raters by chance
cat0 <- (sum(tab[1,]) / sum(tab)) * (sum(tab[,1]) / sum(tab)) * sum(tab);cat0
cat1 <- (sum(tab[2,]) / sum(tab)) * (sum(tab[,2]) / sum(tab)) * sum(tab);cat1
cat2 <- (sum(tab[3,]) / sum(tab)) * (sum(tab[,3]) / sum(tab)) * sum(tab);cat2
cat3 <- (sum(tab[4,]) / sum(tab)) * (sum(tab[,4]) / sum(tab)) * sum(tab);cat3
pr.e <- (cat0 + cat1 + cat2 + cat3) / sum(tab);pr.e

kappa <- (pr.a - pr.e) / (1 - pr.e);kappa

# SAMPLE WEIGHTING --------------------------------------------------------
# appendix e table 1
nrow(df.studies)
nrow(df.out)
nrow(df.reg)
nrow(df.out)/nrow(df.studies)
nrow(df.reg)/nrow(df.studies)

# table 2
sum(datp$typestudy == "r")
sum(datp$typestudy == "o")
sum(datp$typestudy == "r")/500
sum(datp$typestudy == "o")/500

# table 3
sum(df.reg$id == 1)
sum(df.out$id == 1)
sum(df.reg$id == 1)/(sum(df.reg$id == 1)+sum(df.out$id == 1))
sum(df.out$id == 1)/(sum(df.reg$id == 1)+sum(df.out$id == 1))

# Errors
for (i in 1:nrow(datp)) {
  if (datp$info[i] != 0) { 
    datp$error[i] <- 1 
  } else { 
    datp$error[i] <- 0 }
}

# Make list with matrices for all MA seperately
datp$metano <- as.numeric(datp$meta)
reg.ma.err <- reg.ma.nerr <- out.ma.err <- out.ma.nerr <- c()
mats <- list()

for (i in 1:33) {
  
  reg.ma.err <- c(reg.ma.err,sum(datp$error[datp$metano %in% i & datp$typestudy %in% "r"]))
  reg.ma.nerr <- c(reg.ma.nerr,sum(datp$metano %in% i & datp$typestudy %in% "r")-sum(datp$error[datp$metano %in% i & datp$typestudy %in% "r"]))
  out.ma.err <-  c(out.ma.err,sum(datp$error[datp$metano %in% i & datp$typestudy %in% "o"]))
  out.ma.nerr <- c(out.ma.nerr,sum(datp$metano %in% i & datp$typestudy %in% "o")-sum(datp$error[datp$metano %in% i & datp$typestudy %in% "o"]))
  
  mat <- matrix(c(reg.ma.nerr[i],out.ma.nerr[i],reg.ma.err[i],out.ma.err[i]), ncol=2, nrow=2)
  rownames(mat) <- c("reg","out")
  colnames(mat) <- c("no error","error")
  
  mats[[i]] <- mat
  
}

# appendix e table 4
mats[[1]]

# First, we calculate correction weights that correct proportions within each MA sample.
# Within each sample, the proportion of non-outliers and outliers will come to represent the proportion in the 
# individual MA population

p1 <- p2 <- d1 <- d2 <- g1 <- g2 <- w1 <- w2 <- c()

for (i in 1:33) {
  N2 <- out <- sum(as.numeric(df.out$type[df.out$id %in% i]))
  N1 <- reg <- sum(as.numeric(df.reg$type[df.reg$id %in% i]))
  N <- tot <- sum(df.studies$id %in% i)
  n2 <- sel.out <-  sum(c(out.ma.err[i],out.ma.nerr[i]))
  n1 <- sel.reg <- sum(c(reg.ma.err[i],reg.ma.nerr[i]))
  n <- sel.tot <- sum(c(reg.ma.err[i],reg.ma.nerr[i],out.ma.err[i],out.ma.nerr[i]))
  # correction weight for all regs / outliers
  g1[i] <- (N1 / N) / (n1 / n)
  g2[i] <- (N2 / N) / (n2 / n)

    if (g2[i] == "NaN" | g2[i] == "Inf") { g2[i] <- 0 }
}

mats[[1]]

# by multiplying all 33 MA matrices with g1 (for the reg row) and g2 (for the outlier row), the proportions within
# the sample are restored to the population samples.

for (i in 1:33) {
  
  mats[[i]][1,] <- mats[[i]][1,] * g1[i]
  mats[[i]][2,] <- mats[[i]][2,] * g2[i]
  
}

# now MA1, which had a 50/50 sample has a 65/35 sample:
sum(mats[[1]][1,])/sum(mats[[1]])

# appendix e table 5
mats[[1]]

# Probability of finding an error in meta-analysis
p.x <- p.x.a <- p.x.b <- c()
for (i in 1:33) {
  
  p.x[i] <- sum(mats[[i]][,2]) / sum(mats[[i]])
  p.x.a[i] <- mats[[i]][3] / sum(mats[[i]][1,])
  p.x.b[i] <- mats[[i]][4] / sum(mats[[i]][2,])
  if (p.x.b[i] == "NaN") { p.x.b[i] <- 0 
  }
}

# weight each probability with the number of studies in that meta-analysis and divide by the total no. of studies
ss <- as.numeric(table(df.studies$id))
p.x.tot <- sum((p.x * ss)) / nrow(df.studies); p.x.tot
p.x.a.tot <- sum((p.x.a * ss)) / nrow(df.studies)
p.x.b.tot <- sum((p.x.b * ss)) / nrow(df.studies)

