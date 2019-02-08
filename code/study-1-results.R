rm(list = ls())
setwd("C:/Users/s421506/tiu/research/effectsizes/codebooks/")
packages <- c("readxl")
options(scipen=999)
#lapply(packages,install.packages(packages),character.only=T)     # if packages are not yet installed
lapply(packages,library,character.only=T)

datp <- read.table("codebook-primary-studies-final.csv", header=T, sep = '')
datm <- read.table("codebook-meta-analyses-final.csv", header=T, sep = '')
df.studies <- read.table("studies-all.csv", header=T, sep = ';')
df.out <- read.table("studies-outliers.csv", header=T, sep = ';')
df.reg <- read.table("studies-nonoutliers.csv", header=T, sep = ';')
datret <- read.table("nonretrieved-primary-studies.csv", header=T, sep = ';')

# PREREGISTRATION RESULTS -------------------------------------------------

# Discrepancies primary studies -------------------------------------------
# What is the proportion of primary studies for which the reproduced effect size (Fisher's transformed r) differs from the effect size estimate of the authors of the meta-analysis 
# with a small, moderate, or large discrepancy? [It is expected that several primary and meta-analytic effect sizes will be incorrectly reproduced but there is no specific 
# prediction as to how large this proportion will be].
attach(datp)
table(disccat.eff)
# none: 386
# small: 62
# moderate: 20
# large: 32

#table(disccat.z)
# none: 382
# small: 45
# moderate: 17
# large: 56

table(info)   # based on original effect size estimate
# category 0 = 276
# category 1 = 74
# category 2 = 54
# category 3 = 96

# Crosstabs effect type / discrepancy -------------------------------------
# For which types of effect sizes are errors found most frequently?  [It is expected that more data extraction errors will be found in primary studies that report SMDs compared to correlations].
attach(datp)
typeof <- table(efftype,disccat.eff) # original effect sizes
sum(typeof["d",2:4]);sum(typeof["d",1:4]);sum(typeof["d",2:4])/sum(typeof["d",1:4])   # effect size d = 12 discrepancies out of 62 = 19.35484%
sum(typeof["g",2:4]);sum(typeof["g",1:4]);sum(typeof["g",2:4])/sum(typeof["g",1:4])   # effect size g = 49 discrepancies out of 185 = 26.48649%
sum(typeof["r",2:4]);sum(typeof["r",1:4]);sum(typeof["r",2:4])/sum(typeof["r",1:4])   # effect size r = 53 discrepancies out of 245 = 21.63265%
sum(typeof["z",2:4]);sum(typeof["z",1:4]);sum(typeof["z",2:4])/sum(typeof["z",1:4])   # effect size z = 0 discrepancies out of 8 = 0%

# typeof <- table(efftype,disccat.z) # transformed effect sizes
# sum(typeof["d",2:4]);sum(typeof["d",1:4]);sum(typeof["d",2:4])/sum(typeof["d",1:4])   # effect size d = 40 discrepancies out of 141 = 28.36879%
# sum(typeof["g",2:4]);sum(typeof["g",1:4]);sum(typeof["g",2:4])/sum(typeof["g",1:4])   # effect size g = 32 discrepancies out of 106 = 30.18868%
# sum(typeof["r",2:4]);sum(typeof["r",1:4]);sum(typeof["r",2:4])/sum(typeof["r",1:4])   # effect size r = 46 discrepancies out of 245 = 18.77551%
# sum(typeof["z",2:4]);sum(typeof["z",1:4]);sum(typeof["z",2:4])/sum(typeof["z",1:4])   # effect size z = 0 discrepancies out of 8 = 0%

typeof.info <- table(efftype,info)
sum(typeof.info["d",2:4]);sum(typeof.info["d",1:4]);sum(typeof.info["d",2:4])/sum(typeof.info["d",1:4])   # effect size d = 44 either different, omitted, or unclear out of 62 = 70.96774%
sum(typeof.info["g",2:4]);sum(typeof.info["g",1:4]);sum(typeof.info["g",2:4])/sum(typeof.info["g",1:4])   # effect size g = 96 either different, omitted, or unclear out of 185 = 51.89189%
sum(typeof.info["r",2:4]);sum(typeof.info["r",1:4]);sum(typeof.info["r",2:4])/sum(typeof.info["r",1:4])   # effect size r = 84 either different, omitted, or unclear out of 245 = 3428571%
sum(typeof.info["z",2:4]);sum(typeof.info["z",1:4]);sum(typeof.info["z",2:4])/sum(typeof.info["z",1:4])   # effect size z = 0 either different, omitted, or unclear out of 8 = 0%

# Overall, more data extraction errors are expected in unpublished studies compared to published ones, due to more varying reporting standards in the former. [To determine in which types of articles more data extraction errors occur, 
# the frequency of errors will be compared for published and unpublished studies.]
attach(datp)
table(typecat,disccat.eff) # original effect sizes
sum(table(typecat,disccat.eff)[1,2:4]);sum(table(typecat,disccat.eff)[1,]);sum(table(typecat,disccat.eff)[1,2:4])/sum(table(typecat,disccat.eff)[1,])  # 109 out of 455 published papers contained a discrepancy = 23.95604%
sum(table(typecat,disccat.eff)[2,2:4]);sum(table(typecat,disccat.eff)[2,]);sum(table(typecat,disccat.eff)[2,2:4])/sum(table(typecat,disccat.eff)[2,])  # 5 out of 45 non-published papers contained a discrepancy = 11.11%

#table(typecat,disccat.z) # transformed effect sizes
#sum(table(typecat,disccat.z)[1,2:4]);sum(table(typecat,disccat.z)[1,]);sum(table(typecat,disccat.z)[1,2:4])/sum(table(typecat,disccat.z)[1,])  # 116 out of 455 published papers contained a discrepancy = 25.495%
#sum(table(typecat,disccat.z)[2,2:4]);sum(table(typecat,disccat.z)[2,]);sum(table(typecat,disccat.z)[2,2:4])/sum(table(typecat,disccat.z)[2,])  # 2 out of 45 non-published papers contained a discrepancy = 4.4444%

typecat.info <- table(typecat,info)
sum(typecat.info[1,2:4]);sum(typecat.info[1,1:4]);sum(table(typecat,info)[1,2:4])/sum(table(typecat,info)[1,])        # 216 out of 455 published papers either different, omitted, or unclear = 47.47253%
sum(typecat.info[2,2:4]);sum(typecat.info[2,1:4]);sum(table(typecat,info)[2,2:4])/sum(table(typecat,info)[2,])        # 8 out of 45 non-published papers either different, omitted, or unclear = 17.7777%

# Non-retrieval -----------------------------------------------------------
# It will be noted how many studies cannot be retrieved, and what kind of studies these are (e.g., book chapters, dissertations).
attach(datret)
table(type)
length(unique(ma));length(unique(ma))/33
# In total 154 primary papers not found in 26 out of 33 meta-analyses = 79%.

sum(type == "article");sum(type == "article")/nrow(datret)*100                                                  # of which 70 (45%) are articles
sum((type == "dissertation")+(type == "unpublished doctoral dissertation"))
sum((type == "dissertation")+(type == "unpublished doctoral dissertation"))/nrow(datret)*100                    # of which 48 (31%) are dissertations
sum(type == "paper presented at conference");sum(type == "paper presented at conference")/nrow(datret)*100      # of which 8 (5%) are papers presented at conferences
sum((type == "book")+(type == "book chapter"));sum((type == "book")+(type == "book chapter"))/nrow(datret)*100  # of which 7 (5%) are books or book chapters
sum((type == "manuscript submitted for publication")+(type == "unpublished manuscript"))
sum((type == "manuscript submitted for publication")+(type == "unpublished manuscript"))/nrow(datret)*100       # of which 7 (5%) are unpublished manuscripts
sum(type == "poster");sum(type == "poster")/nrow(datret)*100                                                    # of which 7 (5%) are posters
sum(type == "unpublished master thesis");sum(type == "unpublished master thesis")/nrow(datret)*100              # of which 3 (2%) are unpublished master theses
sum(type == "unpublished raw data");sum(type == "unpublished raw data")/nrow(datret)*100                        # of which 3 (2%) are unpublished raw data
sum(type == "report");sum(type == "report")/nrow(datret)*100                                                    # of which 1 (1%) is a report

# NON PREREGISTRATION / EXPLORATORY RESULTS -------------------------------

# Outliers ----------------------------------------------------------------
# How many of the outlying studies contained discrepancies
attach(datp)

table(typestudy,disccat.eff)
sum(table(typestudy,disccat.eff)[1,2:4]);sum(table(typestudy,disccat.eff)[1,]);sum(table(typestudy,disccat.eff)[1,2:4])/sum(table(typestudy,disccat.eff)[1,])  # 36 out of 179 outliers contained a discrepancy = 20.1111%
sum(table(typestudy,disccat.eff)[2,2:4]);sum(table(typestudy,disccat.eff)[2,]);sum(table(typestudy,disccat.eff)[2,2:4])/sum(table(typestudy,disccat.eff)[2,])  # 78 out of 321 non-outliers contained a discrepancy = 24.29907%

#table(typestudy,disccat.z)
#sum(table(typestudy,disccat.z)[1,2:4]);sum(table(typestudy,disccat.z)[1,]);sum(table(typestudy,disccat.z)[1,2:4])/sum(table(typestudy,disccat.z)[1,])  # 29 out of 178 outliers contained a discrepancy = 16.292%
#sum(table(typestudy,disccat.z)[2,2:4]);sum(table(typestudy,disccat.z)[2,]);sum(table(typestudy,disccat.z)[2,2:4])/sum(table(typestudy,disccat.z)[2,])  # 89 out of 322 non-outliers contained a discrepancy = 27.63975%

table(typestudy,info)
sum(table(typestudy,info)[1,2:4]);sum(table(typestudy,info)[1,]);sum(table(typestudy,info)[1,2:4])/sum(table(typestudy,info)[1,])        # 70 out of 179 outliers either different, omitted, or unclear = 0.3910615%
sum(table(typestudy,info)[2,2:4]);sum(table(typestudy,info)[2,]);sum(table(typestudy,info)[2,2:4])/sum(table(typestudy,info)[2,])        # 154 out of 321 non-outliers either different, omitted, or unclear = 0.4797508%


# Discrepancies in sample size --------------------------------------------
attach(datp)
length(which(!abs(datp$disc.n) == 0))    # 67


# Sample weighting --------------------------------------------------------
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


# by multiplying all 33 MA matrices with g1 (for the reg row) and g2 (for the outlier row), the proportions within
# the sample are restored to the population samples.

for (i in 1:33) {
  
  mats[[i]][1,] <- mats[[i]][1,] * g1[i]
  mats[[i]][2,] <- mats[[i]][2,] * g2[i]
  
}

# now MA1, which had a 50/50 sample has a 65/35 sample:
sum(mats[[1]][1,])/sum(mats[[1]])
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
p.x.tot <- sum((p.x * ss)) / nrow(df.studies)
p.x.a.tot <- sum((p.x.a * ss)) / nrow(df.studies)
p.x.b.tot <- sum((p.x.b * ss)) / nrow(df.studies)















a <-mats[[1]]+mats[[2]]+mats[[3]]+mats[[4]]+mats[[5]]+mats[[6]]+mats[[7]]+mats[[8]]+mats[[9]]+
  mats[[10]]+mats[[11]]+mats[[12]]+mats[[13]]+mats[[14]]+mats[[15]]+mats[[16]]+mats[[17]]+
  mats[[18]]+mats[[19]]+mats[[20]]+mats[[21]]+mats[[22]]+mats[[23]]+mats[[24]]+mats[[25]]+
  mats[[26]]+mats[[27]]+mats[[28]]+mats[[29]]+mats[[30]]+mats[[31]]+mats[[32]]+mats[[33]]

(a[3]+a[4])/sum(a)
(a[3]/(a[1]+a[3]))
(a[4]/(a[2]+a[4]))

#

















































































# So now, within each meta-analysis, the proportions are adjusted to the population.
# Now, we want to generalize the sample to the population
# one outlier in the population (for MA1 this is 20 total) is represented by 2 in the sample
1 / (10/20) # outlier
1 / (10/37) # reg

# we again multiply all the matrices with these estimates
for (i in 1:33) {
  
  mats[[i]][1,] <- mats[[i]][1,] * d1[i]
  mats[[i]][2,] <- mats[[i]][2,] * d2[i]
  
}

# probability of finding an error in these 33 meta-analyses:
prob <- c()

for (i in 1:33) {
  
  prob <- c(prob,mats[[i]][1,2] / sum(mats[[i]]))
 
}






datp$metano <- as.numeric(datp$meta)

sum(c(reg.ma.err,reg.ma.nerr,out.ma.err,out.ma.nerr))

p1 <- p2 <- d1 <- d2 <- g1 <- g2 <- w1 <- w2 <- c()

for (i in 1:33) {
  
  out <- sum(as.numeric(df.out$type[df.out$id %in% i]))
  reg <- sum(as.numeric(df.reg$type[df.reg$id %in% i]))
  tot <- sum(df.studies$id %in% i)
  sel.out <- sum(c(out.ma.err[i],out.ma.nerr[i]))
  sel.reg <- sum(c(reg.ma.err[i],reg.ma.nerr[i]))
  sel.tot <- sum(c(reg.ma.err[i],reg.ma.nerr[i],out.ma.err[i],out.ma.nerr[i]))
  err.reg <- reg.ma.err[i]
  err.out <- out.ma.err[i]
  inc.prob.out <- out / tot
  inc.prob.reg <- reg / tot
  d1[i] <- 1 / (sel.reg/reg)
  d2[i] <- 1 / (sel.out/out)
  if (d2[i] == "NaN" | d2[i] == "Inf") { d2[i] <- 0 }
  N <- tot
  N1 <- reg
  N2 <- out
  n <- sel.tot
  n1 <- sel.reg
  n2 <- sel.out
  # correction weight for all regs / outliers
  g1[i] <- (N1 / N) / (n1 / n)
  g2[i] <- (N2 / N) / (n2 / n)
  if (g2[i] == "NaN" | g2[i] == "Inf") { g2[i] <- 0 }
  # weights
  w1[i] <- g1[i] * d1[i]
  w2[i] <- g2[i] * d2[i]
  # population weights
  p1[i] <- sel.reg/reg # number of regular studies in sample vs regular studies in population of MA
  p2[i] <- sel.out/tot # number of outliers in sample vs outliers in population of MA
}


# weights
for (i in 1:nrow(datp)) {
  for (q in 1:33) {
  
  if (datp$metano[i] %in% q & datp$typestudy[i] %in% "r") { 
    datp$d.weight[i] <- d1[q]                                      # design weights
    datp$g.weight[i] <- as.numeric(g1[q]) # first one in w1 list is a function, so +1
    datp$w.weight[i] <- as.numeric(w1[q])
  } 
  if (datp$metano[i] %in% q & datp$typestudy[i] %in% "o") { 
    datp$d.weight[i] <- d2[q]
    datp$g.weight[i] <- as.numeric(g2[q]) # first one in w1 list is a function, so +1
    datp$w.weight[i] <- as.numeric(w2[q])
  } 
  }
}




N <- nrow(df.studies)
N1 <- nrow(df.reg)
N2 <- nrow(df.out)
n <- nrow(datp)
n1 <- sum(datp$typestudy == "r")
n2 <- sum(datp$typestudy == "o")
# proportion regulars population P(A) / sample P(a)
N1/N;n1/n 
# proportion outliers population P(B) / sample P(b)
N2/N;n2/n
# probabilitiy of being included in the sample (inclusion probability):
p1 <- n1/N1
p2 <- n2/N2
# inclusion weights
d1 <- N1/n1
d2 <- N2/n2
d <- N/n
# correction weight for all regs / outliers
g1 <-(N1 / N) / (n1 / n)
g2 <- (N2 / N) / (n2 / n) 
# weights
w1 <- g1 * d
w2 <- g2 * d

# make summary matrix
tab <- table(typestudy,info)
out.noerror <- tab[1,1]
out.error <- sum(tab[1,2:4])
reg.noerror <- tab[2,1]
reg.error <- sum(tab[2,2:4])
mat <- matrix(c(reg.noerror,reg.error,out.noerror,out.error), ncol=2, nrow=2)
rownames(mat) <- c("reg","out")
colnames(mat) <- c("no error","error")
mat
p.x <- sum(mat[,2])/sum(mat);p.x
p.x.a <- mat[1,2] / sum(mat[1,]);p.x.a
p.x.b <- mat[2,2] / sum(mat[2,]);p.x.b
p.a.x <- mat[1,2] / sum(mat[,2]);p.a.x
p.b.x <- mat[2,2] / sum(mat[,2]);p.b.x
# multiply all outlier values by g2, regulars by g1
mat[1,] <- mat[1,] * w1
mat[2,] <- mat[2,] * w2
mat
p.x <- sum(mat[,2])/sum(mat);p.x
p.x.a <- mat[1,2] / sum(mat[1,]);p.x.a
p.x.b <- mat[2,2] / sum(mat[2,]);p.x.b
p.a.x <- mat[1,2] / sum(mat[,2]);p.a.x
p.b.x <- mat[2,2] / sum(mat[,2]);p.b.x
# probability error in both P(x|ab) in the sample / with correction
p.x.ab <- p.x.a * p.x.b / p.x            # 25.04%
# count
out.noerror <- mat[2,1]
out.error <- mat[2,2]
reg.noerror <- mat[1,1]
reg.error <- mat[1,2]
#
p.a.x <- mat[1,2] / sum(mat[,2]);p.a.x
p.b.x <- mat[2,2] / sum(mat[,2]);p.b.x


reg.ma.err <- c(reg.ma.err,sum(datp$error[datp$metano %in% i & datp$typestudy %in% "r"]))
reg.ma.nerr <- c(reg.ma.nerr,sum(datp$metano %in% i & datp$typestudy %in% "r")-sum(datp$error[datp$metano %in% i & datp$typestudy %in% "r"]))
out.ma.err <-  c(out.ma.err,sum(datp$error[datp$metano %in% i & datp$typestudy %in% "o"]))
out.ma.nerr <- c(out.ma.nerr,sum(datp$metano %in% i & datp$typestudy %in% "o")-sum(datp$error[datp$metano %in% i & datp$typestudy %in% "o"]))



test2 <- rbind(reg.ma.err,out.ma.err)
barplot(test2,beside=T)

barplot(ma.err, names.arg=unique(datp$metano))
barplot(ma.nerr, names.arg=unique(datp$metano))
library(ggplot2)
ggplot(test2)







plot(density(ma.nerr))
plot(ma.err,unique(datp$meta))

vegLengths <- rbind(carrots, cukes)


sum(datp[i,""])

df <- data.frame(id = rep(LETTERS[1:8], 10), weight = as.integer(rnorm(80, 80, 10)))

datp
df <- data.frame(as.numeric(datp$meta),datp$error,)


library(ggplot2)
ggplot(df, aes(x=weight)) + 
  geom_density(alpha=.2, fill="#FF6666") +
  facet_wrap( ~ id, nrow=2)

plot((datp$info),as.numeric(datp$meta))
plot(density(words2))

# inclusion weights for all regs / outliers
d1 <- N1 / n1
d2 <- N2 / n2
# errors
for (i in 1:nrow(datp)) {
  if (datp$info[i] != 0) { 
    datp$error[i] <- 1 
    } else { 
      datp$error[i] <- 0 }
}

# weights
for (i in 1:nrow(datp)) {
  if (datp$typestudy[i] == "r") { 
    datp$weight[i] <- g1 
  } else { 
    datp$weight[i] <- g2 }
}

# weighted estimate
datp$est.weight <- datp$error * datp$weight

#poststratification estimator
sub.reg <- datp[datp$typestudy == "r",]
sub.out <- datp[datp$typestudy == "o",]
y <- (1 / N) * sum(c(N1 * mean(sub.reg$est.weight)),(N2 * mean(sub.out$est.weight))); y


# responding elements in all regulars
datp$error <- datp$info[datp$typestudy %in% "0"]
datp$weight <- datp$info[datp$typestudy %in% "r"]
res.out <- datp$info[datp$typestudy %in% "o"]  
# weight the estimates per cluster
res.reg <- res.reg * g1
res.out <- res.out * g2
  
  
  
  
  
# propability error in regular P(x|a)
table(typestudy,info)
p.x.a <- sum(table(typestudy,info)[2,2:4])/sum(table(typestudy,info)[2,])  # 48.5%
# propability error in outlier P(x|b)
p.x.b <- sum(table(typestudy,info)[1,2:4])/sum(table(typestudy,info)[1,])  # 39.1%
# probability error in total P(x)
p.x <- sum(table(typestudy,info)[1:2,2:4])/sum(table(typestudy,info)[1:2,]) # 44.8%
# probability error in both P(x|ab) in the sample / without correction      # 42.32%   
p.x.ab <- p.x.a * p.x.b / p.x;p.x.ab
# make summary matrix
tab <- table(typestudy,info)
out.noerror <- tab[1,1]
out.error <- sum(tab[1,2:4])
reg.noerror <- tab[2,1]
reg.error <- sum(tab[2,2:4])
mat <- matrix(c(reg.noerror,reg.error,out.noerror,out.error), ncol=2, nrow=2)
rownames(mat) <- c("reg","out")
colnames(mat) <- c("no error","error")
# multiply all outlier values by g2, regulars by g1
mat[1,] <- mat[1,] * g1 * d1
mat[2,] <- mat[2,] * g2 * d2
# propability error in regular P(x|a)    # 43.47%
p.x.a <- mat[1,2] / sum(mat[1,])
# propability error in outlier P(x|b)    # 34.375
p.x.b <- mat[2,2] / sum(mat[2,])
# probability error in total P(x)        # 59.6%
p.x <- sum(mat[,1])/sum(mat)
# probability error in both P(x|ab) in the sample / with correction
p.x.ab <- p.x.a * p.x.b / p.x            # 25.04%




# inclusion weight
d <- N/n
# weighted estimates
n1 * d * g1
n2 * d * g2
# inverse probaility sampling weight - how many regs in the population are represented by a single reg in the sample
1 / (n1/N1)
# inverse probabilitie sampling weight - how many regs in the population are represented by a single reg in the sample
1 / (n2/N2)
# correction weights
# no of outliers population vs sample




# no of errors in reg
reg <- sum(table(typestudy,info)[2,2:4]);reg
out <- sum(table(typestudy,info)[1,2:4]);out
err.reg.corr <- err.reg * d * g1
err.out.corr <- err.out * d * g3
# 218
w.err.reg = g1 * d
w.err.
err.corr <- err * w1


# design weight * adjustment weight





# propability error in regular P(x|a)
attach(datp)
table(typestudy,info)
p.x.a <- sum(table(typestudy,info)[2,2:4])/sum(table(typestudy,info)[2,])  # 48.5%
# propability error in outlier P(x|b)
p.x.b <- sum(table(typestudy,info)[1,2:4])/sum(table(typestudy,info)[1,])  # 39.1%
# probability error in total P(x)
p.x <- sum(table(typestudy,info)[1:2,2:4])/sum(table(typestudy,info)[1:2,]) # 44.8%
# probability error in both P(x|ab) in the sample / without correction
p.x.ab <- p.x.a * p.x.b / p.x

# correction
# multiply probabilities with weight
p.x.a
p.x.a.corr <- p.x.a * g1
p.x.b
p.x.b.corr <- p.x.b * g2
# probability error in both P(x|ab) in the sample / without correction
p.x.ab
p.x.ab.corr <- p.x.a.corr * p.x.b.corr / p.x; p.x.ab.corr


p.x.a.corr <- p.x.a * g1 * d
# doesnt work because probabilities exceed 1, convert to odds first
# odds
odds.reg <- p.x.a * (1 - p.x.a)
odds.out <- p.x.b * (1 - p.x.b)
# oversampled odds
org_frac <- n1/n # original fraction: 60.1% of 500
over_frac <- N1/N # oversample fraction: 70.2% of 1951 



































# expected no of outliers
pi1 <- N2/N 

# sampling probability - expected no of reg / outlier
a1 <- n1/N1 # from the total no of regs, 22% is selected
p1 <- N1/N

1/n1 * n1







# inverse probaility sampling weight - how many regs in the population are represented by a single reg in the sample
1 / (n1/N1)
# inverse probabilitie sampling weight - how many regs in the population are represented by a single reg in the sample
1 / (n2/N2)
# propability error in regular P(x|a)
attach(datp)
table(typestudy,info)
p.x.a <- sum(table(typestudy,info)[2,2:4])/sum(table(typestudy,info)[2,])  # 48.5%
# propability error in outlier P(x|b)
p.x.b <- sum(table(typestudy,info)[1,2:4])/sum(table(typestudy,info)[1,])  # 39.1%
# probability error in total P(x)
p.x <- sum(table(typestudy,info)[1:2,2:4])/sum(table(typestudy,info)[1:2,]) # 44.8%
# probability error in both P(x|ab) in the sample / without correction
p.x.ab <- p.x.a * p.x.b / p.x

# correction
# multiply probabilities with weight
p.x.a
p.x.a.corr <- p.x.a * g1
p.x.b
p.x.b.corr <- p.x.b * g2
# probability error in both P(x|ab) in the sample / without correction
p.x.ab
p.x.ab.corr <- p.x.a.corr * p.x.b.corr / p.x; p.x.ab.corr











# sample proportion error and non-error in regular
p.x # unadjusted probability
p.x.ab
r0 <- n1/n   # sample proportion reg
r1 <- n2/n
p0 <- N1/N   # population proportion reg
p1 <- N2/N #   population proportion reg

  
p.x.corrected <- (p.x * r0 * p1) / ((1 - p.x) * r1 * p0 + (p.x * r0 * p1))
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
  # inclusion weights
  d1[i] <- 1 / (n1 / N1)
  d2[i] <- 1 / (n2 / N2)
  if (d2[i] == "NaN" | d2[i] == "Inf") { d2[i] <- 0 }
  if (g2[i] == "NaN" | g2[i] == "Inf") { g2[i] <- 0 }
  # weights
  w1[i] <- g1[i] * d1[i]
  w2[i] <- g2[i] * d2[i]
  # population weights
  p1[i] <- n1/N1 # number of regular studies in sample vs regular studies in population of MA
  p2[i] <- n2/N2 # number of outliers in sample vs outliers in population of MA
}


