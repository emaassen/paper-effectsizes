rm(list = ls())
options(scipen=999,warn=0) # show only warning notification 
#options(scipen=999,warn=1) # show all warnings
options(scipen=999,warn=0) # supress warnings (caused by coerced NAs) 
setwd("C:/Users/s421506/tiu/research/effectsizes/codebooks/")
x <- c("MAc","readxl","MAd","writexl") # packages
#lapply(x,install.packages(x),character.only=T)
lapply(x,library,character.only=T)

# Open codebook
df <- read_excel("codebook-primary-studies-final-empty.xlsx",2)
attach(df)
table(input,output)


# Input: 0 / Output: all --------------------------------------------------
df$effestnew[df$input %in% "0"] <- 0


# Input: berry2.agg.r / Output: all ---------------------------------------
dftemp <- df[df$input=="berry2.agg.r",]
detach("package:MAd", unload=TRUE)
library(MAc)

berry2.est <- c()

for (k in 1:nrow(dftemp)) {
  
  r <- c()
  
  for (l in 1:14) {
    
    r <- c(r,as.numeric(eval(parse(text = paste0("dftemp$cor.",l))))[k])
    r <- r[!is.na(r)]
    
    int.cor <- as.numeric(dftemp$int.cor[k])

    n <- as.numeric(dftemp$nnew)[k]
    #n <- c(rep(n,length(r)))
    
    id <- c(rep(1,length(r)))
    
    MAdf <- data.frame(r, n, id)
    
    aggregate <- agg(id = id, r = r, n = n, data=MAdf, cor = int.cor)
    

  }
  
  berry2.est <- c(berry2.est,aggregate$r)

}

df$effestnew[df$input %in% "berry2.agg.r"] <- berry2.est



# Input: beta / Output: d -------------------------------------------------
dftemp <- df[df$input=="beta" & df$output=="d",]
r <- as.numeric(dftemp$cor.1)
d <- (2*r)/sqrt(1-r^2)
d <- d*-1    # reverse effect in correct direction
df$effestnew[df$input %in% "beta" & df$output=="d"] <- d


# Input: beta / Output: r -----------------------------------------------
df$effestnew[df$input %in% "beta" & df$output=="r"] <- df$cor.1[df$input %in% "beta" & df$output=="r"]


# Input: chisquare / Output: r --------------------------------------------
dftemp <- df[df$input=="chisquare" & df$output=="r",]

chisquare.est <- c()

for (k in 1:nrow(dftemp)) {
  
  teststat <- c()
  
  for (l in 1:5) {
    
    teststat <- c(teststat,as.numeric(eval(parse(text = paste0("dftemp$teststat.",l))))[k])
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    n <- nc+nt

    r <- sqrt((teststat/n))
    
  }
  
  chisquare.est <- c(chisquare.est,r)
  
}

df$effestnew[df$input %in% "chisquare" & df$output %in% "r"] <- chisquare.est


# Input: chisquare / Output: d --------------------------------------------
dftemp <- df[df$input=="chisquare" & df$output=="d",]

chisquare.est <- c()

for (k in 1:nrow(dftemp)) {
  
  teststat <- c()
  
  for (l in 1:5) {
    
    teststat <- c(teststat,as.numeric(eval(parse(text = paste0("dftemp$teststat.",l))))[k])
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    n <- nc+nt
    
    r <- sqrt((teststat/n))
    
    d <- (2*r)/sqrt(1-r^2)
    
  }
  
  chisquare.est <- c(chisquare.est,d)
  
}

chisquare.est[1] <- chisquare.est[1]*-1 # first effect reversed to be in right direction

df$effestnew[df$input %in% "chisquare" & df$output %in% "d"] <- chisquare.est


# Input: etasquare / Output: all ------------------------------------------
dftemp <- df[df$input=="etasquare",]
attach(dftemp)

d <- (2*as.numeric(cor.1)/sqrt(1-(as.numeric(cor.1)^2)))*-1 # effect reversed to be in right direction

df$effestnew[df$input %in% "etasquare"] <- d


# Input: F / Output: d ----------------------------------------------------
dftemp <- df[df$input=="F" & df$output=="d",]

F.est <- c()

for (k in 1:nrow(dftemp)) {
  
  teststat <- cohensd.temp <- c()

  for (l in 1:5) {
    
    teststat <- c(teststat,as.numeric(eval(parse(text = paste0("dftemp$teststat.",l))))[k])
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]

    d <- sqrt((teststat*(nc+nt))/(nc*nt))
    
    
  }
  

  # take mean 
  F.est <- c(F.est,mean(d, na.rm=T)) # save primary study's mean cohen's d for later
  
}

F.est[1] <- F.est[1]*-1  # effect reversed to be in right direction
F.est[7] <- F.est[7]*-1 
F.est[8] <- F.est[8]*-1

df$effestnew[df$input %in% "F" & df$output %in% "d"] <- F.est


# Input: F / Output: g ----------------------------------------------------
dftemp <- df[df$input=="F" & df$output=="g",]

F.est <- c()

for (k in 1:nrow(dftemp)) {
  
  teststat <- hedgesg.temp <- c()
  
  for (l in 1:5) {
    
    teststat <- c(teststat,as.numeric(eval(parse(text = paste0("dftemp$teststat.",l))))[k])
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    dfs <- nc+nt-2
    
    d <- sqrt((teststat*(nc+nt))/(nc*nt))
    J <- (1 - 3 / (4 * dfs - 1))
    g <- J * d
    
    hedgesg.temp <- c(hedgesg.temp,g)   # vector with all effect sizes within one primary study 
    
  }
  
  
  # take mean 
  F.est <- c(F.est,mean(hedgesg.temp, na.rm=T)) # save primary study's mean cohen's d for later
  
}

df$effestnew[df$input %in% "F" & df$output %in% "g"] <- F.est


# Input: F / Output: r ----------------------------------------------------
dftemp <- df[df$input=="F" & df$output=="r",]

F.est <- c()

for (k in 1:nrow(dftemp)) {
  
  teststat <- corr.temp <- c()
  
  for (l in 1:5) {
    
    teststat <- c(teststat,as.numeric(eval(parse(text = paste0("dftemp$teststat.",l))))[k])
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    
    d <- sqrt((teststat*(nc+nt))/(nc*nt))
    a <- (nc+nt)^2/(nc*nt)
    r <- d / (sqrt(d^2 + a))
    
    corr.temp <- c(corr.temp,r)   # vector with all effect sizes within one primary study 
    
  }
  
  
  # take mean 
  F.est <- c(F.est,mean(corr.temp, na.rm=T)) # save primary study's mean cohen's d for later
  
}

F.est[2] <- F.est[2]*-1  # effect reversed to be in right direction

df$effestnew[df$input %in% "F" & df$output %in% "r"] <- F.est


# Input: farber.agg.r / Output: all ---------------------------------------
#detach("package:MAd", unload=TRUE)
library(MAc)

dftemp <- df[df$input=="farber.agg.r",]

farber.est <- c()

for (k in 1:nrow(dftemp)) {
  
  r <- c()
  
  for (l in 1:14) {
    
    r <- c(r,as.numeric(eval(parse(text = paste0("dftemp$cor.",l))))[k])
    r <- r[!is.na(r)]
    
    int.cor <- as.numeric(dftemp$int.cor[k])
    
    n <- as.numeric(dftemp$nnew)[k]

    id <- c(rep(1,length(r)))
    
    MAdf <- data.frame(r, n, id)
    
    aggregate <- agg(id = id, r = r, n = n, data=MAdf, cor = int.cor)
    
    
  }
  
  farber.est <- c(farber.est,aggregate$r)
  
}

df$effestnew[df$input %in% "farber.agg.r"] <- farber.est


# Input: kolden.agg.r / Output: all ---------------------------------------
dftemp <- df[df$input == "kolden.agg.r", ]
detach("package:MAc", unload = TRUE)
library(MAd)

kolden.est.d <- kolden.est.r <- c()

for (k in 1:nrow(dftemp)) {
  
  r <- c()
  
  for (l in 1:14) {
    
    r <- c(r, as.numeric(eval(parse(text = paste0("dftemp$cor.", l))))[k])
    r <- r[!is.na(r)]
    
    d <- (2*r)/sqrt(1-r^2)

  }
  
  int.cor <- as.numeric(dftemp$int.cor[k])
  
  n <- as.numeric(dftemp$nnew)[k]
  n1 <- n/2
  n2 <- n/2
  
  id <- c(rep(1, length(d)))
  MAdf <- data.frame(d, n1, n2, id)
  
  aggregate <- agg(id = id, es = d, n.1 = n1, n.2 = n2, data = MAdf, cor = int.cor, method = "GO1")
  
  kolden.est.d <- c(kolden.est.d, aggregate$es1)
  
  a <- (n1+n2)^2/(n1*n2) 
  
  kolden.est.r <- c(kolden.est.r,aggregate$es1 / (sqrt(aggregate$es1^2 + a)))
  
}

df$effestnew[df$input %in% "kolden.agg.r"] <- kolden.est.r
    


# Input: Mann-Whitney U / Output: All -------------------------------------
dftemp <- df[df$input=="mann-whitney-u",]

nc <- as.numeric(dftemp$ncnew)
nt <- as.numeric(dftemp$ntnew)
U <- as.numeric(dftemp$teststat.1)
std.dev <- sqrt(((nc*nt)*(nc+nt+1))/12)

Z <- (U - ((nc*nt)/2))/(std.dev)
r <- Z / (sqrt(nc+nt))

# correct for unequal sample size per Fox
a <- sqrt(0.25/((nc*0.01)*(nt*0.01)))

r <- (a*r)/(sqrt(((a^2-1)*(r^2))+1)) 
r <- r*-1

df$effestnew[df$input %in% "mann-whitney-u"] <- r



# Input: NA / Output: all -------------------------------------------------
df$effestnew[df$input %in% "NA"] <- df$effest[df$input %in% "NA"]


# Input: OR / Output: d ---------------------------------------------------
dftemp <- df[df$input=="OR" & df$output=="d",]

d.est <- c()

for (k in 1:nrow(dftemp)) {
  
   d.temp <- a <- b <- c <- d <- c()
  
  for (l in 1:4) {
    
    a <- c(a,as.numeric(eval(parse(text = paste0("dftemp$or.a.",l))))[k])
    b <- c(b,as.numeric(eval(parse(text = paste0("dftemp$or.b.",l))))[k])
    c <- c(c,as.numeric(eval(parse(text = paste0("dftemp$or.c.",l))))[k])
    d <- c(d,as.numeric(eval(parse(text = paste0("dftemp$or.d.",l))))[k])
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    
    or <- (a*d)/(b*c)
    cohensd <- log(or) * (sqrt(3)/pi)
    
    d.temp <- c(d.temp,cohensd)   # vector with all effect sizes within one primary study 
    
  }
  
  # take mean 
  d.est <- c(d.est,mean(d.temp, na.rm=T)) # save primary study's mean cohen's d for later
  
}

d.est[1] <- d.est[1]*-1  # effect reversed to be in right direction

df$effestnew[df$input %in% "OR" & df$output %in% "d"] <- d.est


# Input: OR / Output: g ---------------------------------------------------
dftemp <- df[df$input=="OR" & df$output=="g",]

g.est <- c()

for (k in 1:nrow(dftemp)) {
  
  or.temp <- a <- b <- c <- d <- c()
  
  for (l in 1:4) {
    
    a <- c(a,as.numeric(eval(parse(text = paste0("dftemp$or.a.",l))))[k])
    b <- c(b,as.numeric(eval(parse(text = paste0("dftemp$or.b.",l))))[k])
    c <- c(c,as.numeric(eval(parse(text = paste0("dftemp$or.c.",l))))[k])
    d <- c(d,as.numeric(eval(parse(text = paste0("dftemp$or.d.",l))))[k])
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    dfs <- nc+nt-2
    
    or <- (a*d)/(b*c)

    
  }
  
  or.mean <- mean(or ,na.rm=T)
  
  cohensd <- log(or.mean) * (sqrt(3)/pi)
  
  J <- (1 - 3 / (4 * dfs - 1))
  g <- J * cohensd

  # take mean 
  g.est <- c(g.est,g) # save primary study's mean hedges g for later

}

g.est[9] <- g.est[9]*-1 # reverse effect to be in right direction
df$effestnew[df$input %in% "OR" & df$output %in% "g"] <- g.est



# Input: OR / Output: r ---------------------------------------------------
dftemp <- df[df$input=="OR" & df$output=="r",]

r.est <- c()

l <- 1

a <- as.numeric(eval(parse(text = paste0("dftemp$or.a.",l))))
b <- as.numeric(eval(parse(text = paste0("dftemp$or.b.",l))))
c <- as.numeric(eval(parse(text = paste0("dftemp$or.c.",l))))
d <- as.numeric(eval(parse(text = paste0("dftemp$or.d.",l))))

nc <- as.numeric(dftemp$ncnew)
nt <- as.numeric(dftemp$ntnew)

or <- (a*d)/(b*c)
cohensd <- log(or) * (sqrt(3)/pi)
a <- (nc+nt)^2/(nc*nt)
r <- cohensd / (sqrt(cohensd^2 + a))

df$effestnew[df$input %in% "OR" & df$output %in% "r"] <- r



# Input: r / Output: d ----------------------------------------------------
dftemp <- df[df$input=="r" & df$output=="d",]

d.est <- c()

for (k in 1:nrow(dftemp)) {
  
  r <- c()
  
  for (l in 1:14) {
    
    n <- as.numeric(dftemp$nnew)[k]

    r <- c(r,as.numeric(eval(parse(text = paste0("dftemp$cor.",l))))[k])
    r <- r[!is.na(r)]
    
    r <- mean(r)
    
  }
  
  d <- (2*r)/sqrt(1-r^2)
  
  # take mean 
  d.est <- c(d.est,d) # save primary study's mean hedges g for later
  
}

df$effestnew[df$input %in% "r" & df$output %in% "d"] <- d.est



# Input: r / Output: r ----------------------------------------------------
dftemp <- df[df$input=="r" & df$output=="r",]

r.est <- c()

for (k in 1:nrow(dftemp)) {
  
  r <- c()
  
  for (l in 1:14) {
    
    r <- c(r,as.numeric(eval(parse(text = paste0("dftemp$cor.",l))))[k])
  }
  
  # take mean 
  r.est <- c(r.est,mean(r, na.rm=T)) # save primary study's mean correlation for later
  
}

df$effestnew[df$input %in% "r" & df$output %in% "r"] <- r.est



# Input: r / Output: z ----------------------------------------------------
dftemp <- df[df$input=="r" & df$output=="z",]

z.est <- c()

for (k in 1:nrow(dftemp)) {
  
  r <- z.temp <- c()
  
  for (l in 1:14) {
    
    r <- c(r,as.numeric(eval(parse(text = paste0("dftemp$cor.",l))))[k])
    r <- r[!is.na(r)]
    z <- 0.5 * log((1 + r)/(1-r))
  
    z.temp <- c(z.temp,z)   # vector with all effect sizes within one primary study 
    
  }
  
  # take mean 
  z.est <- c(z.est,mean(z.temp, na.rm=T)) # save primary study's mean correlation for later
  
}

df$effestnew[df$input %in% "r" & df$output %in% "z"] <- z.est



# Input: smd / Output: d --------------------------------------------------
dftemp <- df[df$input=="smd" & df$output=="d",]

cohensd <- c()

for (k in 1:nrow(dftemp)) {
  
  cohensd.temp <- c()
  
  for (l in 1:14) {
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    n <- nc+nt
    mc <-  as.numeric(eval(parse(text = paste0("dftemp$mc.",l))))[k]
    mt <-  as.numeric(eval(parse(text = paste0("dftemp$mt.",l))))[k]
    sdc <- as.numeric(eval(parse(text = paste0("dftemp$sdc.",l))))[k]
    sdt <- as.numeric(eval(parse(text = paste0("dftemp$sdt.",l))))[k]
    
    d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))

    cohensd.temp <- c(cohensd.temp,d)   # vector with all effect sizes within one primary study
    
  }
  
  # take mean 
  cohensd[k] <- mean(cohensd.temp, na.rm=T) # mean of all effects from one primary study
  
}

cohensd[16] <- cohensd[16]*-1 # effect size reversed to be in right direction
cohensd[17] <- cohensd[17]*-1
#cohensd[18] <- cohensd[18]*-1
cohensd[20] <- cohensd[20]*-1
#cohensd[21] <- cohensd[21]*-1
cohensd[36] <- cohensd[36]*-1
cohensd[37] <- cohensd[37]*-1
#cohensd[38] <- cohensd[38]*-1
cohensd[39] <- cohensd[39]*-1
#cohensd[40] <- cohensd[40]*-1
cohensd[41] <- cohensd[41]*-1
cohensd[42] <- cohensd[42]*-1
cohensd[43] <- cohensd[43]*-1
cohensd[44] <- cohensd[44]*-1
cohensd[45] <- cohensd[45]*-1



df$effestnew[df$input %in% "smd" & df$output %in% "d"] <- cohensd


# Input: smd / Output: g --------------------------------------------------
dftemp <- df[df$input=="smd" & df$output=="g",]

hedgesg <- c()

for (k in 1:nrow(dftemp)) {
  
  hedgesg.temp <- c()
  
  for (l in 1:14) {
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    n <- nc+nt
    mc <-  as.numeric(eval(parse(text = paste0("dftemp$mc.",l))))[k]
    mt <-  as.numeric(eval(parse(text = paste0("dftemp$mt.",l))))[k]
    sdc <- as.numeric(eval(parse(text = paste0("dftemp$sdc.",l))))[k]
    sdt <- as.numeric(eval(parse(text = paste0("dftemp$sdt.",l))))[k]
    
    d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
    J <- 1 - (3 / (4 * (nc + nt - 2) - 1))
    g <- J * d
    
    hedgesg.temp <- c(hedgesg.temp,g)   # vector with all effect sizes within one primary study 
    
  }
  
  # take mean 
  hedgesg[k] <- mean(hedgesg.temp, na.rm=T) # mean of all effects from one primary study
  
}

hedgesg[30] <- hedgesg[30]*-1
hedgesg[31] <- hedgesg[31]*-1 # effect size reversed to be in right direction
hedgesg[32] <- hedgesg[32]*-1
hedgesg[33] <- hedgesg[33]*-1
hedgesg[35] <- hedgesg[35]*-1
hedgesg[37] <- hedgesg[37]*-1
hedgesg[38] <- hedgesg[38]*-1
hedgesg[51] <- hedgesg[51]*-1
hedgesg[52] <- hedgesg[52]*-1
hedgesg[53] <- hedgesg[53]*-1
hedgesg[54] <- hedgesg[54]*-1
hedgesg[55] <- hedgesg[55]*-1
hedgesg[57] <- hedgesg[57]*-1
hedgesg[58] <- hedgesg[58]*-1
hedgesg[59] <- hedgesg[59]*-1
hedgesg[60] <- hedgesg[60]*-1
hedgesg[61] <- hedgesg[61]*-1
hedgesg[62] <- hedgesg[62]*-1
hedgesg[63] <- hedgesg[63]*-1
hedgesg[64] <- hedgesg[64]*-1
hedgesg[65] <- hedgesg[65]*-1
hedgesg[68] <- hedgesg[68]*-1
hedgesg[69] <- hedgesg[69]*-1
hedgesg[70] <- hedgesg[70]*-1


df$effestnew[df$input %in% "smd" & df$output %in% "g"] <- hedgesg


# Input: smd / Output: r --------------------------------------------------
dftemp <- df[df$input=="smd" & df$output=="r",]

r.est <- c()

for (k in 1:nrow(dftemp)) {
  
  d.temp <- c()
  
  for (l in 1:14) {
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    n <- nc+nt
    mc <-  as.numeric(eval(parse(text = paste0("dftemp$mc.",l))))[k]
    mt <-  as.numeric(eval(parse(text = paste0("dftemp$mt.",l))))[k]
    sdc <- as.numeric(eval(parse(text = paste0("dftemp$sdc.",l))))[k]
    sdt <- as.numeric(eval(parse(text = paste0("dftemp$sdt.",l))))[k]
    
    d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
    
    d.temp <- c(d.temp,d)   # vector with all effect sizes within one primary study

  }
  
  # take mean 
  d <- mean(d.temp, na.rm=T)
  a <- (nc+nt)^2/(nc*nt)
  r <- d / (sqrt(d^2 + a))
  r.est <- c(r.est,r)
}

r.est[14] <- r.est[14]*-1
r.est[18] <- r.est[18]*-1
df$effestnew[df$input %in% "smd" & df$output %in% "r"] <- r.est



# Input: t / Output: d ----------------------------------------------------
dftemp <- df[df$input=="t" & df$output=="d",]

d.est <- c()

for (k in 1:nrow(dftemp)) {
  
  d.temp <- teststat <- c()
  
  for (l in 1:5) {
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    n <- nc+nt
    
    teststat <- c(teststat,as.numeric(eval(parse(text = paste0("dftemp$teststat.",l))))[k])
    teststat <- teststat[!is.na(teststat)]
    
    d <- teststat * sqrt((nc+nt)/(nc*nt))
    
    d.temp <- c(d.temp,d)   # vector with all effect sizes within one primary study 
    
  }
  
  # take mean 
  d.est <- c(d.est,mean(d.temp, na.rm=T)) # mean of all effects from one primary study
  
}

d.est[1] <- d.est[1]*-1
d.est[2] <- d.est[2]*-1
d.est[3] <- d.est[3]*-1

df$effestnew[df$input %in% "t" & df$output %in% "d"] <- d.est


# Input: t / Output: g ----------------------------------------------------
dftemp <- df[df$input=="t" & df$output=="g",]

g.est <- c()

for (k in 1:nrow(dftemp)) {
  
  g.temp <- teststat <- c()
  
  for (l in 1:5) {
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    n <- nc+nt
    
    teststat <- c(teststat,as.numeric(eval(parse(text = paste0("dftemp$teststat.",l))))[k])
    teststat <- teststat[!is.na(teststat)]
    
    d <- teststat * sqrt((nc+nt)/(nc*nt))
    J <- 1 - (3 / (4 * (nc + nt - 2) - 1))
    g <- J * d
    
    g.temp <- c(g.temp,g)   # vector with all effect sizes within one primary study 
    
  }
  
  # take mean 
  g.est <- c(g.est,mean(g.temp, na.rm=T)) # mean of all effects from one primary study
  
}

df$effestnew[df$input %in% "t" & df$output %in% "g"] <- g.est


# Input: t / Output: r ----------------------------------------------------
dftemp <- df[df$input=="t" & df$output=="r",]

r.est <- c()

for (k in 1:nrow(dftemp)) {
  
  r.temp <- teststat <- c()
  
  for (l in 1:5) {
    
    nc <- as.numeric(dftemp$ncnew)[k]
    nt <- as.numeric(dftemp$ntnew)[k]
    n <- nc+nt
    
    teststat <- c(teststat,as.numeric(eval(parse(text = paste0("dftemp$teststat.",l))))[k])
    teststat <- teststat[!is.na(teststat)]
    
    d <- teststat * sqrt((nc+nt)/(nc*nt))
    a <- (nc+nt)^2/(nc*nt)
    r <- d / (sqrt(d^2 + a))
    
    r.temp <- c(r.temp,r)   # vector with all effect sizes within one primary study 
    
  }
  
  # take mean 
  r.est <- c(r.est,mean(r.temp, na.rm=T)) # mean of all effects from one primary study
  
}

df$effestnew[df$input %in% "t" & df$output %in% "r"] <- r.est



# Input: Z / Output: g ----------------------------------------------------
dftemp <- df[df$input=="Zstatistic" & df$output=="g",]

nc <- as.numeric(dftemp$ncnew)
nt <- as.numeric(dftemp$ntnew)
n <- nc+nt

r <- as.numeric(dftemp$teststat.1) / sqrt(n)
d <- (2*r)/sqrt(1-r^2)
J <- 1 - (3 / (4 * (nc + nt - 2) - 1))
g <- J * d
g <- g*-1   # effect size reversed to be in right direction

df$effestnew[df$input %in% "Zstatistic" & df$output %in% "g"] <- g


# Input: Z / Output: r ----------------------------------------------------
dftemp <- df[df$input=="Zstatistic" & df$output=="r",]

nc <- as.numeric(dftemp$ncnew)
nt <- as.numeric(dftemp$ntnew)
n <- nc+nt

r <- as.numeric(dftemp$teststat.1) / sqrt(n)

df$effestnew[df$input %in% "Zstatistic" & df$output %in% "r"] <- r



# OUTPUT: D ---------------------------------------------------------------
dftemp <- df[df$output=="d",]

z.temp <- c()

for (k in 1:nrow(dftemp)) {
  
  if (dftemp$ncnew[k] == "NA") {
    
    A <- 4
    
  } else {
    
    A <- (as.numeric(dftemp$ncnew[k]) + as.numeric(dftemp$ntnew[k]))^2 / (as.numeric(dftemp$ncnew[k]) * as.numeric(dftemp$ntnew[k]))
    
  }
  
  d <- as.numeric(dftemp$effestnew[k])
  
  r <- d / (sqrt((d^2) + A))
  
  z <- 0.5 * log((1 + r) / (1 - r))
  
  z.temp <- c(z.temp,z)
  
}

df$znew[df$output=="d"] <- z.temp


#dftemp <- df[471,]

# OUTPUT: G ---------------------------------------------------------------
dftemp <- df[df$output=="g",]

A <- z.temp <- c()

for (k in 1:nrow(dftemp)) {
  
  if (dftemp$ncnew[k] == "NA") {
    
    A <- 4
    
  } else {
    
    A <- (as.numeric(dftemp$ncnew[k]) + as.numeric(dftemp$ntnew[k]))^2 / (as.numeric(dftemp$ncnew[k]) * as.numeric(dftemp$ntnew[k]))
    
  }
  
  g <- as.numeric(dftemp$effestnew[k])
  
  J <- (1 - 3 / (4 * (dftemp$nnew[k] - 2) - 1))
  
  d <- g / J                                    
  
  r <- d / (sqrt((d^2) + A))
  
  z <- 0.5 * log((1 + r) / (1 - r))
  
  z.temp <- c(z.temp,z)
  
}

df$znew[df$output=="g"] <- z.temp



# OUTPUT: R ---------------------------------------------------------------
dftemp <- df[df$output=="r",]

z.temp <- c()

for (k in 1:nrow(dftemp)) {
  
  r <- as.numeric(dftemp$effestnew[k])
  
  z <- 0.5 * log((1 + r) / (1 - r))
  
  z.temp <- c(z.temp,z)
  
}

df$znew[df$output=="r"] <- z.temp



# OUTPUT: z ---------------------------------------------------------------
df$znew[df$output=="z"] <- df$effestnew[df$output %in% "z"]


# Other transformations ---------------------------------------------------
df$effestnew <- as.numeric(df$effestnew)
df$znew <- as.numeric(df$znew)


# REPORTED EFFECT SIZES ---------------------------------------------------
#
# efftype: D ---------------------------------------------------------------
dftemp <- df[df$efftype == "d", ]

z.temp <- c()

for (k in 1:nrow(dftemp)) {
 if (dftemp$nc[k] == "NA") {
   A <- 4

 } else {
   A <-
     (as.numeric(dftemp$nc[k]) + as.numeric(dftemp$nt[k])) ^ 2 / (as.numeric(dftemp$nc[k]) * as.numeric(dftemp$nt[k]))

 }

 d <- as.numeric(dftemp$effest[k])

 r <- d / (sqrt((d ^ 2) + A))

 z <- 0.5 * log((1 + r) / (1 - r))

 z.temp <- c(z.temp, z)

}

df$z[df$efftype == "d"] <- z.temp


#dftemp <- df[471,]

# efftype: G ---------------------------------------------------------------
dftemp <- df[df$efftype == "g", ]

A <- z.temp <- c()

for (k in 1:nrow(dftemp)) {
 if (dftemp$nc[k] == "NA") {
   A <- 4

 } else {
   A <-
     (as.numeric(dftemp$nc[k]) + as.numeric(dftemp$nt[k])) ^ 2 / (as.numeric(dftemp$nc[k]) * as.numeric(dftemp$nt[k]))

 }

 g <- as.numeric(dftemp$effest[k])

 J <- (1 - 3 / (4 * (dftemp$n[k] - 2) - 1))

 d <- g / J

 r <- d / (sqrt((d ^ 2) + A))

 z <- 0.5 * log((1 + r) / (1 - r))

 z.temp <- c(z.temp, z)

}

df$z[df$efftype == "g"] <- z.temp



# efftype: R ---------------------------------------------------------------
dftemp <- df[df$efftype == "r", ]

z.temp <- c()

for (k in 1:nrow(dftemp)) {
 r <- as.numeric(dftemp$effest[k])

 z <- 0.5 * log((1 + r) / (1 - r))

 z.temp <- c(z.temp, z)

}

df$z[df$efftype == "r"] <- z.temp



# efftype: z ---------------------------------------------------------------
df$z[df$efftype == "z"] <- df$effest[df$efftype %in% "z"]



# Discrepancies -----------------------------------------------------------

# Transform discrepancies
g <- c(.05,.149,.150,.249,.250)
J <- (1 - 3 / (4 * (64 - 2) - 1)) # Assuming N = 64
d <- g / J;d 
A <- (31 + 31)^2 / (31*31);A
r <- d / (sqrt((d^2) + A));r
z <- 0.5 * log((1 + r) / (1 - r));z


# Hedges g - small [.050 - .149], moderate [.150 - .249] and large [.250 - inf] 
# Cohens d - small [.051 - .151], moderate [.152 - .252] and large [.253 - inf].
# Correlation - small [.025 - .075], moderate [.076 - .125], and large [.126 - inf].
# Fisher's z - small [.025 - .075], moderate [.076 - .125], and large [.126 - inf].


# Fill in discrepancies ---------------------------------------------------
df$disc.eff <- df$effest - df$effestnew
df$disc.z <- df$z - df$znew
df$disc.n <- df$n - df$nnew


# Fill in discrepancy category 
for (i in 1:nrow(df)) {

if (df$efftype[i] == "g") {
  
  if (abs(df$disc.eff[i]) < 0.05) {
    df$disccat.eff[i] = 0
  } else if (abs(df$disc.eff[i]) >= 0.050 & abs(df$disc.eff[i]) <= 0.149) {
    df$disccat.eff[i] = 1
  } else if (abs(df$disc.eff[i]) > 0.149 & abs(df$disc.eff[i]) <= 0.249) {
    df$disccat.eff[i] = 2
  } else if (abs(df$disc.eff[i]) > 0.249) {
    df$disccat.eff[i] = 3
  } else {
    df$disccat.eff[i] = "check"
  }
  
}

if (df$efftype[i] == "d") {
  
  if (abs(df$disc.eff[i]) < 0.051) {
    df$disccat.eff[i] = 0
  } else if (abs(df$disc.eff[i]) >= 0.051 & abs(df$disc.eff[i]) <= 0.151) {
    df$disccat.eff[i] = 1
  } else if (abs(df$disc.eff[i]) > 0.151 & abs(df$disc.eff[i]) <= 0.252) {
    df$disccat.eff[i] = 2
  } else if (abs(df$disc.eff[i]) > 0.252) {
    df$disccat.eff[i] = 3
  } else {
    df$disccat.eff[i] = "check"
  }
  
}

if (df$efftype[i] == "r" | df$efftype[i] == "z") {
  
  if (abs(df$disc.eff[i]) < 0.025) {
    df$disccat.eff[i] = 0
  } else if (abs(df$disc.eff[i]) >= 0.025 & abs(df$disc.eff[i]) <= 0.075) {
    df$disccat.eff[i] = 1
  } else if (abs(df$disc.eff[i]) > 0.075 & abs(df$disc.eff[i]) <= 0.125) {
    df$disccat.eff[i] = 2
  } else if (abs(df$disc.eff[i]) > 0.125) {
    df$disccat.eff[i] = 3
  } else {
    df$disccat.eff[i] = "check"
  }
  
}
  
  if (abs(df$disc.z[i]) < 0.025) {
    df$disccat.z[i] = 0
  } else if (abs(df$disc.z[i]) >= 0.025 & abs(df$disc.z[i]) <= 0.075) {
    df$disccat.z[i] = 1
  } else if (abs(df$disc.z[i]) > 0.075 & abs(df$disc.z[i]) <= 0.125) {
    df$disccat.z[i] = 2
  } else if (abs(df$disc.z[i]) > 0.125) {
    df$disccat.z[i] = 3
  } else {
    df$disccat.z[i] = "check"
  }
  
  
}



# Reverse effect sizes so all effects are in hypothesized direction -------
df$effest.exp <- df$effest
df$effestnew.exp <- df$effestnew
df$z.exp <- df$z
df$znew.exp <- df$znew

reverse <- c("Babbage","Balliet","deWit","Fischer","Woodin","Woodley","Yoon")

df$effest.exp[which(df$meta %in% reverse)] <- df$effest.exp[which(df$meta %in% reverse)]*-1
df$effestnew.exp[which(df$meta %in% reverse)] <- df$effestnew.exp[which(df$meta %in% reverse)]*-1
df$z.exp[which(df$meta %in% reverse)] <- df$z.exp[which(df$meta %in% reverse)]*-1
df$znew.exp[which(df$meta %in% reverse)] <- df$znew.exp[which(df$meta %in% reverse)]*-1

# Discrepancy categories not coinciding between original effect and after transformation to z.
check <- df[which(!df$disccat.eff == df$disccat.z),]


# Save file in xlsx -------------------------------------------------------
write_xlsx(df,"codebook-primary-studies-final-complete.xlsx",col_names=T)

# Save file in csv --------------------------------------------------------
write.table(df, file = "codebook-primary-studies-final.csv", row.names=F, col.names=T, sep=' ')
