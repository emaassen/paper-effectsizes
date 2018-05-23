rm(list = ls())
setwd("C:/Users/s421506/tiu/research/effectsizes/data-per-ma")
x <- c("metafor","MAc","MAd")
#lapply(x,install.packages(x),character.only=T)
lapply(x,library,character.only=T)

# Adesope -----------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("adesope_complete.csv", header=T, sep = ';')
df$d <- df$g / (1-3/(4*(df$n - 2) - 1))   
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r))
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "adesope_complete.csv", row.names=FALSE, sep=";")

df <- read.table("adesope_subset.csv", header=T, sep=';')
df$d <- df$g / (1-3/(4*(df$n - 2) - 1))                                
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                      
write.table(df, file = "adesope_subset.csv", row.names=FALSE, sep=";")

# Dependency check with g
df <- read.table("adesope_complete.csv", header=T, sep=';')
df <- df[order(df$g),]                                # order df by effect size

J <- (1-3/(4*(df$n - 2) - 1))
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d
df$vg <- J^2 * df$vd                                  # variance hedges' g        

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 18) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s2a1","s2a2","s2a3","s2a4","s2a5") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(g, vg, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4)                          
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value 
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies               
               df[x,]$g,df[x+1,]$g,                   # effect size of two studies               
               as.character(df[x,]$a1),               # authors study 1            
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:18], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Alfieri -----------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("alfieri_complete.csv", header=T, sep=';')
df$n <- df$n1+df$n2
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2) # Unequal groups
df$r <- df$d / (sqrt(((df$d^2)+A)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r))
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "alfieri_complete.csv", row.names=FALSE, sep=";")

df <- read.table("alfieri_subset.csv", header=T, sep=';')
df$n <- df$n1+df$n2
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)
df$r <- df$d / (sqrt(((df$d^2)+A)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z 
write.table(df, file = "alfieri_subset.csv", row.names=FALSE, sep=";")

# Dependency check with d
df <- read.table("alfieri_complete.csv", header=T, sep=',')
df <- df[order(df$d),]                                # order df by effect size
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 18) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s2a1","s2a2","s2a3","s2a4","s2a5") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(d, vd, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$d,df[x+1,]$d,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:18], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Babbage -----------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("babbage_complete.csv", header=T, sep=',')
df$n <- df$n1+df$n2
df$d <- df$g / (1-3/(4*(df$n - 2) - 1))                                
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)                   # Unequal groups
df$r <- df$d / (sqrt(((df$d^2)+A)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r))
df$vz <- 1 / (df$n-3)                                  # variance fischer's z 
df$zexp <- df$z*-1                                     # Effect size direction is not in hypothesized direction
write.table(df, file = "babbage_complete.csv", row.names=FALSE, sep=",")

df <- read.table("babbage_subset.csv", header=T, sep=',')
df$n <- df$n1+df$n2
df$d <- df$g / (1-3/(4*(df$n - 2) - 1))                                
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)
df$r <- df$d / (sqrt(((df$d^2)+A)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r))    
df$vz <- 1 / (df$n-3)                                  # variance fischer's z 
df$zexp <- df$z*-1                                     
write.table(df, file = "babbage_subset.csv", row.names=FALSE, sep=",")

# Dependency check with g
df <- read.table("babbage_complete.csv", header=T, sep=',')
df <- df[order(df$g),]                                 # order df by effect size

J <- 1-3/(4*(df$n - 2) - 1)
df$d <- df$g/J
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)               # variance cohen's d
df$vg <- J^2 * df$vd                                   # variance hedges' g

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 18) # empty matrices 

colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s2a1","s2a2","s2a3","s2a4","s2a5") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(g, vg, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$g,df[x+1,]$g,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:18], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Balliet -----------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("balliet_complete.csv", header=T, sep=',')
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r))    
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1                                    # Effect size direction is not in hypothesized direction
write.table(df, file = "balliet_complete.csv", row.names=FALSE, sep=",")

df <- read.table("balliet_subset.csv", header=T, sep=',')
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r))  
df$vz <- 1 / (df$n-3)                                  # variance fischer's z 
df$zexp <- df$z*-1   
write.table(df, file = "balliet_subset.csv", row.names=FALSE, sep=",")

# Dependency check with d
df <- read.table("balliet_complete.csv", header=T, sep=',')
df <- df[order(df$d),]                                # order df by effect size

df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d
dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 

colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(d, vd, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$d,df[x+1,]$d,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Benish ------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("benish_complete.csv", header=T, sep=';')
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r))     
df$vz <- 1 / (df$n-3)                                  # variance fischer's z 
df$zexp <- df$z                                        # Effect size direction is in hypothesized direction
write.table(df, file = "benish_complete.csv", row.names=FALSE, sep=";")

df <- read.table("benish_subset.csv", header=T, sep=',')
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r))     
df$vz <- 1 / (df$n-3)                                  # variance fischer's z 
df$zexp <- df$z  
write.table(df, file = "benish_subset.csv", row.names=FALSE, sep=",")

# Effect size calculation - aggregating
MA <- data.frame(id=c(1,1,1,1,1),
                 es <- c(0.821697, 0.812402, 0.879328, 0.980306, 0.66471),
                 var.es <- c(0.1399, 0.1396, 0.1415, 0.14452, 0.1362),
                 nT <- c(15,15,15,15,15),
                 nC <- c(16,16,16,16,16))

MA1 <- agg(id=id, es=es, var=var.es, n.1=nT, n.2=nC, method="BHHR",data=MA);MA1

# Formula from Wampold et al.
es <- c(0.821697, 0.812402, 0.879328, 0.980306, 0.66471)
var.es <- c(0.1399, 0.1396, 0.1415, 0.14452, 0.1362)
sum(es/var.es)/sum(1/var.es)

# Dependency check with d
df <- read.table("benish_complete.csv", header=T, sep=',')
df <- df[order(df$d),]                               # order df by effect size
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)             # variance cohen's d

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(d, vd, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$d,df[x+1,]$d,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Berry1 ------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("berry1_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))     
df$vz <- 1 / (df$n-3)                                # variance fischer's z 
df$zexp <- df$z                                      # Effect size direction is in hypothesized direction
write.table(df, file = "berry1_complete.csv", row.names=FALSE, sep=";")

df <- read.table("berry1_subset.csv", header=T, sep=',')
df$z <- 0.5*log((1+df$r)/(1-df$r))     
df$vz <- 1 / (df$n-3)                                # variance fischer's z 
df$zexp <- df$z  
write.table(df, file = "berry1_subset.csv", row.names=FALSE, sep=",")

# Dependency check with r
df <- read.table("berry1_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                               # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                   # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 16) # empty matrices 

colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s2a1","s2a2","s2a3","s2a4") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:16], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Berry2 ------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("berry2_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))  
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "berry 2_complete.csv", row.names=FALSE, sep=";")

df <- read.table("berry2_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))     
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z   
write.table(df, file = "berry2_subset.csv", row.names=FALSE, sep=";")

# Effect size calculation - aggregating
MA <- data.frame(id=c(1,1),
                 r <- c(0.53,0.30),
                 n <- c(583,583))

MA1 <- agg(id = id, r = r, n = n, data=MA, cor = 0.61);MA1

# Dependency check with r
df <- read.table("berry2_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                 # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                     # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 18) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s2a1","s2a2","s2a3","s2a4","s2a5") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:18], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Card --------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("card_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))     
df$vz <- 1 / (df$n-3)                                  # variance fischer's z 
df$zexp <- df$z                                        # Effect size direction is in hypothesized direction
write.table(df, file = "card_complete.csv", row.names=FALSE, sep=";")

df <- read.table("card_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))    
df$vz <- 1 / (df$n-3)                                  # variance fischer's z 
df$zexp <- df$z  
write.table(df, file = "card_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("card_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                 # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                     # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Crook -------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("crook_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))     
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "crook_complete.csv", row.names=FALSE, sep=";")

df <- read.table("crook_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))     
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "crook_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("crook_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# DeWit -------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("dewit_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))  
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1                                    # Effect size direction is not in hypothesized direction
write.table(df, file = "dewit_complete.csv", row.names=FALSE, sep=";")

df <- read.table("dewit_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))     
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1  
write.table(df, file = "dewit_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("dewit_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Elsequest ---------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("elsequest_complete.csv", header=T, sep=';')
df$n <- df$n1+df$n2
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)
df$r <- df$d / (sqrt(((df$d^2)+A)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "elsequest_complete.csv", row.names=FALSE, sep=";")

df <- read.table("elsequest_subset.csv", header=T, sep=';')
df$n <- df$n1+df$n2
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)
df$r <- df$d / (sqrt(((df$d^2)+A)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "elsequest_subset.csv", row.names=FALSE, sep=";")

# Dependency check with d
df <- read.table("elsequest_complete.csv", header=T, sep=',')
df <- df[order(df$d),]                                # order df by effect size
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(d, vd, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$d,df[x+1,]$d,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Farber ------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("farber_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "farber_complete.csv", row.names=FALSE, sep=";")

df <- read.table("farber_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "farber_subset.csv", row.names=FALSE, sep=";")

# Farber, aggregating correlated effect sizes.
id <- rep(1,2)
n <- rep(31,2)
r <- c(0.30,0.05)
df <- data.frame(id,n,r)
agg(id = id, r = r, n = n, data=df, cor = 0.5)
agg(id = id, r = r, n = n, data=df)

# Dependency check with r
df <- read.table("farber_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Fischer -----------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("fischer_complete.csv", header=T, sep=';')
df$n <- df$n1+df$n2
df$d <- df$g / (1-3/(4*(df$n - 2) - 1))  
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)
df$r <- df$d / (sqrt(((df$d^2)+A)))       
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1                                    # Effect size direction is not in hypothesized direction
write.table(df, file = "fischer_complete.csv", row.names=FALSE, sep=";")

df <- read.table("fischer_subset.csv", header=T, sep=';')
df$n <- df$n1+df$n2
df$d <- df$g / (1-3/(4*(df$n - 2) - 1))                                 
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)
df$r <- df$d / (sqrt(((df$d^2)+A)))                                         
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1  
write.table(df, file = "fischer_subset.csv", row.names=FALSE, sep=";")

# Dependency check with g
df <- read.table("fischer_complete.csv", header=T, sep=',')
df <- df[order(df$g),]                                # order df by effect size
J <- 1-3/(4*(df$n - 2) - 1)
df$d <- df$g/J
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d
df$vg <- J^2 * df$vd

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s2a5","s2a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(g, vg, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$g,df[x+1,]$g,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Fox ---------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("fox_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "fox_complete.csv", row.names=FALSE, sep=";")

df <- read.table("fox_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                      
write.table(df, file = "fox_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("fox_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Freund ------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("freund_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "freund_complete.csv", row.names=FALSE, sep=";")

df <- read.table("freund_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                      
write.table(df, file = "freund_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("freund_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Green -------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("green_complete.csv", header=T, sep=';')
df$n <- df$n1+df$n2  
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)
df$r <- df$d / (sqrt(((df$d^2)+A))) 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "green_complete.csv", row.names=FALSE, sep=";")

df <- read.table("green_subset.csv", header=T, sep=';')
df$n <- df$n1+df$n2  
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)
df$r <- df$d / (sqrt(((df$d^2)+A))) 
df$z <- 0.5*log((1+df$r)/(1-df$r))  
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "green_subset.csv", row.names=FALSE, sep=";")

# Dependency check with d
df <- read.table("green_complete.csv", header=T, sep=',')
df <- df[order(df$d),]                                # order df by effect size
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(d, vd, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$d,df[x+1,]$d,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:18], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Hallion -----------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("hallion_complete.csv", header=T, sep=';')
df$d <- df$g / (1-3/(4*(df$n - 2) - 1))  
df$r <- df$d / (sqrt((df$d^2)+4))
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "hallion_complete.csv", row.names=FALSE, sep=";")

df <- read.table("hallion_subset.csv", header=T, sep=';')
df$d <- df$g / (1-3/(4*(df$n - 2) - 1))  
df$r <- df$d / (sqrt((df$d^2)+4))
df$z <- 0.5*log((1+df$r)/(1-df$r))   
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "hallion_subset.csv", row.names=FALSE, sep=";")

# Dependency check with g
df <- read.table("hallion_complete.csv", header=T, sep=',')
df <- df[order(df$g),]                                # order df by effect size
J <- 1-3/(4*(df$n - 2) - 1)
df$d <- df$g/J
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d
df$vg <- J^2 * df$vd                                  # variance hedges' g

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(g, vg, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$g,df[x+1,]$g,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Ihle --------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("ihle_complete.csv", header=T, sep=';')
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)
df$r <- df$d / (sqrt(((df$d^2)+A)))
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "ihle_complete.csv", row.names=FALSE, sep=";")

df <- read.table("ihle_subset.csv", header=T, sep=';')
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2)
df$r <- df$d / (sqrt(((df$d^2)+A)))
df$z <- 0.5*log((1+df$r)/(1-df$r))  
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "ihle_subset.csv", row.names=FALSE, sep=";")

# Dependency check with d
df <- read.table("ihle_complete.csv", header=T, sep=',')
df <- df[order(df$d),]                                # order df by effect size
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(d, vd, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$d,df[x+1,]$d,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Koenig ------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("koenig_complete.csv", header=T, sep=';')
df$d <- df$g / (1-3/(4*(df$n - 2) - 1)) 
df$r <- df$d / (sqrt(((df$d^2)+4)))
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "koenig_complete.csv", row.names=FALSE, sep=";")

df <- read.table("koenig_subset.csv", header=T, sep=';')
df$d <- df$g / (1-3/(4*(df$n - 2) - 1)) 
df$r <- df$d / (sqrt(((df$d^2)+4)))
df$z <- 0.5*log((1+df$r)/(1-df$r))  
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "koenig_subset.csv", row.names=FALSE, sep=";")

# Dependency check with g
df <- read.table("koenig_complete.csv", header=T, sep=',')
df <- df[order(df$g),]                                # order df by effect size
J <- 1-3/(4*(df$n - 2) - 1)
df$d <- df$g/J
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d
df$vg <- J^2 * df$vd                                  # variance hedges' g

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(g, vg, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$g,df[x+1,]$g,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Kolden ------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("kolden_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "kolden_complete.csv", row.names=FALSE, sep=";")

df <- read.table("kolden_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))  
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "kolden_subset.csv", row.names=FALSE, sep=";")

# Aggregating correlated effect sizes.
# transform outcomes to Cohen's d as per Re (2015)- https://www.researchgate.net/profile/A_C_Del_Re/publication/271766790_A_Practical_Tutorial_on_Conducting_Meta-Analysis_in_R/links/54d119e30cf25ba0f040b2f3/A-Practical-Tutorial-on-Conducting-Meta-Analysis-in-R.pdf
# Hoyt  and  Del  Re  found  that  using  Hedges' g (and Vg) is  preferred  when using Borenstein's method (BHHR)
# Cohen's d (and Vd) is preferred with Gleser and Olkin's methods (GO1 or GO2).  

n1 <- 9  # treatment group
n2 <- 9  # control group
r <- c(.40,.47,.61,.25,.4,.41,.37,.36,.33)
d <- es <- (2*r)/(sqrt(1-r^2))  
d.var <- var.es <- (1/n1)+(1/(n2))+(d^2/(2*(n1+n2)))
cor <- .25

MA <- data.frame(id=rep(1, length(r)),
                 es = d,
                 var.es = var.es,
                 nT = rep(n1, length(r)),
                 nC = rep(n2, length(r))); MA

#MA1 uses pooled standard deviation
MA1 <- agg(id=id, es=es, var=var.es, n.1=nT, n.2=nC, cor = cor, method="GO1",data=MA);MA1

#MA2 uses second group SD (control group) - we didn't use this one
#MA2 <- agg(id=id, es=es, var=var.es, n.1=nT, n.2=nC, cor = cor, method="GO2",data=MA);MA2

#MA3 uses second group SD (treatment group) - reverse n1, r, etc. first - we didn;t use this one
#n0 <- n1
#n1 <- n2
#n2 <- n0

#MA <- data.frame(id=rep(1, length(r)),
#es = d,
#var.es = var.es,
#nT = rep(n1, length(r)),
#nC = rep(n2, length(r))); MA  

#MA3 <- agg(id=id, es=es, var=var.es, n.1=nT, n.2=nC, cor = cor, method="GO2",data=MA);MA3

# Lucassen ----------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("lucassen_complete.csv", header=T, sep = ";")
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "lucassen_complete.csv", row.names=FALSE, sep=";")

df <- read.table("lucassen_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r))  
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "lucassen_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("lucassen_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation 

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Mol ---------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("mol_complete.csv", header=T, sep=';')
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "mol_complete.csv", row.names=FALSE, sep=";")

df <- read.table("mol_subset.csv", header=T, sep=';')
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "mol_subset.csv", row.names=FALSE, sep=";")

# Dependency check with g
df <- read.table("mol_complete.csv", header=T, sep=',')
df <- df[order(df$z),]                                # order df by effect size

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(z, vz, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$z,df[x+1,]$z,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Morgan ------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("morgan_complete.csv", header=T, sep=';')
df$r <- df$d / (sqrt(((df$d^2)+4))) 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "morgan_complete.csv", row.names=FALSE, sep=";")

df <- read.table("morgan_subset.csv", header=T, sep=';')
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "morgan_subset.csv", row.names=FALSE, sep=";")

# Dependency check with d
df <- read.table("morgan_complete.csv", header=T, sep=',')
df <- df[order(df$d),]                                # order df by effect size
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(d, vd, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$d,df[x+1,]$d,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Munder ------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("munder_complete.csv", header=T, sep=';')
df$d <- df$g / (1-3/(4*(df$n - 2) - 1))   
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "munder_complete.csv", row.names=FALSE, sep=";")

df <- read.table("munder_subset.csv", header=T, sep=';')
df$d <- df$g / (1-3/(4*(df$n - 2) - 1)) 
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                      
write.table(df, file = "munder_subset.csv", row.names=FALSE, sep=";")

# Dependency check with d
df <- read.table("munder_complete.csv", header=T, sep=',')
df <- df[order(df$d),]                                # order df by effect size
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(d, vd, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$d,df[x+1,]$d,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Piet --------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("piet_complete.csv", header=T, sep=';')
df$d <- df$g / (1-3/(4*(df$n - 2) - 1)) 
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "piet_complete.csv", row.names=FALSE, sep=";")

df <- read.table("piet_subset.csv", header=T, sep=';')
df$d <- df$g / (1-3/(4*(df$n - 2) - 1)) 
df$r <- df$d / (sqrt(((df$d^2)+4)))                                 
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z   
write.table(df, file = "piet_subset.csv", row.names=FALSE, sep=";")

# Dependency check with g
df <- read.table("piet_complete.csv", header=T, sep=',')
df <- df[order(df$g),]                                # order df by effect size
J <- 1-3/(4*(df$n - 2) - 1)
df$d <- df$g/J
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d
df$vg <- J^2 * df$vd                                  # variance hedges'g

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s1a6","s1a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(g, vg, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$g,df[x+1,]$g,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Smith -------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("smith_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "smith_complete.csv", row.names=FALSE, sep=";")

df <- read.table("smith_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "smith_subset.csv", row.names=FALSE, sep=";")

# Dependency check with 
df <- read.table("smith_complete.csv", header=T, sep=',')
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Tillman -----------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("tillman_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "tillman_complete.csv", row.names=FALSE, sep=";")

df <- read.table("tillman_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z   
write.table(df, file = "tillman_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("tillman_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size

df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Toosi -------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("toosi_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "toosi_complete.csv", row.names=FALSE, sep=";")

df <- read.table("toosi_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       
write.table(df, file = "toosi_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("toosi_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# VanIddekinge ------------------------------------------------------------

# Transform to Fisher's z for total and subset
df <- read.table("vaniddekinge_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "vaniddekinge_complete.csv", row.names=FALSE, sep=";")

df <- read.table("vaniddekinge_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z   
write.table(df, file = "vaniddekinge_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("vaniddekinge_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Webb --------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("webb_complete.csv", header=T, sep=';')
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2) # Unequal groups
df$r <- df$d / (sqrt(((df$d^2)+A)))
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z                                       # Effect size direction is in hypothesized direction
write.table(df, file = "webb_complete.csv", row.names=FALSE, sep=";")

df <- read.table("webb_subset.csv", header=T, sep=';')
A <- ((df$n1+df$n2)^2)/(df$n1*df$n2) # Unequal groups
df$r <- df$d / (sqrt(((df$d^2)+A)))
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z  
write.table(df, file = "webb_subset.csv", row.names=FALSE, sep=";")

# Dependency check with d
df <- read.table("webb_complete.csv", header=T, sep=',')
df <- df[order(df$d),]                                # order df by effect size
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 22) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s1a7","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6","s2a7") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(d, vd, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$d,df[x+1,]$d,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x,]$a7),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6),
               as.character(df[x+1,]$a7))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:22], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Woodin ------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("woodin_complete.csv", header=T, sep=';')
df$r <- df$d / (sqrt(((df$d^2)+4)))
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1                                    # Effect size direction is not in hypothesized direction
write.table(df, file = "woodin_complete.csv", row.names=FALSE, sep=";")

df <- read.table("woodin_subset.csv", header=T, sep=';')
df$r <- df$d / (sqrt(((df$d^2)+4)))
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1   
write.table(df, file = "woodin_subset.csv", row.names=FALSE, sep=";")

# Dependency check with d
df <- read.table("woodin_complete.csv", header=T, sep=',')
df <- df[order(df$d),]                                # order df by effect size
df$vd <- df$n/(df$n/2)^2+df$d^2/(2*df$n)              # variance cohen's d

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(d, vd, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$d,df[x+1,]$d,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Woodley -----------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("woodley_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1                                    # Effect size direction is not in hypothesized direction
write.table(df, file = "woodley_complete.csv", row.names=FALSE, sep=";")

df <- read.table("woodley_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1   
write.table(df, file = "woodley_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("woodley_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]


# Yoon --------------------------------------------------------------------
# Transform to Fisher's z for total and subset
df <- read.table("yoon_complete.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1                                    # Effect size direction is not in hypothesized direction
write.table(df, file = "yoon_complete.csv", row.names=FALSE, sep=";")

df <- read.table("yoon_subset.csv", header=T, sep=';')
df$z <- 0.5*log((1+df$r)/(1-df$r)) 
df$vz <- 1 / (df$n-3)                                 # variance fischer's z 
df$zexp <- df$z*-1  
write.table(df, file = "yoon_subset.csv", row.names=FALSE, sep=";")

# Dependency check with r
df <- read.table("yoon_complete.csv", header=T, sep=',')
df <- df[order(df$r),]                                # order df by effect size
df$vr <- ((1-(df$r^2))^2)/(df$n-1)                    # variance correlation

dep <- depfinal <- matrix(NA, nrow = nrow(df), ncol = 20) # empty matrices 
colnames(dep) <- colnames(depfinal) <- c("s1","s2","QE","QEp","n1","n2","eff1","eff2","s1a1","s1a2","s1a3","s1a4","s1a5","s1a6","s2a1","s2a2","s2a3","s2a4","s2a5","s2a6") #Empty frame to store results
x <- 1

for (i in nrow(df)) {  
  
  subset <- df[x:(x+1),]                              # subset of two studies
  res <- rma(r, vr, data = subset, method="FE")       # FE MA
  
  dep[i,] <- c(x, (x+1),                              # two studies that are compared (e.g., 1 with 2, 2 with 3, 3 with 4) 
               round(res$QE,3), round(res$QEp,3),     # QE statistic and p-value     
               df[x,]$n,df[x+1,]$n,                   # sample size of two studies
               df[x,]$r,df[x+1,]$r,                   # effect size of two studies  
               as.character(df[x,]$a1),               # authors study 1
               as.character(df[x,]$a2),
               as.character(df[x,]$a3),
               as.character(df[x,]$a4),
               as.character(df[x,]$a5),
               as.character(df[x,]$a6),
               as.character(df[x+1,]$a1),             # authors study 2
               as.character(df[x+1,]$a2),
               as.character(df[x+1,]$a3),
               as.character(df[x+1,]$a4),
               as.character(df[x+1,]$a5),
               as.character(df[x+1,]$a6))
  x <- x + 1
} 

dep <- dep[-nrow(dep),]                               # removes last row, which is comparison of last study with nothing. 
dep[dep == ""] <- NA
dep2 <- dep[dep[, "QEp"] >= 0.80,]                    # only keep rows where QE p-val > .80 (i.e., homogeneous studies)

for (i in nrow(df)) {                                 # if any authors are duplicates, keep the row
  if (any(duplicated(dep2[i,9:20], incomparables=NA) == TRUE)) {
    depfinal[i,] <- dep2[i,]
  } 
}
# remove emptry rows
depfinal <- depfinal[!apply(is.na(depfinal) | depfinal == "", 1, all),]
