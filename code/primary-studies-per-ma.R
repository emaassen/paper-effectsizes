rm(list = ls())
setwd("C:/Users/s421506/tiu/research/effectsizes/data-per-ma")
packages <- c("metafor","dplyr")
#lapply(packages,install.packages(packages),character.only=T)
lapply(packages,library,character.only=T)

outlier.total <- regular.total <- c() # empty vector to store all studies that are outliers or non-outliers

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

# Adesope -----------------------------------------------------------------
df <- read.table("adesope_complete.csv", header=T, sep = ';')

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

# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg))
write.table(df, file = "adesope_complete.csv", row.names=FALSE, sep=";")

df <- read.table("adesope_subset.csv", header=T, sep=";")

J <- 1 - (3 / (4 * df$n - 9))
df$d <- df$g / J                                        # cohen's d
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd                                    # variance hedges' g 

write.table(df, file = "adesope_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg))
write.table(df, file = "alfieri_complete.csv", row.names=FALSE, sep=";")

df <- read.table("alfieri_subset.csv", header=T, sep=";")

df$n <- df$n1 + df$n2
dfs <- df$n1 + df$n2 - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J                                        # cohen's d
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd                                    # variance hedges' g 

write.table(df, file = "alfieri_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg))        # no possible dependent studies
write.table(df, file = "babbage_complete.csv", row.names=FALSE, sep=";")

df <- read.table("babbage_subset.csv", header=T, sep=";")

J <- 1 - (3 / (4 * df$n - 9))
df$d <- df$g / J                                        # cohen's d
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd                                    # variance hedges' g 

write.table(df, file = "babbage_subset.csv", row.names=FALSE, sep=";")


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


# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg))
write.table(df, file = "balliet_complete.csv", row.names=FALSE, sep=";")


df <- read.table("balliet_subset.csv", header=T, sep=";")
J <- 1 - (3 / (4 * df$n - 9))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
write.table(df, file = "balliet_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg))
write.table(df, file = "benish_complete.csv", row.names=FALSE, sep=";")

df <- read.table("benish_subset.csv", header=T, sep=";")

J <- 1 - (3 / (4 * df$n - 9))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 

write.table(df, file = "benish_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))               
write.table(df, file = "berry1_complete.csv", row.names=FALSE, sep=";")

df <- read.table("berry1_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)             # variance correlation r
write.table(df, file = "berry1_subset.csv", row.names=FALSE, sep=";")



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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))   
write.table(df, file = "berry2_complete.csv", row.names=FALSE, sep=";")

df <- read.table("berry2_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)             # variance correlation r
write.table(df, file = "berry2_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))   
write.table(df, file = "card_complete.csv", row.names=FALSE, sep=";")

df <- read.table("card_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "card_subset.csv", row.names=FALSE, sep=";")


# Crook -------------------------------------------------------------------

df <- read.table("crook_complete.csv", header=T, sep=";")

df$vr <- ((1 - (df$r^2))^2) / (df$n-1)             # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id,  study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))   
write.table(df, file = "crook_complete.csv", row.names=FALSE, sep=";")

df <- read.table("crook_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "crook_subset.csv", row.names=FALSE, sep=";")


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
# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))  
write.table(df, file = "dewit_complete.csv", row.names=FALSE, sep=";")

df <- read.table("dewit_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "dewit_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with d
depdf <- as.data.frame(dependency(df$d,df$vd))
write.table(df, file = "elsequest_complete.csv", row.names=FALSE, sep=";")

df <- read.table("elsequest_subset.csv", header=T, sep=";")
df$n <- df$n1 + df$n2
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
write.table(df, file = "elsequest_subset.csv", row.names=FALSE, sep=";")

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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))        # no possible dependent studies
write.table(df, file = "farber_complete.csv", row.names=FALSE, sep=";")

df <- read.table("farber_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "farber_subset.csv", row.names=FALSE, sep=";")



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

# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "fischer_complete.csv", row.names=FALSE, sep=";")

df <- read.table("fischer_subset.csv", header=T, sep=";")
df$n <- df$n1 + df$n2
dfs <- df$n1 + df$n2 - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
write.table(df, file = "fischer_subset.csv", row.names=FALSE, sep=";")



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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "fox_complete.csv", row.names=FALSE, sep=";")

df <- read.table("fox_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "fox_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "freund_complete.csv", row.names=FALSE, sep=";")

df <- read.table("freund_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "freund_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with d
depdf <- as.data.frame(dependency(df$d,df$vd))
write.table(df, file = "green_complete.csv", row.names=FALSE, sep=";")

df <- read.table("green_subset.csv", header=T, sep=";")
df$n <- df$n1 + df$n2
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
write.table(df, file = "green_subset.csv", row.names=FALSE, sep=";")


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


# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "hallion_complete.csv", row.names=FALSE, sep=";")

df <- read.table("hallion_subset.csv", header=T, sep=";")
J <- 1 - (3 / (4 * (df$n - 1)))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd
write.table(df, file = "hallion_subset.csv", row.names=FALSE, sep=";")


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


# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "ihle_complete.csv", row.names=FALSE, sep=";")

df <- read.table("ihle_subset.csv", header=T, sep=";")
df$n <- df$n1 + df$n2
dfs <- df$n1 + df$n2 - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd
write.table(df, file = "ihle_subset.csv", row.names=FALSE, sep=";")


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


# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "koenig_complete.csv", row.names=FALSE, sep=";")

df <- read.table("koenig_subset.csv", header=T, sep=";")
dfs <- df$n - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
write.table(df, file = "koenig_subset.csv", row.names=FALSE, sep=";")



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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "kolden_complete.csv", row.names=FALSE, sep=";")

# Subset
df <- read.table("kolden_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "kolden_subset.csv", row.names=FALSE, sep=";")

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


# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "lucassen_complete.csv", row.names=FALSE, sep=";")

df <- read.table("lucassen_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "lucassen_subset.csv", row.names=FALSE, sep=";")


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


# Dependency check with r
depdf <- as.data.frame(dependency(df$z,df$vz))
write.table(df, file = "mol_complete.csv", row.names=FALSE, sep=";")

df <- read.table("mol_subset.csv", header=T, sep=";")
df$vz <- 1 / (df$n - 3)                              # variance fisher's z 
write.table(df, file = "mol_subset.csv", row.names=FALSE, sep=";")


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


# Dependency check with d
depdf <- as.data.frame(dependency(df$d,df$vd))
write.table(df, file = "morgan_complete.csv", row.names=FALSE, sep=";")

df <- read.table("morgan_subset.csv", header=T, sep=";")
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
write.table(df, file = "morgan_subset.csv", row.names=FALSE, sep=";")


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


# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "munder_complete.csv", row.names=FALSE, sep=";")

df <- read.table("munder_subset.csv", header=T, sep=";")
dfs <- df$n - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
write.table(df, file = "munder_subset.csv", row.names=FALSE, sep=";")


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


# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "piet_complete.csv", row.names=FALSE, sep=";")

df <- read.table("piet_subset.csv", header=T, sep=";")
dfs <- df$n - 2
J <- 1 - (3 / (4 * dfs - 1))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
write.table(df, file = "piet_subset.csv", row.names=FALSE, sep=";")


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


# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "smith_complete.csv", row.names=FALSE, sep=";")

df <- read.table("smith_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "smith_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "tillman_complete.csv", row.names=FALSE, sep=";")

df <- read.table("tillman_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "tillman_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "toosi_complete.csv", row.names=FALSE, sep=";")

df <- read.table("toosi_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "toosi_subset.csv", row.names=FALSE, sep=";")

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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "vaniddekinge_complete.csv", row.names=FALSE, sep=";")

df <- read.table("vaniddekinge_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "vaniddekinge_subset.csv", row.names=FALSE, sep=";")


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

# Dependency check with d
depdf <- as.data.frame(dependency(df$d,df$vd))
write.table(df, file = "webb_complete.csv", row.names=FALSE, sep=";")

df <- read.table("webb_subset.csv", header=T, sep=";")
df$n <- df$n1 + df$n2
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
write.table(df, file = "webb_subset.csv", row.names=FALSE, sep=";")



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

# Dependency check with g
depdf <- as.data.frame(dependency(df$g,df$vg)) 
write.table(df, file = "woodin_complete.csv", row.names=FALSE, sep=";")

df <- read.table("woodin_subset.csv", header=T, sep=";")
J <- 1 - (3 / (4 * df$n - 9))
df$d <- df$g / J   
df$vd <- df$n / ((df$n / 2)^2 + (df$d^2 / (2 * df$n)))  # variance cohen's d
df$vg <- J^2 * df$vd 
write.table(df, file = "woodin_subset.csv", row.names=FALSE, sep=";")

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

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "woodley_complete.csv", row.names=FALSE, sep=";")

df <- read.table("woodley_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "woodley_subset.csv", row.names=FALSE, sep=";")


# Yoon --------------------------------------------------------------------
df <- read.table("yoon_complete.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
res <- rma(r, vr, data=df)                         # random effects meta-analysis
l1o <- leave1out(res)                              # leave1out analysis
q <- res$QE - l1o$Q                                # change in Q with l1o
outlier <- df[which(q >= 3.84), ]                 # subset of studies that are included as outliers
regular <- df[which(q < 3.84), ]                   # subset of studies that are included as non-outliers
outlier.total <- bind_rows(outlier.total,select(outlier, id, study, r))
regular.total <- bind_rows(regular.total,select(regular, id, study, r))

# Dependency check with r
depdf <- as.data.frame(dependency(df$r,df$vr))
write.table(df, file = "yoon_complete.csv", row.names=FALSE, sep=";")

df <- read.table("yoon_subset.csv", header=T, sep=";")
df$vr <- ((1 - (df$r^2))^2) / (df$n-1)               # variance correlation r
write.table(df, file = "yoon_subset.csv", row.names=FALSE, sep=";")

# dataset with all studies
regular.total[,"type"] <- "r"
outlier.total[,"type"] <- "o"
dftot <- rbind(regular.total,outlier.total)
dftot <- dftot[order(dftot$id),]

# save dataframes
setwd("C:/Users/s421506/tiu/research/effectsizes/codebooks")
write.table(outlier.total, file = "studies-outliers.csv", row.names=FALSE, sep=";")
write.table(regular.total, file = "studies-nonoutliers.csv", row.names=FALSE, sep=";")
write.table(dftot, file = "studies-all.csv", row.names=FALSE, sep=";")

# Further calculations ----------------------------------------------------
k.tot <- nrow(outlier.total)+nrow(regular.total);k.tot     # number of total studies 1951
nrow(outlier.total);nrow(regular.total)                    # 581 outliers 1370 regulars
k.out <- nrow(outlier.total)/k.tot*100; k.out              # perc of outlier studies of total k 29.8%
k.reg <- nrow(regular.total)/k.tot*100; k.reg              # perc of regular studies of total k 70.2%

