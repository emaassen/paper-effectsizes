rm(list = ls())
setwd("C:/Users/s421506/tiu/research/effectsizes/codebooks/")
packages <- c("readxl","metafor","ggplot2")
options(scipen=999)
#lapply(packages,install.packages(packages),character.only=T)     # if packages are not yet installed
lapply(packages,library,character.only=T)


# Study 2 Figures ---------------------------------------------------------
eff.so <- eff.sc <- cilb.so <- cilb.sc <- ciub.so <- ciub.sc <- tau2.so <- tau2.sc <- c() # empty vectors to store results
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# First, we need to get the estimates from all meta-analyses that have Hedges' g and Fisher's z and do not need to be transformed
datm <- read.table("codebook-meta-analyses-final-complete.csv", header=T, sep = '')
setwd("C:/Users/s421506/tiu/research/effectsizes/data-per-ma/")

for (i in 1:nrow(datm)) {
  if (datm$efftype[i] == "g" | datm$efftype[i] == "z") {
    
    eff.so[i] <- datm$eff.so[i]
    eff.sc[i] <- datm$eff.sc[i]
    cilb.so[i] <- datm$cilb.so[i]
    cilb.sc[i] <- datm$cilb.sc[i]
    ciub.so[i] <- datm$ciub.so[i]
    ciub.sc[i] <- datm$ciub.sc[i]
    tau2.so[i] <- datm$tau2.so[i]
    tau2.sc[i] <- datm$tau2.sc[i]

  }
}

# Then, we need to get the estimates from all meta-analyses that have Cohen's d or Correlation r, and transform them first

# Cohen's d to Hedges' g

# Else-quest
df <- read.table("elsequest_subset.csv", header=T, sep=';') 
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$g.o <- df$effest.exp * J.o                              # cohen's d
df$g.c <- df$effestnew.exp * J.c                           # cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 
write.table(df, file = "elsequest_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(g.o, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(g.c, vg.c, data=df, method="DL")               

eff.so[11] <- res.so$b       # MA subset original effect size estimate
cilb.so[11] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[11] <- res.so$ci.ub # upperbound
tau2.so[11] <- res.so$tau2  # tau2 estimate

eff.sc[11] <- res.sc$b 
cilb.sc[11] <- res.sc$ci.lb
ciub.sc[11] <- res.sc$ci.ub
tau2.sc[11] <- res.sc$tau2

# Green
df <- read.table("green_subset.csv", header=T, sep=';') 
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$g.o <- df$effest.exp * J.o                              # cohen's d
df$g.c <- df$effestnew.exp * J.c                           # cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 
write.table(df, file = "green_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(g.o, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(g.c, vg.c, data=df, method="DL") 

eff.so[16] <- res.so$b       # MA subset original effect size estimate
cilb.so[16] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[16] <- res.so$ci.ub # upperbound
tau2.so[16] <- res.so$tau2  # tau2 estimate

eff.sc[16] <- res.sc$b 
cilb.sc[16] <- res.sc$ci.lb
ciub.sc[16] <- res.sc$ci.ub
tau2.sc[16] <- res.sc$tau2

# Morgan
df <- read.table("morgan_subset.csv", header=T, sep=';') 
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$g.o <- df$effest.exp * J.o                              # cohen's d
df$g.c <- df$effestnew.exp * J.c                           # cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 
write.table(df, file = "morgan_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(g.o, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(g.c, vg.c, data=df, method="DL") 

eff.so[23] <- res.so$b       # MA subset original effect size estimate
cilb.so[23] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[23] <- res.so$ci.ub # upperbound
tau2.so[23] <- res.so$tau2  # tau2 estimate

eff.sc[23] <- res.sc$b 
cilb.sc[23] <- res.sc$ci.lb
ciub.sc[23] <- res.sc$ci.ub
tau2.sc[23] <- res.sc$tau2

# Webb
df <- read.table("webb_subset.csv", header=T, sep=';') 
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$g.o <- df$effest.exp * J.o                              # cohen's d
df$g.c <- df$effestnew.exp * J.c                           # cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 
write.table(df, file = "webb_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(g.o, vg.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(g.c, vg.c, data=df, method="DL") 

eff.so[30] <- res.so$b       # MA subset original effect size estimate
cilb.so[30] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[30] <- res.so$ci.ub # upperbound
tau2.so[30] <- res.so$tau2  # tau2 estimate

eff.sc[30] <- res.sc$b 
cilb.sc[30] <- res.sc$ci.lb
ciub.sc[30] <- res.sc$ci.ub
tau2.sc[30] <- res.sc$tau2


# Correlation r to Fisher's z

# Berry 1
df <- read.table("berry1_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "berry1_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="HS")  

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="HS")               

eff.so[6] <- res.so$b # MA subset original effect size estimate
cilb.so[6] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[6] <- res.so$ci.ub # upperbound
tau2.so[6] <- res.so$tau2 # tau2 estimate

eff.sc[6] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[6] <- res.sc$ci.lb
ciub.sc[6] <- res.sc$ci.ub
tau2.sc[6] <- res.sc$tau2

# Berry 2
df <- read.table("berry2_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "berry2_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="HS")  

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="HS")               

eff.so[7] <- res.so$b # MA subset original effect size estimate
cilb.so[7] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[7] <- res.so$ci.ub # upperbound
tau2.so[7] <- res.so$tau2 # tau2 estimate

eff.sc[7] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[7] <- res.sc$ci.lb
ciub.sc[7] <- res.sc$ci.ub
tau2.sc[7] <- res.sc$tau2

# Card
df <- read.table("card_subset.csv", header=T, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="DL")               

eff.so[8] <- res.so$b # MA subset original effect size estimate
cilb.so[8] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[8] <- res.so$ci.ub # upperbound
tau2.so[8] <- res.so$tau2 # tau2 estimate

eff.sc[8] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[8] <- res.sc$ci.lb
ciub.sc[8] <- res.sc$ci.ub
tau2.sc[8] <- res.sc$tau2

# Crook
df <- read.table("crook_subset.csv", header=T, sep=";")
# effest size
rbar.so <- sum(df$n * df$effest.exp) / sum(df$n)
rbar.sc <- sum(df$nnew * df$effestnew.exp) / sum(df$nnew)

zbar.so <- 0.5 * log((1 + rbar.so) / (1 - rbar.so))
zbar.sc <- 0.5 * log((1 + rbar.sc) / (1 - rbar.sc))

# confidence interval
vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
vz.c <- 1 / (df$nnew - 3)                              

r.so <- df$effest.exp
r.sc <- df$effestnew.exp
z.so <- 0.5 * log((1 + r.so) / (1 - r.so))
z.sc <- 0.5 * log((1 + r.sc) / (1 - r.sc))
n.so <- df$n
n.sc <- df$nnew
var.z.o <- sum(n.so * (z.so - zbar.so)^2) / sum(df$n)
var.z.c <- sum(n.sc * (z.sc - zbar.sc)^2) / sum(df$nnew)
var.e.o <- sum(n.so * (1 - zbar.so^2)^2 / (n.so - 1)) / sum(df$n)  
var.e.c <- sum(n.sc * (1 - zbar.sc^2)^2 / (n.sc - 1)) / sum(df$nnew) 
res.var.o <- var.z.o - var.e.o
res.var.c <- var.z.c - var.e.c
k <- nrow(df)
se.o <- ((1 - zbar.so^2)^2 / (sum(df$n) - k)) + (res.var.o / k)^1/2
se.c <- ((1 - zbar.sc^2)^2 / (sum(df$nnew) - k)) + (res.var.c / k)^1/2
ci.lb.so <- zbar.so - (1.96 * se.o)
ci.ub.so <- zbar.so + (1.96 * se.o)
ci.lb.sc <- zbar.sc - (1.96 * se.c)
ci.ub.sc <- zbar.sc + (1.96 * se.c)

#tau2
df$vi.so <- ((1 - (z.so^2))^2) / (df$n - 1) 
df$vi.sc <- ((1 - (z.sc^2))^2) / (df$nnew - 1) 
wi.so <- 1 / df$vi.so
wi.sc <- 1 / df$vi.sc
yi.so <- z.so
yi.sc <- z.sc
theta.so <- sum(wi.so * yi.so) / sum(wi.so)
theta.sc <- sum(wi.sc * yi.sc) / sum(wi.sc)
tausq.so <- (sum(wi.so * (yi.so - theta.so)^2)/sum(wi.so)) - (sum(wi.so * df$vi.so) / sum(wi.so))
tausq.sc <- (sum(wi.sc * (yi.sc - theta.sc)^2)/sum(wi.sc)) - (sum(wi.sc * df$vi.sc) / sum(wi.sc))

eff.so[9] <- zbar.so # MA subset original effect size estimate
cilb.so[9] <- ci.lb.so # MA subset original effect size CI lowerbound
ciub.so[9] <- ci.ub.so # upperbound
tau2.so[9] <- tausq.so # tau2 estimate

eff.sc[9] <- zbar.sc # MA subset checked effect size estimate
cilb.sc[9] <- ci.lb.sc
ciub.sc[9] <- ci.ub.sc
tau2.sc[9] <- tausq.sc

# De Wit
df <- read.table("dewit_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "dewit_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="HS")  

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="HS")               

eff.so[10] <- res.so$b # MA subset original effect size estimate
cilb.so[10] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[10] <- res.so$ci.ub # upperbound
tau2.so[10] <- res.so$tau2 # tau2 estimate

eff.sc[10] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[10] <- res.sc$ci.lb
ciub.sc[10] <- res.sc$ci.ub
tau2.sc[10] <- res.sc$tau2

# Farber
df <- read.table("farber_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "farber_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="DL")  

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="DL")               

eff.so[12] <- res.so$b # MA subset original effect size estimate
cilb.so[12] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[12] <- res.so$ci.ub # upperbound
tau2.so[12] <- res.so$tau2 # tau2 estimate

eff.sc[12] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[12] <- res.sc$ci.lb
ciub.sc[12] <- res.sc$ci.ub
tau2.sc[12] <- res.sc$tau2

# Fox
df <- read.table("card_subset.csv", header=T, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="ML")  

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="ML")               

eff.so[14] <- res.so$b # MA subset original effect size estimate
cilb.so[14] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[14] <- res.so$ci.ub # upperbound
tau2.so[14] <- res.so$tau2 # tau2 estimate

eff.sc[14] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[14] <- res.sc$ci.lb
ciub.sc[14] <- res.sc$ci.ub
tau2.sc[14] <- res.sc$tau2

# Freund
df <- read.table("freund_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "freund_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma.mv(z.o, vz.o, data=df)   

# subset reproduced
res.sc <- rma.mv(z.c, vz.c, data=df)               

eff.so[15] <- res.so$b # MA subset original effect size estimate
cilb.so[15] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[15] <- res.so$ci.ub # upperbound
tau2.so[15] <- res.so$tau2 # tau2 estimate

eff.sc[15] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[15] <- res.sc$ci.lb
ciub.sc[15] <- res.sc$ci.ub
tau2.sc[15] <- res.sc$tau2

# Kolden
df <- read.table("kolden_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "kolden_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="HE")   

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="HE")               

eff.so[20] <- res.so$b # MA subset original effect size estimate
cilb.so[20] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[20] <- res.so$ci.ub # upperbound
tau2.so[20] <- res.so$tau2 # tau2 estimate

eff.sc[20] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[20] <- res.sc$ci.lb
ciub.sc[20] <- res.sc$ci.ub
tau2.sc[20] <- res.sc$tau2

# Lucassen
df <- read.table("lucassen_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "lucassen_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="DL")   

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="DL")               

eff.so[21] <- res.so$b # MA subset original effect size estimate
cilb.so[21] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[21] <- res.so$ci.ub # upperbound
tau2.so[21] <- res.so$tau2 # tau2 estimate

eff.sc[21] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[21] <- res.sc$ci.lb
ciub.sc[21] <- res.sc$ci.ub
tau2.sc[21] <- res.sc$tau2

# Smith
df <- read.table("smith_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "smith_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="DL")   

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="DL")               

eff.so[26] <- res.so$b # MA subset original effect size estimate
cilb.so[26] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[26] <- res.so$ci.ub # upperbound
tau2.so[26] <- res.so$tau2 # tau2 estimate

eff.sc[26] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[26] <- res.sc$ci.lb
ciub.sc[26] <- res.sc$ci.ub
tau2.sc[26] <- res.sc$tau2

# Tillman
df <- read.table("tillman_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                            
write.table(df, file = "tillman_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="DL")   

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="DL")               

eff.so[27] <- res.so$b # MA subset original effect size estimate
cilb.so[27] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[27] <- res.so$ci.ub # upperbound
tau2.so[27] <- res.so$tau2 # tau2 estimate

eff.sc[27] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[27] <- res.sc$ci.lb
ciub.sc[27] <- res.sc$ci.ub
tau2.sc[27] <- res.sc$tau2

# Toosi
df <- read.table("toosi_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                            
write.table(df, file = "toosi_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="DL")   

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="DL")               

eff.so[28] <- res.so$b # MA subset original effect size estimate
cilb.so[28] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[28] <- res.so$ci.ub # upperbound
tau2.so[28] <- res.so$tau2 # tau2 estimate

eff.sc[28] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[28] <- res.sc$ci.lb
ciub.sc[28] <- res.sc$ci.ub
tau2.sc[28] <- res.sc$tau2

# Van Iddekinge
df <- read.table("vaniddekinge_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                            
write.table(df, file = "vaniddekinge_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="HS")   

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="HS")               

eff.so[29] <- res.so$b # MA subset original effect size estimate
cilb.so[29] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[29] <- res.so$ci.ub # upperbound
tau2.so[29] <- res.so$tau2 # tau2 estimate

eff.sc[29] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[29] <- res.sc$ci.lb
ciub.sc[29] <- res.sc$ci.ub
tau2.sc[29] <- res.sc$tau2

# Woodley
df <- read.table("woodley_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                            
write.table(df, file = "woodley_subset.csv", row.names=FALSE, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="DL")   

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="DL")               

eff.so[32] <- res.so$b # MA subset original effect size estimate
cilb.so[32] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[32] <- res.so$ci.ub # upperbound
tau2.so[32] <- res.so$tau2 # tau2 estimate

eff.sc[32] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[32] <- res.sc$ci.lb
ciub.sc[32] <- res.sc$ci.ub
tau2.sc[32] <- res.sc$tau2

# Yoon
df <- read.table("yoon_subset.csv", header=T, sep=";")

# subset original
res.so <- rma(z.o, vz.o, data=df, method="DL")   

# subset reproduced
res.sc <- rma(z.c, vz.c, data=df, method="DL")               

eff.so[33] <- res.so$b # MA subset original effect size estimate
cilb.so[33] <- res.so$ci.lb # MA subset original effect size CI lowerbound
ciub.so[33] <- res.so$ci.ub # upperbound
tau2.so[33] <- res.so$tau2 # tau2 estimate

eff.sc[33] <- res.sc$b # MA subset checked effect size estimate
cilb.sc[33] <- res.sc$ci.lb
ciub.sc[33] <- res.sc$ci.ub
tau2.sc[33] <- res.sc$tau2


# Figure X: average meta-analytic effect  ---------------------------------
subg <- which(datm$efftype == "g" | datm$efftype == "d")
subz <- which(datm$efftype == "r" | datm$efftype == "z")
eff.so.g <- eff.so[subg]
eff.so.z <- eff.so[subz]
eff.sc.g <- eff.sc[subg]
eff.sc.z <- eff.sc[subz]

datplot <- subset(datm,efftype=="g" | efftype=="d")
ggplot(datplot, aes(eff.so.g, eff.sc.g, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=8, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_y_continuous("Reproduced pooled MA effect size (Hedges' g)", breaks=number_ticks(6)) +
  scale_x_continuous("Original pooled MA effect size (Hedges' g)", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.6, 1.55), y = c(-0.6, 1.55)) +
  scale_shape_manual(values=rep(c(21:25), times=4)) +
  theme(legend.position="none") +
  theme(axis.title=element_text(size=14)) 

datplot <- subset(datm,efftype=="r" | efftype=="z")
ggplot(datplot, aes(eff.so.z, eff.sc.z, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=8, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_y_continuous("Reproduced pooled MA effect size (Fisher's z)", breaks=number_ticks(6)) +
  scale_x_continuous("Original pooled MA effect size (Fisher's z)", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.1, 0.5), y = c(-0.1, 0.5)) +
  scale_shape_manual(values=rep(c(21:25), times=4)) +
  theme(legend.position="none")+
  theme(axis.title=element_text(size=14)) 

# Figure X: meta-analytic effect CI ---------------------------------------
ci.so <- ciub.so - cilb.so
ci.sc <- ciub.sc - cilb.sc
ci.so.g <- ci.so[subg]
ci.sc.g <- ci.sc[subg]
ci.so.z <- ci.so[subz]
ci.sc.z <- ci.sc[subz]

datplot <- subset(datm,efftype=="g" | efftype=="d")
ggplot(datplot, aes(ci.so.g, ci.sc.g, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=8, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_y_continuous("Reproduced pooled MA effect size confidence interval (Hedges' g)", breaks=number_ticks(6)) +
  scale_x_continuous("Original pooled MA effect size confidence interval (Hedges' g)", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.6, 1.55), y = c(-0.6, 1.55)) +
  scale_shape_manual(values=rep(c(21:25), times=4)) +
  theme(legend.position="none") +
  theme(axis.title=element_text(size=14)) 

datplot <- subset(datm,efftype=="r" | efftype=="z")
ggplot(datplot, aes(ci.so.z, ci.sc.z, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=8, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_y_continuous("Reproduced pooled MA effect size confidence interval (Fisher's z)", breaks=number_ticks(6)) +
  scale_x_continuous("Original pooled MA effect size confidence interval (Fisher's z)", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.1, 0.5), y = c(-0.1, 0.5)) +
  scale_shape_manual(values=rep(c(21:25), times=4)) +
  theme(legend.position="none")+
  theme(axis.title=element_text(size=14)) 


# Figure X: tau2 ----------------------------------------------------------
tau2.so.g <- tau2.so[subg]
tau2.so.z <- tau2.so[subz]
tau2.sc.g <- tau2.sc[subg]
tau2.sc.z <- tau2.sc[subz]

datplot <- subset(datm,efftype=="g" | efftype=="d")
ggplot(datplot, aes(tau2.so.g, tau2.sc.g, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=8, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_y_continuous("Reproduced MA tau2 estimate (Hedges' g)", breaks=number_ticks(6)) +
  scale_x_continuous("Original MA tau2 estimate (Hedges' g)", breaks=number_ticks(6)) +
  expand_limits(x = c(0, 1.55), y = c(0, 1.55)) +
  scale_shape_manual(values=rep(c(21:25), times=4)) +
  theme(legend.position="none") +
  theme(axis.title=element_text(size=14)) 

datplot <- subset(datm,efftype=="r" | efftype=="z")
ggplot(datplot, aes(tau2.so.z, tau2.sc.z, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=8, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_y_continuous("Reproduced pooled MA effect size confidence interval (Fisher's z)", breaks=number_ticks(6)) +
  scale_x_continuous("Original pooled MA effect size confidence interval (Fisher's z)", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.1, 0.5), y = c(-0.1, 0.5)) +
  scale_shape_manual(values=rep(c(21:25), times=4)) +
  theme(legend.position="none")+
  theme(axis.title=element_text(size=14)) 

