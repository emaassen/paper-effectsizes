library(lme4)       # library packages
library(statmod)
options(scipen=999)

datp <- read.table("../codebooks/codebook-primary-studies-complete.csv", header=T, sep = '')

# number of effect sizes in total per MA
k.tot <- c(57,108,13,270,21,21,283,11,67,80,93,18,105,94,154,16,43,20,48,16,16,8,15,29,9,184,30,37,17,31,40,12,12)

# number of effect sizes sampled per MA
k.sample <- c(20,20,11,21,11,11,23,9,22,21,20,11,20,20,20,13,11,14,20,10,11,8,9,16,8,21,19,20,10,20,13,6,11)

# number of studies reproducible per MA
k.rep <- c(12,5,10,14,0,1,21,5,17,17,5,1,15,13,20,5,3,7,8,1,7,8,1,5,4,14,14,6,9,7,6,4,11)

# number of studies irreproducible per MA
k.irrep <- k.sample - k.rep

# number of studies not sampled
k.nsample <- k.tot - k.sample
 
# number of effect sizes not sampled in three MAs that were fully reproducible
k.nsample[c(15,22,33)]

# probability of succes (= reproducible)
p.success <- k.rep / k.sample
success <- weighted.mean(p.success, k.sample)
success

### 'Stupid' model with same reproducibility probabilities for each effect size
p <- success             # probability an effect size is reproducible
k <- 1:max(k.nsample)    # the number of effect sizes in a meta-analysis with uncertain reproducibility
rep <- p^k               # probability that all k effect sizes are reproducible
cbind(k,rep)

# But is the reproducibility for individual effect sizes equal for all meta-analyses? No:
tab <- data.frame(k.rep,k.irrep)
chisq.test(tab) 

# multilevel logistic model
df <- data.frame(datp$id,datp$meta,datp$info)
colnames(df) <- c("id","meta","correct")
df$correct[df$correct == 0] <- "c" 
df$correct[df$correct != "c"] <- 0 
df$correct[df$correct == "c"] <- 1 
df$correct <- as.numeric(df$correct) # reproducible = 1 // irreproducible = 0

fit <- glmer(correct ~ 1 + (1 | meta), data=df, family=binomial)
summary(fit) 

### Less stupid model, assuming heterogeneity in reproducibility
int <- coef(summary(fit))[ , "Estimate"]    # intercept of multilevel logistic regression model
var <- VarCorr(fit)$meta[1]                 # variance of random effects

n <- 1000          # number of points in approximation
n.approx <- gauss.quad.prob(n,dist="normal",mu=int,sigma=sqrt(var)) 
# approximation of normal distribution by n discrete points
p <- 1/(1+exp(-n.approx$nodes))
# approximation of prob reproduced, based on the normal approximation

rep <- 1:max(k.nsample)

for (i in 1:max(k.nsample)) { # for each number of effect sizes, compute the average probability of prob reproduced
  rep[i] <- sum(p^i*n.approx$weights)
}

rep[c(k.nsample[c(15,22,33)])] # probability of succesfully reproducing remaining effect sizes from MA 15, 22, and 33 (no. 22 is excluded bc k.nsample = 0)
rep[c(1,5,10,15,16)] # probability of succesfully reproducing a set of 1, 5, 10, 15, 16 individual effect sizes

