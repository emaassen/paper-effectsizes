library(psych);library(irr);library(readxl)
dfe <- read_excel("C:/Users/s421506/tiu/research/effectsizes/codebooks/codebook-primary-studies-em.xlsx", 2)
dfa <- read_excel("C:/Users/s421506/tiu/research/effectsizes/codebooks/codebook-primary-studies-aoc.xlsx", 2)

ze <- dfe$znew
za <- dfa$znew
dfz <- cbind(ze,za)

table(ze == za)
cor(ze,za)                           # 0.97
cohen.kappa(x=cbind(ze,za))          # unweighted 0.74, weighted 0.97
icc(dfz, "twoway")                   # 0.964 (CI 0.97 ; 0.975)

# by category
infoe <- dfe$info
infoa <- dfa$info
dfinfo <- cbind(infoe,infoa)

tab <- table(infoe,infoa)
# simple agreement among raters
pr.a <- (tab[1,1] + tab[2,2] + tab[3,3] + tab[4,4])  / sum(tab)
# expected agreement among raters by chance
cat0 <- (sum(tab[1,]) / sum(tab)) * (sum(tab[,1]) / sum(tab)) * sum(tab)
cat1 <- (sum(tab[2,]) / sum(tab)) * (sum(tab[,2]) / sum(tab)) * sum(tab)
cat2 <- (sum(tab[3,]) / sum(tab)) * (sum(tab[,3]) / sum(tab)) * sum(tab)
cat3 <- (sum(tab[4,]) / sum(tab)) * (sum(tab[,4]) / sum(tab)) * sum(tab)
pr.e <- (cat0 + cat1 + cat2 + cat3) / sum(tab)

kappa <- (pr.a - pr.e) / (1 - pr.e)
