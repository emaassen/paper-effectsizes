rm(list = ls())
setwd("C:/Users/s421506/tiu/research/effectsizes/study-1/")
packages <- c("metafor","readxl","knitr","kableExtra")
options(scipen=999)
#lapply(packages,install.packages(packages),character.only=T)     # if packages are not yet installed
lapply(packages,library,character.only=T)

datp <- read.table("codebook-primary-studies-final.csv", header=T, sep = '') 
datret <- read.table("nonretrieved-primary-studies.csv", header=T, sep = ';')

# Table 1 -----------------------------------------------------------------

# Type of effect
attach(datp)
typeof.info <- table(efftype,info);typeof.info

smds.disc <- sum(typeof.info[1:2,2:4])
smds.nodisc <- sum(typeof.info[1:2,1])
cors.disc <- sum(typeof.info[3:4,2:4])
cors.nodisc <- sum(typeof.info[3:4,1])

# Outlier studies
out.info <- table(typestudy,info)

out.disc <- sum(out.info[1,2:4])
out.nodisc <- sum(out.info[1,1])
nout.disc <- sum(out.info[2,2:4])
nout.nodisc <- sum(out.info[2,1])

# Published studies
pub.info <- table(typecat,info)

pub.disc <- sum(pub.info[1,2:4])
pub.nodisc <- sum(pub.info[1,1])
npub.disc <- sum(pub.info[2,2:4])
npub.nodisc <- sum(pub.info[2,1])

# Table 1
table.1 <- data.frame(matrix(ncol = 6, nrow = 3))
col <- c("SMD", "correlation", "outlier", "non-outlier", "published", "non-published")
row <- c("Discrepancy","No discrepancy","Total")
colnames(table.1) <- col
rownames(table.1) <- row

table.1[1,1] <- smds.disc; table.1[2,1] <- smds.nodisc; table.1[3,1] <- sum(table.1[1,1]+table.1[2,1])
table.1[1,2] <- cors.disc; table.1[2,2] <- cors.nodisc; table.1[3,2] <- sum(table.1[1,2]+table.1[2,2])
table.1[1,3] <- out.disc; table.1[2,3] <- out.nodisc; table.1[3,3] <- sum(table.1[1,3]+table.1[2,3])
table.1[1,4] <- nout.disc; table.1[2,4] <- nout.nodisc; table.1[3,4] <- sum(table.1[1,4]+table.1[2,4])
table.1[1,5] <- pub.disc; table.1[2,5] <- pub.nodisc; table.1[3,5] <- sum(table.1[1,5]+table.1[2,5])
table.1[1,6] <- npub.disc; table.1[2,6] <- npub.nodisc; table.1[3,6] <- sum(table.1[1,6]+table.1[2,6])

kable(table.1) %>%
  kable_styling(bootstrap_options = c("striped", full_width = F)) %>%
  column_spec(3, border_right = T) %>% 
  column_spec(5, border_right = T) %>% 
  column_spec(6, border_right = F) 

