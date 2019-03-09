rm(list=ls())

# Scatterplot 500 effect sizes.
library(readxl);library(ggplot2);library(gridExtra)

# Study 2
dat <- read.table("C:/Users/s421506/tiu/research/effectsizes/codebooks/codebook-meta-analyses-final-complete.csv", header=T, sep = '')
number_ticks <- function(n) {function(limits) pretty(limits, n)}

## EFFECT SIZE ESIMATES
dat.g <- subset(dat,efftype=="g")
dat.d <- subset(dat,efftype=="d")
dat.r <- subset(dat,efftype=="r")
dat.z <- subset(dat,efftype=="z")

ggplot(dat.g, aes(eff.so, eff.sc, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced pooled MA effect size (Hedges' g)", breaks=number_ticks(6)) +
  scale_y_continuous("Original pooled MA effect size (Hedges' g)", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.6, 1.55), y = c(-0.6, 1.55)) +
  scale_shape_manual(values=rep(c(21:24), times=3)) +
  theme(legend.position="none") 

ggplot(dat.d, aes(eff.so, eff.sc, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced pooled MA effect size (Cohen's d)", breaks=number_ticks(6)) +
  scale_y_continuous("Original pooled MA effect size (Cohen's d)", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.1, 2.5), y = c(-0.1, 2.5)) +
  scale_shape_manual(values=rep(c(21:24), times=1)) +
  theme(legend.position="none") 

ggplot(dat.r, aes(eff.so, eff.sc, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced pooled MA effect size (Pearson r)", breaks=number_ticks(6)) +
  scale_y_continuous("Original pooled MA effect size (Pearson r)", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.05, 0.45), y = c(-0.05, 0.45)) +
  scale_shape_manual(values=rep(c(21:24), times=4)) +
  theme(legend.position="none") 

ggplot(dat.z, aes(eff.so, eff.sc, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced pooled MA effect size (Fisher's z)", breaks=number_ticks(6)) +
  scale_y_continuous("Original pooled MA effect size (Fisher's z)", breaks=number_ticks(6)) +
  expand_limits(x = c(0.2, 0.3), y = c(0.2, 0.3)) +
  scale_shape_manual(values=21) +
  theme(legend.position="none") 

## CONFIDENCE INTERVAL
ci.so <- dat.g$ciub.so-dat.g$cilb.so
ci.sc <- dat.g$ciub.sc-dat.g$cilb.sc

ggplot(dat.g, aes(ci.sc, ci.so, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced CI range (Hedges' g)", breaks=number_ticks(6)) +
  scale_y_continuous("Original CI range (Hedges' g)", breaks=number_ticks(6)) + 
  expand_limits(x = c(0, 0.9), y = c(0, 0.9)) +
  scale_shape_manual(values=rep(c(21:24),3)) +
  theme(legend.position="none") 

ci.so <- dat.d$ciub.so-dat.d$cilb.so
ci.sc <- dat.d$ciub.sc-dat.d$cilb.sc

ggplot(dat.d, aes(ci.sc, ci.so, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced CI range (Hedges' g)", breaks=number_ticks(6)) +
  scale_y_continuous("Original CI range (Hedges' g)", breaks=number_ticks(6)) + 
  expand_limits(x = c(0, 2), y = c(0, 2)) +
  scale_shape_manual(values=rep(c(21:24),1)) +
  theme(legend.position="none") 

ci.so <- dat.r$ciub.so-dat.r$cilb.so
ci.sc <- dat.r$ciub.sc-dat.r$cilb.sc

ggplot(dat.r, aes(ci.sc, ci.so, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced CI range (Pearson r)", breaks=number_ticks(6)) +
  scale_y_continuous("Original CI range (Pearson r)", breaks=number_ticks(6)) + 
  expand_limits(x = c(0,0.35), y = c(0,0.35)) +
  scale_shape_manual(values=rep(c(21:24),4)) +
  theme(legend.position="none") 

ci.so <- dat.z$ciub.so-dat.z$cilb.so
ci.sc <- dat.z$ciub.sc-dat.z$cilb.sc

ggplot(dat.z, aes(ci.sc, ci.so, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced CI range (Pearson r)", breaks=number_ticks(6)) +
  scale_y_continuous("Original CI range (Pearson r)", breaks=number_ticks(6)) + 
  expand_limits(x = c(0.1,0.2), y = c(0.1,0.2)) +
  scale_shape_manual(values=21) +
  theme(legend.position="none") 

## TAU2 ESIMATES
ggplot(dat.g, aes(tau2.sc, tau2.so, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced tau2 (Hedges' g)", breaks=number_ticks(6)) +
  scale_y_continuous("Original tau2 (Hedges' g)", breaks=number_ticks(6)) + 
  expand_limits(x = c(-0.05, 1), y = c(-0.05, 1)) +
  scale_shape_manual(values=rep(c(21:24), times=3)) +
  theme(legend.position="none") 

ggplot(dat.d, aes(tau2.sc, tau2.so, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced tau2 (Cohen's d)", breaks=number_ticks(6)) +
  scale_y_continuous("Original tau2 (Cohen's d)", breaks=number_ticks(6)) + 
  expand_limits(x = c(-0.1, 1), y = c(-0.1, 1)) +
  scale_shape_manual(values=rep(c(21:24), times=1)) +
  theme(legend.position="none") 

ggplot(dat.r, aes(tau2.sc, tau2.so, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced tau2 (Pearson r)", breaks=number_ticks(6)) +
  scale_y_continuous("Original tau2 (Pearson r)", breaks=number_ticks(6)) + 
  expand_limits(x = c(-0.005, 0.09), y = c(-0.005, 0.09)) +
  scale_shape_manual(values=rep(c(21:24), times=4)) +
  theme(legend.position="none") 

ggplot(dat.z, aes(tau2.sc, tau2.so, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced tau2 (Fisher's z)", breaks=number_ticks(2)) +
  scale_y_continuous("Original tau2 (Fisher's z)", breaks=number_ticks(2)) + 
  expand_limits(x = c(-0.005, 0.01), y = c(-0.005, 0.01)) +
  scale_shape_manual(values=rep(c(21:24), times=1)) +
  theme(legend.position="none") 








df <- read_excel("C:/Users/s421506/tiu/research/effectsizes/codebooks/codebook-primary-studies-final-complete.xlsx", 1)

# Study 1 scatterplots - SMD
df.smd <- df[df$efftype == "g" | df$efftype == "d",]
df.smd$discrepancy <- as.character(df.smd$info)
sum(df.smd$discrepancy == "0");sum(df.smd$discrepancy == "1");sum(df.smd$discrepancy == "2");sum(df.smd$discrepancy == "3")
df.smd$discrepancy[df.smd$discrepancy=="0"]<-a<-"No discrepancy, k = 107, 43%"
df.smd$discrepancy[df.smd$discrepancy=="2"]<-b<-"Not enough information, k = 40, 16%"
df.smd$discrepancy[df.smd$discrepancy=="1"]<-c<-"Different effect, k = 40, 16%"
df.smd$discrepancy[df.smd$discrepancy=="3"]<-d<-"Ambiguous effect, k =60, 24%"

df.smd$discrepancy <- factor(df.smd$discrepancy, levels = c(a,b,c,d))
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

# transform cohen's d to hedges g
g.reported <- g.reproduced <- c()

for (i in 1:nrow(df.smd)){
  
  if (df.smd$efftype[i] == "d" ) {
    
    J.reported <- (1 - 3 / (4 * df.smd$n[i] - 9)) 
    g.reported[i] <- J.reported * df.smd$effest.exp[i]
    
    J.reproduced <- (1 - 3 / (4 * df.smd$nnew[i] - 9)) 
    g.reproduced[i] <- J.reproduced * df.smd$effestnew.exp[i]
    
  } else {
    
    g.reported[i] <- df.smd$effest.exp[i]
    g.reproduced[i] <- df.smd$effestnew.exp[i]
    
  }
  
}

ggplot(df.smd, aes(g.reproduced, g.reported, col=discrepancy)) + 
  geom_point(alpha = 0.9, size = 3, stroke=0.2) +
  theme_dark() +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes (Hedges' g)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Hedges' g)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "") +
  expand_limits(x = c(-1.2, 1.5), y = c(-0.75, 1.5)) +
  theme(legend.text=element_text(size=11)) +
  theme(legend.position=c(.85, .1)) +
  theme(legend.title=element_blank()) 


# Study 1 scatterplots - COR
df.cor <- df[df$efftype == "r" | df$efftype == "z",]
df.cor$discrepancy <- as.character(df.cor$info)
sum(df.cor$discrepancy == "0");sum(df.cor$discrepancy == "1");sum(df.cor$discrepancy == "2");sum(df.cor$discrepancy == "3")
df.cor$discrepancy[df.cor$discrepancy=="0"]<-a<-"No discrepancy, k = 169, 67%"
df.cor$discrepancy[df.cor$discrepancy=="2"]<-b<-"Not enough information, k = 14, 6%"
df.cor$discrepancy[df.cor$discrepancy=="1"]<-c<-"Different effect, k = 34, 13%"
df.cor$discrepancy[df.cor$discrepancy=="3"]<-d<-"Ambiguous effect, k = 36, 14%"

df.cor$discrepancy <- factor(df.cor$discrepancy, levels = c(a,b,c,d))
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

# transform correlation r to fisher's z
z.reported <- z.reproduced <- c()

for (i in 1:nrow(df.cor)){
  
  if (df.cor$efftype[i] == "r" ) {
    
    z.reported[i] <- 0.5 * log((1 + df.cor$effest.exp[i]) / (1 - df.cor$effest.exp[i]))
    z.reproduced[i] <- 0.5 * log((1 + df.cor$effestnew.exp[i]) / (1 - df.cor$effestnew.exp[i]))
    
  } else {
    
    z.reported[i] <- df.cor$effest.exp[i]
    z.reproduced[i] <- df.cor$effestnew.exp[i]
    
  }
  
}

ggplot(df.cor, aes(z.reproduced, z.reported, col=discrepancy)) + 
  geom_point(alpha = 0.9, size = 3, stroke=0.2) +
  theme_dark() +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "") +
  expand_limits(x = c(-1, 1), y = c(-0.75, 1)) +
  theme(legend.text=element_text(size=11)) +
  theme(legend.position=c(.85, .1)) +
  theme(legend.title=element_blank()) 



# Figure 4. ---------------------------------------------------------------
df$discrepancy <- as.character(df$info)
table(df$discrepancy)

df$discrepancy[df$discrepancy=="0"]<-a<-"No discrepancy, k = 291, 58%"
df$discrepancy[df$discrepancy=="2"]<-b<-"Not enough information, k = 49, 10%"
df$discrepancy[df$discrepancy=="1"]<-c<-"Different effect, k = 62, 12%"
df$discrepancy[df$discrepancy=="3"]<-d<-"Ambiguous effect, k = 98, 20%"

df$discrepancy <- factor(df$discrepancy, levels = c(a,b,c,d))
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

#df$discrepancy <- factor(df$discrepancy, levels = c("No discrepancy, k = 280, 56%", "Not enough information, k = 49, 10%", "Different effect, k = 73, 15%", "Ambiguous effect, k = 98, 20%"))
#mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(df, aes(znew.exp, z.exp, col=discrepancy)) + 
  geom_point(alpha = 0.9, size = 3, stroke=0.2) +
  theme_dark() +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "") +
  expand_limits(x = c(-1.2, 1.5), y = c(-0.75, 1.5)) +
  theme(legend.text=element_text(size=11)) +
  theme(legend.position=c(.85, .1)) +
  theme(legend.title=element_blank()) 

ggplot(df, aes(znew.exp, z.exp, col=discrepancy)) + 
  geom_point(alpha = 0.9, size = 3, stroke=0.2) +
  theme_dark() +
  theme(legend.key=element_blank()) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "") +
  expand_limits(x = c(-1.2, 1.5), y = c(-0.75, 1.5)) +
  theme(legend.text=element_text(size=11)) +
  theme(legend.position=c(.85, .1)) +
  theme(legend.title=element_blank())  +
  theme(legend.background = element_rect(colour = 'black', fill = 'lightgrey', linetype='solid'))

ggplot(df, aes(znew.exp, z.exp, col=discrepancy)) + 
  geom_point(alpha = 0.9, size = 3, stroke=0.2) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_fill_manual(values=cbbPalette) +
  expand_limits(x = c(-1.2, 1.5), y = c(-0.75, 1.5)) +
  theme(legend.text=element_text(size=11)) +
  theme(legend.position=c(.85, .1)) +
  theme(legend.title=element_blank())




# Scatterplot 500 effect sizes.
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/s421506/tiu/research/effectsizes/codebooks/codebook-primary-studies-final-complete.xlsx", 1)


library(ggplot2)
# Figure 4. ---------------------------------------------------------------
df$discrepancy <- as.character(df$info)
table(df$discrepancy)

df$discrepancy[df$discrepancy=="0"]<-a<-"No discrepancy, k = 288, 58%"
df$discrepancy[df$discrepancy=="2"]<-b<-"Not enough information, k = 49, 10%"
df$discrepancy[df$discrepancy=="1"]<-c<-"Different effect, k = 65, 13%"
df$discrepancy[df$discrepancy=="3"]<-d<-"Ambiguous effect, k = 98, 20%"

df$discrepancy <- factor(df$discrepancy, levels = c(a,b,c,d))
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

#df$discrepancy <- factor(df$discrepancy, levels = c("No discrepancy, k = 280, 56%", "Not enough information, k = 49, 10%", "Different effect, k = 73, 15%", "Ambiguous effect, k = 98, 20%"))
#mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(df, aes(znew.exp, z.exp, col=discrepancy)) + 
  geom_point(shape = 19, alpha = 0.9, size = 3, stroke=0.2) +
  geom_abline() + 
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  expand_limits(x = c(-1.2, 1.5), y = c(-0.75, 1.5)) +
  theme(legend.text=element_text(size=11)) +
  theme(legend.position=c(.85, .1)) +
  scale_fill_manual(values=cbbPalette, name = "") 
scale_color_manual(values=mycolors, name = "") 





# APPENDIX ??

# HEDGES G
df.g <- subset(df,efftype=="g")
df.g$discrepancy <- as.character(df.g$info)

df.g$discrepancy[df.g$discrepancy=="0"]<-"No discrepancy, k = 50, 47%"
df.g$discrepancy[df.g$discrepancy=="2"]<-"Not enough information, k = 11, 10%"
df.g$discrepancy[df.g$discrepancy=="1"]<-"Different effect, k = 23, 22%"
df.g$discrepancy[df.g$discrepancy=="3"]<-"Ambiguous effect, k = 22, 21%"

df.g$discrepancy <- factor(df.g$discrepancy, levels = c("No discrepancy, k = 50, 47%", "Not enough information, k = 11, 10%", "Different effect, k = 23, 22%", "Ambiguous effect, k = 22, 21%"))

mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(df.g, aes(effestnew.exp, effest.exp, col=discrepancy)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "Effect size: Hedges' g") + 
  expand_limits(x = c(-1.2, 2), y = c(-1, 2)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position=c(.80, .15))

# COHENS D
df.d <- subset(df,efftype=="d")
df.d$discrepancy <- as.character(df.d$info)

df.d$discrepancy[df.d$discrepancy=="0"]<-"No discrepancy, k = 52, 37%"
df.d$discrepancy[df.d$discrepancy=="2"]<-"Not enough information, k = 24, 17%"
df.d$discrepancy[df.d$discrepancy=="1"]<-"Different effect, k = 25, 18%"
df.d$discrepancy[df.d$discrepancy=="3"]<-"Ambiguous effect, k = 40, 28%"

df.d$discrepancy <- factor(df.d$discrepancy, levels = c("No discrepancy, k = 52, 37%", "Not enough information, k = 24, 17%", "Different effect, k = 25, 18%", "Ambiguous effect, k = 40, 28%"))

mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(df.d, aes(effestnew.exp, effest.exp, col=discrepancy)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "Effect size: Cohen's d") + 
  expand_limits(x = c(-1.2, 2), y = c(-1, 2)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position=c(.80, .15))

# CORRELATION R
df.r <- subset(df,efftype=="r")
df.r$discrepancy <- as.character(df.r$info)

df.r$discrepancy[df.r$discrepancy=="0"]<-"No discrepancy, k = 170, 69%"
df.r$discrepancy[df.r$discrepancy=="2"]<-"Not enough information, k = 14, 6%"
df.r$discrepancy[df.r$discrepancy=="1"]<-"Different effect, k = 25, 10%"
df.r$discrepancy[df.r$discrepancy=="3"]<-"Ambiguous effect, k = 36, 15%"

df.r$discrepancy <- factor(df.r$discrepancy, levels = c("No discrepancy, k = 170, 69%", "Not enough information, k = 14, 6%", "Different effect, k = 25, 10%", "Ambiguous effect, k = 36, 15%"))

mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(df.r, aes(effestnew.exp, effest.exp, col=discrepancy)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "Effect size: correlation r") + 
  expand_limits(x = c(-0.7, 1.2), y = c(-0.8, 1.2)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.position=c(.80, .15))

# CORRELATION Z
df.z <- subset(df,efftype=="z")
df.z$discrepancy <- as.character(df.z$info)

df.z$discrepancy[df.z$discrepancy=="0"]<-"No discrepancy, k = 8, 100%"
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(df.z, aes(effestnew.exp, effest.exp, col=discrepancy)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "Effect size: correlation z") + 
  expand_limits(x = c(-0.3, 0.8), y = c(-0.3, 0.8)) +
  theme(legend.text=element_text(size=8)) +
  theme(legend.position=c(.80, .15))

```





# Figure 5. ---------------------------------------------------------------
# HEDGES G
df.g <- subset(df,efftype=="g")
df.g$discrepancy <- as.character(df.g$info)
table(df.g$discrepancy)

df.g$discrepancy[df.g$discrepancy=="0"]<-"No discrepancy, k = 50, 47%"
df.g$discrepancy[df.g$discrepancy=="2"]<-"Not enough information, k = 11, 10%"
df.g$discrepancy[df.g$discrepancy=="1"]<-"Different effect, k = 23, 22%"
df.g$discrepancy[df.g$discrepancy=="3"]<-"Ambiguous effect, k = 22, 21%"

df.g$discrepancy <- factor(df.g$discrepancy, levels = c("No discrepancy, k = 50, 47%", "Not enough information, k = 11, 10%", "Different effect, k = 23, 22%", "Ambiguous effect, k = 22, 21%"))

mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(df.g, aes(effestnew.exp, effest.exp, col=discrepancy)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "") + 
  expand_limits(x = c(-1.2, 2), y = c(-1, 2)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position=c(.87, .15))

# COHENS D
df.d <- subset(df,efftype=="d")
df.d$discrepancy <- as.character(df.d$info)
table(df.d$discrepancy)

df.d$discrepancy[df.d$discrepancy=="0"]<-"No discrepancy, k = 52, 37%"
df.d$discrepancy[df.d$discrepancy=="2"]<-"Not enough information, k = 24, 17%"
df.d$discrepancy[df.d$discrepancy=="1"]<-"Different effect, k = 25, 18%"
df.d$discrepancy[df.d$discrepancy=="3"]<-"Ambiguous effect, k = 40, 28%"

df.d$discrepancy <- factor(df.d$discrepancy, levels = c("No discrepancy, k = 52, 37%", "Not enough information, k = 24, 17%", "Different effect, k = 25, 18%", "Ambiguous effect, k = 40, 28%"))

mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(df.d, aes(effestnew.exp, effest.exp, col=discrepancy)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "") + 
  expand_limits(x = c(-1.2, 2), y = c(-1, 2)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position=c(.87, .15))


# CORRELATION R
df.r <- subset(df,efftype=="r")
df.r$discrepancy <- as.character(df.r$info)
table(df.r$discrepancy)

df.r$discrepancy[df.r$discrepancy=="0"]<-"No discrepancy, k = 170, 69%"
df.r$discrepancy[df.r$discrepancy=="2"]<-"Not enough information, k = 14, 6%"
df.r$discrepancy[df.r$discrepancy=="1"]<-"Different effect, k = 25, 10%"
df.r$discrepancy[df.r$discrepancy=="3"]<-"Ambiguous effect, k = 36, 15%"

df.r$discrepancy <- factor(df.r$discrepancy, levels = c("No discrepancy, k = 170, 69%", "Not enough information, k = 14, 6%", "Different effect, k = 25, 10%", "Ambiguous effect, k = 36, 15%"))

mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(df.r, aes(effestnew.exp, effest.exp, col=discrepancy)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "") + 
  expand_limits(x = c(-0.7, 1.2), y = c(-0.8, 1.2)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position=c(.87, .15))


# CORRELATION Z
df.z <- subset(df,efftype=="z")
df.z$discrepancy <- as.character(df.z$info)
table(df.z$discrepancy)

df.z$discrepancy[df.z$discrepancy=="0"]<-"No discrepancy, k = 8, 100%"

df.z$discrepancy <- factor(df.z$discrepancy, levels = c("No discrepancy, k = 8, 100"))

mycolors <- c("#FFFFFF")

ggplot(df.z, aes(effestnew.exp, effest.exp, col=discrepancy)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors, name = "") + 
  expand_limits(x = c(-0.3, 0.8), y = c(-0.3, 0.8)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.position=c(.87, .15))



# MA Barplot --------------------------------------------------------------




# Bar graph
library(ggplot2); library(plyr); library(readxl); library(reshape2); library(dplyr);library(scales)
df <- read_excel("C:/Users/s421506/tiu/research/effectsizes/codebooks/codebook-primary-studies-final-complete.xlsx", 1)
df$discrepancy <- as.character(df$info)

metaauthors <- levels(factor(df$meta))

# discrepancy variable counted per meta-analysis
zero <- df %>% group_by(meta) %>% summarise(discrepancy = sum(discrepancy == 0))
one <- df %>% group_by(meta) %>% summarise(discrepancy = sum(discrepancy == 1))
two <- df %>% group_by(meta) %>% summarise(discrepancy = sum(discrepancy == 2))
three <- df %>% group_by(meta) %>% summarise(discrepancy = sum(discrepancy == 3))
df <- cbind(metaauthors,three[,2],one[,2],two[,2],zero[,2])
colnames(df) <- c("MA","Ambiguous effect","Different effect","Not enough information","No discrepancy")

dfma <- melt(df, id.var="MA")
dfma$MA <- factor(dfma$MA, levels=rev(sort(dfma$MA)))
dfma$MAnum <- as.character(as.numeric(dfma$MA))
dfma$variable <- factor(dfma$variable)
colnames(dfma)[colnames(dfma)=="value"] <- "Frequency"
mycolors <- c("#FF0000", "#FFa500", "#FFFF00", "#7cfa57")

dfma %>%
  arrange(MAnum) %>%
  mutate(MAnum=factor(MAnum, levels=33:1)) %>% 
  ggplot(aes(x = MAnum, y = Frequency, fill = variable)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(legend.title=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank()) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.position="bottom") +
  scale_fill_manual(values=mycolors, guide = guide_legend(reverse=TRUE))














ggplot(dat, aes(z.sc, z.so)) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(6)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(6)) + 
  expand_limits(x = c(-0.2, 0.8), y = c(-0.2, 0.8)) +
  scale_shape_manual(values=rep(c(21:25), times=33)) +
  theme(legend.position="none") 

## CONFIDENCE INTERVAL
ci.co <- dat$ci.ul.co-dat$ci.li.co
ci.cc <- dat$ci.ul.cc-dat$ci.li.cc
ci.so <- dat$ci.ul.so-dat$ci.li.so
ci.sc <- dat$ci.ul.sc-dat$ci.li.sc

ggplot(dat, aes(ci.cc, ci.co, shape = factor(author), fill = factor(author), col = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced tau2", breaks=number_ticks(6)) +
  scale_y_continuous("Original tau2", breaks=number_ticks(6)) + 
  expand_limits(x = c(0, 0.2), y = c(0, 0.2)) +
  scale_shape_manual(values=rep(c(21), times=33)) +
  theme(legend.position="none") 

ggplot(dat, aes(ci.sc, ci.so, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced CI range", breaks=number_ticks(6)) +
  scale_y_continuous("Original CI range", breaks=number_ticks(6)) + 
  expand_limits(x = c(0, 0.2), y = c(0, 0.2)) +
  scale_shape_manual(values=rep(c(21:25), times=33)) +
  theme(legend.position="none") 






n <- 33
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)

