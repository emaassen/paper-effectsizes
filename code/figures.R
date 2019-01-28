rm(list=ls())

# Scatterplot 500 effect sizes.
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/s421506/tiu/research/effectsizes/codebooks/codebook-primary-studies-final-complete.xlsx", 1)

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




#  scale_color_manual(values=mycolors, name = "") + 


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
df <- read_excel("C:/Users/s421506/tiu/research/effectsizes/xlsxfiles/codebook-primary-studies-final.xlsx", 2)
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
  theme(legend.text=element_text(size=14)) +
  theme(legend.position="bottom") +
  scale_fill_manual(values=mycolors, guide = guide_legend(reverse=TRUE))




# MA Scatterplot ----------------------------------------------------------

dat <- read.table("C:/Users/s421506/tiu/research/effectsizes/codebook-meta-analyses-finalbackup.csv", header=T, sep = '')


## EFFECT SIZE ESIMATES
library(randomcoloR)

ggplot(dat, aes(z.cc, z.co, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=4, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(6)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(6)) + 
  expand_limits(x = c(-0.2, 0.8), y = c(-0.2, 0.8)) +
  scale_shape_manual(values=rep(c(21:25), times=33)) +
  theme(legend.position="none") 


ggplot(dat, aes(z.sc, z.so, shape = factor(author), fill = factor(author))) + 
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




## TAU2 ESIMATES

ggplot(dat, aes(tau.cc, tau.co, shape = factor(author), fill = factor(author), col = factor(author))) + 
  geom_point(alpha = 0.6, size=4, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced tau2", breaks=number_ticks(6)) +
  scale_y_continuous("Original tau2", breaks=number_ticks(6)) + 
  expand_limits(x = c(0, 0.2), y = c(0, 0.2)) +
  scale_shape_manual(values=rep(c(21), times=33)) +
  theme(legend.position="none") 

ggplot(dat, aes(tau.sc, tau.so, shape = factor(author), fill = factor(author))) + 
  geom_point(alpha = 0.6, size=5, stroke=0.7) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  scale_x_continuous("Reproduced tau2", breaks=number_ticks(6)) +
  scale_y_continuous("Original tau2", breaks=number_ticks(6)) + 
  expand_limits(x = c(0, 0.2), y = c(0, 0.2)) +
  scale_shape_manual(values=rep(c(21:25), times=7)) +
  theme(legend.position="none") 
  


n <- 33
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)



# Appendices - this still needs to be fixed

# Adesope
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Adesope")

dfma$Discr <- factor(dfma$info)

levels(dfma$Discr) = c("No discrepancy, k = 16, 80%", "Different effect, k = 2, 10%", "Ambiguous effect, k = 2, 10%")

subheader <- paste("Adesope (k = 20)")
mycolors <- c("#FFFFFF", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.65, y = 0.8, label = " Original MA effect = 0.077 \n Reproduced MA effect = 0.088", hjust = 0) 


# Alfieri
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Alfieri")

dfma$infonew

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 11, 55%", "Not enough information, k = 1, 5%", "Different effect, k = 1, 5%", "Ambiguous effect, k = 7, 35%")

subheader <- paste("Alfieri (k = 20)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.6, y = 0.95, label = " Original MA effect = 0.237 \n Reproduced MA effect = 0.222", hjust = 0) 


# Babbage
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Babbage")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 10, 91%", "Not enough information, k = 1, 9%")

subheader <- paste("Babbage (k = 11)")
mycolors <- c("#FFFFFF", "#FFFF00")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.2, y = 0.9, label = " Original MA effect = 0.509 \n Reproduced MA effect = 0.509", hjust = 0) 


# Balliet 
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Balliet")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 18, 86%", "Not enough information, k = 1, 5%", "Ambiguous effect, k = 2, 10%")

subheader <- paste("Balliet (k = 21)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.5, y = 0.8, label = " Original MA effect = 0.040 \n Reproduced MA effect = 0.039", hjust = 0) 


# Benish
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Benish")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 2, 18%", "Not enough information, k = 3, 27%", "Ambiguous effect, k = 6, 54%")

subheader <- paste("Benish (k = 11)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.05, y = 0.44, label = " Original MA effect = 0.229 \n Reproduced MA effect = 0.252", hjust = 0) 


# Berry1
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Berry 1")

dfma$Discr <- factor(dfma$info)

levels(dfma$Discr) = c("No discrepancy, k = 10, 91%", "Different effect, k = 1, 9%")

subheader <- paste("Berry 1 (k = 11)")
mycolors <- c("#FFFFFF", "#FFa500")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.15, y = 0.52, label = " Original MA effect = 0.328 \n Reproduced MA effect = 0.323", hjust = 0) 


# Berry2
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Berry 2")

dfma$Discr <- factor(dfma$info)

levels(dfma$Discr) = c("No discrepancy, k = 23, 100%")

subheader <- paste("Berry 2 (k = 23)")
mycolors <- c("#FFFFFF")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.2, y = 0.75, label = " Original MA effect = 0.337 \n Reproduced MA effect = 0.337", hjust = 0) 


# Card
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Card")

dfma$Discr <- factor(dfma$info)

levels(dfma$Discr) = c("No discrepancy, k = 7, 78%", "Not enough information, k = 1, 11%", "Ambiguous effect, k = 1, 11%")

subheader <- paste("Card (k = 9)")
mycolors <- c("#FFFFFF", "#FFff00", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.05, y = 0.06, label = " Original MA effect = -0.003 \n Reproduced MA effect = 0.005", hjust = 0) 


# Crook
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Crook")

dfma$Discr <- factor(dfma$info)

levels(dfma$Discr) = c("No discrepancy, k = 21, 95%", "Different effect, k = 1, 5%")

subheader <- paste("Crook (k = 22)")
mycolors <- c("#FFFFFF", "#FFa500")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.05, y = 0.61, label = " Original MA effect = 0.200 \n Reproduced MA effect = 0.193", hjust = 0) 


# De wit
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="de Wit")

dfma$Discr <- factor(dfma$info)

levels(dfma$Discr) = c("No discrepancy, k = 20, 95%", "Ambiguous effect, k = 1, 5%")

subheader <- paste("de Wit (k = 21)")
mycolors <- c("#FFFFFF", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.2, y = 0.75, label = " Original MA effect = 0.231 \n Reproduced MA effect = 0.226", hjust = 0) 


# Else-quest
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Else-quest")

dfma$Discr <- factor(dfma$info)

levels(dfma$Discr) = c("No discrepancy, k = 6, 30%", "Not enough information, k = 14, 70%")

subheader <- paste("Else-Quest (k = 20)")
mycolors <- c("#FFFFFF", "#FFFF00")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.15, y = 0.23, label = " Original MA effect = 0.018 \n Reproduced MA effect = 0.019", hjust = 0) 


# Farber
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Farber")

dfma$Discr <- factor(dfma$info)

levels(dfma$Discr) = c("No discrepancy, k = 4, 36%", "Ambiguous effect, k = 7, 64%")

subheader <- paste("Farber (k = 11)")
mycolors <- c("#FFFFFF","#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.15, y = 0.95, label = " Original MA effect = 0.285 \n Reproduced MA effect = 0.303", hjust = 0) 


# Fischer
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Fischer")

dfma$Discr <- factor(dfma$info)

levels(dfma$Discr) = c("No discrepancy, k = 15, 75%", "Different effect, k = 4, 20%", "Ambiguous effect, k = 1, 5%")

subheader <- paste("Fischer (k = 20)")
mycolors <- c("#FFFFFF", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.39, y = 0.9, label = " Original MA effect = 0.233 \n Reproduced MA effect = 0.235", hjust = 0) 


# Fox
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Fox")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 16, 80%", "Not enough information, k = 1, 5%", "Different effect, k = 2, 10%", "Ambiguous effect, k = 1, 5%")

subheader <- paste("Fox (k = 20)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.65, y = 0.53, label = " Original MA effect = -0.084 \n Reproduced MA effect = 0.083", hjust = 0) 


# Freund
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Freund")

dfma$Discr <- factor(dfma$info)

levels(dfma$Discr) = c("No discrepancy, k = 20, 100%")

subheader <- paste("Freund (k = 20)")
mycolors <- c("#FFFFFF")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.1, y = 0.77, label = " Original MA effect = 0.359 \n Reproduced MA effect = 0.358", hjust = 0) 


# Green
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Green")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 7, 50%", "Not enough information, k = 5, 36%", "Different effect, k = 2, 14%")

subheader <- paste("Green (k = 14)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0, y = 1.7, label = " Original MA effect = 0.867 \n Reproduced MA effect = 0.865", hjust = 0) 


# Hallion
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Hallion")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 6, 55%", "Not enough information, k = 1, 9%", "Different effect, k = 1, 9%", "Ambiguous effect, k = 3, 27%")

subheader <- paste("Hallion (k = 11)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = TRUE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.982, 0.73), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.01, y = 0.6, label = " Original MA effect = 0.161 \n Reproduced MA effect = 0.015", hjust = 0) 


# Ihle
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Ihle")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 7, 50%", "Not enough information, k = 4, 29%", "Different effect, k = 2, 14%", "Ambiguous effect, k = 1, 7%")

subheader <- paste("Ihle (k = 14)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.57), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.1, y = 1.2, label = " Original MA effect = 0.020 \n Reproduced MA effect = 0.110", hjust = 0) 


# Koenig
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Koenig")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 9, 45%", "Not enough information, k = 3, 15%", "Different effect, k = 6, 30%", "Ambiguous effect, k = 2, 10%")

subheader <- paste("Koenig (k = 20)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.5, y = 1.1, label = " Original MA effect = 0.641 \n Reproduced MA effect = 0.688", hjust = 0) 


# Kolden
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Kolden")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 1, 10%", "Not enough information, k = 3, 30%", "Different effect, k = 1, 10%", "Ambiguous effect, k = 5, 50%")

subheader <- paste("Kolden (k = 10)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.25, y = 0.83, label = " Original MA effect = 0.252 \n Reproduced MA effect = 0.220", hjust = 0) 


#Lucassen
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Lucassen")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 9, 82%", "Different effect, k = 1, 9%", "Ambiguous effect, k = 1, 9%")

subheader <- paste("Lucassen (k = 11)")
mycolors <- c("#FFFFFF", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0, y = 0.45, label = " Original MA effect = 0.109 \n Reproduced MA effect = 0.101", hjust = 0) 


#Mol
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Mol")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 7, 88%", "Different effect, k = 1, 12%")

subheader <- paste("Mol (k = 8)")
mycolors <- c("#FFFFFF", "#FFa500")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.05, y = 0.43, label = " Original MA effect = 0.271 \n Reproduced MA effect = 0.264", hjust = 0) 


# Morgan
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Morgan")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 4, 44%", "Not enough information, k = 1, 11%", "Different effect, k = 2, 22%", "Ambiguous effect, k = 2, 22%")

subheader <- paste("Morgan (k = 9)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.2, y = 0.78, label = " Original MA effect = 0.428 \n Reproduced MA effect = 0.420", hjust = 0) 


# Munder
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Munder")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 13, 81%", "Not enough information, k = 1, 6%", "Different effect, k = 2, 13%", "Ambiguous effect, k = 2, 22%")

subheader <- paste("Munder (k = 16)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.5, y = 0.55, label = " Original MA effect = -0.009 \n Reproduced MA effect = -0.007", hjust = 0) 


# Piet
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Piet")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 5, 63%", "Not enough information, k = 2, 25%", "Different effect, k = 1, 13%")

subheader <- paste("Piet (k = 8)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.1, y = 0.31, label = " Original MA effect = 0.182 \n Reproduced MA effect = 0.186", hjust = 0) 


# Smith
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Smith")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 15, 71%", "Different effect, k = 1, 5%", "Ambiguous effect, k = 5, 24%")

subheader <- paste("Smith (k = 21)")
mycolors <- c("#FFFFFF", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.02, y = 0.51, label = " Original MA effect = 0.250 \n Reproduced MA effect = 0.248", hjust = 0) 


# Tillman
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Tilman")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 19, 100%")

subheader <- paste("Tillman (k = 19)")
mycolors <- c("#FFFFFF")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.1, y = 0.75, label = " Original MA effect = 0.421 \n Reproduced MA effect = 0.420", hjust = 0) 


# Toosi
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Toosi")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 8, 40%", "Not enough information, k = 6, 30%", "Different effect, k = 3, 15%", "Ambiguous effect, k = 3, 15%")

subheader <- paste("Toosi (k = 20)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.4, y = 0.70, label = " Original MA effect = 0.093 \n Reproduced MA effect = 0.086", hjust = 0) 


# Van Iddekinge
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="van Iddekinge")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 9, 90%", "Different effect, k = 1, 10%")

subheader <- paste("van Iddekinge (k = 10)")
mycolors <- c("#FFFFFF", "#FFa500")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.15, y = 0.31, label = " Original MA effect = 0.137 \n Reproduced MA effect = 0.141", hjust = 0) 


# Webb
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Webb")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 8, 40%", "Not enough information, k = 9, 45%", "Different effect, k = 2, 10%", "Ambiguous effect, k = 1, 5%")

subheader <- paste("Webb (k = 20)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.05, y = 0.8, label = " Original MA effect = 0.338 \n Reproduced MA effect = 0.343", hjust = 0) 


# Woodin
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Woodin")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 7, 58%", "Different effect, k = 2, 17%", "Ambiguous effect, k = 3, 25%")

subheader <- paste("Woodin (k = 12)")
mycolors <- c("#FFFFFF", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 0.05, y = 0.53, label = " Original MA effect = 0.306 \n Reproduced MA effect = 0.313", hjust = 0) 


# Woodley
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Woodley")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 3, 50%", "Not enough information, k = 2, 33%", "Different effect, k = 1, 17%")

subheader <- paste("Woodley (k = 6)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.05, y = 0.28, label = " Original MA effect = 0.104 \n Reproduced MA effect = 0.103", hjust = 0) 


# Yoon
library(readxl);library(ggplot2);library(gridExtra)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

df <- read_excel("C:/Users/Esther/Dropbox/Thesis/Codebook primary studies.xlsx", 2)
dfma <- subset(df, meta=="Yoon")

dfma$Discr <- factor(dfma$infonew)

levels(dfma$Discr) = c("No discrepancy, k = 11, 100%")

subheader <- paste("Yoon (k = 11)")
mycolors <- c("#FFFFFF", "#FFFF00", "#FFa500", "#FF0000")

ggplot(dfma, aes(x = znewexp, y = zexp, col=Discr)) + 
  geom_point(size = 3) +
  geom_abline() + 
  scale_shape(solid = FALSE) +
  theme(legend.justification = c(0.9, 0), legend.position = c(0.985, 0.02), legend.title=element_blank()) +
  scale_x_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_y_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=mycolors) +
  ggtitle(subheader) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.3, y = 0.47, label = " Original MA effect = 0.069 \n Reproduced MA effect = 0.069", hjust = 0) 