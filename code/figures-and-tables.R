# Configuring settings and loading packages and datasets ------------------
# rm(list = ls()) #clean workspace
oldw <- getOption("warn")                                         # temporarily surpress warning messages so it doesnt print in Rmd file
options(warn = -1)                                                # temporarily surpress warning messages so it doesnt print in Rmd file
options(warn = oldw, scipen=999)                                  # turn warning messages back on, no scientifc notiation
number_ticks <- function(n) {function(limits) pretty(limits, n)}  # Number of ticks in figures

# # Study 1 Scatterplot - SMD [Figure 2] ----------------------------------
df.smd <- df[df$efftype == "g" | df$efftype == "d",]
df.smd$discrepancy <- as.character(df.smd$info)
# table(df.smd$info) # check frequencies

df.smd$discrepancy[df.smd$discrepancy=="0"] <- a <- "Reproducible, k = 107, 43%"
df.smd$discrepancy[df.smd$discrepancy=="2"] <- b <- "Incomplete, k = 40, 16%"
df.smd$discrepancy[df.smd$discrepancy=="1"] <- c <- "Incorrect, k = 40, 16%"
df.smd$discrepancy[df.smd$discrepancy=="3"] <- d <- "Ambiguous, k = 60, 24%"

df.smd$discrepancy <- factor(df.smd$discrepancy, levels = c(a,b,c,d))

col.magma <- rev(magma(4)) # These are colorblind friendly

# Transform Cohen's d to Hedges' g
g.reported <- g.reproduced <- c()    # empty vectors to store results

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

scatter.smd <- ggplot(df.smd, aes(g.reproduced, g.reported, col=discrepancy, alpha = discrepancy)) + 
  geom_point(aes(shape = discrepancy), size = 6, stroke=0.2) +
  geom_abline() + 
  scale_y_continuous("Reproduced effect sizes (Hedges' g)", breaks=number_ticks(5)) +
  scale_x_continuous("Original effect sizes (Hedges' g)", breaks=number_ticks(5)) + 
  scale_color_manual(values=col.magma, name = "") +
  scale_shape_manual("", values=c(15:18)) +
  scale_alpha_manual("",values=c(0.2, 0.5, 1, 1)) +
  expand_limits(x = c(-2, 6), y = c(-2, 6)) +
  theme(legend.text=element_text(size=14)) +
  theme(legend.position=c(.85, .17)) +
  theme(legend.title=element_blank()) +
  theme(axis.title=element_text(size=18)) +
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) +
  theme(legend.key = element_rect(fill = "#C0C0C0")) +
  theme(panel.background = element_rect(fill = "#C0C0C0"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))



# Study 1 Scatterplot - Correlation [Figure 3] ----------------------------
df.cor <- df[df$efftype == "r" | df$efftype == "z",]
df.cor$discrepancy <- as.character(df.cor$info)
# table(df.cor$info) # check frequencies

df.cor$discrepancy[df.cor$discrepancy=="0"] <- a <- "Reproducible, k = 169, 67%"
df.cor$discrepancy[df.cor$discrepancy=="2"] <- b <- "Incomplete, k = 14, 6%"
df.cor$discrepancy[df.cor$discrepancy=="1"] <- c <- "Incorrect, k = 34, 13%"
df.cor$discrepancy[df.cor$discrepancy=="3"] <- d <- "Ambiguous, k = 36, 14%"

df.cor$discrepancy <- factor(df.cor$discrepancy, levels = c(a,b,c,d))

# Transform correlation r to Fisher's z
z.reported <- z.reproduced <- c()  # empty vectors to store results

for (i in 1:nrow(df.cor)){
  if (df.cor$efftype[i] == "r" ) {
    z.reported[i] <- 0.5 * log((1 + df.cor$effest.exp[i]) / (1 - df.cor$effest.exp[i]))
    z.reproduced[i] <- 0.5 * log((1 + df.cor$effestnew.exp[i]) / (1 - df.cor$effestnew.exp[i]))
  } else {
    z.reported[i] <- df.cor$effest.exp[i]
    z.reproduced[i] <- df.cor$effestnew.exp[i]
  }
}

scatter.cor <- ggplot(df.cor, aes(z.reported, z.reproduced, col=discrepancy, alpha = discrepancy)) +
  geom_point(aes(shape = discrepancy),size = 6, stroke=0.2) +
  geom_abline() + 
  scale_y_continuous("Reproduced effect sizes (Fisher's z)", breaks=number_ticks(5)) +
  scale_x_continuous("Original effect sizes (Fisher's z)", breaks=number_ticks(5)) + 
  scale_color_manual(values=col.magma, name = "") +
  scale_shape_manual("",values=c(15:18)) + 
  scale_alpha_manual("",values=c(0.2, 0.5, 1, 1)) +
  expand_limits(x = c(-0.75, 1.25), y = c(-0.75, 1.25)) +
  theme(legend.text=element_text(size=14)) +
  theme(legend.position=c(.85, .17)) +
  theme(legend.title=element_blank()) +
  theme(axis.title=element_text(size=18))  +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  theme(legend.key = element_rect(fill = "#C0C0C0")) +
  theme(panel.background = element_rect(fill = "#C0C0C0"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))


# Study 1 Frequency table -------------------------------------------------
# Type of effect
typeof.info <- table(df$efftype,df$info)

smds.disc <- sum(typeof.info[1:2,2:4])
smds.nodisc <- sum(typeof.info[1:2,1])
cors.disc <- sum(typeof.info[3:4,2:4])
cors.nodisc <- sum(typeof.info[3:4,1])

# Outlier studies
out.info <- table(df$typestudy,df$info)

out.disc <- sum(out.info[1,2:4])
out.nodisc <- sum(out.info[1,1])
nout.disc <- sum(out.info[2,2:4])
nout.nodisc <- sum(out.info[2,1])

# Published studies
pub.info <- table(df$typecat,df$info)

pub.disc <- sum(pub.info[1,2:4])
pub.nodisc <- sum(pub.info[1,1])
npub.disc <- sum(pub.info[2,2:4])
npub.nodisc <- sum(pub.info[2,1])

# Table 1
table.1 <- data.frame(matrix(ncol = 6, nrow = 6))
row <- c("SMD", "correlation", "outlier", "non-outlier", "published","unpublished")
col <- c("Irreproducible","","Reproducible","","Total","")
colnames(table.1) <- col
rownames(table.1) <- row

table.1[1,1] <- smds.disc;
table.1[1,2] <- c(paste0("(",round((smds.disc/(smds.disc+smds.nodisc))*100,0),"%)"))
table.1[1,3] <- smds.nodisc; 
table.1[1,4] <- c(paste0("(",round((smds.nodisc/(smds.disc+smds.nodisc))*100,0),"%)"))
table.1[1,5] <- smds.disc+smds.nodisc;
table.1[1,6] <- c("(100%)")
table.1[2,1] <- cors.disc; 
table.1[2,2] <- c(paste0("(",round((cors.disc/(cors.disc+cors.nodisc))*100,0),"%)"))
table.1[2,3] <- cors.nodisc; 
table.1[2,4] <- c(paste0("(",round((cors.nodisc/(cors.disc+cors.nodisc))*100,0),"%)"))
table.1[2,5] <- cors.disc+cors.nodisc;
table.1[2,6] <- c("(100%)")
table.1[3,1] <- out.disc;  
table.1[3,2] <- c(paste0("(",round((out.disc/(out.disc+out.nodisc))*100,0),"%)"))
table.1[3,3] <- out.nodisc;
table.1[3,4] <- c(paste0("(",round((out.nodisc/(out.disc+out.nodisc))*100,0),"%)"))
table.1[3,5] <- out.disc+out.nodisc;
table.1[3,6] <- c("(100%)")
table.1[4,1] <- nout.disc; 
table.1[4,2] <- c(paste0("(",round((nout.disc/(nout.disc+nout.nodisc))*100,0),"%)"))
table.1[4,3] <- nout.nodisc; 
table.1[4,4] <- c(paste0("(",round((nout.nodisc/(nout.disc+nout.nodisc))*100,0),"%)"))
table.1[4,5] <- nout.disc+nout.nodisc;
table.1[4,6] <- c("(100%)")
table.1[5,1] <- pub.disc;  
table.1[5,2] <- c(paste0("(",round((pub.disc/(pub.disc+pub.nodisc))*100,0),"%)"))
table.1[5,3] <- pub.nodisc;  
table.1[5,4] <- c(paste0("(",round((pub.nodisc/(pub.disc+pub.nodisc))*100,0),"%)"))
table.1[5,5] <- pub.disc+pub.nodisc;
table.1[5,6] <- c("(100%)")
table.1[6,1] <- npub.disc; 
table.1[6,2] <- c(paste0("(",round((npub.disc/(npub.disc+npub.nodisc))*100,0),"%)"))
table.1[6,3] <- npub.nodisc; 
table.1[6,4] <- c(paste0("(",round((npub.nodisc/(npub.disc+npub.nodisc))*100,0),"%)"))
table.1[6,5] <- npub.disc+npub.nodisc;
table.1[6,6] <- c("(100%)")


# Study 1 Barplot ---------------------------------------------------------
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
dfma$MA <- factor(dfma$MA, levels=unique(sort(dfma$MA)))
dfma$MAnum <- as.character(as.numeric(dfma$MA))
dfma$variable <- factor(dfma$variable)
colnames(dfma)[colnames(dfma)=="value"] <- "Frequency"

barplot.ma <- dfma %>%
  arrange(MAnum) %>%
  mutate(MAnum=factor(MAnum, levels=1:33)) %>% 
  ggplot(aes(x = MAnum, y = Frequency, fill = variable)) +
  geom_bar(stat = "identity") +
  theme(legend.title=element_blank()) +
  ylab("Primary study frequency") + 
  xlab("Meta-analysis") +
  theme(legend.text=element_text(size=14)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  theme(axis.title=element_text(size=18))  +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  scale_fill_manual(values=rev(col.magma), guide = guide_legend(reverse=TRUE)) +
  theme(legend.key = element_rect(fill = "#C0C0C0")) +
  theme(panel.background = element_rect(fill = "#C0C0C0"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))

# To make it horizontal, add + coord_flip()


# Study 2 SMD figures -----------------------------------------------------
eff.so <- eff.sc <- cilb.so <- cilb.sc <- ciub.so <- ciub.sc <- tau2.so <- tau2.sc <- c() # empty vectors to store results

# First, we need to get the estimates from all meta-analyses that have Hedges' g and Fisher's z and do not need to be transformed

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
df <- read.table("../data-per-ma/elsequest_subset.csv", header=T, sep=';') 
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$g.o <- df$effest.exp * J.o                              # cohen's d
df$g.c <- df$effestnew.exp * J.c                           # cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 
write.table(df, file = "../data-per-ma/elsequest_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/green_subset.csv", header=T, sep=';') 
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$g.o <- df$effest.exp * J.o                              # cohen's d
df$g.c <- df$effestnew.exp * J.c                           # cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 
write.table(df, file = "../data-per-ma/green_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/morgan_subset.csv", header=T, sep=';') 
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$g.o <- df$effest.exp * J.o                              # cohen's d
df$g.c <- df$effestnew.exp * J.c                           # cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 
write.table(df, file = "../data-per-ma/morgan_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/webb_subset.csv", header=T, sep=';') 
J.o <- 1 - (3 / (4 * df$n - 9))
J.c <- 1 - (3 / (4 * df$nnew - 9)) 
df$g.o <- df$effest.exp * J.o                              # cohen's d
df$g.c <- df$effestnew.exp * J.c                           # cohen's d
df$vg.o <- J.o^2 * df$vd.o                                 # variance hedges' g 
df$vg.c <- J.c^2 * df$vd.c                                 # variance hedges' g 
write.table(df, file = "../data-per-ma/webb_subset.csv", row.names=FALSE, sep=";")

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

# Study 2 correlation figures ---------------------------------------------

# Correlation r to Fisher's z

# Berry 1
df <- read.table("../data-per-ma/berry1_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "../data-per-ma/berry1_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/berry2_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "../data-per-ma/berry2_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/card_subset.csv", header=T, sep=";")

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
df <- read.table("../data-per-ma/crook_subset.csv", header=T, sep=";")
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
df <- read.table("../data-per-ma/dewit_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "../data-per-ma/dewit_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/farber_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "../data-per-ma/farber_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/card_subset.csv", header=T, sep=";")

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
df <- read.table("../data-per-ma/freund_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "../data-per-ma/freund_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/kolden_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "../data-per-ma/kolden_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/lucassen_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "../data-per-ma/lucassen_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/smith_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                              # variance fisher's z 
write.table(df, file = "../data-per-ma/smith_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/tillman_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                            
write.table(df, file = "../data-per-ma/tillman_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/toosi_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                            
write.table(df, file = "../data-per-ma/toosi_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/vaniddekinge_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                            
write.table(df, file = "../data-per-ma/vaniddekinge_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/woodley_subset.csv", header=T, sep=";")
df$z.o <- 0.5 * log((1 + df$effest.exp) / (1 - df$effest.exp))
df$z.c <- 0.5 * log((1 + df$effestnew.exp) / (1 - df$effestnew.exp))
df$vz.o <- 1 / (df$n - 3)                              # variance fisher's z 
df$vz.c <- 1 / (df$nnew - 3)                            
write.table(df, file = "../data-per-ma/woodley_subset.csv", row.names=FALSE, sep=";")

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
df <- read.table("../data-per-ma/yoon_subset.csv", header=T, sep=";")

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

subg <- which(datm$efftype == "g" | datm$efftype == "d")
subz <- which(datm$efftype == "r" | datm$efftype == "z")
eff.so.g <- eff.so[subg]
eff.so.z <- eff.so[subz]
eff.sc.g <- eff.sc[subg]
eff.sc.z <- eff.sc[subz]

# How many average effects I estimate larger than reported?
sum(datm$eff.so-datm$eff.sc < 0)

# How many average effects I estimate smaller than reported?
sum(datm$eff.so-datm$eff.sc > 0)

ci.so <- ciub.so - cilb.so
ci.sc <- ciub.sc - cilb.sc
ci.so.g <- ci.so[subg]
ci.sc.g <- ci.sc[subg]
ci.so.z <- ci.so[subz]
ci.sc.z <- ci.sc[subz]

# How many confidence intervals I estimated to be wider?
sum(ci.so - ci.sc < 0)
# How many confidence intervals I estimated smaller?
sum(ci.so - ci.sc > 0)

tau2.so.g <- tau2.so[subg]
tau2.so.z <- tau2.so[subz]
tau2.sc.g <- tau2.sc[subg]
tau2.sc.z <- tau2.sc[subz]

# How many tau2 parameters I estimated to be higher?
sum(tau2.so - tau2.sc < 0)
# How many tau2 parameters I estimated smaller?
sum(tau2.so - tau2.sc > 0)
# How many tau2 parameters show no difference?
sum(tau2.so - tau2.sc == 0)

col.magma.long.g <- magma(17)
col.magma.long.z <- magma(17)

# Figure 5 panel a: average meta-analytic effect --------------------------
datplot <- subset(datm,efftype=="g" | efftype=="d")

# Data points are different shapes and colors

#ma.smd.eff <- ggplot(datplot, aes(eff.so.g, eff.sc.g, col = factor(author))) +
#  geom_point(aes(shape = factor(author)), alpha = 0.8, size = 9, stroke=0.9) +
#  geom_abline() + 
#  scale_y_continuous("Reproduced pooled MA effect size", breaks=number_ticks(6)) +
#  scale_x_continuous("Original pooled MA effect size (Hedges' g)", breaks=number_ticks(6)) +
#  expand_limits(x = c(-0.5, 2.5), y = c(-0.5, 2.5)) +
#  scale_shape_manual(values=rep(c(15:19), times=4)) +
#  scale_color_manual("", values=col.magma.long.g) +
#  theme(legend.position="none") +
#  theme(axis.title=element_text(size=16)) +
#  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
#  theme(panel.background = element_rect(fill = "#C0C0C0"),
#        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
#        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))

ma.smd.eff <- ggplot(datplot, aes(eff.so.g, eff.sc.g, label = as.numeric(author))) +
  geom_abline() +
  geom_point() +
  geom_label_repel() +
  scale_y_continuous("Reproduced MA effect size ", breaks=number_ticks(6)) +
  scale_x_continuous("Original pooled MA effect size (Hedges' g)", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.5, 2.5), y = c(-0.5, 2.5)) +
  theme(legend.position="none") +
  theme(axis.title=element_text(size=16)) +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  theme(panel.background = element_rect(fill = "#C0C0C0"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))

# Figure 5 panel b: confidence interval of MA effect ----------------------

# Data points are different shapes and colors

#ma.smd.ci <- ggplot(datplot, aes(ci.so.g, ci.sc.g, col = factor(author))) + 
#  geom_point(aes(shape = factor(author)), alpha = 0.8, size = 9, stroke=0.9) +
#  geom_abline() + 
#  scale_y_continuous("Reproduced MA effect size CI", breaks=number_ticks(6)) +
#  scale_x_continuous("Original pooled MA effect size CI, difference lower and upper bound (Hedges' g)", breaks=number_ticks(6)) +
#  expand_limits(x = c(-0.5, 2.5), y = c(-0.5, 2.5)) +
#  scale_shape_manual(values=rep(c(15:19), times=4)) +
#  scale_color_manual("", values=col.magma.long.g) +
#  theme(legend.position="none") +
#  theme(axis.title=element_text(size=16)) +
#  theme(axis.text.x = element_text(size=14), 
#        axis.text.y = element_text(size=14)) + 
#  theme(panel.background = element_rect(fill = "#C0C0C0"),
#        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
#        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))

ma.smd.ci <- ggplot(datplot, aes(ci.so.g, ci.sc.g, label = as.numeric(author))) + 
  geom_abline() +
  geom_point() +
  geom_label_repel() +
  scale_y_continuous("Reproduced MA effect size CI", breaks=number_ticks(6)) +
  scale_x_continuous("Original pooled MA effect size CI, difference lower and upper bound (Hedges' g)", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.5, 2.5), y = c(-0.5, 2.5)) +
  theme(legend.position="none") +
  theme(axis.title=element_text(size=16)) +
  theme(axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14)) + 
  theme(panel.background = element_rect(fill = "#C0C0C0"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))

# Figure 5 panel c: tau2 smd ----------------------------------------------

# Data points are different shapes and colors

#ma.smd.tau <- ggplot(datplot, aes(tau2.so.g, tau2.sc.g, col = factor(author))) + 
#  geom_point(aes(shape = factor(author)), alpha = 0.8, size = 9, stroke=0.9) +
#  geom_abline() + 
#  scale_y_continuous("Reproduced MA tau2 estimate", breaks=number_ticks(6)) +
#  scale_x_continuous("Original MA tau2 estimate", breaks=number_ticks(6)) +
#  expand_limits(x = c(-0.5, 0.12), y = c(-0.5, 0.12)) +
#  scale_shape_manual(values=rep(c(15:19), times=4)) +
#  scale_color_manual("", values=col.magma.long.g) +
#  theme(legend.position="none") +
#  theme(axis.title=element_text(size=16)) +
#  theme(axis.text.x = element_text(size=14), 
#        axis.text.y = element_text(size=14)) +   
#  theme(panel.background = element_rect(fill = "#C0C0C0"),
#        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
#        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))
#

ma.smd.tau <- ggplot(datplot, aes(tau2.so.g, tau2.sc.g, label = as.numeric(author))) + 
  geom_abline() +
  geom_point() +
  geom_label_repel() +
  scale_y_continuous("Reproduced MA tau2 estimate", breaks=number_ticks(6)) +
  scale_x_continuous("Original MA tau2 estimate", breaks=number_ticks(6)) +
  expand_limits(x = c(-0.5, 0.12), y = c(-0.5, 0.12)) +
  theme(legend.position="none") +
  theme(axis.title=element_text(size=16)) +
  theme(axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14)) +   
  theme(panel.background = element_rect(fill = "#C0C0C0"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))

# Figure 6 panel a: average meta-analytic effect --------------------------
datplot <- subset(datm,efftype=="r" | efftype=="z")

# Data points are different shapes and colors

#ma.cor.eff <- ggplot(datplot, aes(eff.so.z, eff.sc.z, col = factor(author))) + 
#  geom_point(aes(shape = factor(author)), alpha = 0.8, size = 8, stroke=0.9) +
#  geom_abline() + 
#  scale_y_continuous("Reproduced pooled MA effect size", breaks=number_ticks(6)) +
#  scale_x_continuous("Original pooled MA effect size (Fisher's z)", breaks=number_ticks(6)) +
#  expand_limits(x = c(0, 0.45), y = c(0, 0.45)) +
#  scale_color_manual("", values=col.magma.long.z) +
#  scale_shape_manual("", values=rep(c(15:19), times=4)) +
#  theme(legend.position="none")+
#  theme(axis.title=element_text(size=16)) +
#  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
#  theme(panel.background = element_rect(fill = "#C0C0C0"),
#        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
#        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))

ma.cor.eff <- ggplot(datplot, aes(eff.so.z, eff.sc.z, label = as.numeric(author))) + 
  geom_abline() +
  geom_point() +
  geom_label_repel() +
  scale_y_continuous("Reproduced MA effect size", breaks=number_ticks(6)) +
  scale_x_continuous("Original pooled MA effect size (Fisher's z)", breaks=number_ticks(6)) +
  expand_limits(x = c(0, 0.45), y = c(0, 0.45)) +
  theme(legend.position="none")+
  theme(axis.title=element_text(size=16)) +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  theme(panel.background = element_rect(fill = "#C0C0C0"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))


# Figure 6 panel b: confidence interval of MA effect ----------------------

# Data points are different shapes and colors

#ma.cor.ci <- ggplot(datplot, aes(ci.so.z, ci.sc.z, col = factor(author))) + 
#  geom_point(aes(shape = factor(author)), alpha = 0.8, size = 9, stroke=0.9) +
#  geom_abline() + 
#  scale_y_continuous("Reproduced confidence interval", breaks=number_ticks(6)) +
#  scale_x_continuous("Original pooled MA effect size confidence interval (Fisher's z)", breaks=number_ticks(6)) +
#  expand_limits(x = c(0, 0.45), y = c(0, 0.45)) +
#  scale_color_manual("", values=col.magma.long.z) +
#  scale_shape_manual(values=rep(c(15:19), times=4)) +
#  theme(legend.position="none")+
#  theme(axis.title=element_text(size=16)) +
#  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
#  theme(panel.background = element_rect(fill = "#C0C0C0"),
#        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
#        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))

ma.cor.ci <- ggplot(datplot, aes(ci.so.z, ci.sc.z, label = as.numeric(author))) + 
  geom_abline() +
  geom_point() +
  geom_label_repel() +
  scale_y_continuous("Reproduced confidence interval", breaks=number_ticks(6)) +
  scale_x_continuous("Original pooled MA effect size CI, difference lower and upper bound (Fisher's z)", breaks=number_ticks(6)) +
  expand_limits(x = c(0, 0.45), y = c(0, 0.45)) +
  theme(legend.position="none")+
  theme(axis.title=element_text(size=16)) +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  theme(panel.background = element_rect(fill = "#C0C0C0"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))


# Figure 6 panel c: tau2 cor ----------------------------------------------

# Data points are different shapes and colors

#ma.cor.tau <- ggplot(datplot, aes(tau2.so.z, tau2.sc.z, col = factor(author))) + 
#  geom_point(aes(shape = factor(author)), alpha = 0.8, size = 9, stroke=0.6) +
#  geom_abline() + 
#  scale_y_continuous("Reproduced MA tau2", breaks=number_ticks(6)) +
#  scale_x_continuous("Original MA tau2 estimate", breaks=number_ticks(6)) +
#  expand_limits(x = c(0, 0.11), y = c(0, 0.11)) +
#  scale_color_manual("", values=col.magma.long.z) +
#  scale_shape_manual(values=rep(c(15:19), times=4)) +
#  theme(legend.position="none")+
#  theme(axis.title=element_text(size=16)) +
#  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
#  theme(panel.background = element_rect(fill = "#C0C0C0"),
#        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
#        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))

ma.cor.tau <- ggplot(datplot, aes(tau2.so.z, tau2.sc.z, label = as.numeric(author))) + 
  geom_abline() +
  geom_point() +
  geom_label_repel() +
  scale_y_continuous("Reproduced MA tau2", breaks=number_ticks(6)) +
  scale_x_continuous("Original MA tau2 estimate", breaks=number_ticks(6)) +
  expand_limits(x = c(0, 0.06), y = c(0, 0.06)) +
  theme(legend.position="none")+
  theme(axis.title=element_text(size=16)) +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  theme(panel.background = element_rect(fill = "#C0C0C0"),
        panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"),
        panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "#a0a0a0"))
