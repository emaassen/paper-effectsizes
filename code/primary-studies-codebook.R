# rm(list = ls()) # clear workspace
options(scipen=999,warn=0) # show only warning notification 
#options(scipen=999,warn=1) # show all warnings
#options(scipen=999,warn=-1) # supress warnings (caused by coerced NAs) 
packages <- c("readxl","writexl","MAd")
sapply(packages,install.packages(packages),character.only=T)
sapply(packages,library,character.only=T)

# Open codebook
df <- read_excel("../codebooks/codebook-primary-studies-empty.xlsx")

# Functions effect size transformations -----------------------------------
agg_f_to_g1_benish <- function(x) {
  
  hedgesg <- c()       
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
      
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    # transform to cohens d
    d <- sqrt((teststat * n) / (nc * nt))
    var.cohensd <- ((nc + nt) / (nc * nt)) + ((d^2) / (2 * (nc + nt)))
    
    benish.agg <- sum(d / var.cohensd) / sum(1 / var.cohensd)
    
    # transform to hedges g
    J <- 1 - (3 / ((4 * n) - 9))        # this correction factor is different in g2 and g3
    hedgesg[k] <- J * benish.agg
    
  }
  return(hedgesg)
}

agg_smd_to_g1_benish <- function(x) {
  
  hedgesg <- c()  # empty vector to store results
  
  for (k in 1:nrow(x)) {
    
    cohensd <- var.cohensd <- c()      # empty vector to store results
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    for (l in 1:14) {
      
      mc <-  as.numeric(eval(parse(text = paste0("x$mc.",l))))[k]
      mt <-  as.numeric(eval(parse(text = paste0("x$mt.",l))))[k]
      sdc <- as.numeric(eval(parse(text = paste0("x$sdc.",l))))[k]
      sdt <- as.numeric(eval(parse(text = paste0("x$sdt.",l))))[k]
      
      d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
      cohensd <- c(cohensd,d)   # vector with all effect sizes within one primary study 
      
      var.d <- ((nc + nt) / (nc * nt)) + ((d^2) / (2 * (nc + nt)))
      var.cohensd <- c(var.cohensd,var.d)
    }
    
    cohensd <- cohensd[!is.na(cohensd)]
    var.cohensd <- var.cohensd[!is.na(var.cohensd)]
   
    benish.agg <- sum(cohensd / var.cohensd) / sum(1 / var.cohensd)
    
    # transform to hedges g
    J <- 1 - (3 / ((4 * n) - 9))        # this correction factor is different in g2 and g3
    hedgesg[k] <- J * benish.agg
  }
  return(hedgesg)
}
      
agg_r_to_r <- function(x) {
  
  r <- c()
  
  for (k in 1:nrow(x)) {
    
    r.temp <- c()
    
    for (l in 1:14) {
      
      r.est <- as.numeric(eval(parse(text = paste0("x$cor.",l))))[k]
      r.temp <- c(r.temp,r.est)                   # vector with all effect sizes within one primary study 
      
    }
    
    r.temp <- r.temp[!is.na(r.temp)]
    int.cor <- as.numeric(x$int.cor[k])
    c <- length(r.temp)

    r[k] <- sum(r.temp) / sqrt(c + c*(c-1) * int.cor)
  }
  return(r)
}  

agg_d_to_r <- function(x) {
  
 library(MAd)
  
  r <- c()
  
  for (k in 1:nrow(x)) {
    
    r.temp <- c()
    
    for (l in 1:14) {
      
      r.est <- as.numeric(eval(parse(text = paste0("x$cor.",l))))[k]
      r.temp <- c(r.temp,r.est)                   # vector with all effect sizes within one primary study 
    
    }
    
    r.temp <- r.temp[!is.na(r.temp)]
    int.cor <- as.numeric(x$int.cor[k])
    
    #transform to cohens d
    d <- (2 * r.temp) / sqrt(1 - r.temp^2)
    
    n <- as.numeric(x$nnew)[k]
    n1 <- n/2
    n2 <- n/2
    
    id <- c(rep(1, length(d)))
    MAdf <- data.frame(d, n1, n2, id)
    
    aggregate <- agg(id = id, es = d, n.1 = n1, n.2 = n2, data = MAdf, cor = int.cor, method = "GO1")
    
    cohensd <- aggregate$es1
    a <- (n1 + n2)^2 / (n1 * n2)
    r[k] <- cohensd / (sqrt(cohensd^2 + a))
    
  }
  return(r)
} 

etasq_to_g1 <- function(x) {
  
  hedgesg <- c()                  # empty vector to store results                               
  
  for (k in 1:nrow(x)) {
    
    etasq.temp <- c()                 # empty vectors to store results 
    
    for (l in 1:14) {
      
      etasquare <- as.numeric(eval(parse(text = paste0("x$cor.",l))))[k]
      etasq.temp <- c(etasq.temp,etasquare)                   # vector with all effect sizes within one primary study 
    }
    
    n <- as.numeric(x$nnew)[k]
    
    # take mean, transform to r
    etasq <- mean(etasq.temp, na.rm=T) # mean of all effects from one primary study
    r <- sqrt(etasq)
    
    # transform to cohens d and hedges g
    cohensd <- (2*r)/sqrt(1-r^2)
    J <- 1 - (3 / ((4 * n) - 9))        
    hedgesg[k] <- J * cohensd
  }
  return(hedgesg)
}

f_to_d <- function(x) {
  
  cohensd <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
  
    }
    
    teststat <- teststat[!is.na(teststat)]
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    # transform to cohens d
    d <- sqrt((teststat * n) / (nc * nt))
    
    # take mean 
    d <- mean(d, na.rm=T) # mean of all effects from one primary study
    cohensd[k] <- d
    
  }
  return(cohensd)
}

f_to_d_elseq <- function(x) {
  
  cohensd <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
      
    }
    
    teststat <- teststat[!is.na(teststat)]
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    # transform to cohens d
    d <- 2 * sqrt(teststat / n)
    
    # take mean 
    d <- mean(d, na.rm=T) # mean of all effects from one primary study
    cohensd[k] <- d
    
  }
  return(cohensd)
}
    
f_to_g1 <- function(x) {
  
  hedgesg <- c()                              # empty vector to store results                                # empty vector to store result

  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])

    }
    
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    # transform to cohens d
    d <- sqrt((teststat*n)/(nc*nt))
    
    # take mean 
    cohensd <- mean(d, na.rm=T) # mean of all effects from one primary study
    
    # transform to hedges g
    J <- 1 - (3 / ((4 * n) - 9))        # this correction factor is different in f_to_g2
    hedgesg[k] <- J * cohensd
    
  }
  return(hedgesg)
}
      
f_to_g2 <- function(x) {
  
  hedgesg <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
      
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    dfs <- nc+nt-2
    
    # transform to cohens d
    d <- sqrt((teststat*n)/(nc*nt))
    
    # take mean 
    cohensd <- mean(d, na.rm=T) # mean of all effects from one primary study
    
    # transform to hedges g
    J <- 1 - (3 / ((4 * dfs) - 1))       # this correction factor is different in f_to_g1
    hedgesg[k] <- J * cohensd
    
  }
  return(hedgesg)
}

f_to_r1 <- function(x) {
  
  r <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])

    }
    
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    # transform to cohens d
    d <- sqrt((teststat * n) / (nc * nt))
    
    # take mean 
    cohensd <- mean(d, na.rm=T) # mean of all effects from one primary study
    
    # transform to r
    a <- (nc + nt)^2 / (nc * nt)
    r.est <- cohensd / (sqrt(cohensd^2 + a))
    r[k] <- r.est
    
  }
  return(r)
}

f_to_r2 <- function(x) {
  
  r <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]

    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
      
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    #transform to cohens r
    r.est <- sqrt((sqrt(teststat)) / ((sqrt(teststat)) + (nc + nt - 2))) 
    
    # take mean 
    r.est <- mean(r.est, na.rm=T) # mean of all effects from one primary study
    r[k] <- r.est
    
  }
  return(r)
}

chi_to_d <- function(x) {
  
  cohensd <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()                 # empty vector to store results 
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
     
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    # transform to r
    r <- sqrt((teststat/n))
    
    # take mean 
    r <- mean(r, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohens d to hedges g
    d <- (2 * r) / sqrt(1 - r^2)
    cohensd[k] <- d

  }
  return(d)
}

chi_to_g2 <- function(x) {
  
  hedgesg <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()                 # empty vector to store results 
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
      
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    dfs <- nc+nt-2
    
    # transform to r
    r <- sqrt((teststat/n))
    
    # take mean 
    r <- mean(r, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohens d to hedges g
    cohensd <- (2*r)/sqrt(1-r^2)
    J <- 1 - (3 / ((4 * dfs) - 1))       # this correction factor is different in f_to_g1
    hedgesg[k] <- J * cohensd
    
  }
  return(hedgesg)
}

chi_to_r <- function(x) {
  
  r <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()                 # empty vector to store results 
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
      
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    dfs <- nc+nt-2
    
    # transform to r
    r <- sqrt((teststat/n))
    
    # take mean 
    r <- mean(r, na.rm=T) # mean of all effects from one primary study
  }
  return(r)
}

mann_to_r <- function(x) {
  
  rb <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
     
    }
  
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]

    teststat <- teststat[!is.na(teststat)]
    
    #transform to r-pb
    rpb.est <- 1 - (teststat*2 / (nc * nt))
    
    # take mean 
    rpb <- mean(rpb.est, na.rm=T) # mean of all effects from one primary study
    
    # find pvalue based on sample size calculation
    pval <- nt / (nt + nc)
    z <- qnorm(pval) # zscore where p-val on normal distribution
    z.df <- dnorm(z) # density at this specific pvalue
    rb <- (rpb * (sqrt(nt * nc))) / (abs(z.df) * (nt * nc))
    
  }
  return(rb)
}

or_to_d <- function(x) {
  
  cohensd <- c()                                                 # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    a <- b <- c <- d <- c()    # empty vectors to store results
    
    for (l in 1:4) {
      
      a <- c(a,as.numeric(eval(parse(text = paste0("x$or.a.",l))))[k])
      b <- c(b,as.numeric(eval(parse(text = paste0("x$or.b.",l))))[k])
      c <- c(c,as.numeric(eval(parse(text = paste0("x$or.c.",l))))[k])
      d <- c(d,as.numeric(eval(parse(text = paste0("x$or.d.",l))))[k])
    
    }
    or <- (a*d)/(b*c)

    # take mean 
    or <- mean(or, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohen's d and hedges g
    cohensd[k] <- log(or) * (sqrt(3)/pi)
    
  }
  return(cohensd)
}

or_to_g1 <- function(x) {
  
  hedgesg <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    a <- b <- c <- d <- c()    # empty vectors to store results
    
    for (l in 1:4) {
      
      a <- c(a,as.numeric(eval(parse(text = paste0("x$or.a.",l))))[k])
      b <- c(b,as.numeric(eval(parse(text = paste0("x$or.b.",l))))[k])
      c <- c(c,as.numeric(eval(parse(text = paste0("x$or.c.",l))))[k])
      d <- c(d,as.numeric(eval(parse(text = paste0("x$or.d.",l))))[k])
      
   }
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    or <- (a*d)/(b*c)

    # take mean 
    or <- mean(or, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohen's d and hedges g
    cohensd <- log(or) * (sqrt(3)/pi)
    J <- 1 - (3 / ((4 * n) - 9))        # this correction factor is different in or_to_g2
    hedgesg[k] <- J * cohensd
  }
  return(hedgesg)
}

or_to_g2 <- function(x) {
  
  hedgesg <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    a <- b <- c <- d <- c()    # empty vectors to store results
    
    for (l in 1:4) {

      a <- c(a,as.numeric(eval(parse(text = paste0("x$or.a.",l))))[k])
      b <- c(b,as.numeric(eval(parse(text = paste0("x$or.b.",l))))[k])
      c <- c(c,as.numeric(eval(parse(text = paste0("x$or.c.",l))))[k])
      d <- c(d,as.numeric(eval(parse(text = paste0("x$or.d.",l))))[k])
      
    }
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    dfs <- nc+nt-2
    
    or <- (a*d)/(b*c)
    
    # take mean 
    or <- mean(or, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohen's d and hedges g
    cohensd <- log(or) * (sqrt(3)/pi)
    J <- 1 - (3 / ((4 * dfs) - 1))          # this correction factor is different in or_to_g1   
    hedgesg[k] <- J * cohensd
  }
  return(hedgesg)
}

or_to_g2_koenig <- function(x) {
  
  hedgesg <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    a <- b <- c <- d <- c()    # empty vectors to store results
    
    for (l in 1:4) {
      
      a <- c(a,as.numeric(eval(parse(text = paste0("x$or.a.",l))))[k])
      b <- c(b,as.numeric(eval(parse(text = paste0("x$or.b.",l))))[k])
      c <- c(c,as.numeric(eval(parse(text = paste0("x$or.c.",l))))[k])
      d <- c(d,as.numeric(eval(parse(text = paste0("x$or.d.",l))))[k])
      
    }
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    dfs <- nc+nt-2
    
    or <- (a*d)/(b*c)
    
    # take mean 
    or <- mean(or, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohen's d and hedges g
    cohensd <- log(or) / 1.65               # this calculation is different from or_to_g2 general- only used by koenig
    J <- 1 - (3 / ((4 * dfs) - 1))          # this correction factor is different in or_to_g1   
    hedgesg[k] <- J * cohensd
  }
  return(hedgesg)
}

or_to_r <- function(x) {
  
  r <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    a <- b <- c <- d <- c()    # empty vectors to store results
    
    for (l in 1:4) {
      
      a <- c(a,as.numeric(eval(parse(text = paste0("x$or.a.",l))))[k])
      b <- c(b,as.numeric(eval(parse(text = paste0("x$or.b.",l))))[k])
      c <- c(c,as.numeric(eval(parse(text = paste0("x$or.c.",l))))[k])
      d <- c(d,as.numeric(eval(parse(text = paste0("x$or.d.",l))))[k])
      
    }
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    dfs <- nc + nt - 2
    
    or <- (a*d)/(b*c)
    
    # take mean 
    or <- mean(or, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohen's d and r
    cohensd <- log(or) * (sqrt(3)/pi)
    a <- (nc + nt)^2 / (nc * nt)
    r.est <- cohensd / (sqrt(cohensd^2 + a))
    r[k] <- r.est
  }
  return(r)
}

r_corr <- function(x) {
  
  r <- c()
  
  for (k in 1:nrow(x)) {
    
    rs <- rxx <- ryy <- c()
    
    for (l in 1:14) {
      
      rs <- c(rs,as.numeric(eval(parse(text = paste0("x$cor.",l))))[k])

      
    }
    
    for (m in 1:2) {     # only max 2 * 2 reliabilities to extract
      
      rxx <- c(rxx,as.numeric(eval(parse(text = paste0("x$rxx.",m))))[k])             
      ryy <- c(ryy,as.numeric(eval(parse(text = paste0("x$ryy.",m))))[k])             
      
    }
    
    rs <- rs[!is.na(rs)]
    rxx <- rxx[!is.na(rxx)]
    ryy <- ryy[!is.na(ryy)]
    
    corr.r <- rs / (sqrt(rxx) * sqrt(ryy))
    corr.r <- mean(corr.r)
    
    r[k] <- corr.r
  }
  return(r)
}  

r_to_d <- function(x) {
  
  d <- c()
  
  for (k in 1:nrow(x)) {
    
    r.temp <- c()
    
    for (l in 1:14) {
      
      r.est <- as.numeric(eval(parse(text = paste0("x$cor.",l))))[k]
      r.temp <- c(r.temp,r.est) # vector with all effect sizes within one primary study 
    }
    
    # take mean 
    r <- mean(r.temp, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohens d
    d.est <- (2 * r) / sqrt(1 - r^2)
    d[k] <- d.est
    
  }
  
  return(d)
  
  
  
}

r_to_d_elseq <- function(x) {
  
  d <- c()
  
  for (k in 1:nrow(x)) {
    
    r.temp <- c()
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc + nt
    
    for (l in 1:14) {
      
      r.est <- as.numeric(eval(parse(text = paste0("x$cor.",l))))[k]
      r.temp <- c(r.temp,r.est) # vector with all effect sizes within one primary study 
    }
    
    # take mean 
    r <- mean(r.temp, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohens d
    p = nc / n
    d.est <- r / sqrt(p * (1 - p) * (1 - r^2))
    d[k] <- d.est
    
  }
  
  return(d)
  
}

r_to_g1 <- function(x) {
  
  hedgesg <- c()                  # empty vector to store results                              
  
  for (k in 1:nrow(x)) {
    
    r.temp <- c()                 # empty vectors to store results 
    n <- as.numeric(x$nnew)[k]
    
    for (l in 1:14) {
      
      r <- as.numeric(eval(parse(text = paste0("x$cor.",l))))[k]
      r.temp <- c(r.temp,r)                   # vector with all effect sizes within one primary study 
    }
    
    # take mean
    r <- mean(r.temp, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohens d and hedges g
    cohensd <- (2 * r) / sqrt(1 - r^2)
    J <- 1 - (3 / ((4 * n) - 9))        
    hedgesg[k] <- J * cohensd
  }
  return(hedgesg)
}

r_to_r <- function(x) {
  
  r <- c()
  
  for (k in 1:nrow(x)) {
    
    r.temp <- c()
    
    for (l in 1:14) {
      
      r.est <- as.numeric(eval(parse(text = paste0("x$cor.",l))))[k]
      r.temp <- c(r.temp,r.est) # vector with all effect sizes within one primary study 
    }
    
    # take mean 
    r[k] <- mean(r.temp, na.rm=T) # mean of all effects from one primary study
  }
  
  return(r)
}   

r_to_d <- function(x) {
  
  d <- c()
  
  for (k in 1:nrow(x)) {
    
    r.temp <- c()
    
    for (l in 1:14) {
      
      r.est <- as.numeric(eval(parse(text = paste0("x$cor.",l))))[k]
      r.temp <- c(r.temp,r.est) # vector with all effect sizes within one primary study 
    }
    
    # take mean 
    r <- mean(r.temp, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohens d
    d.est <- (2 * r) / sqrt(1 - r^2)
    d[k] <- d.est
    
  }
  
  return(d)
  
  
  
}

r_to_z <- function(x) {
  
  z <- c()
  
  for (k in 1:nrow(x)) {
    
    r.temp <- c()
    
    for (l in 1:14) {
      
      r.est <- as.numeric(eval(parse(text = paste0("x$cor.",l))))[k]
      r.temp <- c(r.temp,r.est) # vector with all effect sizes within one primary study 
      
    }
    
    # take mean 
    r <- mean(r.temp, na.rm=T) # mean of all effects from one primary study
    
    #transform to fishers z
    z[k] <- 0.5 * log((1 + r) / (1 - r))
  }
  
  return(z)
}

smd_to_d <- function(x) {
  
  cohensd <- c()  # empty vector to store results
  
  for (k in 1:nrow(x)) {
    
    cohensd.temp <- c()      # empty vector to store results
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]

    for (l in 1:14) {
      
      mc <-  as.numeric(eval(parse(text = paste0("x$mc.",l))))[k]
      mt <-  as.numeric(eval(parse(text = paste0("x$mt.",l))))[k]
      sdc <- as.numeric(eval(parse(text = paste0("x$sdc.",l))))[k]
      sdt <- as.numeric(eval(parse(text = paste0("x$sdt.",l))))[k]
      
      d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
      cohensd.temp <- c(cohensd.temp,d)   # vector with all effect sizes within one primary study 
      
    }
    
    # take mean 
    d <- mean(cohensd.temp, na.rm=T) # mean of all effects from one primary study
    cohensd[k] <- d
    
  }
  return(cohensd)
}

smd_to_g1 <- function(x) {
  
  hedgesg <- c()  # empty vector to store results
  
  for (k in 1:nrow(x)) {
    
    cohensd.temp <- c()      # empty vector to store results
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    for (l in 1:14) {
      
      mc <-  as.numeric(eval(parse(text = paste0("x$mc.",l))))[k]
      mt <-  as.numeric(eval(parse(text = paste0("x$mt.",l))))[k]
      sdc <- as.numeric(eval(parse(text = paste0("x$sdc.",l))))[k]
      sdt <- as.numeric(eval(parse(text = paste0("x$sdt.",l))))[k]
      
      d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
      cohensd.temp <- c(cohensd.temp,d)   # vector with all effect sizes within one primary study 
      
    }
    
    # take mean 
    cohensd <- mean(cohensd.temp, na.rm=T) # mean of all effects from one primary study
    
    # transform to hedges g
    J <- 1 - (3 / ((4 * n) - 9))        # this correction factor is different in smd_to_g2
    hedgesg[k] <- J * cohensd
    
  }
  return(hedgesg)
}

smd_to_g2 <- function(x) {
  
  hedgesg <- c()            # empty vector to store results
  
  for (k in 1:nrow(x)) {
    
    cohensd.temp <- c()     # empty vector to store results
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    dfs <- nc+nt-2
    
    for (l in 1:14) {
      
      mc <-  as.numeric(eval(parse(text = paste0("x$mc.",l))))[k]
      mt <-  as.numeric(eval(parse(text = paste0("x$mt.",l))))[k]
      sdc <- as.numeric(eval(parse(text = paste0("x$sdc.",l))))[k]
      sdt <- as.numeric(eval(parse(text = paste0("x$sdt.",l))))[k]
      
      d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
      cohensd.temp <- c(cohensd.temp,d)   # vector with all effect sizes within one primary study 
    }
    
    # take mean 
    cohensd <- mean(cohensd.temp, na.rm=T) # mean of all effects from one primary study
    
    # transform to hedges g
    J <- 1 - (3 / ((4 * dfs) - 1))        # this correction factor is different in smd_to_g1
    hedgesg[k] <- J * cohensd
  }
  return(hedgesg)
}

smd_to_g3 <- function(x) {
  
  hedgesg <- c()            # empty vector to store results
  
  for (k in 1:nrow(x)) {
    
    cohensd.temp <- c()     # empty vector to store results
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    for (l in 1:14) {
      
      mc <-  as.numeric(eval(parse(text = paste0("x$mc.",l))))[k]
      mt <-  as.numeric(eval(parse(text = paste0("x$mt.",l))))[k]
      sdc <- as.numeric(eval(parse(text = paste0("x$sdc.",l))))[k]
      sdt <- as.numeric(eval(parse(text = paste0("x$sdt.",l))))[k]
      
      d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
      cohensd.temp <- c(cohensd.temp,d)   # vector with all effect sizes within one primary study 
    }
    
    # take mean 
    cohensd <- mean(cohensd.temp, na.rm=T) # mean of all effects from one primary study
    
    # transform to hedges g
    J <- 1 - (3 / (4 * (n - 1)))       # this correction factor is different in smd_to_g1 and smd_to_g2
    hedgesg[k] <- J * cohensd
  }
  return(hedgesg)
}

smd_to_r <- function(x) {
  
  r <- c()            # empty vector to store results
  
  for (k in 1:nrow(x)) {
    
    cohensd.temp <- c()     # empty vector to store results
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]

    for (l in 1:14) {
      
      mc <-  as.numeric(eval(parse(text = paste0("x$mc.",l))))[k]
      mt <-  as.numeric(eval(parse(text = paste0("x$mt.",l))))[k]
      sdc <- as.numeric(eval(parse(text = paste0("x$sdc.",l))))[k]
      sdt <- as.numeric(eval(parse(text = paste0("x$sdt.",l))))[k]
      
      d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
      cohensd.temp <- c(cohensd.temp,d)   # vector with all effect sizes within one primary study 
    }
    
    # take mean 
    cohensd <- mean(cohensd.temp, na.rm=T) # mean of all effects from one primary study
    
    # transform to r
    a <- (nc + nt)^2 / (nc * nt)
    r.est <- cohensd / (sqrt(cohensd^2 + a))
    r[k] <- r.est
  }
  return(r)
}

smd_to_r_fox <- function(x) {
  
  r <- c()            # empty vector to store results
  
  for (k in 1:nrow(x)) {
    
    cohensd.temp <- c()     # empty vector to store results
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    for (l in 1:14) {
      
      mc <-  as.numeric(eval(parse(text = paste0("x$mc.",l))))[k]
      mt <-  as.numeric(eval(parse(text = paste0("x$mt.",l))))[k]
      sdc <- as.numeric(eval(parse(text = paste0("x$sdc.",l))))[k]
      sdt <- as.numeric(eval(parse(text = paste0("x$sdt.",l))))[k]
      
      d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
      cohensd.temp <- c(cohensd.temp,d)   # vector with all effect sizes within one primary study 
    }
    
    # take mean 
    cohensd <- mean(cohensd.temp, na.rm=T) # mean of all effects from one primary study
    
    # transform to r
    a.1 <- (nc + nt)^2 / (nc * nt)
    r.est <- cohensd / (sqrt(cohensd^2 + a.1))
    
    #correct for unequal sampling per Fox
    a.2 <- sqrt(0.25 / ((nc/n) * (nt/n)))
    r.est.corr <- (a.2 * r.est) / (sqrt(((a.2^2 - 1) * r.est^2 + 1)))
    
    r[k] <- r.est.corr
  }
  return(r)
}

smdpaired_to_d <- function(x) {
  
  cohensd <- c()  # empty vector to store results
  
  for (k in 1:nrow(x)) {
    
    cohensd.temp <- c()      # empty vector to store results
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    int.cor <- as.numeric(x$int.cor[k])
    
    for (l in 1:14) {
      
      mc <-  as.numeric(eval(parse(text = paste0("x$mc.",l))))[k]
      mt <-  as.numeric(eval(parse(text = paste0("x$mt.",l))))[k]
      sdc <- as.numeric(eval(parse(text = paste0("x$sdc.",l))))[k]
      sdt <- as.numeric(eval(parse(text = paste0("x$sdt.",l))))[k]
      
      sdiff <- sqrt(sdc^2 + sdt^2 - (2 * int.cor * sdc * sdt)) 
      d.est <- (mt - mc) / (sdiff / sqrt(2 * (1 - int.cor)))
      
      cohensd.temp <- c(cohensd.temp,d.est)   # vector with all effect sizes within one primary study 
      
    }
    
    # take mean 
    d <- mean(cohensd.temp, na.rm=T) # mean of all effects from one primary study
    cohensd[k] <- d
    
  }
  return(cohensd)
}

smd_to_d_webb <- function(x) {
  
  cohensd <- c()  # empty vector to store results
  
  for (k in 1:nrow(x)) {
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    df <- nc+nt-2
    N = nc+nt
    
    pval <- 0.05                      # one-tailed test
    t.value <- abs(qt(pval,df)) ;t.value
    #2*pt(-abs(t.value),df) # corresponding p-value calculated from two-tailed t-test
    #getAnywhere(t.test.default)
    
    # calculate standardized mean difference for this test-statistic
    es <- (2 * t.value) / sqrt(N);es
    cohensd[k] <- es
    
  }
  return(cohensd)
}

t_to_d <- function(x) {
  
  cohensd <- c()                              # empty vector to store results                           
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]

    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    # transform to cohens d
    d <- teststat * sqrt((nc+nt)/(nc*nt))
    
    # take mean 
    cohensd[k] <- mean(d, na.rm=T) # mean of all effects from one primary study
    
  }
  return(cohensd)
}

t_to_g1 <- function(x) {
  
  hedgesg <- c()              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc+nt
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])

    }

    teststat <- teststat[!is.na(teststat)]
    
    # transform to cohens d
    d <- teststat * sqrt((nc+nt)/(nc*nt))
    
    # take mean 
    cohensd <- mean(d, na.rm=T) # mean of all effects from one primary study
    
    # transform to hedges g
    J <- 1 - (3 / ((4 * n) - 9))        # this correction factor is different in t_to_g2
    hedgesg[k] <- J * cohensd
    
  }
  return(hedgesg)
}

t_to_g2 <- function(x) {
  
  hedgesg <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    dfs <- nc+nt-2
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])

    }
    
    teststat <- teststat[!is.na(teststat)]
    
    #transform to cohens d
    d <- teststat * sqrt((nc + nt) / (nc * nt))
   
    # take mean 
    cohensd <- mean(d, na.rm=T) # mean of all effects from one primary study
    
    # transform to hedges g
    J <- 1 - (3 / ((4 * dfs) - 1))       # this correction factor is different in t_to_g1
    hedgesg[k] <- J * cohensd
    
  }
  return(hedgesg)
}

t_to_r1 <- function(x) {
  
  r <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    dfs <- nc+nt-2
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
      
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    #transform to cohens d
    d <- teststat * sqrt((nc + nt) / (nc * nt))

    # take mean 
    cohensd <- mean(d, na.rm=T) # mean of all effects from one primary study
    
    # transform to r
    a <- (nc + nt)^2 / (nc * nt)
    r.est <- cohensd / (sqrt(cohensd^2 + a))
    r[k] <- r.est
    
  }
  return(r)
}

t_to_r2 <- function(x) {
  
  r <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]

    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
      
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    #transform to cohens r
    r.est <- sqrt((teststat ^ 2) / ((teststat^2) + (nc + nt - 2))) 
    
    # take mean 
    r.est <- mean(r.est, na.rm=T) # mean of all effects from one primary study
    r[k] <- r.est
    
  }
  return(r)
  

  
}

zstat_to_g2 <- function(x) {
  
  hedgesg <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc + nt
    dfs <- nc+nt-2
    
    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
      
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    #transform to r
    r.est <- teststat / sqrt(n)

    # take mean 
    r <- mean(r.est, na.rm=T) # mean of all effects from one primary study
    
    # transform to cohens d to hedges g
    cohensd <- (2*r)/sqrt(1-r^2)
    J <- 1 - (3 / ((4 * dfs) - 1))       # this correction factor is different in f_to_g1
    hedgesg[k] <- J * cohensd
    
  }
  return(hedgesg)
}

zstat_to_r_fox <- function(x) {
  
  r <- c()                              # empty vector to store results                                # empty vector to store result
  
  for (k in 1:nrow(x)) {
    
    teststat <- c()           # empty vectors to store results 
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc + nt

    for (l in 1:5) {
      
      teststat <- c(teststat,as.numeric(eval(parse(text = paste0("x$teststat.",l))))[k])
      
    }
    
    teststat <- teststat[!is.na(teststat)]
    
    #transform to r
    rs <- teststat / sqrt(n)
    
    # take mean 
    r.est <- mean(rs, na.rm=T) # mean of all effects from one primary study
    
    #correct for unequal sampling per Fox
    a <- sqrt(0.25 / ((nc/n) * (nt/n)))
    r.est.corr <- (a * r.est) / (sqrt(((a^2 - 1) * r.est^2 + 1)))
    
    r[k] <- r.est.corr
    
  }
  return(r)
}

final_d_to_z <- function(x) {
  
  z <- c()
  
  for (k in 1:nrow(x)) {
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc + nt
    d <- as.numeric(x$effest[k])
    
    # transform d or g to point biserial correlation
    m <- nc + nt - 2
    h <- (m / nc) + (m / nt)
    rpb <- d / sqrt(d^2 + h)
    
    # transform point biserial correlation to biserial correlation
    p <- nc/n
    q <- 1 - p
    z <- qnorm(p)
    fzp <- dnorm(z)
    rb <- (sqrt(p * q) / fzp) * rpb
    
    # truncate all estimates of rb > 1 to 1, and < -1 to -1 as per Jacobs & Viechtbauer (2016).
    if (rb > 1) { rb <- 1 }
    if (rb < -1) { rb <- -1 }
    
    # transform biserial to fishers r-to-Z
    a <- sqrt(fzp) / ((p*q)^(1/4))
    z[k] <- (a/2) * log((1 + a*rb) / (1 - a*rb))
    
  }
  return(z)
 }

final_r_to_z <- function(x) {
  
  z <- c()
  
  for (k in 1:nrow(x)) {
    
    nc <- as.numeric(x$ncnew)[k]
    nt <- as.numeric(x$ntnew)[k]
    n <- nc + nt
    rb <- as.numeric(x$effest[k])
    
    # transform point biserial correlation to biserial correlation
    p <- nc/n
    q <- 1 - p
    z <- qnorm(p)
    fzp <- dnorm(z)

    # truncate all estimates of rb > 1 to 1, and < -1 to -1 as per Jacobs & Viechtbauer (2016).
    if (rb > 1) { rb <- 1 }
    if (rb < -1) { rb <- -1 }
    
    # transform biserial to fishers r-to-Z
    a <- sqrt(fzp) / ((p*q)^(1/4))
    z[k] <- (a/2) * log((1 + a*rb) / (1 - a*rb))
    
  }
  return(z)
}


# Manual coding -----------------------------------------------------------
# Some studies had multiple different effect sizes that had to be combined, e.g. a SMD and an OR. 
# In these specific cases, we entered the effect size estimate by hand. The code we used for calculation 
# can be found here, sorted by study number (ID, first column of the codebook primary studies):

# STUDY 25 - Alfieri - Hendrix
OR.a <- 57
OR.b <- 8
OR.c <- 47
OR.d <-  21
OR <- (OR.a * OR.d) / (OR.b * OR.c)
nc <- 14
nt <- 6.5 
dfs <- nc + nt - 2

d1 <- log(OR) * ((sqrt(3)) / pi)
J <- 1 - (3 / ((4 * dfs) - 1))
g1 <- d1 * J

OR.a <- 57
OR.b <- 8
OR.c <- 49
OR.d <-  16
OR <- (OR.a * OR.d) / (OR.b * OR.c)
nc <- 13
nt <- 6.5 
dfs <- nc + nt - 2

d2 <- log(OR) * ((sqrt(3)) / pi)
J <- 1 - (3 / ((4 * dfs) - 1))
g2 <- d2 * J

eff25 <- mean(g1,g2)
df$effestnew[df$id %in% "25"] <- eff25


# STUDY 33 - Alfieri - Quilici
nc <- 27
nt <- 27
dfs = nc + nt - 2

mc <- 0.323	
mt <- 0.049	
sdc <- 0.162	
sdt <- 0.113	
d1 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))

mc <- 0.327	
mt <- 0.049	
sdc <- 0.135	
sdt <- 0.113	
d2 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))

mc <- 0.488	
mt <- 0.873	
sdc <- 0.319	
sdt <- 0.27	
d3 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))

mc <- 0.441	
mt <- 0.873	
sdc <- 0.282	
sdt <- 0.27
d4 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))

chisq <- 8.358
phi <- sqrt(chisq/(nc+nt))
d5 <- (2*phi) / sqrt(1-phi^2)
d5 <- d5*-1 # effect size reversed to right direction

d <- c(d1,d2,d3,d4,d5)
d <- mean(d)
J <- 1 - (3 / ((4 * dfs) - 1))
g <- d * J

eff33 <- g <- d * J
df$effestnew[df$id %in% "33"] <- eff33


# STUDY 74 - Benish - Bradley
nc <- 25
nt <- 25
mc <- 46	
mt <- 44.12	
sdc <- 10.44	
sdt <- 8.82
d1 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
var.d1 <- ((nc + nt) / (nc * nt)) + ((d1^2) / (2 * (nc + nt)))

nc <- 23
nt <- 25
mc <- 26.94	
mt <- 32.78	
sdc <- 15.59	
sdt <- 15.31
d2 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
var.d2 <- ((nc + nt) / (nc * nt)) + ((d2^2) / (2 * (nc + nt)))

nc <- 21
nt <- 25
mc <- 63.22	
mt <- 58.17	
sdc <- 22.34	
sdt <- 14.25
d3 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
d3 <- d3*-1 # reverse effect size to right direction
var.d3 <- ((nc + nt) / (nc * nt)) + ((d3^2) / (2 * (nc + nt)))

nc <- 24
nt <- 25
mc <- 7.66	
mt <- 9.6	
sdc <- 4.85	
sdt <- 4.63
d4 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
var.d4 <- ((nc + nt) / (nc * nt)) + ((d4^2) / (2 * (nc + nt)))

nc <- 22
nt <- 24
mc <- 9.38	
mt <- 18.95	
sdc <- 8.1	
sdt <- 15.39
d5 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
var.d5 <- ((nc + nt) / (nc * nt)) + ((d5^2) / (2 * (nc + nt)))

cohensd <- c(d1,d2,d3,d4,d5)
var.cohensd <- c(var.d1,var.d2,var.d3,var.d4,var.d5)
benish.agg <- sum(cohensd / var.cohensd) / sum(1 / var.cohensd)

# transform to hedges g
n <- 50                       # we use n=50 to calculate, even though all effects make use of a different sample size, of which only one adds up to 25
J <- 1 - (3 / ((4 * n) - 9))       
eff74 <- J * benish.agg
df$effestnew[df$id %in% "74"] <- eff74


# STUDY 75 - Benish - Costantino
eff75 <- .63                                # as reported in the primary paper
df$effestnew[df$id %in% "75"] <- eff75


# STUDY 77 - Benish - Huey
nc <- 4
nt <- 5

d1 <- -0.4
var.d1 <- ((nc + nt) / (nc * nt)) + ((d1^2) / (2 * (nc + nt)))

d2 <- 0.26
var.d2 <- ((nc + nt) / (nc * nt)) + ((d2^2) / (2 * (nc + nt)))

d3 <- 1.13
var.d3 <- ((nc + nt) / (nc * nt)) + ((d3^2) / (2 * (nc + nt)))

d4 <- 1.02
var.d4 <- ((nc + nt) / (nc * nt)) + ((d4^2) / (2 * (nc + nt)))

d5 <- 1.29
var.d5 <- ((nc + nt) / (nc * nt)) + ((d5^2) / (2 * (nc + nt)))

d6 <- 1.63
var.d6 <- ((nc + nt) / (nc * nt)) + ((d6^2) / (2 * (nc + nt)))

cohensd <- c(d1,d2,d3,d4,d5,d6)
var.cohensd <- c(var.d1,var.d2,var.d3,var.d4,var.d5,var.d6)
benish.agg <- sum(cohensd / var.cohensd) / sum(1 / var.cohensd)

# transform to hedges g
n <- nc+nt                   
J <- 1 - (3 / ((4 * n) - 9))       
eff77 <- J * benish.agg
df$effestnew[df$id %in% "77"] <- eff77


# STUDY 121 - Card - Jensen 1988
nc <- 119
nt <- 56
mc <- 55.5	
mt <- 55.4	
sdc <- 9.1	
sdt <- 10.1	
d1 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
A <- (nc+nt)^2 / (nc*nt)
r1 <- d1 / (sqrt(d1^2 + A)) # transform to r seperately, because the sample size is different in each case

nc <- 109
nt <- 52
mc <- 52.5	
mt <- 52.1	
sdc <- 9.7	
sdt <- 9.6	
d2 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
A <- (nc+nt)^2 / (nc*nt)
r2 <- d2 / (sqrt(d2^2 + A)) # transform to r seperately, because the sample size is different in each case

nc <- 114
nt <- 51
mc <- 50.4	
mt <- 49.9	
sdc <- 9.7	
sdt <- 10.2
d3 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
A <- (nc+nt)^2 / (nc*nt)
r3 <- d3 / (sqrt(d3^2 + A)) # transform to r seperately, because the sample size is different in each case

r <- c(r1,r2,r3)
eff121 <- mean(r)
df$effestnew[df$id %in% "121"] <- eff121


# STUDY 125 - Card - Lester
nc <- 104
nt <- 45
mc <- 46.81
mt <- 44.43	
sdc <- 10.64	
sdt <- 9.81
d1 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
A <- (nc+nt)^2 / (nc*nt)
r1 <- d1 / (sqrt(d1^2 + A)) # transform to r seperately, because the sample size is different in each case

nc <- 83
nt <- 40
mc <- 44.08	
mt <- 48.45	
sdc <- 8.77	
sdt <- 11.59
d2 <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
A <- (nc+nt)^2 / (nc*nt)
r2 <- d2 / (sqrt(d2^2 + A)) # transform to r seperately, because the sample size is different in each case

r <- c(r1,r2)
eff125 <- mean(r)
df$effestnew[df$id %in% "125"] <- eff125


# STUDY 140 - Crook - Ployhart
r <- c(0,-.01,-.03,-.04,-.04,-.02,-.02,-.02,.02,.02,0,.01,.05,.05,.07,.06,.06,.08,.07,.04,.01,.02,0,.04,.05,.02,.07)
eff140 <- mean(r)
df$effestnew[df$id %in% "140"] <- eff140


# STUDY 146 - Crook - Van Iddekinge
r <- c(.12,.13,.06,.19,.17,.19,.05,.16,.15,.21,.10,.19,.19,.20,.12,.21,.18,.19,.08,.23,.16,.23,.15,.20)
eff146 <- mean(r)
df$effestnew[df$id %in% "146"] <- eff146

# STUDY 211 - Fischer - Howard
g <- 0    # because OR = 1 so no difference
eff211 <- g
df$effestnew[df$id %in% "211"] <- eff211

# STUDY 214 - Fischer - Karakashian
OR <- 2.86
nc <- nt <- 41.5
dfs <- nc + nt -2
d <- log(OR) * (sqrt(3)/pi)
J <- 1 - (3 / ((4 * dfs) - 1))
eff214 <- d * J
df$effestnew[df$id %in% "214"] <- eff214

# STUDY 216 - Fischer - Lewis
nc <- 7
nt <- 8
dfs <- nc + nt - 2

mc <- 3.93
mt <- 38.92	
sdc <- 10.13	
sdt <- 45.82
d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
J <- 1 - (3 / ((4 * dfs) - 1))
g1 <- d * J # calculate effects separately because different Ns

nc <- 49
nt <- 134
dfs <- nc + nt - 2

OR.a <- 0.2985
OR.b <- 0.7015
OR.c <- 0.2276
OR.d <- 0.7724

OR <- (OR.a * OR.d) / (OR.b * OR.c)
d <- log(OR) * ((sqrt(3)) / pi)
J <- 1 - (3 / ((4 * dfs) - 1))
g2 <- d * J # calculate effects separately because different Ns

g <- c(g1,g2)
eff216 <- mean(g)
df$effestnew[df$id %in% "216"] <- eff216
			

# STUDY 230 - Fox - Karahasanovic
mc <- 120
mt <- 138
sdc <- 25
sdt <- 51
nc <- 7
nt <- 9 
n <- nc+nt

A <- (nc + nt)^2 / (nc*nt)
d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
r <- d / (sqrt((d^2) + A))
a <- sqrt(0.25 / ((nc/n) * (nt/n)))

r1 <- (a * r) / (sqrt(((a^2 - 1) * r^2 + 1)))

OR.a <- 0.67
OR.b <- 0.33
OR.c <- 0.57
OR.d <- 0.43

OR <- (OR.a * OR.d) / (OR.b * OR.c)
d <- log(OR) * ((sqrt(3)) / pi)
r <- d / (sqrt((d^2) + A))
a <- sqrt(0.25 / ((nc/n) * (nt/n)))
r2 <- (a * r) / (sqrt(((a^2 - 1) * r^2 + 1))) 

eff230 <- mean(c(r1,r2)) # final estimate (r)
df$effestnew[df$id %in% "230"] <- eff230

# STUDY 353 - Morgan - Nelson
eff353 <- .75   # in text
df$effestnew[df$id %in% "353"] <- eff353

# STUDY 384 - Smith - Carlson
eff384 <- ((462 * 0.27) + (223 * 0.39)) / 685
df$effestnew[df$id %in% "384"] <- eff384


# Study 394 - Smith - Roberts
r1 <- c(.21,.13,.14,.14,.04,.07) # inverting last two correlations, loneliness and depression
r2 <- c(.2,.12,.14,.1,.08,.01) # inverting last two correlations, loneliness and depression
n1 <- 1237
n2 <- 755
eff394 <- ((mean(r1) * n1) + (mean(r2) * n2)) / (n1+n2)
df$effestnew[df$id %in% "394"] <- eff394


# Study 397 - Smith - Umana-Taylor
r1 <- .15
n1 <- 196
r2 <- .19
n2 <- 137
r3 <- .15
n3 <- 729

eff397 <- ((r1 * n1) + (r2 * n2) + (r3 * n3))  / (n1+n2+n3)
df$effestnew[df$id %in% "397"] <- eff397


# STUDY 426 - Toosi - Dovidio
t <- c(0.4,1.65,1.07,2.42)
n <- c(15,20,18,17)
r <- sqrt((t^2) / ((t^2) + (n - 2))) 

eff426 <- mean(r)
df$effestnew[df$id %in% "426"] <- eff426


# MA1: Adesope ------------------------------------------------------------
df_adesope <- df[df$meta %in% "Adesope" & df$input %in% "smd",]
adesope_smd <- smd_to_g1(x = df_adesope)

df_adesope <- df[df$meta %in% "Adesope" & df$input %in% "OR",]
adesope_or <- or_to_g1(x = df_adesope)

df_adesope <- df[df$meta %in% "Adesope" & df$input %in% "F",]
adesope_f <- f_to_g1(x = df_adesope)

df$effestnew[df$meta %in% "Adesope" & df$input %in% "smd"] <- adesope_smd
df$effestnew[df$meta %in% "Adesope" & df$input %in% "OR"] <- adesope_or
df$effestnew[df$meta %in% "Adesope" & df$input %in% "F"] <- adesope_f


# MA2: Alfieri ------------------------------------------------------------
df_alfieri <- df[df$meta %in% "Alfieri" & df$input %in% "smd",]
alfieri_smd <- smd_to_g2(x = df_alfieri)

df_alfieri <- df[df$meta %in% "Alfieri" & df$input %in% "F",]
alfieri_f <- f_to_g2(x = df_alfieri)
alfieri_f[1] <- alfieri_f[1]*-1 # effect size reversed to be in right direction

df_alfieri <- df[df$meta %in% "Alfieri" & df$input %in% "chisquare",]
alfieri_chi <- chi_to_g2(x = df_alfieri)
alfieri_chi[1] <- alfieri_chi[1]*-1 # effect size reversed to right direction 

df_alfieri <- df[df$meta %in% "Alfieri" & df$input %in% "OR",]
alfieri_or <- or_to_g2(x = df_alfieri)
alfieri_or[1] <- alfieri_or[1]*-1 # effect size reversed to right direction 

df_alfieri <- df[df$meta %in% "Alfieri" & df$input %in% "t",]
alfieri_t <- t_to_g2(x = df_alfieri)
alfieri_t[1] <- alfieri_t[1]*-1 # effect size reversed to right direction 

df$effestnew[df$meta %in% "Alfieri" & df$input %in% "0"] <- 0
df$effestnew[df$meta %in% "Alfieri" & df$input %in% "NA"] <- df$effest[df$meta %in% "Alfieri" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Alfieri" & df$input %in% "smd"] <- alfieri_smd
df$effestnew[df$meta %in% "Alfieri" & df$input %in% "F"] <- alfieri_f
df$effestnew[df$meta %in% "Alfieri" & df$input %in% "chisquare"] <- alfieri_chi
df$effestnew[df$meta %in% "Alfieri" & df$input %in% "OR"] <- alfieri_or
df$effestnew[df$meta %in% "Alfieri" & df$input %in% "t"] <- alfieri_t


# MA3: Babbage ------------------------------------------------------------
df_babbage <- df[df$meta %in% "Babbage" & df$input %in% "smd",]
babbage_smd <- smd_to_g1(x = df_babbage)

df$effestnew[df$meta %in% "Babbage" & df$input %in% "NA"] <- df$effest[df$meta %in% "Babbage" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Babbage" & df$input %in% "smd"] <- babbage_smd


# MA4: Balliet ------------------------------------------------------------
df_balliet <- df[df$meta %in% "Balliet" & df$input %in% "r",]
balliet_r <- r_to_g1(x = df_balliet)

df_balliet <- df[df$meta %in% "Balliet" & df$input %in% "t",]
balliet_t <- t_to_g1(x = df_balliet)
balliet_t[1] <- balliet_t[1]*-1 # effect size reversed to right direction 
balliet_t[2] <- balliet_t[2]*-1 # effect size reversed to right direction 

df_balliet <- df[df$meta %in% "Balliet" & df$input %in% "OR",]
balliet_or <- or_to_g1(x = df_balliet)

df_balliet <- df[df$meta %in% "Balliet" & df$input %in% "F",]
balliet_f <- f_to_g1(x = df_balliet)
balliet_f[3] <- balliet_f[3]*-1 # effect size reversed to right direction 
balliet_f[4] <- balliet_f[4]*-1 # effect size reversed to right direction 

df_balliet <- df[df$meta %in% "Balliet" & df$input %in% "smd",]
balliet_smd <- smd_to_g1(x = df_balliet)

df$effestnew[df$meta %in% "Balliet" & df$input %in% "0"] <- 0
df$effestnew[df$meta %in% "Balliet" & df$input %in% "NA"] <- df$effest[df$meta %in% "Balliet" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Balliet" & df$input %in% "r"] <- balliet_r
df$effestnew[df$meta %in% "Balliet" & df$input %in% "t"] <- balliet_t
df$effestnew[df$meta %in% "Balliet" & df$input %in% "OR"] <- balliet_or
df$effestnew[df$meta %in% "Balliet" & df$input %in% "F"] <- balliet_f
df$effestnew[df$meta %in% "Balliet" & df$input %in% "smd"] <- balliet_smd


# MA5: Benish -------------------------------------------------------------
df_benish <- df[df$meta %in% "Benish" & df$input %in% "smd",]
benish_smd <- smd_to_g1(x = df_benish)
benish_smd[2] <- benish_smd[2]*-1 # effect size reversed to right direction 
benish_smd[3] <- benish_smd[3]*-1 # effect size reversed to right direction 

df_benish <- df[df$meta %in% "Benish" & df$input %in% "benish.agg.smd",]
benish_agg_smd <- agg_smd_to_g1_benish(x = df_benish)
benish_agg_smd <- benish_agg_smd*-1 # effect size reversed to right direction

df_benish <- df[df$meta %in% "Benish" & df$input %in% "benish.agg.f",]
benish_agg_f <- agg_f_to_g1_benish(x = df_benish)

df$effestnew[df$meta %in% "Benish" & df$input %in% "NA"] <- df$effest[df$meta %in% "Benish" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Benish" & df$input %in% "smd"] <- benish_smd
df$effestnew[df$meta %in% "Benish" & df$input %in% "benish.agg.smd"] <- benish_agg_smd
df$effestnew[df$meta %in% "Benish" & df$input %in% "benish.agg.f"] <- benish_agg_f


# MA6: Berry1 -------------------------------------------------------------
df_berry1 <- df[df$meta %in% "Berry 1" & df$input %in% "berry1.correction",]
berry1_corr <- r_corr(x = df_berry1)

df$effestnew[df$meta %in% "Berry 1" & df$input %in% "berry1.correction"] <- berry1_corr



# MA7: Berry2 -------------------------------------------------------------
df_berry2 <- df[df$meta %in% "Berry 2" & df$input %in% "r",]
berry2_r <- r_to_r(x = df_berry2)

df_berry2 <- df[df$meta %in% "Berry 2" & df$input %in% "berry2.agg.r",]
berry2_agg <- agg_r_to_r(x = df_berry2)

df$effestnew[df$meta %in% "Berry 2" & df$input %in% "r"] <- berry2_r
df$effestnew[df$meta %in% "Berry 2" & df$input %in% "berry2.agg.r"] <- berry2_agg


# MA8: Card ---------------------------------------------------------------
df_card <- df[df$meta %in% "Card" & df$input %in% "smd",]
card_smd <- smd_to_r(x = df_card)

df$effestnew[df$meta %in% "Card" & df$input %in% "NA"] <- df$effest[df$meta %in% "Card" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Card" & df$input %in% "smd"] <- card_smd



# MA9: Crook --------------------------------------------------------------
df_crook <- df[df$meta %in% "Crook" & df$input %in% "r",]
crook_r <- r_to_r(x = df_crook)

df$effestnew[df$meta %in% "Crook" & df$input %in% "r"] <- crook_r



# MA10: de Wit ------------------------------------------------------------
df_dewit <- df[df$meta %in% "de Wit" & df$input %in% "r",]
dewit_r <- r_to_r(x = df_dewit)

df_dewit <- df[df$meta %in% "de Wit" & df$input %in% "dewit.agg.r",]
dewit_agg <- agg_r_to_r(x = df_dewit)

df_dewit <- df[df$meta %in% "de Wit" & df$input %in% "dewit.correction",]
dewit_corr <- r_corr(x = df_dewit)


df$effestnew[df$meta %in% "de Wit" & df$input %in% "0"] <- 0
df$effestnew[df$meta %in% "de Wit" & df$input %in% "r"] <- dewit_r
df$effestnew[df$meta %in% "de Wit" & df$input %in% "dewit.agg.r"] <- dewit_agg
df$effestnew[df$meta %in% "de Wit" & df$input %in% "dewit.correction"] <- dewit_corr




# MA11: Else-quest --------------------------------------------------------
df_elseq <- df[df$meta %in% "Else-quest" & df$input %in% "r.elseq",]
elseq_r <- r_to_d_elseq(x = df_elseq)

df_elseq <- df[df$meta %in% "Else-quest" & df$input %in% "chisquare",]
elseq_chi <- chi_to_d(x = df_elseq)

df_elseq <- df[df$meta %in% "Else-quest" & df$input %in% "smd",]
elseq_smd <- smd_to_d(x = df_elseq)

df_elseq <- df[df$meta %in% "Else-quest" & df$input %in% "F.elseq",]
elseq_f <- f_to_d_elseq(x = df_elseq)

df$effestnew[df$meta %in% "Else-quest" & df$input %in% "NA"] <- df$effest[df$meta %in% "Else-quest" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Else-quest" & df$input %in% "r.elseq"] <- elseq_r
df$effestnew[df$meta %in% "Else-quest" & df$input %in% "chisquare"] <- elseq_chi
df$effestnew[df$meta %in% "Else-quest" & df$input %in% "smd"] <- elseq_smd
df$effestnew[df$meta %in% "Else-quest" & df$input %in% "F.elseq"] <- elseq_f


# MA12: Farber ------------------------------------------------------------
df_farber <- df[df$meta %in% "Farber" & df$input %in% "r",]
farber_r <- r_to_r(x = df_farber)

df_farber <- df[df$meta %in% "Farber" & df$input %in% "farber.agg.r",]
farber_agg <- agg_r_to_r(x = df_farber)

df_farber <- df[df$meta %in% "Farber" & df$input %in% "chisquare",]
farber_chi <- chi_to_r(x = df_farber)

df_farber <- df[df$meta %in% "Farber" & df$input %in% "t",]
farber_t<- t_to_r2(x = df_farber)

df$effestnew[df$meta %in% "Farber" & df$input %in% "r"] <- farber_r
df$effestnew[df$meta %in% "Farber" & df$input %in% "farber.agg.r"] <- farber_agg
df$effestnew[df$meta %in% "Farber" & df$input %in% "chisquare"] <- farber_chi
df$effestnew[df$meta %in% "Farber" & df$input %in% "t"] <- farber_t

# MA13: Fischer -----------------------------------------------------------
df_fischer <- df[df$meta %in% "Fischer" & df$input %in% "OR",]
fischer_or <- or_to_g2(x = df_fischer)
fischer_or[8] <- fischer_or[8]*-1 # effect reversed to right direction

df_fischer <- df[df$meta %in% "Fischer" & df$input %in% "smd",]
fischer_smd <- smd_to_g2(x = df_fischer)

df_fischer <- df[df$meta %in% "Fischer" & df$input %in% "Zstatistic",]
fischer_zstat <- zstat_to_g2(x = df_fischer)
fischer_zstat[1] <- fischer_zstat[1]*-1 # effect reversed to right direction

df$effestnew[df$meta %in% "Fischer" & df$input %in% "OR"] <- fischer_or
df$effestnew[df$meta %in% "Fischer" & df$input %in% "smd"] <- fischer_smd
df$effestnew[df$meta %in% "Fischer" & df$input %in% "Zstatistic"] <- fischer_zstat


# MA14: Fox ---------------------------------------------------------------
df_fox <- df[df$meta %in% "Fox" & df$input %in% "F",]
fox_f <- f_to_r1(x = df_fox)
fox_f[2] <- fox_f[2]*-1 # effect size reversed to right direction 

df_fox <- df[df$meta %in% "Fox" & df$input %in% "fox.correction.smd",]
# this one corrects for unequal sampling
fox_smd_corr <- smd_to_r_fox(x = df_fox) 
fox_smd_corr[3] <- fox_smd_corr[3]*-1 # effect reversed to right direction

df_fox <- df[df$meta %in% "Fox" & df$input %in% "smd",]
# this one does not correct for unequal sampling (because the samples are equal)
fox_smd <- smd_to_r(x = df_fox)
fox_smd[9] <- fox_smd[9]*-1 # effect size reversed to right direction 

df_fox <- df[df$meta %in% "Fox" & df$input %in% "fox.correction.z",]
fox_z_corr <- zstat_to_r_fox(x = df_fox) 

df_fox <- df[df$meta %in% "Fox" & df$input %in% "OR",]
fox_or <- or_to_r(x = df_fox)

df_fox <- df[df$meta %in% "Fox" & df$input %in% "mann-whitney-u",]
fox_mann <- mann_to_r(x = df_fox)

df$effestnew[df$meta %in% "Fox" & df$input %in% "NA"] <- df$effest[df$meta %in% "Fox" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Fox" & df$input %in% "F"] <- fox_f
df$effestnew[df$meta %in% "Fox" & df$input %in% "OR"] <- fox_or
df$effestnew[df$meta %in% "Fox" & df$input %in% "fox.correction.smd"] <- fox_smd_corr
df$effestnew[df$meta %in% "Fox" & df$input %in% "smd"] <- fox_smd
df$effestnew[df$meta %in% "Fox" & df$input %in% "fox.correction.z"] <- fox_z_corr
df$effestnew[df$meta %in% "Fox" & df$input %in% "mann-whitney-u"] <- fox_mann

# MA15: Freund ------------------------------------------------------------
df_freund <- df[df$meta %in% "Freund" & df$input %in% "r",]
freund_r <- r_to_r(x = df_freund)

df$effestnew[df$meta %in% "Freund" & df$input %in% "r"] <- freund_r


# MA16: Green -------------------------------------------------------------
df_green <- df[df$meta %in% "Green" & df$input %in% "smd",]
green_smd <- smd_to_d(x = df_green)
green_smd[1] <- green_smd[1]*-1 # effect size reversed to right direction 

df$effestnew[df$meta %in% "Green" & df$input %in% "NA"] <- df$effest[df$meta %in% "Green" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Green" & df$input %in% "smd"] <- green_smd


# MA17: Hallion -----------------------------------------------------------
df_hallion <- df[df$meta %in% "Hallion" & df$input %in% "smd",]
hallion_smd <- smd_to_g3(x = df_hallion)
hallion_smd[1] <- hallion_smd[1]*-1 # effect size reversed to right direction 
hallion_smd[2] <- hallion_smd[2]*-1 # effect size reversed to right direction 
hallion_smd[3] <- hallion_smd[3]*-1 # effect size reversed to right direction 
hallion_smd[5] <- hallion_smd[5]*-1 # effect size reversed to right direction 
hallion_smd[7] <- hallion_smd[7]*-1 # effect size reversed to right direction 
hallion_smd[8] <- hallion_smd[8]*-1 # effect size reversed to right direction 

df$effestnew[df$meta %in% "Hallion" & df$input %in% "NA"] <- df$effest[df$meta %in% "Hallion" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Hallion" & df$input %in% "smd"] <- hallion_smd


# MA18: Ihle --------------------------------------------------------------
df_ihle <- df[df$meta %in% "Ihle" & df$input %in% "smd",]
ihle_smd <- smd_to_g2(x = df_ihle)

df_ihle <- df[df$meta %in% "Ihle" & df$input %in% "OR",]
ihle_or <- or_to_g2(x = df_ihle)

df$effestnew[df$meta %in% "Ihle" & df$input %in% "NA"] <- df$effest[df$meta %in% "Ihle" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Ihle" & df$input %in% "smd"] <- ihle_smd
df$effestnew[df$meta %in% "Ihle" & df$input %in% "OR"] <- ihle_or


# MA19: Koenig ------------------------------------------------------------
df_koenig <- df[df$meta %in% "Koenig" & df$input %in% "smd",]
koenig_smd <- smd_to_g2(x = df_koenig)

df_koenig <- df[df$meta %in% "Koenig" & df$input %in% "OR",]
koenig_or <- or_to_g2_koenig(x = df_koenig)

df_koenig <- df[df$meta %in% "Koenig" & df$input %in% "t",]
koenig_t <- t_to_g2(x = df_koenig)

df$effestnew[df$meta %in% "Koenig" & df$input %in% "NA"] <- df$effest[df$meta %in% "Koenig" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Koenig" & df$input %in% "smd"] <- koenig_smd
df$effestnew[df$meta %in% "Koenig" & df$input %in% "OR"] <- koenig_or
df$effestnew[df$meta %in% "Koenig" & df$input %in% "t"] <- koenig_t


# MA20: Kolden ------------------------------------------------------------
df_kolden <- df[df$meta %in% "Kolden" & df$input %in% "kolden.agg.r",]
kolden_agg <- agg_d_to_r(x = df_kolden)

df_kolden <- df[df$meta %in% "Kolden" & df$input %in% "r",]
kolden_r <- r_to_r(x = df_kolden)

df$effestnew[df$meta %in% "Kolden" & df$input %in% "kolden.agg.r"] <- kolden_agg
df$effestnew[df$meta %in% "Kolden" & df$input %in% "r"] <- kolden_r


# MA21: Lucassen ----------------------------------------------------------
df_lucassen <- df[df$meta %in% "Lucassen" & df$input %in% "r",]
lucassen_r <- r_to_r(x = df_lucassen)

df$effestnew[df$meta %in% "Lucassen" & df$input %in% "0"] <- 0
df$effestnew[df$meta %in% "Lucassen" & df$input %in% "NA"] <- df$effest[df$meta %in% "Lucassen" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Lucassen" & df$input %in% "r"] <- lucassen_r


# MA22: Mol ---------------------------------------------------------------
df_mol <- df[df$meta %in% "Mol" & df$input %in% "r",]
mol_r <- r_to_z(x = df_mol)

df$effestnew[df$meta %in% "Mol" & df$input %in% "r"] <- mol_r


# MA23: Morgan ------------------------------------------------------------
df_morgan <- df[df$meta %in% "Morgan" & df$input %in% "smd",]
morgan_smd <- smd_to_d(x = df_morgan)
morgan_smd <- morgan_smd*-1 # effect size reversed to right direction 

df_morgan <- df[df$meta %in% "Morgan" & df$input %in% "smd.paired",]
morgan_smdpaired <- smdpaired_to_d(x = df_morgan)

df_morgan <- df[df$meta %in% "Morgan" & df$input %in% "OR",]
morgan_or <- or_to_d(x = df_morgan)

df$effestnew[df$meta %in% "Morgan" & df$input %in% "NA"] <- df$effest[df$meta %in% "Morgan" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Morgan" & df$input %in% "smd"] <- morgan_smd
df$effestnew[df$meta %in% "Morgan" & df$input %in% "smd.paired"] <- morgan_smdpaired
df$effestnew[df$meta %in% "Morgan" & df$input %in% "OR"] <- morgan_or

# MA24: Munder ------------------------------------------------------------
df_munder <- df[df$meta %in% "Munder" & df$input %in% "smd",]
munder_smd <- smd_to_g1(x = df_munder)
munder_smd <- munder_smd*-1 # all effects reversed to be in right direction

df$effestnew[df$meta %in% "Munder" & df$input %in% "NA"] <- df$effest[df$meta %in% "Munder" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Munder" & df$input %in% "smd"] <- munder_smd



# MA25: Piet --------------------------------------------------------------
df_piet <- df[df$meta %in% "Piet" & df$input %in% "smd",]
piet_smd <- smd_to_g2(x = df_piet)
piet_smd[3] <- piet_smd[3]*-1 # effect size reversed to right direction 
piet_smd[4] <- piet_smd[4]*-1 # effect size reversed to right direction 
piet_smd[5] <- piet_smd[5]*-1 # effect size reversed to right direction 

df$effestnew[df$meta %in% "Piet" & df$input %in% "NA"] <- df$effest[df$meta %in% "Piet" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Piet" & df$input %in% "smd"] <- piet_smd


# MA26: Smith -------------------------------------------------------------
df_smith <- df[df$meta %in% "Smith" & df$input %in% "r",]
smith_r <- r_to_r(x = df_smith)

df_smith <- df[df$meta %in% "Smith" & df$input %in% "smd",]
smith_smd <- smd_to_r(x = df_smith)

df$effestnew[df$meta %in% "Smith" & df$input %in% "r"] <- smith_r
df$effestnew[df$meta %in% "Smith" & df$input %in% "smd"] <- smith_smd


# MA27: Tillman -----------------------------------------------------------
df_tillman <- df[df$meta %in% "Tillman" & df$input %in% "r",]
tillman_r <- r_to_r(x = df_tillman)

df$effestnew[df$meta %in% "Tillman" & df$input %in% "NA"] <- df$effest[df$meta %in% "Tillman" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Tillman" & df$input %in% "r"] <- tillman_r


# MA28: Toosi -------------------------------------------------------------
df_toosi <- df[df$meta %in% "Toosi" & df$input %in% "smd",]
toosi_smd <- smd_to_r(x = df_toosi)
toosi_smd[3] <- toosi_smd[3]*-1 # effect size reversed to right direction 

df_toosi <- df[df$meta %in% "Toosi" & df$input %in% "F",]
toosi_f <- f_to_r2(x = df_toosi)

df_toosi <- df[df$meta %in% "Toosi" & df$input %in% "t",]
toosi_t <- t_to_r2(x = df_toosi)

df$effestnew[df$meta %in% "Toosi" & df$input %in% "0"] <- 0
df$effestnew[df$meta %in% "Toosi" & df$input %in% "NA"] <- df$effest[df$meta %in% "Toosi" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Toosi" & df$input %in% "smd"] <- toosi_smd
df$effestnew[df$meta %in% "Toosi" & df$input %in% "F"] <- toosi_f
df$effestnew[df$meta %in% "Toosi" & df$input %in% "t"] <- toosi_t

# MA29: van Iddekinge -----------------------------------------------------
df_vanidde <- df[df$meta %in% "van Iddekinge" & df$input %in% "r",]
vanidde_r <- r_to_r(x = df_vanidde)

df$effestnew[df$meta %in% "van Iddekinge" & df$input %in% "r"] <- vanidde_r


# MA30: Webb --------------------------------------------------------------
df_webb <- df[df$meta %in% "Webb" & df$input %in% "smd",]
webb_smd <- smd_to_d(x = df_webb)
webb_smd[1] <- webb_smd[1]*-1 # effect size reversed to right direction 
webb_smd[2] <- webb_smd[2]*-1 # effect size reversed to right direction 
webb_smd[3] <- webb_smd[3]*-1 # effect size reversed to right direction 
webb_smd[4] <- webb_smd[4]*-1 # effect size reversed to right direction 
webb_smd[5] <- webb_smd[5]*-1 # effect size reversed to right direction 

df_webb <- df[df$meta %in% "Webb" & df$input %in% "smd.webb",]
webb_smd_webb <- smd_to_d_webb(x = df_webb)

df_webb <- df[df$meta %in% "Webb" & df$input %in% "t",]
webb_t <- t_to_d(x = df_webb)

df$effestnew[df$meta %in% "Webb" & df$input %in% "0"] <- 0
df$effestnew[df$meta %in% "Webb" & df$input %in% "NA"] <- df$effest[df$meta %in% "Webb" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Webb" & df$input %in% "smd"] <- webb_smd
df$effestnew[df$meta %in% "Webb" & df$input %in% "smd.webb"] <- webb_smd_webb
df$effestnew[df$meta %in% "Webb" & df$input %in% "t"] <- webb_t



# MA31: Woodin ------------------------------------------------------------
df_woodin <- df[df$meta %in% "Woodin" & df$input %in% "r",]
woodin_r <- r_to_g1(x = df_woodin)

df_woodin <- df[df$meta %in% "Woodin" & df$input %in% "etasquare",]
woodin_etasq <- etasq_to_g1(x = df_woodin)
woodin_etasq <- woodin_etasq*-1 # effect size reversed to right direction 

df_woodin <- df[df$meta %in% "Woodin" & df$input %in% "smd",]
woodin_smd <- smd_to_g1(x = df_woodin)

df$effestnew[df$meta %in% "Woodin" & df$input %in% "r"] <- woodin_r
df$effestnew[df$meta %in% "Woodin" & df$input %in% "etasquare"] <- woodin_etasq
df$effestnew[df$meta %in% "Woodin" & df$input %in% "smd"] <- woodin_smd

# MA32: Woodley -----------------------------------------------------------
df_woodley <- df[df$meta %in% "Woodley" & df$input %in% "r",]
woodley_r <- r_to_r(x = df_woodley)

df$effestnew[df$meta %in% "Woodley" & df$input %in% "NA"] <- df$effest[df$meta %in% "Woodley" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Woodley" & df$input %in% "r"] <- woodley_r


# MA33: Yoon --------------------------------------------------------------
df_yoon <- df[df$meta %in% "Yoon" & df$input %in% "r",]
yoon_r <- r_to_r(x = df_yoon)

df$effestnew[df$meta %in% "Yoon" & df$input %in% "NA"] <- df$effest[df$meta %in% "Yoon" & df$input %in% "NA"]
df$effestnew[df$meta %in% "Yoon" & df$input %in% "r"] <- yoon_r

# Discrepancies -----------------------------------------------------------

# Effect size discrepancies -----------------------------------------------
r.disc <- c(.025, .076, .126)
z.disc <- 0.5 * log((1 + r.disc) / (1 - r.disc))
d.disc <- (2 * r.disc) / sqrt(1 - r.disc^2)
J.disc <- (1 - 3 / (4 * (64 - 2) - 1)) # Assuming N = 64
g.disc <- d.disc * J.disc


# Fill in discrepancies ---------------------------------------------------
df$disc.eff <- df$effest - df$effestnew
df$disc.n <- df$n - df$nnew


# Fill in discrepancy category 
for (i in 1:nrow(df)) {
  
  if (df$efftype[i] == "g") {
    
    if (abs(df$disc.eff[i]) < g.disc[1]) {
      df$disccat.eff[i] = 0
    } else if (abs(df$disc.eff[i]) >= g.disc[1] & abs(df$disc.eff[i]) < g.disc[2]) {
      df$disccat.eff[i] = 1
    } else if (abs(df$disc.eff[i]) >= g.disc[2] & abs(df$disc.eff[i]) < g.disc[3]) {
      df$disccat.eff[i] = 2
    } else if (abs(df$disc.eff[i]) >= g.disc[3]) {
      df$disccat.eff[i] = 3
    } else {
      df$disccat.eff[i] = "check"
      
    }
  }
  
  if (df$efftype[i] == "d") {
  
  if (abs(df$disc.eff[i]) < d.disc[1]) {
    df$disccat.eff[i] = 0
  } else if (abs(df$disc.eff[i]) >= d.disc[1] & abs(df$disc.eff[i]) < d.disc[2]) {
    df$disccat.eff[i] = 1
  } else if (abs(df$disc.eff[i]) >= d.disc[2]& abs(df$disc.eff[i]) < d.disc[3]) {
    df$disccat.eff[i] = 2
  } else if (abs(df$disc.eff[i]) >= d.disc[3]) {
    df$disccat.eff[i] = 3
  } else {
    df$disccat.eff[i] = "check"
  }
 }

if (df$efftype[i] == "r" | df$efftype[i] == "z") {
    
    if (abs(df$disc.eff[i]) < r.disc[1]) {
      df$disccat.eff[i] = 0
    } else if (abs(df$disc.eff[i]) >= r.disc[1] & abs(df$disc.eff[i]) < r.disc[2]) {
      df$disccat.eff[i] = 1
    } else if (abs(df$disc.eff[i]) >= r.disc[2] & abs(df$disc.eff[i]) < r.disc[3]) {
      df$disccat.eff[i] = 2
    } else if (abs(df$disc.eff[i]) >= r.disc[3]) {
      df$disccat.eff[i] = 3
    } else {
      df$disccat.eff[i] = "check"
    }
  }
} 

# Reverse effect sizes so all effects are in hypothesized direction -------
df$effest.exp <- df$effest
df$effestnew.exp <- df$effestnew

reverse <- c("Babbage","de Wit","Fischer","Woodin","Yoon")

df$effest.exp[which(df$meta %in% reverse)] <- df$effest.exp[which(df$meta %in% reverse)]*-1
df$effestnew.exp[which(df$meta %in% reverse)] <- df$effestnew.exp[which(df$meta %in% reverse)]*-1

# Save file in xlsx -------------------------------------------------------
write_xlsx(df,"../codebooks/codebook-primary-studies-complete.xlsx",col_names=T)

# Save file in csv --------------------------------------------------------
write.table(df, file = "../codebooks/codebook-primary-studies-complete.csv", row.names=F, col.names=T, sep=' ')
