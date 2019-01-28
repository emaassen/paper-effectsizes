rm(list = ls())
library(metafor)

# Run all meta-analyses 
setwd("C:/Users/s421506/tiu/research/effectsizes/data-per-ma")
authors_complete <- list.files("C:/Users/s421506/tiu/research/effectsizes/data-per-ma", pattern="_complete.csv")
effect_size <- c("g","d","g","d","d","r","r","r","r","r","d","r","g","r","r","d","g","d","g","r","r","z","d","g","g","r","r","r","r","d","r","r","r")

class_primary_out <- matrix(NA, ncol=5)
class_primary_nout <- matrix(NA, ncol=5)

for (i in 1:length(authors_complete)) {
  
  # Complete original MAs ---------------------------------------------------
  df <- read.table(paste(authors_complete[i]), header=T, sep=';')
  
  if (effect_size[i] == "g") {
    
    res <- rma(g, vg, data=df)                                            # random effects meta-analysis
    l1o <- leave1out(res)                                                 # leave1out analysis
    df$q <- res$QE - l1o$Q                                                # change in Q with l1o
                                   
    nout <- df[which(df$q < 3.84), ]                                      # subset of studies that are included as non-outliers
    out = tryCatch ({                                                     # subset of studies that are included as outliers
      out <- df[which(df$q >= 3.84), ] },
                       error = function() { 
                         out <- c() })
    
    
  } else if (effect_size[i] == "d") {
    
    res <- rma(d, vd, data=df)                                            # random effects meta-analysis
    l1o <- leave1out(res)                                                 # leave1out analysis
    df$q <- res$QE - l1o$Q                                                # change in Q with l1o

    nout <- df[which(df$q < 3.84), ]                                      # subset of studies that are included as non-outliers
    out = tryCatch ({                                                     # subset of studies that are included as outliers
      out <- df[which(df$q >= 3.84), ] },
      error = function() { 
        out <- c() })
    
    
  } else if (effect_size[i] == "r") {
    
    z <- 0.5 * log((1 + df$r) / (1 - df$r)) 
    vz <- 1 / (df$n - 3) 
    
    res <- rma(z, vz, data=df)                                            # random effects meta-analysis
    l1o <- leave1out(res)                                                 # leave1out analysis
    df$q <- res$QE - l1o$Q                                                # change in Q with l1o

    nout <- df[which(df$q < 3.84), ]                                      # subset of studies that are included as non-outliers
    out = tryCatch ({                                                     # subset of studies that are included as outliers
      out <- df[which(df$q >= 3.84), ] },
      error = function() { 
        out <- c() })
    
  } else if (effect_size[i] == "z") {
    
    res <- rma(z, vz, data=df)                                            # random effects meta-analysis
    l1o <- leave1out(res)                                                 # leave1out analysis
    df$q <- res$QE - l1o$Q                                                # change in Q with l1o

    nout <- df[which(df$q < 3.84), ]                                      # subset of studies that are included as non-outliers
    out = tryCatch ({                                                     # subset of studies that are included as outliers
      out <- df[which(df$q >= 3.84), ] },
      error = function() { 
        out <- c() })
    
  }
  
  out <- tryCatch ({ 
    out <- df[which(df$q >= 3.84), ] },
    error = function() { 
      out <- c() })
  
  if (typeof(out) != "list") {
    
    out <- c()

  } else {
      
    out <- out[,c("study","n",effect_size[i])]
    out$meta <- rep(authors_complete[i],nrow(out))
    out$type <- rep("outlier",nrow(out))
    colnames(class_primary_out) <- colnames(out)
    class_primary_out <- rbind(class_primary_out,out)
    
    }
  
  nout <- nout[,c("study","n",effect_size[i])]
  nout$meta <- rep(authors_complete[i],nrow(nout))
  nout$type <- rep("non-outlier",nrow(nout))
  colnames(class_primary_nout) <- colnames(nout)
  class_primary_nout <- rbind(class_primary_nout,nout)
  
}

# remove first row
class_primary_nout <- class_primary_nout[-1,]
class_primary_out <- class_primary_out[-1,]

setwd("C:/Users/s421506/tiu/research/effectsizes/codebooks/")
write.table(class_primary_nout, file = "primary-study-non-outlier.csv", row.names=F, col.names=T, sep=';')
write.table(class_primary_out, file = "primary-study-outlier.csv", row.names=F, col.names=T, sep=';')


  
  