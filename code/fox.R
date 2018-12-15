rm(list = ls())

# Study 222
mc <- 1.7
mt <- 2.45
sdc <- 0.73
sdt <- 1.3

nc <- 20
nt <- 29

A <- (nc + nt)^2 / (nc*nt)
d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))

r <- d / (sqrt((d^2) + A));r

a <- sqrt(0.25/((nc*0.01)*(nt*0.01)))

(a*r)/(sqrt(((a^2-1)*(r^2))+1)) # final estimate (r)

# Study 225
mc <- 20.31
mt <- 17.12
sdc <- 5.79
sdt <- 5.57

nc <- 35
nt <- 42

A <- (nc + nt)^2 / (nc*nt)
d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))

r <- d / (sqrt((d^2) + A));r

a <- sqrt(0.25/((nc*0.01)*(nt*0.01)))

(a*r)/(sqrt(((a^2-1)*(r^2))+1)) # final estimate (r)


# Study 229
mc <- 368
mt <- 331
sdc <- 249
sdt <- 206

nc <- 24
nt <- 44

A <- (nc + nt)^2 / (nc*nt)
d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))

r <- d / (sqrt((d^2) + A));r

a <- sqrt(0.25/((nc*0.01)*(nt*0.01)))

(a*r)/(sqrt(((a^2-1)*(r^2))+1)) # final estimate (r)


#Study 230
mc <- 120
mt <- 138
sdc <- 25
sdt <- 51

nc <- 7
nt <- 9 

A <- (nc + nt)^2 / (nc*nt)
d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
r <- d / (sqrt((d^2) + A))
a <- sqrt(0.25/((nc*0.01)*(nt*0.01)))

r1 <- (a*r)/(sqrt(((a^2-1)*(r^2))+1)) 

OR.a <- 0.67
OR.b <- 0.33
OR.c <- 0.57
OR.d <- 0.43

OR <- (OR.a*OR.d)/(OR.b*OR.c)
d <- log(OR) * ((sqrt(3))/pi)
r <- d / (sqrt((d^2) + A))
a <- sqrt(0.25/((nc*0.01)*(nt*0.01)))

r2 <- (a*r)/(sqrt(((a^2-1)*(r^2))+1)) 

mean(c(r1,r2)) # final estimate (r)


# Study 234
Z <- -1.988
nc <- 10
nt <- 5
n <- nt+nc

r <- Z / sqrt(n)
a <- sqrt(0.25/((nc*0.01)*(nt*0.01)))

(a*r)/(sqrt(((a^2-1)*(r^2))+1)) # Final estimate r

#Study 235
mc <- 7.17
mt <- 9.69
sdc <- 3.01
sdt <- 2.56

nc <- 13
nt <- 11

A <- (nc + nt)^2 / (nc*nt)
d <- (mt - mc) / sqrt((((nc - 1) * (sdc^2)) + ((nt - 1) * (sdt^2))) / (nc + nt - 2))
r <- d / (sqrt((d^2) + A))
a <- sqrt(0.25/((nc*0.01)*(nt*0.01)))

(a*r)/(sqrt(((a^2-1)*(r^2))+1)) # Final estimate r
