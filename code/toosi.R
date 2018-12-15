rm(list = ls())

# Study 422
t <- 2.15
nc <- 7.5
nt <- 7.5

r1 <- sqrt((t^2)/((t^2)+(nc+nt-2))) # Rosenthal formula
r2 <- 0

r <- mean(c(r1,r2));r # final estimate

# Study 423
t <- 3.31
nc <- 9
nt <- 9

r1 <- sqrt((t^2)/((t^2)+(nc+nt-2))) # Rosenthal formula
r2 <- 0

r <- mean(c(r1,r2));r # final estimate

# Study 426
t <- c(0.4,1.65,1.07,2.42)
n <- c(15,20,18,17)
nc <- nt <- n/2

r <- sqrt((t^2)/((t^2)+(nc+nt-2))) # Rosenthal formula
mean(r) # final estimate


# Study 440
t <- c(2.79)
n <- c(15)
nc <- nt <- n/2

r <- sqrt((t^2)/((t^2)+(nc+nt-2)));r # Rosenthal formula
