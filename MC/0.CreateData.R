# Monte Carlo Experiment ----
setwd("~/Dropbox/work/Packages/Testing-ivhetLc")
# Load packages
rm(list = ls(all.names = TRUE))

# Globals
set.seed(26051986)
N <- c(100, 500, 1000, 5000)
S <- 1000                       # MC samples

## Experiment 1 ----
data     <- list(S)
Data     <- list(length(N))
prop     <- c(0.3, 0.7) # Proportion of individuals in each class
rho      <- c(0.5, 0.5) # Degree of the endogeneity 
sigma_e  <- c(1, 1)     # Variance of epsilon
sigma_v  <- c(1, 1)     # Variance of upsilon
sigma_z  <- c(3, 3)
Ez       <- 0
# Structural parameters
b0       <- c(-1, 1)
gamma    <- c(-1, 2)
# First-Equation parameters
d0       <- c(-1, 1)
d1       <- c(-1, 2)


for (n in N){
  cat("=================", "\n")
  cat("Sample size: ", n, "\n")
  cat("=================", "\n")
  
  # True coefficients in the linear model
  b0.v    <- c(rep(b0[1],    n * prop[1]),  rep(b0[2],    n * prop[2]))
  gamma.v <- c(rep(gamma[1], n * prop[1]),  rep(gamma[2], n * prop[2]))
  
  # True coefficients in the first equation
  d0.v    <- c(rep(d0[1], n * prop[1]), rep(d0[2], n * prop[2]))
  d1.v    <- c(rep(d1[1], n * prop[1]), rep(d1[2], n * prop[2]))
  
  count.n <- which(N %in% n)
  print(count.n)
  s <- 1
  while(s <= S){
    cat( "Monte Carlo Sample: ", s, "\n")
    
    # Generate upsilon for each class
    v1 <- rnorm(n * prop[1], 0, sigma_v[1])
    v2 <- rnorm(n * prop[2], 0, sigma_v[2])
    v  <- c(v1, v2)
    
    # Generate epsilon for each class
    e1 <- ((rho[1] * sigma_e[1]) / sigma_v[1]) * v1 + rnorm(n * prop[1], 0, sqrt((1 - rho[1]^2)* sigma_e[1]^2))
    e2 <- ((rho[2] * sigma_e[2]) / sigma_v[2]) * v2 + rnorm(n * prop[2], 0, sqrt((1 - rho[2]^2)* sigma_e[2]^2))
    e  <- c(e1, e2)
    
    # Generate instrument
    z1 <- rnorm(n* prop[1], Ez, sigma_z[1])
    z2 <- rnorm(n* prop[2], Ez, sigma_z[2])
    z <- c(z1, z2)
  
    # DGP
    y2   <- d0.v + d1.v * z     + v # Reduced equation
    y1   <- b0.v + gamma.v * y2 + e # Structural equation
    
    # Estimands
    ate.p <- prop[1]* gamma[1] +  prop[2] * gamma[2]
    cat(" - ATE  --> E(gamma_1)", " Population:", ate.p, "\n")
    aie.p <- prop[1]* d1[1] +  prop[2] * d1[2]
    cat(" - AIE  --> E(delta_1)", " Population:", aie.p, ", Sample:",  round(cov(y2, z) / var(z), 3), "\n")
    itt.p <- gamma[1] * prop[1] * d1[1] + gamma[2] * prop[2] * d1[2]
    cat(" - ITT  --> E(delta_1 * gamma) ", " Population:", itt.p, ", Sample:", round(cov(y1, z) / var(z), 3),"\n")
    late.p <- itt.p / aie.p
    cat(" - LATE --> E(delta_1 * gamma) / E(gamma)", " Population:", round(late.p, 3), ", Sample:", round(cov(z, y1) / cov(y2, z), 3), "\n")
    cat(" - OLS", ", Sample:", round(cov(y1, y2) / var(y2), 3), "\n")

    
    Ed0    <- prop[1]* d0[1] +  prop[2] * d0[2]
    Eb0    <- prop[1]* b0[1] +  prop[2] * b0[2]
    vd0    <- prop[1]*(d0[1]- Ed0)^2   + prop[2]*(d0[2]- Ed0)^2
    vd1    <- prop[1]*(d1[1]- aie.p)^2 + prop[2]*(d1[2]- aie.p)^2
    vgamma <- prop[1]*(gamma[1]- ate.p)^2 + prop[2]*(gamma[2]- ate.p)^2
    vv     <- prop[1]*(sigma_v[1])^2   + prop[2]*(sigma_v[2])^2
    ve     <- prop[1]*(sigma_e[1])^2   + prop[2]*(sigma_e[2])^2
    
    cov.gb <- gamma[1] * prop[1] * b0[1] + gamma[2] * prop[2] * b0[2] - Eb0*ate.p
    cat(" - Cov(b0, gamma): ", "Population", cov.gb, ", Sample", round(cov(b0.v, gamma.v),3), "\n")
    cov.db <- d0[1] * prop[1] * b0[1] + d0[2] * prop[2] * b0[2] - Eb0*Ed0
    cat(" - Cov(b0, d0): ", "Population", cov.db, ", Sample", round(cov(b0.v, d0.v),3), "\n")
    cov.gd <- itt.p - ate.p * aie.p
    cat(" - Cov(delta, gamma): ", "Population", cov.gd, ", Sample", round(cov(d1.v, gamma.v),3), "\n")
    vy2 <- vd0 + (vd1*sigma_z[1]^2 + vd1*Ez^2 + sigma_z[1]^2*aie.p^2) + vv
    cor.zy2 <- aie.p * sigma_z[1]^2 / (sigma_z[1] * sqrt(vy2))
    cat(" - Corr(y2, z): ", "Population", round(cor.zy2, 3), ", Sample", round(cor(y2, z),3), "\n")
    cat(" - Var(y2)    : ", "Population", vy2, ", Sample", round(var(y2),3), "\n")
    data[[s]] <- as.data.frame(cbind(y1, y2, z))
    cat(" - Corr(v, e) : ", "Population", rho[1]*prop[1] + rho[2]*prop[2], ",Sample", round(cor(e,v),3), "\n")
    cov.y2e <- (rho[1]*prop[1] + rho[2]*prop[2])
    cat(" - Corr(y2, e): ","Population", round(cov.y2e / (sqrt(vy2) * sqrt(ve)), 3),"Sample", round(cor(y2,e),3), "\n")
    cat(" - Cov(d * y2, y2): ", "Sample", round(cov(d1.v* y2,y2),3), "\n")
    cat(" - E(y1): ", "Sample", round(mean(y1),3), "\n")
    cat(" - E(y2): ", "Sample", round(mean(y2),3), "\n")
    s <- s + 1
  }
  Data[[count.n]] <- data
}

All_exp1 <- list(data = Data)
save(All_exp1, file  = "MC1_datasets.Rdata")


## Experiment 2 ----
data     <- list(S)
Data     <- list(length(N))
prop     <- c(0.3, 0.7) # Proportion of individuals in each class
rho      <- c(0,   0.5) # Degree of the endogeneity 
sigma_e  <- c(1, 1)     # Variance of epsilon
sigma_v  <- c(1, 1)     # Variance of upsilon
sigma_z  <- c(3, 3)
Ez       <- 0
# Structural parameters
b0       <- c(-1, 1)
gamma    <- c(-1, 2)
# First-Equation parameters
d0       <- c(-1, 1)
d1       <- c(-1, 2)


for (n in N){
  cat("=================", "\n")
  cat("Sample size: ", n, "\n")
  cat("=================", "\n")
  
  # True coefficients in the linear model
  b0.v    <- c(rep(b0[1],    n * prop[1]),  rep(b0[2],    n * prop[2]))
  gamma.v <- c(rep(gamma[1], n * prop[1]),  rep(gamma[2], n * prop[2]))
  
  # True coefficients in the first equation
  d0.v    <- c(rep(d0[1], n * prop[1]), rep(d0[2], n * prop[2]))
  d1.v    <- c(rep(d1[1], n * prop[1]), rep(d1[2], n * prop[2]))
  
  count.n <- which(N %in% n)
  print(count.n)
  s <- 1
  while(s <= S){
    cat( "Monte Carlo Sample: ", s, "\n")
    
    # Generate upsilon for each class
    v1 <- rnorm(n * prop[1], 0, sigma_v[1])
    v2 <- rnorm(n * prop[2], 0, sigma_v[2])
    v  <- c(v1, v2)
    
    # Generate epsilon for each class
    e1 <- ((rho[1] * sigma_e[1]) / sigma_v[1]) * v1 + rnorm(n * prop[1], 0, sqrt((1 - rho[1]^2)* sigma_e[1]^2))
    e2 <- ((rho[2] * sigma_e[2]) / sigma_v[2]) * v2 + rnorm(n * prop[2], 0, sqrt((1 - rho[2]^2)* sigma_e[2]^2))
    e  <- c(e1, e2)
    
    # Generate instrument
    z1 <- rnorm(n* prop[1], Ez, sigma_z[1])
    z2 <- rnorm(n* prop[2], Ez, sigma_z[2])
    z <- c(z1, z2)
    
    # DGP
    y2   <- d0.v + d1.v * z     + v # Reduced equation
    y1   <- b0.v + gamma.v * y2 + e # Structural equation
    
    # Estimands
    ate.p <- prop[1]* gamma[1] +  prop[2] * gamma[2]
    cat(" - ATE  --> E(gamma_1)", " Population:", ate.p, "\n")
    aie.p <- prop[1]* d1[1] +  prop[2] * d1[2]
    cat(" - AIE  --> E(delta_1)", " Population:", aie.p, ", Sample:",  round(cov(y2, z) / var(z), 3), "\n")
    itt.p <- gamma[1] * prop[1] * d1[1] + gamma[2] * prop[2] * d1[2]
    cat(" - ITT  --> E(delta_1 * gamma) ", " Population:", itt.p, ", Sample:", round(cov(y1, z) / var(z), 3),"\n")
    late.p <- itt.p / aie.p
    cat(" - LATE --> E(delta_1 * gamma) / E(gamma)", " Population:", round(late.p, 3), ", Sample:", round(cov(z, y1) / cov(y2, z), 3), "\n")
    cat(" - OLS", ", Sample:", round(cov(y1, y2) / var(y2), 3), "\n")
    
    
    Ed0    <- prop[1]* d0[1] +  prop[2] * d0[2]
    Eb0    <- prop[1]* b0[1] +  prop[2] * b0[2]
    vd0    <- prop[1]*(d0[1]- Ed0)^2   + prop[2]*(d0[2]- Ed0)^2
    vd1    <- prop[1]*(d1[1]- aie.p)^2 + prop[2]*(d1[2]- aie.p)^2
    vgamma <- prop[1]*(gamma[1]- ate.p)^2 + prop[2]*(gamma[2]- ate.p)^2
    vv     <- prop[1]*(sigma_v[1])^2   + prop[2]*(sigma_v[2])^2
    ve     <- prop[1]*(sigma_e[1])^2   + prop[2]*(sigma_e[2])^2
    
    cov.gb <- gamma[1] * prop[1] * b0[1] + gamma[2] * prop[2] * b0[2] - Eb0*ate.p
    cat(" - Cov(b0, gamma): ", "Population", cov.gb, ", Sample", round(cov(b0.v, gamma.v),3), "\n")
    cov.db <- d0[1] * prop[1] * b0[1] + d0[2] * prop[2] * b0[2] - Eb0*Ed0
    cat(" - Cov(b0, d0): ", "Population", cov.db, ", Sample", round(cov(b0.v, d0.v),3), "\n")
    cov.gd <- itt.p - ate.p * aie.p
    cat(" - Cov(delta, gamma): ", "Population", cov.gd, ", Sample", round(cov(d1.v, gamma.v),3), "\n")
    vy2 <- vd0 + (vd1*sigma_z[1]^2 + vd1*Ez^2 + sigma_z[1]^2*aie.p^2) + vv
    cor.zy2 <- aie.p * sigma_z[1]^2 / (sigma_z[1] * sqrt(vy2))
    cat(" - Corr(y2, z): ", "Population", round(cor.zy2, 3), ", Sample", round(cor(y2, z),3), "\n")
    cat(" - Var(y2)    : ", "Population", vy2, ", Sample", round(var(y2),3), "\n")
    data[[s]] <- as.data.frame(cbind(y1, y2, z))
    cat(" - Corr(v, e) : ", "Population", rho[1]*prop[1] + rho[2]*prop[2], ",Sample", round(cor(e,v),3), "\n")
    cov.y2e <- (rho[1]*prop[1] + rho[2]*prop[2])
    cat(" - Corr(y2, e): ","Population", round(cov.y2e / (sqrt(vy2) * sqrt(ve)), 3),"Sample", round(cor(y2,e),3), "\n")
    cat(" - Cov(d * y2, y2): ", "Sample", round(cov(d1.v* y2,y2),3), "\n")
    cat(" - E(y1): ", "Sample", round(mean(y1),3), "\n")
    cat(" - E(y2): ", "Sample", round(mean(y2),3), "\n")
    s <- s + 1
  }
  Data[[count.n]] <- data
}

All_exp2 <- list(data = Data)
save(All_exp2, file  = "MC2_datasets.Rdata")

## Experiment 3 ----
data     <- list(S)
Data     <- list(length(N))
prop     <- c(0.3, 0.7) # Proportion of individuals in each class
rho      <- c(0.5, 0.5) # Degree of the endogeneity 
sigma_e  <- c(1, 1)     # Variance of epsilon
sigma_v  <- c(1, 1)     # Variance of upsilon
sigma_z  <- c(3, 3)
Ez       <- 0
# Structural parameters
b0       <- c(-1, 1)
gamma    <- c(-1, 2)
# First-Equation parameters
d0       <- c(-1, 1)
d1       <- c(0, 2)


for (n in N){
  cat("=================", "\n")
  cat("Sample size: ", n, "\n")
  cat("=================", "\n")
  
  # True coefficients in the linear model
  b0.v    <- c(rep(b0[1],    n * prop[1]),  rep(b0[2],    n * prop[2]))
  gamma.v <- c(rep(gamma[1], n * prop[1]),  rep(gamma[2], n * prop[2]))
  
  # True coefficients in the first equation
  d0.v    <- c(rep(d0[1], n * prop[1]), rep(d0[2], n * prop[2]))
  d1.v    <- c(rep(d1[1], n * prop[1]), rep(d1[2], n * prop[2]))
  
  count.n <- which(N %in% n)
  print(count.n)
  s <- 1
  while(s <= S){
    cat( "Monte Carlo Sample: ", s, "\n")
    
    # Generate upsilon for each class
    v1 <- rnorm(n * prop[1], 0, sigma_v[1])
    v2 <- rnorm(n * prop[2], 0, sigma_v[2])
    v  <- c(v1, v2)
    
    # Generate epsilon for each class
    e1 <- ((rho[1] * sigma_e[1]) / sigma_v[1]) * v1 + rnorm(n * prop[1], 0, sqrt((1 - rho[1]^2)* sigma_e[1]^2))
    e2 <- ((rho[2] * sigma_e[2]) / sigma_v[2]) * v2 + rnorm(n * prop[2], 0, sqrt((1 - rho[2]^2)* sigma_e[2]^2))
    e  <- c(e1, e2)
    
    # Generate instrument
    z1 <- rnorm(n* prop[1], Ez, sigma_z[1])
    z2 <- rnorm(n* prop[2], Ez, sigma_z[2])
    z <- c(z1, z2)
    
    # DGP
    y2   <- d0.v + d1.v * z     + v # Reduced equation
    y1   <- b0.v + gamma.v * y2 + e # Structural equation
    
    # Estimands
    ate.p <- prop[1]* gamma[1] +  prop[2] * gamma[2]
    cat(" - ATE  --> E(gamma_1)", " Population:", ate.p, "\n")
    aie.p <- prop[1]* d1[1] +  prop[2] * d1[2]
    cat(" - AIE  --> E(delta_1)", " Population:", aie.p, ", Sample:",  round(cov(y2, z) / var(z), 3), "\n")
    itt.p <- gamma[1] * prop[1] * d1[1] + gamma[2] * prop[2] * d1[2]
    cat(" - ITT  --> E(delta_1 * gamma) ", " Population:", itt.p, ", Sample:", round(cov(y1, z) / var(z), 3),"\n")
    late.p <- itt.p / aie.p
    cat(" - LATE --> E(delta_1 * gamma) / E(gamma)", " Population:", round(late.p, 3), ", Sample:", round(cov(z, y1) / cov(y2, z), 3), "\n")
    cat(" - OLS", ", Sample:", round(cov(y1, y2) / var(y2), 3), "\n")
    
    
    Ed0    <- prop[1]* d0[1] +  prop[2] * d0[2]
    Eb0    <- prop[1]* b0[1] +  prop[2] * b0[2]
    vd0    <- prop[1]*(d0[1]- Ed0)^2   + prop[2]*(d0[2]- Ed0)^2
    vd1    <- prop[1]*(d1[1]- aie.p)^2 + prop[2]*(d1[2]- aie.p)^2
    vgamma <- prop[1]*(gamma[1]- ate.p)^2 + prop[2]*(gamma[2]- ate.p)^2
    vv     <- prop[1]*(sigma_v[1])^2   + prop[2]*(sigma_v[2])^2
    ve     <- prop[1]*(sigma_e[1])^2   + prop[2]*(sigma_e[2])^2
    
    cov.gb <- gamma[1] * prop[1] * b0[1] + gamma[2] * prop[2] * b0[2] - Eb0*ate.p
    cat(" - Cov(b0, gamma): ", "Population", cov.gb, ", Sample", round(cov(b0.v, gamma.v),3), "\n")
    cov.db <- d0[1] * prop[1] * b0[1] + d0[2] * prop[2] * b0[2] - Eb0*Ed0
    cat(" - Cov(b0, d0): ", "Population", cov.db, ", Sample", round(cov(b0.v, d0.v),3), "\n")
    cov.gd <- itt.p - ate.p * aie.p
    cat(" - Cov(delta, gamma): ", "Population", cov.gd, ", Sample", round(cov(d1.v, gamma.v),3), "\n")
    vy2 <- vd0 + (vd1*sigma_z[1]^2 + vd1*Ez^2 + sigma_z[1]^2*aie.p^2) + vv
    cor.zy2 <- aie.p * sigma_z[1]^2 / (sigma_z[1] * sqrt(vy2))
    cat(" - Corr(y2, z): ", "Population", round(cor.zy2, 3), ", Sample", round(cor(y2, z),3), "\n")
    cat(" - Var(y2)    : ", "Population", vy2, ", Sample", round(var(y2),3), "\n")
    data[[s]] <- as.data.frame(cbind(y1, y2, z))
    cat(" - Corr(v, e) : ", "Population", rho[1]*prop[1] + rho[2]*prop[2], ",Sample", round(cor(e,v),3), "\n")
    cov.y2e <- (rho[1]*prop[1] + rho[2]*prop[2])
    cat(" - Corr(y2, e): ","Population", round(cov.y2e / (sqrt(vy2) * sqrt(ve)), 3),"Sample", round(cor(y2,e),3), "\n")
    cat(" - Cov(d * y2, y2): ", "Sample", round(cov(d1.v* y2,y2),3), "\n")
    cat(" - E(y1): ", "Sample", round(mean(y1),3), "\n")
    cat(" - E(y2): ", "Sample", round(mean(y2),3), "\n")
    s <- s + 1
  }
  Data[[count.n]] <- data
}

All_exp3 <- list(data = Data)
save(All_exp3, file  = "MC3_datasets.Rdata")



