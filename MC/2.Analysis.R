# Analysis

# Directory and packages
setwd("~/Dropbox/work/Packages/Testing-ivhetLc")
# Load packages
rm(list = ls(all.names = TRUE))
library("ivhetLc")
library("AER")
library("memisc")

# Change bad values
load("results_MC1.Rdata")


sanitize <- function(x, par, class){
  #if(class == "class 1") q <-1 else q <- 2
  share1 <- x[["lciv"]]["share", "est", "class 1"]
  share2 <- x[["lciv"]]["share", "est", "class 2"]
  if (class == "class 1"){
    if(share1 > share2) {
      return(x[["lciv"]][par, "est", "class 2"])
    } else {
      return(x[["lciv"]][par, "est", "class 1"])
    }
  } else {
    if(share1 < share2) {
      return(x[["lciv"]][par, "est", "class 2"])
    } else {
      return(x[["lciv"]][par, "est", "class 1"])
    }
  }
}

# Plot Gamma ----
gamma_q2_n1 <- as.numeric(lapply(FL1_n1, sanitize, par = "eq.1.y2", class = "class 2")) - 2
gamma_q2_n2 <- as.numeric(lapply(FL1_n2, sanitize, par = "eq.1.y2", class = "class 2")) - 2
gamma_q2_n3 <- as.numeric(lapply(FL1_n3, sanitize, par = "eq.1.y2", class = "class 2")) - 2
gamma_q2_n4 <- as.numeric(lapply(FL1_n4, sanitize, par = "eq.1.y2", class = "class 2")) - 2

gamma_q1_n1 <- as.numeric(lapply(FL1_n1, sanitize, par = "eq.1.y2", class = "class 1")) + 1
gamma_q1_n2 <- as.numeric(lapply(FL1_n2, sanitize, par = "eq.1.y2", class = "class 1")) + 1
gamma_q1_n3 <- as.numeric(lapply(FL1_n3, sanitize, par = "eq.1.y2", class = "class 1")) + 1
gamma_q1_n4 <- as.numeric(lapply(FL1_n4, sanitize, par = "eq.1.y2", class = "class 1")) + 1

par(mfrow = c(1,2))
boxplot(cbind(gamma_q1_n1, gamma_q1_n2, gamma_q1_n3, gamma_q1_n4), 
        main = "Class 1",
        names = c("n=100", "n = 500", "n = 1000", "n=5000"))
abline(h = 0, col = "red")
boxplot(cbind(gamma_q2_n1, gamma_q2_n2, gamma_q2_n3, gamma_q2_n4), 
        main = "Class 2",
        names = c("n=100", "n = 500", "n = 1000", "n=5000"))
abline(h = 0, col = "red")

# Plot Delta ----
delta_q1_n1 <- as.numeric(lapply(FL1_n1, sanitize, par = "eq.2.z", class = "class 1"))
delta_q1_n2 <- as.numeric(lapply(FL1_n2, sanitize, par = "eq.2.z", class = "class 1"))
delta_q1_n3 <- as.numeric(lapply(FL1_n3, sanitize, par = "eq.2.z", class = "class 1"))
delta_q1_n4 <- as.numeric(lapply(FL1_n4, sanitize, par = "eq.2.z", class = "class 1"))

delta_q2_n1 <- as.numeric(lapply(FL1_n1, sanitize, par = "eq.2.z", class = "class 2"))
delta_q2_n2 <- as.numeric(lapply(FL1_n2, sanitize, par = "eq.2.z", class = "class 2"))
delta_q2_n3 <- as.numeric(lapply(FL1_n3, sanitize, par = "eq.2.z", class = "class 2"))
delta_q2_n4 <- as.numeric(lapply(FL1_n4, sanitize, par = "eq.2.z", class = "class 2"))

par(mfrow = c(1,2))
boxplot(cbind(delta_q1_n1, delta_q1_n2, delta_q1_n3, delta_q1_n4), 
        main = "Class 1",
        names = c("n=100", "n = 500", "n = 1000", "n=5000"))
abline(h = -1, col = "red")
boxplot(cbind(delta_q2_n1, delta_q2_n2, delta_q2_n3, delta_q2_n4), 
        main = "Class 2",
        names = c("n=100", "n = 500", "n = 1000", "n=5000"))
abline(h = 2, col = "red")

# Plot rho ----
rho_q1_n1 <- as.numeric(lapply(FL1_n1, sanitize, par = "rho", class = "class 1"))
rho_q1_n2 <- as.numeric(lapply(FL1_n2, sanitize, par = "rho", class = "class 1"))
rho_q1_n3 <- as.numeric(lapply(FL1_n3, sanitize, par = "rho", class = "class 1"))
rho_q1_n4 <- as.numeric(lapply(FL1_n4, sanitize, par = "rho", class = "class 1"))

rho_q2_n1 <- as.numeric(lapply(FL1_n1, sanitize, par = "rho", class = "class 2"))
rho_q2_n2 <- as.numeric(lapply(FL1_n2, sanitize, par = "rho", class = "class 2"))
rho_q2_n3 <- as.numeric(lapply(FL1_n3, sanitize, par = "rho", class = "class 2"))
rho_q2_n4 <- as.numeric(lapply(FL1_n4, sanitize, par = "rho", class = "class 2"))

par(mfrow = c(1,2))
boxplot(cbind(rho_q1_n1, rho_q1_n2, rho_q1_n3, rho_q1_n4), 
        main = "Class 1",
        names = c("n=100", "n = 500", "n = 1000", "n=5000"))
abline(h = 0.5, col = "red")
boxplot(cbind(rho_q2_n1, rho_q2_n2, rho_q2_n3, rho_q2_n4), 
        main = "Class 2",
        names = c("n=100", "n = 500", "n = 1000", "n=5000"))
abline(h = 0.5, col = "red")

# Late ----
ate_n1 <- as.numeric(lapply(FL1_n1, function(x) x[["lciv"]]["wate", "est", "class 1"]))
ate_n2 <- as.numeric(lapply(FL1_n2, function(x) x[["lciv"]]["wate", "est", "class 1"]))
ate_n3 <- as.numeric(lapply(FL1_n3, function(x) x[["lciv"]]["wate", "est", "class 1"]))
ate_n4 <- as.numeric(lapply(FL1_n4, function(x) x[["lciv"]]["wate", "est", "class 1"]))

atem_n1 <- unlist(lapply(FL1_n1, function(x) x[["lcivmiss"]]["wate", "est", "class 1"]))
atem_n2 <- unlist(lapply(FL1_n2, function(x) x[["lcivmiss"]]["wate", "est", "class 1"]))
atem_n3 <- unlist(lapply(FL1_n3, function(x) x[["lcivmiss"]]["wate", "est", "class 1"]))
atem_n4 <- unlist(lapply(FL1_n4, function(x) x[["lcivmiss"]]["wate", "est", "class 1"]))

iv_n1 <- as.numeric((lapply(FL1_n1, function(x) x[["iv"]]["eq.1.y2", "est", "y1"])))
iv_n2 <- as.numeric((lapply(FL1_n2, function(x) x[["iv"]]["eq.1.y2", "est", "y1"])))
iv_n3 <- as.numeric((lapply(FL1_n3, function(x) x[["iv"]]["eq.1.y2", "est", "y1"])))
iv_n4 <- as.numeric((lapply(FL1_n4, function(x) x[["iv"]]["eq.1.y2", "est", "y1"])))

ols_n1 <- as.numeric((lapply(FL1_n1, function(x) x[["ols"]]["eq.1.y2", "est", "y1"])))
ols_n2 <- as.numeric((lapply(FL1_n2, function(x) x[["ols"]]["eq.1.y2", "est", "y1"])))
ols_n3 <- as.numeric((lapply(FL1_n3, function(x) x[["ols"]]["eq.1.y2", "est", "y1"])))
ols_n4 <- as.numeric((lapply(FL1_n4, function(x) x[["ols"]]["eq.1.y2", "est", "y1"])))


par(mfrow = c(2,2))
boxplot(cbind(ate_n1, ate_n2, ate_n3, ate_n4), 
        main = "ATE",
        names = c("n=100", "n = 500", "n = 1000", "n=5000"))
abline(h = 0.3 * -1 + 0.7 * 2, col = "red")
boxplot(list(atem_n1, atem_n2, atem_n3, atem_n4), 
        main = "ATE",
        names = c("n=100", "n = 500", "n = 1000", "n=5000"))
abline(h = 0.3 * -1 + 0.7 * 2, col = "red")
boxplot(cbind(iv_n1, iv_n2, iv_n3, iv_n4), 
        main = "Estimated IV",
        names = c("n=100", "n = 500", "n = 1000", "n=5000"))
abline(h = 2.818, col = "red")
boxplot(cbind(ols_n1, ols_n2, ols_n3, ols_n4), 
        main = "Estimated OLS",
        names = c("n=100", "n = 500", "n = 1000", "n=5000"))

