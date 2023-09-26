# Run MC experiments

# Directory and packages
setwd("~/Dropbox/work/Packages/Testing-ivhetLc")
# Load packages
rm(list = ls(all.names = TRUE))
library("ivhetLc")
library("AER")
library("memisc")

# Function to run the model across samples
all_models <- function(MC_data){
  cat("\n", "Estimating OLS model", "\n")
  ols <- lm(y1 ~  y2,
              data = MC_data)
  
  cat("\n", "Estimating IV model", "\n")
  iv <- ivreg(y1 ~  y2 | z ,
                    data = MC_data)
  
  cat("\n", "Estimating LCIV model", "\n")
  lciv <- NULL
  lciv <- tryCatch({ivhetLc(y1 ~  y2 | z | 1,
                   data = MC_data,
                   model = "lciv",
                   method = "bhhh")}, error = function(e) {
                     cat("ERROR :", conditionMessage(e), "\n")})
  
  lcivMiss <- NULL
  lcivMiss <- tryCatch({ivhetLc(y1 ~  y2 | z | 1,
                  data = MC_data,
                  model = "lciv",
                  method = "bhhh",
                  Q = 3)}, error = function(e) {
                    cat("ERROR :", conditionMessage(e), "\n")})
  
  if (is.null(lciv) & is.null(lcivMiss)) {
    table <- mtable("OLS" = ols, "IV" = iv,
                    summary.stats=c("logLik","N", "AIC", "BIC"),
                    signif.symbols = c("***" = .01, "**" = 0.05, "*" = 0.1))
    rownames(table$IV$coef) <- rownames(table$OLS$coef) <- c("eq.1.(Intercept)", "eq.1.y2")
    
    r.ols  <- table$OLS$coef
    r.iv   <- table$IV$coef
    r.lciv <-  NULL
    r.lcivmiss <- NULL
  }  else if (!is.null(lciv) & is.null(lcivMiss)) {
    table <- mtable("OLS" = ols, "IV" = iv, "LCIV" = lciv,
                    summary.stats=c("logLik","N", "AIC", "BIC"),
                    signif.symbols = c("***" = .01, "**" = 0.05, "*" = 0.1))
    rownames(table$IV$coef) <- rownames(table$OLS$coef) <- c("eq.1.(Intercept)", "eq.1.y2")
    
    r.ols  <- table$OLS$coef
    r.iv   <- table$IV$coef
    r.lciv <- table$LCIV$coef
    r.lcivmiss <- NULL
  } else if (is.null(lciv) & !is.null(lcivMiss)){
    table <- mtable("OLS" = ols, "IV" = iv, "LCIVM" = lcivMiss,
                    summary.stats=c("logLik","N", "AIC", "BIC"),
                    signif.symbols = c("***" = .01, "**" = 0.05, "*" = 0.1))
    rownames(table$IV$coef) <- rownames(table$OLS$coef) <- c("eq.1.(Intercept)", "eq.1.y2")
    
    r.ols  <- table$OLS$coef
    r.iv   <- table$IV$coef
    r.lciv <- NULL
    r.lcivmiss <- table$LCIVM$coef
  } else {
    table <- mtable("OLS" = ols, "IV" = iv, "LCIV" = lciv, "LCIVM" = lcivMiss,
                    summary.stats=c("logLik","N", "AIC", "BIC"),
                    signif.symbols = c("***" = .01, "**" = 0.05, "*" = 0.1))
    rownames(table$IV$coef) <- rownames(table$OLS$coef) <- c("eq.1.(Intercept)", "eq.1.y2")
    
    r.ols  <- table$OLS$coef
    r.iv   <- table$IV$coef
    r.lciv <- table$LCIV$coef
    r.lcivmiss <- table$LCIVM$coef
  }
  return(list("ols" = r.ols, "iv" = r.iv, "lciv" = r.lciv, "lcivmiss" = r.lcivmiss))
}

## Run MC E1 ----
load("MC1_datasets.Rdata")
FL1_n1 <- lapply(All_exp1$data[[1]], all_models)
FL1_n2 <- lapply(All_exp1$data[[2]], all_models)
FL1_n3 <- lapply(All_exp1$data[[3]], all_models)
FL1_n4 <- lapply(All_exp1$data[[4]], all_models)
save(FL1_n1, FL1_n2, FL1_n3, FL1_n4, 
     file = "results_MC1.Rdata")

## Run MC E2 ----
load("MC2_datasets.Rdata")
FL2_n1 <- lapply(All_exp2$data[[1]], all_models)
FL2_n2 <- lapply(All_exp2$data[[2]], all_models)
FL2_n3 <- lapply(All_exp2$data[[3]], all_models)
FL2_n4 <- lapply(All_exp2$data[[4]], all_models)
save(FL2_n1, FL2_n2, FL2_n3, FL2_n4, 
     file = "results_MC2.Rdata")

## Run MC E2 ----
load("MC3_datasets.Rdata")
FL3_n1 <- lapply(All_exp3$data[[1]], all_models)
FL3_n2 <- lapply(All_exp3$data[[2]], all_models)
FL3_n3 <- lapply(All_exp3$data[[3]], all_models)
FL3_n4 <- lapply(All_exp3$data[[4]], all_models)
save(FL3_n1, FL3_n2, FL3_n3, FL3_n4, 
     file = "results_MC3.Rdata")



