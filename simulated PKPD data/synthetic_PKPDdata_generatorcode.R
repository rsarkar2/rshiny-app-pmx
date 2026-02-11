set.seed(123)
setwd("/Users/reetikasarkar/Projects/RShiny")
# ---- Population parameters ----
theta <- list(
  ka = 1.2,     # 1/h
  cl = 5,       # L/h
  v  = 50       # L
)

omega <- list(
  cl = 0.25,    # IIV SD on CL
  v  = 0.30     # IIV SD on V
)

sigma <- 0.20   # residual error SD

dose <- 100     # mg
times <- c(0.5, 1, 2, 4, 8)

# ---- PK function (no covariates) ----
pk_conc <- function(time, dose, ka, cl, v) {
  ke <- cl / v
  dose * ka / (v * (ka - ke)) *
    (exp(-ke * time) - exp(-ka * time))
}

# ---- Simulate data ----
library(dplyr)
library(purrr)
library(tibble)

sim_data <- map_dfr(1:10, function(id) {
  
  eta_cl <- rnorm(1, 0, omega$cl)
  eta_v  <- rnorm(1, 0, omega$v)
  
  cl_i <- theta$cl * exp(eta_cl)
  v_i  <- theta$v  * exp(eta_v)
  
  conc <- pk_conc(times, dose, theta$ka, cl_i, v_i)
  conc_obs <- pmax(conc + rnorm(length(times), 0, sigma), 0.01)
  
  tibble(
    id   = id,
    time = times,
    conc = round(conc_obs, 2),
    dose = dose,
    WT   = 70,   # fixed constant
    SEX  = 0    # fixed constant
  )
})

# ---- Save to CSV ----
write.csv(sim_data, "synthetic_pk_data.csv", row.names = FALSE)
