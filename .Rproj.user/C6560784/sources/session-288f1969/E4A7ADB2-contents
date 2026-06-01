# utility_PMB.R
library(rjags)
library(coda)
library(dplyr)
library(tidyr)

# ==============================================================================
# 1. Get the linear model string (fixed the position of dt definition)
# ==============================================================================
get_pmb_model_string <- function(model_type) {
  
  # --- A. General Linear Model ---
  if (model_type == "General") {
    return("
    model {
      # [Fix] Precompute time increments outside the loop
      for (k in 1:m) {
        dt[k] <- t1[k+1] - t1[k]
      }

      # Priors
      mu0 ~ dmnorm(mu_mu, invsigma_mu)
      
      for (j in 1:p) { xii[j] ~ dgamma(0.5, 1.0E-6) }
      Delta[1, 1] <- 4 * xii[1]; Delta[1, 2] <- 0
      Delta[2, 1] <- 0;          Delta[2, 2] <- 4 * xii[2]
      invSIG0 ~ dwish(inverse(Delta), 3)
      SIG0 <- inverse(invSIG0[,])
      
      # csi (xi)
      for (j in 1:p) {
        invcsi[j] ~ dgamma(0.001, 0.001)
        csi[j] <- sqrt(1/invcsi[j])
      }
      
      # Likelihood
      for (i in 1:n) {
        eta[1:p, i] ~ dmnorm(mu0, invSIG0)
        for (j in 1:p) {
          for (k in 1:m) {
            # General Model: Mean = dt / eta
            mu_x[j,k,i]   <- dt[k] / eta[j, i]
            prec_x[j,k,i] <- (eta[j, i]^2) / (csi[j]^2 * dt[k] + 1.0E-10)
            x_diff[j,k,i] ~ dnorm(mu_x[j,k,i], prec_x[j,k,i])
          }
        }
      }
      
      # Derived
      for (j in 1:p) { sigma[j] <- sqrt(SIG0[j,j]) }
      rho_12 <- SIG0[1,2] / (sigma[1] * sigma[2])
    }
    ")
  }
  
  # --- B. Wiener Linear Model ---
  if (model_type == "Wiener") {
    return("
    model {
      for (k in 1:m) {
        dt[k] <- t1[k+1] - t1[k]
      }

      mu0 ~ dmnorm(mu_mu, invsigma_mu) 
      
      for (j in 1:p) { xii[j] ~ dgamma(0.5, 1.0E-6) }
      Delta[1, 1] <- 4 * xii[1]; Delta[1, 2] <- 0
      Delta[2, 1] <- 0;          Delta[2, 2] <- 4 * xii[2]
      Tau_b ~ dwish(inverse(Delta), 3)
      SIG_b <- inverse(Tau_b[,])
      
      # sigma2_err (xi)
      for (j in 1:p) {
        inv_sigma2[j] ~ dgamma(0.001, 0.001)
        sigma2_err[j] <- sqrt(1/inv_sigma2[j]) 
      }
      
      # Likelihood
      for (i in 1:n) {
        b[i, 1:p] ~ dmnorm(mu0, Tau_b)
        for (j in 1:p) {
          for (k in 1:m) {
            # Wiener Model: Mean = dt * b
            mu_x[j,k,i]   <- dt[k] * b[i,j]
            prec_x[j,k,i] <- 1 / (sigma2_err[j]^2 * dt[k] + 1.0E-10)
            x_diff[j,k,i] ~ dnorm(mu_x[j,k,i], prec_x[j,k,i])
          }
        }
      }
      
      for (j in 1:p) { sigma[j] <- sqrt(SIG_b[j,j]) }
      rho_12 <- SIG_b[1,2] / (sigma[1] * sigma[2])
    }
    ")
  }
  
  # --- C. Student-t Linear Model ---
  if (model_type == "Student") {
    return("
    model {
      for (k in 1:m) {
        dt[k] <- t1[k+1] - t1[k]
      }

      eta ~ dmnorm(mu_mu, invsigma_mu)
      
      for (j in 1:p) { xii[j] ~ dgamma(0.5, 1.0E-6) }
      Delta[1, 1] <- 4 * xii[1]; Delta[1, 2] <- 0
      Delta[2, 1] <- 0;          Delta[2, 2] <- 4 * xii[2]
      SIG0_inv ~ dwish(inverse(Delta), 3)
      SIG0 <- inverse(SIG0_inv[,])
      
      # sig.w (xi)
      for (j in 1:p) {
        inv_sig_w[j] ~ dgamma(0.001, 0.001)
        sig.w[j] <- sqrt(1/inv_sig_w[j])
      }
      
      for(j in 1:p) { zeros[j] <- 0 }
      
      # Likelihood
      for (i in 1:n) {
        tau[i] ~ dgamma(5/2, 5/2) 
        z[i, 1:p] ~ dmnorm(zeros[], SIG0_inv[,])
        
        for(j in 1:p){
          theta[i,j] <- eta[j] + z[i,j] / sqrt(tau[i])
        }
        
        for (k in 1:m) {
          for (j in 1:p) {
            # Student Model: Mean = theta * dt
            mu_x[j,k,i]   <- theta[i,j] * dt[k]
            prec_x[j,k,i] <- tau[i] / (sig.w[j]^2 * dt[k] + 1.0E-10)
            x_diff[j,k,i] ~ dnorm(mu_x[j,k,i], prec_x[j,k,i])
          }
        }
      }
      
      for (j in 1:p) { sigma[j] <- sqrt(SIG0[j,j]) }
      rho_12 <- SIG0[1,2] / (sigma[1] * sigma[2])
    }
    ")
  }
}

# ==============================================================================
# 2. Data loading function
# ==============================================================================
load_pmb_data <- function(csv_file) {
  raw_data <- read.csv(csv_file)
  units <- unique(raw_data$Unit)
  n <- length(units)
  times <- sort(unique(raw_data$Time))
  m <- length(times)
  p <- 2
  
  y <- array(NA, dim = c(p, m, n))
  for (i in 1:n) {
    u_data <- filter(raw_data, Unit == units[i]) %>% arrange(Time)
    y[1, , i] <- u_data$PC1
    y[2, , i] <- u_data$PC2
  }
  
  m_diff <- m - 1
  x_diff <- array(NA, dim = c(p, m_diff, n))
  for (j in 1:p) {
    for (i in 1:n) {
      x_diff[j, , i] <- diff(y[j, , i])
    }
  }
  
  return(list(
    x_diff = x_diff,
    t1 = times, 
    n = n, 
    m = m_diff, 
    p = p,
    mu_mu = rep(0, p),
    invsigma_mu = 0.001 * diag(p)
  ))
}

# ==============================================================================
# 3. Inference function
# ==============================================================================
run_pmb_inference <- function(data_list, model_type, inits_params, 
                              n.chains = 2, n.adapt = 10000, n.iter = 50000) {
  
  model_string <- get_pmb_model_string(model_type)
  n <- data_list$n
  p <- data_list$p
  
  gen_inits <- function() {
    base_list <- list()
    
    if (model_type == "General") {
      base_list$eta <- matrix(rnorm(n * p, inits_params$mu0, 0.1), p, n)
      base_list$mu0 <- inits_params$mu0
      base_list$invSIG0 <- diag(p)
      base_list$invcsi <- rep(1, p)
    } else if (model_type == "Wiener") {
      base_list$b <- matrix(rnorm(n * p, inits_params$mu0, 0.1), n, p) 
      base_list$mu0 <- inits_params$mu0
      base_list$Tau_b <- diag(p)
      base_list$inv_sigma2 <- rep(1, p)
    } else { 
      base_list$eta <- inits_params$mu0
      base_list$SIG0_inv <- diag(p)
      base_list$tau <- rep(1, n)
      base_list$z <- matrix(0, n, p)
      base_list$inv_sig_w <- rep(1, p)
    }
    return(base_list)
  }
  
  inits_list <- replicate(n.chains, gen_inits(), simplify = FALSE)
  
  message(paste0(">>> Initializing PMB Model: ", model_type))
  
  jags_model <- jags.model(textConnection(model_string), 
                           data = data_list, inits = inits_list, 
                           n.chains = n.chains, n.adapt = n.adapt)
  
  message(">>> Burning in...")
  update(jags_model, n.adapt)
  
  monitor_vars <- c("mu0", "eta", "sigma", "rho_12", "csi", "sigma2_err", "sig.w")
  
  message(">>> Sampling...")
  samples <- coda.samples(jags_model, variable.names = monitor_vars, n.iter = n.iter)
  
  dic <- dic.samples(jags_model, n.iter = n.iter)
  
  return(list(samples = samples, dic = dic))
}

# ==============================================================================
# 4. Automated initial value estimation function
# ==============================================================================
estimate_initial_values <- function(data_list) {
  n <- data_list$n
  p <- data_list$p
  
  # Store the slope (drift) and residual standard deviation for each unit
  slopes <- matrix(NA, nrow = n, ncol = p)
  res_sds <- matrix(NA, nrow = n, ncol = p)
  
  # 1. Fit each unit separately using a simple linear approximation
  for (i in 1:n) {
    # Extract the time information for the current unit.
    # Note: data_list$x_diff contains increments. We can either reconstruct
    # cumulative observations or directly use the mean increments.
    # Here, the mean-increment approach is used for simplicity:
    # Mean(dX) = Slope * Mean(dt).
    # Assume equally spaced observations with dt = 1, or use the average dt.
    
    for (j in 1:p) {
      # Extract the increment sequence
      dx <- data_list$x_diff[j, , i]
      
      # Compute the observed average drift rate.
      # Here dt depends on the data. For the PMB data, Time is an index
      # such as 1, 2, 3, ..., so dt = 1.
      # If real time is used, divide by mean(diff(times)).
      # Here we simply use the mean increment as the slope estimate.
      slopes[i, j] <- mean(dx)
      
      # Compute volatility, approximated by the standard deviation of increments
      res_sds[i, j] <- sd(dx)
    }
  }
  
  # 2. Compute population-level statistics
  pop_drift_mean <- colMeans(slopes) # Population drift mean [p]
  pop_drift_sd   <- apply(slopes, 2, sd) # Standard deviation of population drift, used for Sigma [p]
  pop_volatility <- colMeans(res_sds) # Population volatility [p]
  
  # Return a list of parameters with their original physical meanings
  return(list(
    drift_mean = pop_drift_mean,
    drift_sd   = pop_drift_sd,
    volatility = pop_volatility
  ))
}

# ==============================================================================
# Function: generate and save trace plots
# ==============================================================================
library(bayesplot)
library(ggplot2)
color_scheme_set("mix-blue-red")

# Define plotting function
save_trace_plots <- function(samples, model_name) {
  # Select parameters automatically according to model type.
  # Matched parameters: mu0, sigma, rho, csi, true_scale, tau_scale.
  # A regular expression is used to exclude a large number of random effects eta[...].
  pars <- c("mu0", "sigma", "rho", "csi", "true_scale", "tau_scale")
  # Find variable names that actually exist in the mcmc.list object
  all_vars <- varnames(samples)
  # Select only core parameters, excluding eta
  target_vars <- grep("^(mu0|sigma|csi|true_scale|tau_scale)", all_vars, value = TRUE)
  # Handle rho separately, since it may appear as rho[1,2]
  rho_vars <- grep("rho", all_vars, value = TRUE)
  target_vars <- c(target_vars, rho_vars)
  
  if(length(target_vars) > 0) {
    p <- mcmc_trace(samples, pars = target_vars) + 
      labs(title = paste("Trace Plots:", model_name)) +
      theme_bw()
    
    # Save
    ggsave(paste0("PMB_Results/Trace_", gsub(" ", "_", model_name), ".pdf"), 
           plot = p, width = 12, height = 8)
  }
}

# Helper function 3: compute the ergodic mean of MCMC samples
# This function computes the cumulative average for each parameter and each chain
# from the first iteration to the current iteration.
# Helper function 3: compute the ergodic mean of MCMC samples
compute_ergodic_mean <- function(samples, par_names) {
  if (!inherits(samples, "mcmc.list")) {
    stop("Input must be an mcmc.list object.")
  }
  # ... Function body unchanged, same as the previously accepted code
  n_iter <- nrow(samples[[1]])
  n_chains <- length(samples)
  
  ergodic_data_list <- list()
  
  for (i in 1:n_chains) {
    chain_df <- as.data.frame(samples[[i]])
    
    for (par in par_names) {
      if (par %in% colnames(chain_df)) {
        # Compute cumulative mean
        cumulative_mean <- cumsum(chain_df[, par]) / 1:n_iter
        
        ergodic_data_list[[length(ergodic_data_list) + 1]] <- data.frame(
          Parameter = par,
          Chain = factor(i),
          Iteration = 1:n_iter,
          Ergodic_Mean = cumulative_mean
        )
      }
    }
  }
  return(bind_rows(ergodic_data_list))
}

# Revised helper function 4: plotmath label mapping function with hat symbols
parameter_plotmath_map <- function(r_name) {
  # Define mapping rules
  mapping <- c(
    # Population mean (mu0) with hat()
    "mu0[1]" = "hat(mu)[0 * ',' * 1]",
    "mu0[2]" = "hat(mu)[0 * ',' * 2]",
    
    # Student model mean (eta) with hat()
    "eta[1]" = "hat(eta)[1]",
    "eta[2]" = "hat(eta)[2]",
    
    # Population standard deviation (sigma) with hat()
    "sigma[1]" = "hat(sigma)[1]",
    "sigma[2]" = "hat(sigma)[2]",
    
    # Correlation coefficient (rho_12) with hat()
    "rho_12" = "hat(rho)[12]",
    
    # Volatility scale (csi/xi) with hat()
    "csi[1]" = "hat(xi)[1]",
    "csi[2]" = "hat(xi)[2]",
    
    # Wiener error term (sigma2_err) with hat()
    "sigma2_err[1]" = "hat(sigma)^2*phantom()[epsilon]*phantom()[1]", 
    "sigma2_err[2]" = "hat(sigma)^2*phantom()[epsilon]*phantom()[2]",
    
    # Student error term (sig.w) with hat()
    "sig.w[1]" = "hat(sigma)[w * ',' * 1]",
    "sig.w[2]" = "hat(sigma)[w * ',' * 2]"
  )
  
  plotmath_name <- mapping[r_name]
  if (is.na(plotmath_name)) return(r_name)
  return(plotmath_name)
}