# ==============================================================================
# Case Study II: Permanent Magnet Brake
# Linear results: proposed model, Wiener model, and Student-t model
# Final version for the manuscript dated 12.09
# Includes Figures 7-8 and Tables 4-5. Figure 9 is provided in PMB_Relab.R
# ==============================================================================
rm(list = ls())
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS) 
source("Real-data/PMB/utility_PMB.R")

# 1. Data preparation
# ----------------------------------------
csv_file <- "Real-data/PMB/PMB_data.csv" # Make sure this file has been generated
# PMB data are equally spaced; therefore, the index is directly used as time
jags_data <- load_pmb_data(csv_file)

# ==============================================================================
# 2. Automated Initialization
# ==============================================================================
# Initial values are estimated from the data rather than manually specified
auto_vals <- estimate_initial_values(jags_data)

cat(">>> Automatically estimated physical parameters:\n")
cat("PC1 (Release Time) drift mean:", round(auto_vals$drift_mean[1], 4), "\n")
cat("PC2 (Torque)       drift mean:", round(auto_vals$drift_mean[2], 4), "\n")

# Construct the list of initial values
inits_list <- list()

# --- General Model ---
# Definition: Mean = dt / eta  =>  eta = dt / Mean = 1 / Slope (assuming dt=1)
# PC2 has a negative slope, so eta can be negative, which is allowed under the Normal prior
inits_list$General <- list(
  # Use 1/drift as the initial mean value of eta
  # To avoid division by zero, a small term may be added if necessary
  mu0 = 1 / auto_vals$drift_mean, 
  
  # csi: volatility coefficient
  # In the General model, Var = eta^2 / (csi^2). SD_data = eta / csi
  # => csi = eta / SD_data = (1/Drift) / Volatility
  # Here invcsi is the precision, while csi is the SD coefficient.
  # A reasonable initial value is sufficient.
  # For simplicity, initializing csi as 1.0 is usually safe.
  # Here we keep the original logic or set it to 1.
  invcsi = rep(1, jags_data$p) 
)

# --- Wiener Model ---
# Definition: Mean = dt * b => b = Slope
# sigma2_err represents volatility
inits_list$Wiener <- list(
  mu0 = auto_vals$drift_mean, # Directly use the estimated slope
  
  # sigma2_err = volatility
  # inv_sigma2 = 1 / var
  inv_sigma2 = 1 / (auto_vals$volatility^2)
)

# --- Student Model ---
# The definition is the same as the Wiener model
inits_list$Student <- list(
  mu0 = auto_vals$drift_mean, # Directly use the estimated slope
  
  # sig.w = volatility
  inv_sig_w = 1 / (auto_vals$volatility^2)
)


# 3. Run scenarios
models <- c("General", "Wiener", "Student")
final_results <- list()
final_results_with_samples <- list() # New: used to store samples and summary statistics
for (mod in models) {
  message(paste0("\n=== Running ", mod, " Linear Model ==="))
  
  # Run inference
  res <- run_pmb_inference(
    data_list = jags_data, 
    model_type = mod, 
    inits_params = inits_list[[mod]],
    n.chains = 3, 
    n.adapt = 10000, 
    n.iter = 40000
  )
  
  # ********** [Add MCMC diagnostics here] **********
  
  # Check whether multiple chains were successfully sampled
  if (inherits(res$samples, "mcmc.list") && length(res$samples) > 1) {
    # 1. Gelman-Rubin Rhat diagnostic
    gelman_diag <- coda::gelman.diag(res$samples, multivariate = FALSE)
    max_rhat <- max(gelman_diag$psrf[, 1], na.rm = TRUE)
    
    # 2. Effective sample size (ESS), used as a proxy for sampling efficiency
    ess_diag <- coda::effectiveSize(res$samples)
    min_ess <- min(ess_diag)
    
    message(">>> MCMC diagnostics:")
    cat(paste0(mod, " model - Max Rhat: "), round(max_rhat, 4), "\n")
    cat(paste0(mod, " model - Min ESS: "), round(min_ess, 0), "\n")
    
    # 3. Extract Rhat and ESS results
    Max_Rhat_Val <- max_rhat
    Min_ESS_Val <- min_ess
  } else {
    # If sampling failed or the number of chains is insufficient
    Max_Rhat_Val <- NA
    Min_ESS_Val <- NA
  }
  
  # ********** [End of MCMC diagnostics] **********
  
  # Extract summary statistics
  stats <- summary(res$samples)$statistics[, "Mean"]
  dic_val <- sum(res$dic$deviance) + sum(res$dic$penalty)
  
  get_v <- function(n) { if(n %in% names(stats)) stats[n] else NA }
  
  # Parameter mapping
  if (mod == "General") {
    mu <- c(get_v("mu0[1]"), get_v("mu0[2]"))
    xi <- c(get_v("csi[1]"), get_v("csi[2]"))
  } else if (mod == "Wiener") {
    mu <- c(get_v("mu0[1]"), get_v("mu0[2]"))
    xi <- c(get_v("sigma2_err[1]"), get_v("sigma2_err[2]"))
  } else { 
    mu <- c(get_v("eta[1]"), get_v("eta[2]"))
    xi <- c(get_v("sig.w[1]"), get_v("sig.w[2]"))
  }
  
  sigma <- c(get_v("sigma[1]"), get_v("sigma[2]"))
  rho   <- get_v("rho_12")
  
  row_df <- data.frame(
    Model = paste(mod, "Linear"),
    mu_1 = mu[1], mu_2 = mu[2],
    sigma_1 = sigma[1], sigma_2 = sigma[2],
    rho_12 = rho,
    xi_1 = xi[1], xi_2 = xi[2],
    DIC = dic_val,
    Max_Rhat = Max_Rhat_Val, # New diagnostic metric
    Min_ESS = Min_ESS_Val,# New diagnostic metric
    stringsAsFactors = FALSE
  )
  
  final_results[[mod]] <- row_df
  # Key modification: store summary statistics together with MCMC samples
  final_results_with_samples[[mod]] <- list(
    stats = row_df,
    samples = res$samples # Store complete MCMC samples
  )
}

# 4. Summarize results

table_pmb <- do.call(rbind, final_results)
num_cols <- names(table_pmb)[-1]
table_pmb[num_cols] <- lapply(table_pmb[num_cols], function(x) sprintf("%.3f", as.numeric(x)))

message("\n>>> PMB Case Study Results (Linear Models):")
print(table_pmb, row.names = FALSE)

# Save
# write.csv(table_pmb, "Real-data/PMB/result/1209-PMB_Linear_Models_Result.csv", row.names = FALSE)


# Assume final_results_with_samples already exists in the environment and contains the result of the "General" model
res_general <- final_results_with_samples[["General"]]
mcmc_samples <- res_general$samples
merged_samples <- as.matrix(mcmc_samples)

report_params <- c(
  "mu0[1]", "mu0[2]",
  "rho_12",
  "sigma[1]", "sigma[2]",
  "csi[1]", "csi[2]",
  "tau_scale[1]", "tau_scale[2]", # Drift time-scale parameter r
  "true_scale[1]", "true_scale[2]" # Diffusion time-scale parameter s
)
cols_to_select <- report_params[report_params %in% colnames(merged_samples)]
merged_report_data <- merged_samples[, cols_to_select]
# 4. Compute quantiles using the standard R apply + quantile functions
quantiles_matrix <- apply(
  merged_report_data, 
  2, # Apply by column, i.e., by parameter
  quantile, 
  probs = c(0.025, 0.5, 0.975)
)
quantiles_df <- as.data.frame((quantiles_matrix))
row.names(quantiles_df) <- c("2.5% quantile", "Posterior mean", "97.5% quantile")
quantiles_df 


# ==============================================================================
# Posterior convergence plots
# ==============================================================================
# List of key parameters: group-level parameters (mu0, sigma, rho_12) and volatility parameters (csi, sigma2_err, sig.w)
my_colors <- c("#E64B35FF", # Red, used for emphasis
               "#45498C", # Blue
               "#00A087FF") # Green

par_names_to_plot <- c(
  "mu0[1]", "mu0[2]", "sigma[1]", "sigma[2]", "rho_12", 
  "csi[1]", "csi[2]", "sigma2_err[1]", "sigma2_err[2]", "sig.w[1]", "sig.w[2]"
)

for (mod in models) {
  
  res_list <- final_results_with_samples[[mod]]
  if (is.null(res_list)) next
  
  samples <- res_list$samples
  
  # Select the parameters that are actually included in the current model
  if (mod == "General") {
    model_pars <- par_names_to_plot[grepl("mu0|sigma|rho_12|csi", par_names_to_plot)]
  } else if (mod == "Wiener") {
    model_pars <- par_names_to_plot[grepl("mu0|sigma|rho_12|sigma2_err", par_names_to_plot)]
  } else { # Student
    model_pars <- par_names_to_plot[grepl("eta|sigma|rho_12|sig.w", par_names_to_plot)]
  }
  
  message(paste0(">>> Generating Ergodic Mean Plot for: ", mod))
  
  # ******** 1. Ergodic mean plots ********
  df_ergodic <- compute_ergodic_mean(samples, model_pars)
  # 1. Determine the expected order of parameters in the current model
  # Assume Chain 1 is the chain to be displayed. Adjust this condition according to your chain names if needed.
  df_ergodic <- df_ergodic[df_ergodic$Chain == 3, ]
  
  ordered_levels <- intersect(par_names_to_plot, model_pars)
  # 2. Reset the factor levels of the Parameter column in df_ergodic
  df_ergodic$Parameter <- factor(df_ergodic$Parameter, levels = ordered_levels)
  # Create a new label column that stores plotmath strings
  df_ergodic$Plotmath_Label <- sapply(as.character(df_ergodic$Parameter), parameter_plotmath_map)
  # Reset factor levels of Plotmath_Label to preserve the order
  df_ergodic$Plotmath_Label <- factor(df_ergodic$Plotmath_Label, levels = unique(df_ergodic$Plotmath_Label[match(ordered_levels, df_ergodic$Parameter)]))
  
  
  p_ergodic <- ggplot(df_ergodic, aes(x = Iteration, y = Ergodic_Mean, color = Chain)) +
    geom_line(linewidth = 0.5) +
    # ******** Key modification: apply customized labeller and color ********
    facet_wrap(~Plotmath_Label, scales = "free_y", ncol = 4, labeller = label_parsed) +
    scale_color_manual(values = my_colors[2]) +
    labs(
      # title = paste("Ergodic Mean Plots for", mod, "Model"),
      x = latex2exp::TeX("Iteration ($\\times 10^4$)"),
      y = "Ergodic mean value"
    ) +
    scale_x_continuous(
      breaks = seq(0, 40000, by = 10000),
      labels = seq(0, 40000, by = 10000) / 10000, # Convert tick values to 0, 1, 2, ..., 5
      limits = c(0, 40000)
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      # strip.background = element_rect(fill = "grey70", color = "black"),
      strip.text = element_text(family = "Times", face = "plain", size = 10),
      axis.title.x = element_text(family = "Times", face = "plain", size = 12),
      axis.title.y = element_text(family = "Times", face = "plain", size = 12),
      legend.title = element_text(family = "Times", face = "plain", size = 12),
      legend.text  = element_text(family = "Times", face = "plain", size = 11)
    ) 
  
}


# ==============================================================================
# Extension 1: Mean-path fitting plot for the best PMB model (General Linear)
# ==============================================================================


# 1. Prepare parameters extracted from previous results
# Assume final_results is still available in the environment.
# If not, manually enter the values from the table.
# Here the "General" model is used by default.
model_res <- final_results[["General"]] 

if(is.null(model_res)) stop("Please run main_PMB.R first to generate final_results")

params_pmb <- list(
  mu    = c(as.numeric(model_res$mu_1), as.numeric(model_res$mu_2)),
  sigma = c(as.numeric(model_res$sigma_1), as.numeric(model_res$sigma_2)),
  rho   = as.numeric(model_res$rho_12),
  xi    = c(as.numeric(model_res$xi_1), as.numeric(model_res$xi_2))
)

# 2. Simulation function optimized for the linear model
# General Linear: X(t) = t / eta + Diffusion
get_pmb_fit <- function(params, time_max, n_sim = 1000) {
  p <- 2
  time_seq <- seq(0, time_max, length.out = 100) # Continuous time grid for plotting
  
  # Construct covariance matrix
  SIG0 <- diag(params$sigma) %*% matrix(c(1, params$rho, params$rho, 1), 2) %*% diag(params$sigma)
  
  res_list <- lapply(1:p, function(j) {
    # 1. Sample random effect eta
    # Note: PC2 shows a decreasing trend, so the mean of eta is negative. No absolute value is needed.
    eta_vec <- mvrnorm(n_sim, params$mu, SIG0)
    eta <- eta_vec[, j]
    
    # 2. Compute the mean path
    # Drift = t / eta
    # Use outer-product form: (1/eta) * t
    # Dimension: [n_sim, 1] * [1, n_time]
    drift_paths <- (1 / eta) %*% t(time_seq)
    
    # 3. Add diffusion term for uncertainty-band simulation
    # Var(X(t)) = xi^2 * t / eta^2
    # SD(X(t))  = |xi * sqrt(t) / eta|
    # Generate a standard normal noise matrix
    Z <- matrix(rnorm(n_sim * length(time_seq)), nrow = n_sim)
    
    # Diffusion component: sigma(t) * Z
    # Note: sqrt(t) is a vector and needs to be broadcast properly
    sd_term <- (abs(params$xi[j]) / abs(eta)) %*% t(sqrt(time_seq))
    diffusion_paths <- sd_term * Z
    
    # Total paths
    total_paths <- drift_paths + diffusion_paths
    
    # 4. Compute summary statistics
    data.frame(
      Time = time_seq,
      Mean = colMeans(total_paths),
      Lower = apply(total_paths, 2, quantile, probs = 0.05),
      Upper = apply(total_paths, 2, quantile, probs = 0.95),
      PC = paste0("PC", j)
    )
  })
  
  do.call(rbind, res_list)
}

# 3. Run simulation
df_fit_pmb <- get_pmb_fit(params_pmb, time_max = 30)

# Update PC labels for display
df_fit_pmb$PC <- factor(df_fit_pmb$PC, levels = c("PC1", "PC2"), 
                        labels = c("PC1", "PC2"))

# 4. Prepare observed data
raw_pmb <- read.csv("Real-data/PMB/PMB_data.csv") %>%
  pivot_longer(cols = c("PC1", "PC2"), names_to = "PC", values_to = "Value") %>%
  mutate(PC = factor(PC, levels = c("PC1", "PC2"), 
                     labels = c("PC1", "PC2")))


init_vals <- raw_pmb %>% 
  filter(Time == 1) %>% 
  group_by(PC) %>% 
  summarise(Init = mean(Value))

df_fit_pmb <- df_fit_pmb %>% 
  left_join(init_vals, by = "PC") %>%
  mutate(
    Mean = Mean + Init,
    Lower = Lower + Init,
    Upper = Upper + Init
  )


# 5. Plot
p_fit <- ggplot() +
  # Fitted uncertainty band
  geom_ribbon(data = df_fit_pmb, aes(x = Time, ymin = Lower, ymax = Upper, fill = PC), alpha = 0.3) +
  # Mean path
  geom_line(data = df_fit_pmb, aes(x = Time, y = Mean, color = PC), linewidth = 0.8) +
  # Observed points
  geom_point(data = raw_pmb, aes(x = Time, y = Value), alpha = 0.5, size = 0.4) +
  
  facet_wrap(~PC, scales = "free_y") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  
  labs(
    # title = "PMB Degradation Fitting (General Linear Model)",
    x = "Time", y = "Degradation value") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    # strip.background = element_rect(fill = "grey70", color = "black"),
    strip.text = element_text(family = "Times", face = "plain", size = 10),
    axis.title.x = element_text(family = "Times", face = "plain", size = 12),
    axis.title.y = element_text(family = "Times", face = "plain", size = 12),
    legend.title = element_text(family = "Times", face = "plain", size = 12),
    legend.text  = element_text(family = "Times", face = "plain", size = 11)
  ) 

print(p_fit)
ggsave("Real-data/PMB/result/PMB_Path_Fitting.pdf", width = 4, height = 3)



# ==============================================================================
# Extension 2: PMB goodness-of-fit assessment (Residual Q-Q Plot)
# ==============================================================================

# 1. Residual calculation function
calculate_pmb_residuals <- function(df, mu_vec, xi_vec) {
  # For the General Linear model, the posterior mean mu_vec can be used as the population-level prior.
  # However, to compute residuals for each specific unit, we need an estimate of the unit-specific eta.
  # The simplest approach is to estimate eta_i using the empirical average rate of each unit.
  
  residuals_list <- list()
  
  for(pc_idx in 1:2) {
    pc_name <- paste0("PC", pc_idx)
    xi <- xi_vec[pc_idx]
    
    # Extract data for the current PC
    pc_data <- df %>% dplyr::select(Unit, Time, all_of(pc_name)) %>% rename(Value = !!pc_name)
    
    res_vec <- c()
    
    for(u in unique(pc_data$Unit)) {
      sub <- pc_data %>% filter(Unit == u) %>% arrange(Time)
      x <- sub$Value
      t <- sub$Time
      
      # Calculate increments
      dx <- diff(x)
      dt <- diff(t) # For PMB data, dt should always be 1
      
      # Estimate eta for the current unit
      # Linear General: X ~ t / eta => eta ~ t / X => eta ~ sum(dt)/sum(dx)
      eta_est <- sum(dt) / sum(dx)
      
      # Calculate theoretical mean and variance
      mean_theo <- dt / eta_est
      sd_theo   <- (xi * sqrt(dt)) / abs(eta_est)
      
      # Standardized residuals
      z <- (dx - mean_theo) / sd_theo
      res_vec <- c(res_vec, z)
    }
    
    residuals_list[[pc_name]] <- res_vec
  }
  
  return(residuals_list)
}

# 2. Run calculation
res_list <- calculate_pmb_residuals(read.csv("Real-data/PMB/PMB_data.csv"), 
                                    params_pmb$mu, params_pmb$xi)

# 3. K-S tests
ks1 <- ks.test(res_list$PC1, "pnorm")
ks2 <- ks.test(res_list$PC2, "pnorm")

# 4. Plot
df_qq <- rbind(
  data.frame(Res = res_list$PC1, PC = "PC1"),
  data.frame(Res = res_list$PC2, PC = "PC2")
)
df_label <- data.frame(
  PC = c("PC1", "PC2"), # Corresponding faceting variable
  label = c(paste0("p-value = ", sprintf("%.3f", ks1$p.value)),
            paste0("p-value = ", sprintf("%.3f", ks2$p.value))),
  x = c(-3.5, -3.5), # x-coordinate of the text label, theoretical quantile
  y = c(2.3, 2)    # y-coordinate of the text label, sample quantile
)

p_qq <- ggplot(df_qq, aes(sample = Res)) +
  stat_qq(color = "#3C5488FF", alpha = 0.7, size = 0.6) + 
  stat_qq_line(color = "#E64B35FF", linewidth = 0.8, linetype = "dashed") + 
  facet_wrap(~PC, scales = "free") +
  labs(
    # title = "Residual Q-Q Plot (PMB Data)", 
    x = "Theoretical Quantiles", y = "Sample Quantiles") +
  # geom_label(data = df_label, 
  #            aes(x = x, y = y, label = label), 
  #            inherit.aes = FALSE,
  #            
  #            fill = "white",      # White background fill
  #            color = "#E64B35FF",     # Font color
  #            family = "Times",    # Font family
  #            size = 3, 
  #            hjust = 0,           # Left alignment
  #            
  #            # Optional parameters: adjust corner radius and padding
  #            label.padding = unit(0.2, "lines"), 
  #            label.r = unit(0.15, "lines")       
  # ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    # strip.background = element_rect(fill = "grey70", color = "black"),
    strip.text = element_text(family = "Times", face = "plain", size = 10),
    axis.title.x = element_text(family = "Times", face = "plain", size = 12),
    axis.title.y = element_text(family = "Times", face = "plain", size = 12),
    legend.title = element_text(family = "Times", face = "plain", size = 12),
    legend.text  = element_text(family = "Times", face = "plain", size = 11)
  ) 

p_qq 
ggsave("Real-data/PMB/result/PMB_GoF_QQ.pdf", width = 4, height = 3)


save.image(file = paste("Real-data/PMB/result/Final_PMB_Model_all.RData", sep=''))