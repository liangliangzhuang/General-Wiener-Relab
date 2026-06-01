# ============================= Case Study =====================================
# ====================== Fatigue Crack Growth Data =============================
# Summary: General (4) + Wiener (2) + Student-t (2) = 8 models
# Workflow: MCMC inference -> diagnostic plots -> reliability simulation ->
#           visualization and model comparison
# ==============================================================================
rm(list = ls())
# Load utility scripts
source("Real-data/utility.R")             # General models: PP, PE, EP, EE
source("Real-data/Wiener_model.R")        # Wiener models: PP, EE
source("Real-data/Student_model.R")       # Student-t models: PP, EE
source("Real-data/Reliability_utils.R")
library(dplyr)
library(ggplot2)
library(tidyr)

# ==============================================================================
# 1. Data preparation
# ==============================================================================

current_case_folder <- "Real-data/crack/result"
csv_file_path <- "Real-data/crack/crack_growth_data.csv"

t_vector <- c(0, 10, 15, 20, 25, 30, 35, 40)
jags_data <- load_jags_data_from_csv(csv_file_path, t_vector)

# Failure thresholds for PC1 and PC2.
# Modify these values according to the application setting if needed.
thresholds <- c(4, 7)

# ==============================================================================
# 2. Initial values
# ==============================================================================

# --- A. Initial values for General models ---
inits_general_power <- list(
  eta_mean = 25,
  eta_sd = 1,
  mu0 = 25,
  invSIG0_scale = 1,
  true_scale = c(1, 1),
  tau_scale = c(1, 1),
  invcsi = c(1, 100)
)

inits_general_exp <- list(
  eta_mean = 0.5,
  eta_sd = 0.1,
  mu0 = 0.5,
  invSIG0_scale = 0.1,
  true_scale = c(0.05, 0.05),
  tau_scale = c(1, 1),
  invcsi = c(40000, 40000)
)

# --- B. Initial values for Wiener models ---
inits_wiener <- list(
  mu0 = 0.6,
  sigma2_err = 0.1,
  Tau_b_scale = 1,
  b_val = 0.6,
  scale = c(1, 1.5)
)

# --- C. Initial values for Student-t models ---
inits_student <- list(
  eta = c(0.03, 0.03),
  SIG0_inv_scale = 1,
  scale = c(1, 1),
  sig_w = c(1, 1)
)

# ==============================================================================
# 3. Batch model configuration
# ==============================================================================

scenarios <- data.frame(
  Model_Name = c(
    "General PP", "General PE", "General EP", "General EE",
    "Wiener PP", "Wiener EE",
    "Student PP", "Student EE"
  ),
  Class = c(
    "General", "General", "General", "General",
    "Wiener", "Wiener",
    "Student", "Student"
  ),
  Drift_Type = c(
    "Power", "Power", "Exponential", "Exponential",
    "Power", "Exponential",
    "Power", "Exponential"
  ),
  Diff_Type = c(
    "Power", "Exponential", "Power", "Exponential",
    NA, NA, NA, NA
  ),
  stringsAsFactors = FALSE
)

final_results_list <- list()

# ==============================================================================
# 4. Model fitting
# ==============================================================================

for (i in 1:nrow(scenarios)) {
  
  model_name  <- scenarios$Model_Name[i]
  model_class <- scenarios$Class[i]
  d_type      <- scenarios$Drift_Type[i]
  diff_type   <- scenarios$Diff_Type[i]
  
  message("\n--------------------------------------------------------")
  message(paste0("Running [", i, "/8]: ", model_name))
  message("--------------------------------------------------------")
  
  # ---------------------------------------------------------------------------
  # A. Run MCMC inference
  # ---------------------------------------------------------------------------
  
  if (model_class == "General") {
    
    # General models
    curr_inits <- if (d_type == "Power") inits_general_power else inits_general_exp
    
    # Model-specific adjustment for the General EP model
    if (model_name == "General EP") curr_inits$tau_scale <- 1
    
    res <- run_inference(
      jags_data,
      d_type,
      diff_type,
      curr_inits,
      n.chains = 2,
      n.adapt = 1000,
      n.iter = 5000
    )
    
    stats <- summary(res$samples)$statistics[, "Mean"]
    dic_val <- sum(res$dic$deviance) + sum(res$dic$penalty)
    
    # Parameter mapping for General models
    mu_1 <- stats["mu0[1]"]
    mu_2 <- stats["mu0[2]"]
    
    r_1 <- stats["true_scale[1]"]
    r_2 <- stats["true_scale[2]"]
    
    s_1 <- stats["tau_scale[1]"]
    s_2 <- stats["tau_scale[2]"]
    
    xi_1 <- stats["csi[1]"]
    xi_2 <- stats["csi[2]"]
    
  } else if (model_class == "Wiener") {
    
    # Wiener models
    res <- run_wiener_inference(
      jags_data,
      d_type,
      inits_wiener,
      n.chains = 2,
      n.adapt = 1000,
      n.iter = 5000
    )
    
    stats <- summary(res$samples)$statistics[, "Mean"]
    dic_val <- sum(res$dic$deviance) + sum(res$dic$penalty)
    
    # Parameter mapping for Wiener models.
    # Here, mu0 represents the drift mean.
    mu_1 <- stats["mu0[1]"]
    mu_2 <- stats["mu0[2]"]
    
    # The Wiener model has one scale parameter, which is used for both r and s.
    r_1 <- stats["scale[1]"]
    r_2 <- stats["scale[2]"]
    
    s_1 <- r_1
    s_2 <- r_2
    
    xi_1 <- stats["sigma2_err[1]"]
    xi_2 <- stats["sigma2_err[2]"]
    
  } else {
    
    # Student-t models
    res <- run_student_inference(
      jags_data,
      d_type,
      inits_student,
      n.chains = 2,
      n.adapt = 1000,
      n.iter = 5000
    )
    
    stats <- summary(res$samples)$statistics[, "Mean"]
    dic_val <- sum(res$dic$deviance) + sum(res$dic$penalty)
    
    # Parameter mapping for Student-t models.
    # Here, eta represents the drift mean.
    mu_1 <- stats["eta[1]"]
    mu_2 <- stats["eta[2]"]
    
    # The Student-t model has one scale parameter.
    r_1 <- stats["scale[1]"]
    r_2 <- stats["scale[2]"]
    
    s_1 <- r_1
    s_2 <- r_2
    
    xi_1 <- stats["sig.w[1]"]
    xi_2 <- stats["sig.w[2]"]
  }
  
  # ---------------------------------------------------------------------------
  # B. Save diagnostic plots
  # ---------------------------------------------------------------------------
  
  # save_diagnostics(res$samples, model_name, output_dir = current_case_folder)
  
  # ---------------------------------------------------------------------------
  # C. Store model summary
  # ---------------------------------------------------------------------------
  # For the Student-t and Wiener models, sigma and rho are the covariance
  # parameters of the random effects.
  
  row_data <- data.frame(
    Model   = model_name,
    mu_1    = mu_1,
    mu_2    = mu_2,
    sigma_1 = stats["sigma[1]"],
    sigma_2 = stats["sigma[2]"],
    rho_12  = stats["rho[1,2]"],
    r_1     = r_1,
    r_2     = r_2,
    s_1     = s_1,
    s_2     = s_2,
    xi_1    = xi_1,
    xi_2    = xi_2,
    DIC     = dic_val,
    stringsAsFactors = FALSE
  )
  
  final_results_list[[i]] <- row_data
}

# ==============================================================================
# 5. Model comparison table
# ==============================================================================

table_final <- do.call(rbind, final_results_list)

message("\n>>> All eight models have been fitted. Final model comparison table:\n")
print(table_final)

rownames(table_final) <- NULL

num_cols <- names(table_final)[-1]

table_final[num_cols] <- lapply(table_final[num_cols], function(x) {
  # Format numerical values to three decimal places.
  sprintf("%.3f", as.numeric(x))
})

table_final

# Save model comparison results if needed.
# write.csv(
#   table_final,
#   "Real-data/crack/result/Final_Model_Comparison_Table.csv",
#   row.names = FALSE
# )
#
# rdata_path <- file.path(
#   "Real-data/crack/result/Final_Model_Comparison_Table.RData"
# )
# save(table_final, file = rdata_path)
#
# save.image(
#   file = "Real-data/crack/result/Final_Model_all.RData"
# )

# ==============================================================================
# 6. Fitted degradation paths for the best model
# ==============================================================================

library(MASS)  # for mvrnorm

# Best model used for path fitting
best_model_name <- "General EE"

# Extract the parameter estimates of the selected model
best_params <- final_results_list[[
  which(sapply(final_results_list, function(x) x$Model) == best_model_name)
]]

# Time grid for plotting fitted degradation paths
plot_time <- seq(0, 45, length.out = 100)

# Compute fitted mean paths and uncertainty bands
df_fit <- get_fitted_paths(best_params, plot_time)

# Load observed degradation data
raw_data <- read.csv(csv_file_path) %>%
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "Value"
  )

# Plot fitted degradation paths
col_vals <- c("#3C5488FF", "#E64B35FF")

p_fit <- ggplot() +
  geom_point(
    data = raw_data,
    aes(x = Time, y = Value),
    alpha = 0.5,
    size = 1
  ) +
  geom_ribbon(
    data = df_fit,
    aes(x = Time, ymin = Lower, ymax = Upper, fill = PC),
    alpha = 0.3
  ) +
  geom_line(
    data = df_fit,
    aes(x = Time, y = Mean, color = PC),
    size = 1.2
  ) +
  facet_wrap(~PC, scales = "free_y") +
  scale_color_manual(values = col_vals) +
  scale_fill_manual(values = col_vals) +
  labs(
    x = "Time",
    y = "Degradation"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(family = "Times", face = "plain", size = 10),
    axis.text.y = element_text(family = "serif", size = 10),
    axis.title.x = element_text(family = "serif", size = 12),
    axis.title.y = element_text(family = "serif", size = 12),
    legend.title = element_text(family = "serif", size = 10),
    legend.text = element_text(family = "serif", size = 10)
  )

print(p_fit)

ggsave(
  file.path(current_case_folder, "Path_Fitting_Best_Model.pdf"),
  plot = p_fit,
  width = 4,
  height = 3
)

# ==============================================================================
# 7. Goodness-of-fit assessment: PC1 and PC2 separately
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# Extract parameter estimates of the selected model
mu_vec <- c(as.numeric(best_params$mu_1), as.numeric(best_params$mu_2))
r_vec  <- c(as.numeric(best_params$r_1),  as.numeric(best_params$r_2))
s_vec  <- c(as.numeric(best_params$s_1),  as.numeric(best_params$s_2))
xi_vec <- c(as.numeric(best_params$xi_1), as.numeric(best_params$xi_2))

df_raw <- read.csv(csv_file_path)

# ------------------------------------------------------------------------------
# Function for standardized residual calculation
# ------------------------------------------------------------------------------

calculate_residuals <- function(df, r, s, xi, pc_col) {
  
  vals <- df[[pc_col]]
  times <- df$Time
  units <- df$Unit
  
  residuals <- c()
  
  for (u in unique(units)) {
    
    idx <- which(units == u)
    
    t_u <- times[idx]
    x_u <- vals[idx]
    
    # Observed degradation increments
    dx <- diff(x_u)
    
    # Corresponding increments of the drift and diffusion time-scale functions
    t_curr <- t_u[-1]
    t_prev <- t_u[-length(t_u)]
    
    d_lambda <- exp(r * t_curr) - exp(r * t_prev)
    d_tau    <- exp(s * t_curr) - exp(s * t_prev)
    
    # Moment estimator of the unit-specific eta.
    # The approximation follows X ~ Lambda / eta, implying eta ~ Lambda / X.
    eta_est <- sum(d_lambda) / sum(dx)
    
    # Standardized residual:
    # Z = (observed increment - theoretical mean) / theoretical standard deviation
    mean_theo <- d_lambda / eta_est
    sd_theo   <- (xi * sqrt(d_tau)) / eta_est
    
    res <- (dx - mean_theo) / sd_theo
    residuals <- c(residuals, res)
  }
  
  return(residuals)
}

# Calculate standardized residuals for PC1 and PC2
res_pc1 <- calculate_residuals(df_raw, r_vec[1], s_vec[1], xi_vec[1], "PC1")
res_pc2 <- calculate_residuals(df_raw, r_vec[2], s_vec[2], xi_vec[2], "PC2")

# Kolmogorov-Smirnov tests against the standard normal distribution
ks_pc1 <- ks.test(res_pc1, "pnorm")
ks_pc2 <- ks.test(res_pc2, "pnorm")

cat("PC1 K-S test p-value:", ks_pc1$p.value, "\n")
cat("PC2 K-S test p-value:", ks_pc2$p.value, "\n")

# Prepare data for Q-Q plots
df_plot <- rbind(
  data.frame(Residuals = res_pc1, PC = "PC1"),
  data.frame(Residuals = res_pc2, PC = "PC2")
)

# Q-Q plot of standardized residuals
p_qq <- ggplot(df_plot, aes(sample = Residuals)) +
  stat_qq(color = "#3C5488FF", alpha = 0.6, size = 0.8) +
  stat_qq_line(color = "#E64B35FF", linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~PC, scales = "free") +
  labs(
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(family = "Times", face = "plain", size = 10),
    axis.text.y = element_text(family = "serif", size = 10),
    axis.title.x = element_text(family = "serif", size = 12),
    axis.title.y = element_text(family = "serif", size = 12),
    legend.title = element_text(family = "serif", size = 10),
    legend.text = element_text(family = "serif", size = 10)
  )

print(p_qq)

ggsave(
  file.path(current_case_folder, "Goodness_of_Fit_QQ.pdf"),
  plot = p_qq,
  width = 4,
  height = 3
)