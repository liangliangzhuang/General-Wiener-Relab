# ==============================================================================
# Exploratory data analysis: Assessing the time-varying mean-variance ratio
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)

# 1. Data preparation using the predefined matrices
pc1_raw <- matrix(c(
  1,0,0.71,1.07,1.49,1.86,2.36,2.96,3.46,
  2,0,0.73,1.03,1.41,1.71,2.30,2.88,3.51,
  3,0,0.63,0.98,1.54,1.91,2.40,2.99,3.67,
  4,0,0.58,0.88,1.34,1.92,2.36,2.96,3.68,
  5,0,0.60,1.02,1.44,1.89,2.45,3.05,3.73,
  6,0,0.77,1.20,1.64,2.10,2.62,3.17,3.80,
  7,0,0.75,1.23,1.67,2.17,2.71,3.46,4.19,
  8,0,0.80,1.28,1.80,2.35,2.96,3.63,4.45,
  9,0,0.91,1.32,1.87,2.44,3.06,3.74,4.60,
  10,0,0.78,1.20,1.81,2.27,3.00,3.84,4.73,
  11,0,0.89,1.37,1.99,2.51,3.25,4.09,4.78,
  12,0,0.81,1.35,1.91,2.60,3.28,4.02,5.11,
  13,0,0.92,1.41,2.07,2.70,3.48,4.44,5.65,
  14,0,0.95,1.63,2.28,3.12,3.89,5.06,6.11,
  15,0,1.32,1.91,2.37,3.24,3.99,5.00,6.58
), nrow = 15, byrow = TRUE)

pc2_raw <- matrix(c(
  1,0,1.23,2.04,2.82,3.87,5.15,6.87,9.46,
  2,0,1.01,1.65,2.41,3.21,3.96,5.14,6.53,
  3,0,0.66,1.34,2.03,2.59,3.55,4.71,6.16,
  4,0,1.25,2.16,3.03,3.74,4.80,6.09,7.83,
  5,0,0.96,1.63,2.28,3.18,4.12,5.31,6.80,
  6,0,1.19,1.77,2.58,3.32,4.46,5.58,7.05,
  7,0,1.05,1.71,2.57,3.29,4.27,5.51,7.16,
  8,0,1.04,1.77,2.58,3.27,4.38,5.64,7.47,
  9,0,1.14,1.94,2.78,3.71,4.66,6.87,9.46,
  10,0,1.14,1.94,2.78,3.71,4.66,6.10,8.05,
  11,0,1.28,1.80,2.84,3.61,4.83,6.36,8.53,
  12,0,1.02,1.53,2.21,2.89,3.92,4.95,6.32,
  13,0,1.18,1.96,2.75,3.85,5.03,6.86,9.53,
  14,0,1.46,2.16,3.12,4.34,5.50,7.33,9.93,
  15,0,1.19,2.00,3.04,3.90,5.39,7.34,10.69
), nrow = 15, byrow = TRUE)

times <- c(0, 10, 15, 20, 25, 30, 35, 40)
col_vals <- c("#3C5488FF", "#E64B35FF")


## ============================ Degradation paths ===============================

process_matrix_to_df <- function(mat, pc_name, time_vec) {
  df <- as.data.frame(mat)
  colnames(df) <- c("Unit", paste0("Time_", seq_along(time_vec)))
  
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = starts_with("Time_"),
      names_to = "TimeIndex",
      values_to = "Degradation"
    ) %>%
    dplyr::mutate(
      TimeIndex = as.numeric(gsub("Time_", "", TimeIndex)),
      Time = time_vec[TimeIndex],
      PC = pc_name,
      Unit = factor(as.integer(Unit))
    ) %>%
    dplyr::select(Unit, Time, PC, Degradation)
  
  return(df_long)
}

# Process PC1 and PC2 data separately
df_p1 <- process_matrix_to_df(pc1_raw, "PC1", times)
df_p2 <- process_matrix_to_df(pc2_raw, "PC2", times)

# Combine the two data sets
plot_data <- rbind(df_p1, df_p2)

library(ggsci)

# 3. Plot degradation paths using ggplot2
p <- ggplot(plot_data, aes(x = Time, y = Degradation, group = Unit, color = PC)) +
  # Draw degradation trajectories
  geom_line(linewidth = 0.5, alpha = 0.6) +
  # Add observation points
  geom_point(size = 0.8) +
  # Faceted display by PC
  facet_wrap(~PC, scales = "free_y") + 
  scale_color_manual(values = col_vals) +
  # Axis labels
  labs(
    x = "Time",
    y = "Degradation"
  ) +
  
  # Theme settings
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    # strip.background = element_rect(fill = "grey70", color = "black"),
    strip.text = element_text(family = "Times", face = "plain", size = 10),
    axis.title.x = element_text(family = "Times", face = "plain", size = 12),
    axis.title.y = element_text(family = "Times", face = "plain", size = 12),
    # The following legend settings are retained for consistency.
    legend.title = element_text(family = "Times", face = "plain", size = 12),
    legend.text  = element_text(family = "Times", face = "plain", size = 11)
  )

# 4. Display the plot
print(p)
# Save the figure
ggsave("Real-data/crack/result/1-crack.pdf", plot = p, width = 4, height = 3)



## ================================ Ratio ======================================

# 2. Function for calculating increment-based statistics at each time interval
calculate_mv_ratio <- function(mat, pc_name, time_vec) {
  # Remove the first column, which contains unit IDs
  data_mat <- mat[, -1]
  n_samples <- nrow(data_mat)
  n_times <- ncol(data_mat)
  
  stats_list <- list()
  
  # Loop over each time interval, starting from the second time point
  for(k in 2:n_times) {
    # Current time point
    t_curr <- time_vec[k]
    # Length of the current time interval
    dt <- time_vec[k] - time_vec[k-1]
    
    # Calculate degradation increments for all units in the current interval
    increments <- data_mat[, k] - data_mat[, k-1]
    
    # Calculate cross-sectional statistics.
    # To remove the effect of unequal time intervals, degradation rates are used.
    rates <- increments / dt
    
    sample_mean <- mean(rates)
    sample_var <- var(rates)
    
    # Calculate the empirical mean-standard deviation ratio.
    # If this ratio remains approximately constant, the mean and variance evolve
    # proportionally over time.
    ratio <- sample_mean / sqrt(sample_var)
    
    stats_list[[k-1]] <- data.frame(
      Time = t_curr,
      Mean_Rate = sample_mean,
      Var_Rate = sample_var,
      Ratio = ratio,
      PC = pc_name
    )
  }
  
  return(do.call(rbind, stats_list))
}

# 3. Calculate empirical ratios
stats_pc1 <- calculate_mv_ratio(pc1_raw, "PC1", times)
stats_pc2 <- calculate_mv_ratio(pc2_raw, "PC2", times)
plot_data <- rbind(stats_pc1, stats_pc2)

# 4. Plot the mean-standard deviation ratio over time.
# Under a linear model or a single-scale model, this ratio is expected to be
# approximately constant. A clear increasing or decreasing trend indicates that
# the mean and variance evolve at different rates, motivating the TMVR model.

p_eda <- ggplot(plot_data, aes(x = Time, y = Ratio, color = PC, group = PC)) +
  geom_point(size = 2) +
  geom_line(size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", alpha = 0.5) +
  
  facet_wrap(~PC, scales = "free_y") +
  
  labs(
    # title = "Empirical Evidence for Time-Varying Mean-Variance Ratio",
    # subtitle = "Ratio of (Mean Increment / Std.Dev Increment) over time",
    x = "Time",
    y = "Ratio"
    # y = expression(Ratio ~ (hat(mu) / hat(sigma)))
  ) +
  theme_bw() +
  scale_color_manual(values = col_vals) +
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

print(p_eda)

ggsave("Real-data/crack/result/2-ratio.pdf", plot = p_eda, width = 4, height = 3)


# ==============================================================================
# Additional exploratory analysis 1: Scatter plot of PC1 and PC2 increments
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Prepare increment data
calc_increments <- function(mat, pc_name) {
  df <- as.data.frame(mat)
  colnames(df) <- c("Unit", paste0("T", 1:8))
  df_long <- df %>%
    pivot_longer(cols = -Unit, names_to = "TimeStep", values_to = "Value") %>%
    group_by(Unit) %>%
    mutate(Increment = c(NA, diff(Value))) %>%
    na.omit() %>%
    mutate(PC = pc_name)
  return(df_long)
}

df_inc1 <- calc_increments(pc1_raw, "PC1")
df_inc2 <- calc_increments(pc2_raw, "PC2")

# Merge and reshape the two increment data sets for plotting
df_corr <- left_join(df_inc1, df_inc2, by = c("Unit", "TimeStep")) %>%
  rename(Inc1 = Increment.x, Inc2 = Increment.y)

# Calculate the Pearson correlation coefficient
cor_val <- cor(df_corr$Inc1, df_corr$Inc2)

# 2. Plot the increment correlation
p_corr <- ggplot(df_corr, aes(x = Inc1, y = Inc2)) +
  geom_point(alpha = 0.6, color = "#3C5488FF", size = 2) +
  geom_smooth(method = "lm", color = "#E64B35FF", se = TRUE, fill = "grey80") +
  
  # Add the correlation coefficient annotation
  annotate("text", x = min(df_corr$Inc1), y = max(df_corr$Inc2), 
           label = paste0("Pearson r = ", round(cor_val, 3)), 
           hjust = 0, vjust = 1, size = 4,  family = "serif") +
  
  labs(
    # title = "Correlation of Degradation Increments",
    x = expression(Delta * PC1),
    y = expression(Delta * PC2)
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

print(p_corr)

ggsave("Real-data/crack/result/3-corr.pdf", plot = p_corr, width = 4, height = 3)


# ==============================================================================
# Additional exploratory analysis 2: Visualization of unit heterogeneity
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Calculate the average degradation rate for each unit
# ------------------------------------------------------------------------------

unit_rates <- rbind(df_p1, df_p2) %>%
  group_by(Unit, PC) %>%
  summarise(
    Avg_Rate = (max(Degradation) - min(Degradation)) / max(Time),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 2. Plot boxplots and jittered points by PC
# ------------------------------------------------------------------------------

p_hetero <- ggplot(unit_rates, aes(x = "", y = Avg_Rate, fill = PC)) +
  
  # Boxplot
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(alpha = 0.5, width = 0.4, outlier.shape = NA) +
  
  # Jittered unit-level points
  geom_jitter(width = 0.1, size = 2, aes(color = PC), alpha = 0.8) +
  
  # Add unit labels
  geom_text(aes(label = Unit), position = position_jitter(width = 0.1), 
            vjust = -0.5, size = 2, check_overlap = TRUE, family = "Times") +
  
  # Color settings
  scale_fill_manual(values = c("PC1" = "#3C5488FF", "PC2" = "#E64B35FF")) +
  scale_color_manual(values = c("PC1" = "#3C5488FF", "PC2" = "#E64B35FF")) +
  
  # Faceted display by PC
  facet_wrap(~PC, scales = "fixed") + 
  
  # Axis labels
  labs(
    x = NULL,
    y = "Average degradation rate"
  ) +
  
  # Theme settings
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

print(p_hetero)

ggsave("Real-data/crack/result/4-hetero.pdf", plot = p_hetero, width = 4, height = 3)