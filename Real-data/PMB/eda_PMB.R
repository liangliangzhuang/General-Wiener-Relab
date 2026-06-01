# ==============================================================================
# Case Study II: Permanent Magnet Brake
# Results for Figure 1 in the main text, dated 12.09
# ==============================================================================

rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(latex2exp) # Load this package for LaTeX-style axis labels

# ==============================================================================
# 1. Load data
# ==============================================================================

# Make sure PMB_data.csv is located in the specified directory
df_pmb <- read.csv("Real-data/PMB/PMB_data.csv")

# Convert Unit to a factor for plotting
df_pmb$Unit <- as.factor(df_pmb$Unit)

# Convert the data into long format for faceted plotting with ggplot2.
# Although the original CSV file is already in long format, PC1 and PC2 are
# combined into a single column here for convenient visualization.
df_plot <- df_pmb %>%
  pivot_longer(
    cols = c("PC1", "PC2"),
    names_to = "PerformanceCharacteristic",
    values_to = "Degradation"
  )

# Specify the display order of the performance characteristics
df_plot$PerformanceCharacteristic <- factor(
  df_plot$PerformanceCharacteristic,
  levels = c("PC1", "PC2")
) # labels = c("Release Time (PC1)", "Braking Torque (PC2)")

# ==============================================================================
# 2. Plot degradation paths
# ==============================================================================

col_vals <- c("#3C5488FF", "#E64B35FF")

p1 <- ggplot(
  df_plot,
  aes(
    x = Time,
    y = Degradation,
    color = PerformanceCharacteristic,
    group = Unit
  )
) +
  geom_line(size = 0.6, alpha = 0.7) +
  geom_point(size = 0.8) +
  # Faceted display by performance characteristic
  facet_wrap(~PerformanceCharacteristic, scales = "free_y") +
  scale_color_manual(values = col_vals) +
  # Axis labels
  labs(
    # title = "Degradation Paths of Permanent Magnet Brake",
    x = "Time",
    y = "Degradation Value",
    color = "Unit"
  ) +
  
  # Theme settings for publication-style figures
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

print(p1)

ggsave(
  "Real-data/PMB/result/1-pmb_path.pdf",
  plot = p1,
  width = 4,
  height = 2.6
)


# ==============================================================================
# Exploratory data analysis for PMB data
# ==============================================================================

# ==============================================================================
# 0. Reconstruct matrix-format data for increment calculation
# ==============================================================================

# Extract basic information
units <- unique(df_pmb$Unit)
n_units <- length(units)
times <- sort(unique(df_pmb$Time))
n_times <- length(times)

# Initialize matrices.
# Rows correspond to units, and columns correspond to Unit ID and observation times.
mat_pc1 <- matrix(NA, nrow = n_units, ncol = n_times + 1) # +1 for ID column
mat_pc2 <- matrix(NA, nrow = n_units, ncol = n_times + 1)

# Fill the first column with unit IDs
mat_pc1[, 1] <- 1:n_units
mat_pc2[, 1] <- 1:n_units

# Fill the degradation measurements for each unit.
# PC1 denotes release time, and PC2 denotes braking torque.
for(i in 1:n_units) {
  u_data <- df_pmb %>% filter(Unit == i) %>% arrange(Time)
  
  mat_pc1[i, 2:(n_times+1)] <- u_data$PC1
  mat_pc2[i, 2:(n_times+1)] <- u_data$PC2
}

# Define color values
col_vals <- c("#3C5488FF", "#E64B35FF")
names(col_vals) <- c("PC1", "PC2")


# ==============================================================================
# 1. Time-varying mean-variance ratio analysis
# ==============================================================================

calculate_mv_ratio <- function(mat, pc_name, time_vec) {
  data_mat <- mat[, -1] # Remove the ID column
  stats_list <- list()
  
  # Calculate increments starting from the second time point
  for(k in 2:length(time_vec)) {
    dt <- time_vec[k] - time_vec[k-1] # For PMB data, dt = 1 when Time is an index
    increments <- data_mat[, k] - data_mat[, k-1]
    
    # For PC2, braking torque decreases over time, so the increments are negative.
    # Absolute increments are used to characterize the magnitude of degradation.
    rates <- abs(increments) / dt
    
    sample_mean <- mean(rates)
    sample_var <- var(rates)
    
    # Ratio = Mean / SD
    ratio <- sample_mean / sqrt(sample_var)
    
    stats_list[[k-1]] <- data.frame(
      Time = time_vec[k],
      Ratio = ratio,
      PC = pc_name
    )
  }
  
  return(do.call(rbind, stats_list))
}

stats_pc1 <- calculate_mv_ratio(mat_pc1, "PC1", times)
stats_pc2 <- calculate_mv_ratio(mat_pc2, "PC2", times)

plot_data_ratio <- rbind(stats_pc1, stats_pc2)

p_ratio <- ggplot(
  plot_data_ratio,
  aes(x = Time, y = Ratio, color = PC, group = PC)
) +
  geom_point(size = 1) +
  geom_line(size = 0.6) +
  # Add a linear trend line to help assess time-varying behavior
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", alpha = 0.5) +
  
  facet_wrap(~PC, scales = "free_y") +
  scale_color_manual(values = col_vals) +
  
  labs(
    x = "Time",
    # y = expression(Ratio ~ (hat(mu) / hat(sigma)))
    y = "Ratio"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(family = "Times", face = "plain", size = 10),
    axis.title.x = element_text(family = "Times", size = 12),
    axis.title.y = element_text(family = "Times", size = 12),
    axis.text = element_text(family = "Times")
  )

print(p_ratio)

ggsave(
  "Real-data/PMB/result/2-pmb_ratio.pdf",
  plot = p_ratio,
  width = 4,
  height = 2.6
)


# ==============================================================================
# 2. Correlation of degradation rates
# ==============================================================================

# In this exploratory analysis, the degradation rate of each unit is estimated
# by fitting a simple linear regression: Value ~ Time.
# The slope represents the average degradation speed of each unit.

unit_slopes <- df_pmb %>%
  group_by(Unit) %>%
  summarise(
    # Slope for PC1, i.e., release time
    Rate_PC1 = coef(lm(PC1 ~ Time))[2],
    # Slope for PC2, i.e., braking torque
    Rate_PC2 = coef(lm(PC2 ~ Time))[2]
  ) %>%
  ungroup() %>%
  # Use absolute values to compare degradation magnitudes
  mutate(
    Abs_Rate_PC1 = abs(Rate_PC1),
    Abs_Rate_PC2 = abs(Rate_PC2)
  )

# Calculate the Pearson correlation coefficient
cor_val <- cor(unit_slopes$Abs_Rate_PC1, unit_slopes$Abs_Rate_PC2)

print(paste("Rate Correlation:", cor_val))

p_corr_rate <- ggplot(unit_slopes, aes(x = Abs_Rate_PC1, y = Abs_Rate_PC2)) +
  # Add a fitted trend line
  geom_smooth(method = "lm", color = "#E64B35FF", fill = "grey85", alpha = 0.5) +
  # Plot unit-level points
  geom_point(size = 3, color = "#3C5488FF", alpha = 0.8) +
  # Add unit labels
  geom_text(aes(label = Unit), vjust = -1, size = 3, family = "Times") +
  
  # Add correlation coefficient annotation
  annotate(
    "text",
    x = min(unit_slopes$Abs_Rate_PC1),
    y = max(unit_slopes$Abs_Rate_PC2),
    label = paste0("Pearson r = ", sprintf("%.3f", cor_val)),
    hjust = 0,
    vjust = 1,
    size = 3,
    family = "Times"
  ) +
  labs(
    x = TeX(r"(|$\Delta$ PC1|)"),
    y = TeX(r"(|$\Delta$ PC2|)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(family = "Times", face = "plain", size = 10),
    axis.title.y = element_text(family = "Times", size = 12),
    axis.text = element_text(family = "Times"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

print(p_corr_rate)

ggsave(
  "Real-data/PMB/result/3-pmb_corr.pdf",
  plot = p_corr_rate,
  width = 4,
  height = 2.6
)


# ==============================================================================
# 3. Unit heterogeneity
# ==============================================================================

# ------------------------------------------------------------------------------
# Prepare increment data for PC1 and PC2
# ------------------------------------------------------------------------------

calc_increments <- function(mat, pc_name, time_vec) {
  df <- as.data.frame(mat)
  colnames(df) <- c("Unit", paste0("T", seq_along(time_vec)))
  
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = -Unit,
      names_to = "TimeStep",
      values_to = "Value"
    ) %>%
    dplyr::group_by(Unit) %>%
    dplyr::mutate(
      Increment = c(NA, diff(Value)),
      TimeIndex = as.numeric(gsub("T", "", TimeStep)),
      Time = time_vec[TimeIndex]
    ) %>%
    tidyr::drop_na(Increment) %>%
    dplyr::mutate(PC = pc_name) %>%
    dplyr::ungroup()
  
  return(df_long)
}

df_inc1 <- calc_increments(mat_pc1, "PC1", times)
df_inc2 <- calc_increments(mat_pc2, "PC2", times)

# Calculate the average absolute degradation rate for each unit
unit_rates <- rbind(df_inc1, df_inc2) %>%
  group_by(Unit, PC) %>%
  summarise(
    # Rate = Mean(|Increment|) / dt. Here dt = 1.
    Avg_Rate = mean(abs(Increment)),
    .groups = "drop"
  )

p_hetero <- ggplot(unit_rates, aes(x = "", y = Avg_Rate, fill = PC)) +
  # Boxplot
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(alpha = 0.5, width = 0.4, outlier.shape = NA) +
  
  # Jittered points
  geom_jitter(width = 0.1, size = 2, aes(color = PC), alpha = 0.8) +
  
  # Add unit labels
  geom_text(
    aes(label = Unit),
    position = position_jitter(width = 0.1),
    vjust = -0.5,
    size = 2,
    check_overlap = TRUE,
    family = "Times"
  ) +
  
  scale_fill_manual(values = col_vals) +
  scale_color_manual(values = col_vals) +
  
  # Use fixed scales for direct comparison between PC1 and PC2
  facet_wrap(~PC, scales = "fixed") +
  
  labs(
    x = NULL,
    y = "Average degradation rate (abs)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(family = "Times", face = "plain", size = 10),
    axis.title.y = element_text(family = "Times", size = 12),
    axis.text = element_text(family = "Times"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

print(p_hetero)

ggsave(
  "Real-data/PMB/result/4-pmb_hetero.pdf",
  plot = p_hetero,
  width = 4,
  height = 2.6
)