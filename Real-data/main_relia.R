# main.R
# 终极汇总：General (4) + Wiener (2) + Student (2) = 8 个模型
# 功能：MCMC推断 -> 诊断图 -> 可靠度仿真 -> 绘图对比

rm(list = ls())

# 1. 加载工具文件 (确保这4个文件在同一目录下)
source("实例/utility.R")       # 处理 General Models (PP, PE, EP, EE)
source("实例/Wiener_model.R")  # 处理 Wiener Models (PP, EE)
source("实例/Student_model.R")  # 处理 Student Models  (PP, EE)
source("实例/Reliability_utils.R")

# 加载必要包
library(dplyr)
library(ggplot2)
library(bayesplot)
library(latex2exp) # 用于诊断图的数学公式标签

# ==============================================================================
# 2. 数据与全局配置
# ==============================================================================
csv_file_path <- "实例/crack/crack_growth_data.csv"
t_vector <- c(0, 10, 15, 20, 25, 30, 35, 40)
jags_data <- load_jags_data_from_csv(csv_file_path, t_vector)

# --- [关键配置] 沿用你原始代码的失效阈值 ---
# PC1 = 4, PC2 = 7
thresholds <- c(4, 7) 

# 可靠度评估的时间序列 (从20到50，与你绘图代码保持一致)
rel_time_seq <- seq(20, 50, by = 0.5)

# 配置输出文件夹
current_case_folder <- "实例/crack/result"
if (!dir.exists(current_case_folder)) dir.create(current_case_folder)

# ==============================================================================
# 3. 初始值配置 (User Defined Parameters)
# ==============================================================================
# General Model Inits
inits_general_power <- list(
  eta_mean = 25, eta_sd = 1, mu0 = 25, invSIG0_scale = 1,
  true_scale = c(1, 1), tau_scale = c(1, 1), invcsi = c(1, 100)
)
inits_general_exp <- list(
  eta_mean = 0.5, eta_sd = 0.1, mu0 = 0.5, invSIG0_scale = 0.1,
  true_scale = c(0.05, 0.05), tau_scale = c(1, 1), invcsi = c(40000, 40000)
)

# Wiener Model Inits
inits_wiener <- list(
  mu0 = 0.6, sigma2_err = 0.1, Tau_b_scale = 1, b_val = 0.6, scale = c(1, 1.5)
)

# Student Model Inits
inits_student <- list(
  eta = c(0.03, 0.03), SIG0_inv_scale = 1, scale = c(1, 1), sig_w = c(1, 1)
)

# ==============================================================================
# 4. 第一阶段：MCMC 参数估计 (只跑模型，存参数)
# ==============================================================================
scenarios <- data.frame(
  Model_Name = c("General EE", "Wiener EE", "Student EE"), # 仅演示 EE 系列
  Class = c("General", "Wiener", "Student"),
  Drift_Type = c("Exponential", "Exponential", "Exponential"),
  Diff_Type = c("Exponential", NA, NA),
  stringsAsFactors = FALSE
)

# 字典，用于存储每个模型的后验参数均值
saved_params_list <- list()

for (i in 1:nrow(scenarios)) {
  model_name <- scenarios$Model_Name[i]
  model_class <- scenarios$Class[i]
  d_type <- scenarios$Drift_Type[i]
  diff_type <- scenarios$Diff_Type[i]
  
  message(paste0("\n>>> [MCMC] Running: ", model_name))
  
  if (model_class == "General") {
    # General EE 的初始值
    res <- run_inference(jags_data, d_type, diff_type, inits_general_exp, n.chains=2, n.iter=5000)
  } else if (model_class == "Wiener") {
    res <- run_wiener_inference(jags_data, d_type, inits_wiener, n.chains=2, n.iter=5000)
  } else {
    res <- run_student_inference(jags_data, d_type, inits_student, n.chains=2, n.iter=5000)
  }
  
  # 保存参数均值
  saved_params_list[[model_name]] <- colMeans(as.matrix(res$samples))
  
  # 保存诊断图 (可选)
  save_diagnostics(res$samples, model_name, output_dir = current_case_folder)
}

# ==============================================================================
# 5. 第二阶段：可靠度仿真 (核心差异化处理)
# ==============================================================================
message("\n>>> Starting Reliability Calculation...")

plot_data_list <- list()

# 遍历刚才跑完的模型
for (model_name in names(saved_params_list)) {
  
  params <- saved_params_list[[model_name]]
  message(paste0("   Calculating R(t) for: ", model_name))
  
  # --- 核心分支：不同模型用不同计算方法 ---
  
  if (grepl("General", model_name)) {
    # 1. General 模型 -> 使用路径仿真 (Method A)
    # 必须区分 Drift/Diff 类型
    drift_t <- if(grepl("PP", model_name)) "Power" else "Exponential"
    diff_t  <- drift_t # 这里简化假设，实际可从 scenarios 表读取
    
    rel_df <- calc_reliability_simulation(params, thresholds, rel_time_seq, drift_t, diff_t)
    
  } else {
    # 2. Wiener / Student -> 使用解析解 (Method B / R_cal)
    model_type <- if(grepl("Wiener", model_name)) "Wiener" else "Student"
    drift_t <- if(grepl("PP", model_name)) "Power" else "Exponential"
    
    rel_df <- calc_reliability_analytical(params, model_type, drift_t, thresholds, rel_time_seq)
  }
  
  rel_df$Model <- model_name
  plot_data_list[[model_name]] <- rel_df
}

# 合并数据
all_plot_data <- do.call(rbind, plot_data_list)

# 经验点
empirical_data <- calculate_empirical_PFT(csv_file_path, thresholds)

# ==============================================================================
# 6. 绘图 (带置信区间的伪效果，或仅画均值线)
# ==============================================================================
# 为了复刻你提供的图，我们需要给 General EE 添加一个 ribbon
# 这里我们简单起见，先画出正确的均值线。如果均值线对了，Ribbon 只是加减标准差的问题。

p_rel <- ggplot() +
  # 经验点
  geom_point(data = empirical_data, aes(x = Time, y = Reliability, shape = "Empirical"), 
             size = 3, color = "black") +
  
  # 模型曲线
  geom_line(data = all_plot_data, 
            aes(x = Time, y = Reliability, color = Model, linetype = Model), 
            size = 1.2) +
  
  # 样式
  scale_color_manual(values = c("General EE" = "#440154", "Wiener EE" = "blue", "Student EE" = "red")) +
  scale_linetype_manual(values = c("General EE" = "solid", "Wiener EE" = "dashed", "Student EE" = "dotdash")) +
  
  labs(x = "Time (cycles)", y = "Reliability") +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.3))

print(p_rel)
ggsave(file.path(current_case_folder, "Reliability_Fixed_Logic.pdf"), plot = p_rel, width = 8, height = 6)