library(rjags)
library(coda)
library(dplyr) # 用于数据处理
library(latex2exp)
# ==============================================================================
# 1. 动态构建 JAGS 模型字符串 (保持不变)
# ==============================================================================
get_model_string <- function(drift_type = "Exponential", diff_type = "Exponential") {
  
  # 定义漂移项公式
  if (drift_type == "Power") {
    drift_formula <- "true_diff[j, k] <- pow(t1[k+1], true_scale[j]) - pow(t1[k], true_scale[j])"
    drift_prior   <- "true_scale[j] ~ dgamma(0.01, 0.01)" 
  } else {
    drift_formula <- "true_diff[j, k] <- exp(t1[k+1] * true_scale[j]) - exp(t1[k] * true_scale[j])"
    drift_prior   <- "true_scale[j] ~ dt(0, 0.01, 1)T(0,)" # 使用你之前的 dt 截断先验
  }
  
  # 定义扩散项公式
  if (diff_type == "Power") {
    diff_formula <- "tau_diff[j, k] <- pow(t1[k+1], tau_scale[j]) - pow(t1[k], tau_scale[j])"
    diff_prior   <- "tau_scale[j] ~ dgamma(0.01, 0.01)"
  } else {
    diff_formula <- "tau_diff[j, k] <- exp(t1[k+1] * tau_scale[j]) - exp(t1[k] * tau_scale[j])"
    diff_prior   <- "tau_scale[j] ~ dgamma(0.01, 0.01)"
  }
  
  model_str <- paste0("
  model {
    # --- Time Scale Priors ---
    for (j in 1:p) {
      ", drift_prior, "
      ", diff_prior, "
    }
    
    # --- Pre-computation ---
    for (j in 1:p) {
      for (k in 1:m) {
        ", drift_formula, "
        ", diff_formula, "
      }
    }
    
    # --- Likelihood ---
    for (i in 1:n) {
      eta[1:p, i] ~ dmnorm(mu0, invSIG0)
      for (j in 1:p) {
        for (k in 1:m) {
          mu_xdiff[j, k, i]  <- true_diff[j, k] / eta[j, i]
          var_xdiff[j, k, i] <- (csi[j]^2 * tau_diff[j, k]) / (eta[j, i]^2)
          prec_xdiff[j, k, i] <- 1 / (var_xdiff[j, k, i] + 1.0E-10)
          x_diff[j, k, i] ~ dnorm(mu_xdiff[j, k, i], prec_xdiff[j, k, i])
        }
      }
    }
    
    # --- Other Priors ---
    for (j in 1:p) {
      invcsi[j] ~ dgamma(0.01, 0.01)
      csi[j] <- sqrt(1 / invcsi[j])
      xii[j] ~ dgamma(0.5, 1.0E-6) 
    }
    
    # Delta Matrix Construction (Assuming p=2)
    Delta[1, 1] <- 4 * xii[1]
    Delta[1, 2] <- 0
    Delta[2, 1] <- 0
    Delta[2, 2] <- 4 * xii[2]
    
    invSIG0 ~ dwish(inverse(Delta), 3)
    SIG0 <- inverse(invSIG0[,])
    mu0 ~ dmnorm(mu_mu, invsigma_mu)
    
    # --- Derived ---
    for (j in 1:p) {
      sigma[j] <- sqrt(SIG0[j,j])
    }
    for (j in 1:(p-1)) {
      for (k in (j+1):p) {
        rho[j,k] <- SIG0[j,k] / (sigma[j] * sigma[k])
      }
    }
  }
  ")
  return(model_str)
}

# ==============================================================================
# 2. 数据读取与处理函数 (新增)
# ==============================================================================
load_jags_data_from_csv <- function(csv_file, t_vec) {
  
  # 1. 读取数据
  raw_data <- read.csv(csv_file)
  
  # 确保按 Unit 和 Time 排序
  raw_data <- raw_data %>% arrange(Unit, Time)
  
  # --- [新增] 2. 准备 ggplot2 绘图数据 (长格式) ---
  # 直接在这里处理好，避免在 main.R 里重复读取和转换
  raw_data_long <- raw_data %>%
    pivot_longer(cols = starts_with("PC"), 
                 names_to = "PC", 
                 values_to = "Value")
  
  # 3. 提取维度信息
  units <- unique(raw_data$Unit)
  n <- length(units)
  m <- length(unique(raw_data$Time)) 
  p <- 2 
  
  # 4. 构建 3D 数组 y[p, m, n] (用于 JAGS 增量计算)
  y <- array(NA, dim = c(p, m, n))
  for (i in 1:n) {
    unit_data <- raw_data %>% filter(Unit == units[i])
    y[1, , i] <- unit_data$PC1
    y[2, , i] <- unit_data$PC2
  }
  
  # 5. 计算增量 y_diff
  y_diff <- array(NA, dim = c(p, m, n))
  for (j in 1:p) {
    for (i in 1:n) {
      y_diff[j, 1, i] <- y[j, 1, i] 
      for (k in 2:m) {
        y_diff[j, k, i] <- y[j, k, i] - y[j, k - 1, i]
      }
    }
  }
  
  # 返回列表 (新增了 plot_data)
  return(list(
    x_diff = y_diff,
    t1 = t_vec,
    n = n,
    m = m,
    p = p,
    mu_mu = rep(0, p),
    invsigma_mu = 0.001 * diag(p),
    plot_data = raw_data_long # <--- 这里直接返回长格式数据
  ))
}
# ==============================================================================
# 3. 运行推断函数 (保持不变)
# ==============================================================================
# utility.R 中的 run_inference 函数更新版本

run_inference <- function(data_list, drift_type, diff_type, 
                          inits_params,  # 接收外部传入的初始值列表
                          n.chains = 3, n.adapt = 1000, n.iter = 5000) {
  
  model_string <- get_model_string(drift_type, diff_type)
  n <- data_list$n
  p <- data_list$p
  
  # 辅助函数：处理参数长度 (兼容标量和向量)
  # 如果传入的是 c(1, 100)，它就保持原样；如果传入 1，它就变成 c(1, 1)
  get_param <- function(val, len) {
    if (length(val) == len) return(val)
    return(rep(val, length.out = len))
  }
  
  # 内部生成函数：直接使用传入的 inits_params 中的数值
  gen_inits <- function() {
    list(
      # 1. 随机效应 eta
      eta = matrix(rnorm(n * p, mean = inits_params$eta_mean, sd = inits_params$eta_sd), nrow = p, ncol = n),
      
      # 2. 均值向量 mu0
      mu0 = rep(inits_params$mu0, p),
      
      # 3. 协方差矩阵 invSIG0
      invSIG0 = inits_params$invSIG0_scale * diag(p),
      
      # 4. 时间尺度参数 (支持向量输入)
      true_scale = get_param(inits_params$true_scale, p),
      tau_scale  = get_param(inits_params$tau_scale, p),
      
      # 5. 波动率精度 (支持向量输入，例如 c(1, 100))
      invcsi = get_param(inits_params$invcsi, p)
    )
  }
  
  inits_list <- replicate(n.chains, gen_inits(), simplify = FALSE)
  
  message(paste0(">>> Initializing JAGS Model: Drift=", drift_type, ", Diff=", diff_type))
  
  jags_model <- jags.model(textConnection(model_string), 
                           data = data_list, inits = inits_list, 
                           n.chains = n.chains, n.adapt = n.adapt)
  
  message(">>> Sampling Posterior...")
  samples <- coda.samples(jags_model, 
                          variable.names = c("mu0", "csi", "sigma", "rho", "true_scale", "tau_scale"), 
                          n.iter = n.iter)
  
  dic <- dic.samples(jags_model, n.iter = n.iter)
  
  return(list(samples = samples, dic = dic))
}
# ==============================================================================
# 0. 定义通用的诊断绘图函数 (新增功能)
# ==============================================================================
# ==============================================================================
# 0. 定义通用的诊断绘图函数 (Ergodic Mean + LaTeX Labels)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(latex2exp) # 确保加载了包

# ==============================================================================
# 0. 定义通用的诊断绘图函数 (Ergodic Mean + latex2exp Labels)
# ==============================================================================
save_diagnostics <- function(mcmc_samples, model_name, output_dir = "Diagnostics_Output") {
  
  # --- 1. 路径处理 ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # --- 2. 提取数据并转换为长格式 ---
  mcmc_df <- do.call(rbind, lapply(1:length(mcmc_samples), function(i) {
    df <- as.data.frame(mcmc_samples[[i]])
    df$Iteration <- 1:nrow(df)
    df$Chain <- as.factor(i)
    return(df)
  }))
  
  mcmc_long <- mcmc_df %>%
    pivot_longer(cols = -c(Iteration, Chain), names_to = "Parameter", values_to = "Value")
  
  # --- 3. 筛选核心参数 ---
  target_params <- grep("^(mu0|eta|sigma|rho|scale|true_scale|tau_scale)\\[", 
                        unique(mcmc_long$Parameter), value = TRUE)
  
  plot_data <- mcmc_long %>% 
    filter(Parameter %in% target_params)
  
  # --- 4. 计算遍历均值 (Ergodic Mean) ---
  plot_data <- plot_data %>%
    arrange(Chain, Parameter, Iteration) %>%
    group_by(Chain, Parameter) %>%
    mutate(ErgodicMean = cumsum(Value) / Iteration) %>%
    ungroup()
  
  # --- 5. 构建 LaTeX 标签 (使用 latex2exp 语法) ---
  # 注意：在 R 字符串中，LaTeX 的反斜杠 \ 需要转义为 \\
  # 例如: $\hat{\mu}_1$ 写作 "$\\hat{\\mu}_1$"
  
  plot_data <- plot_data %>%
    mutate(Label = case_when(
      # --- General Model ---
      grepl("mu0\\[1\\]", Parameter) ~ "$\\hat{\\mu}_1$",
      grepl("mu0\\[2\\]", Parameter) ~ "$\\hat{\\mu}_2$",
      grepl("true_scale\\[1\\]", Parameter) ~ "$\\hat{r}_1$",
      grepl("true_scale\\[2\\]", Parameter) ~ "$\\hat{r}_2$",
      grepl("tau_scale\\[1\\]", Parameter)  ~ "$\\hat{s}_1$",
      grepl("tau_scale\\[2\\]", Parameter)  ~ "$\\hat{s}_2$",
      
      # --- Wiener Model (scale 同时代表 r 和 s) ---
      grepl("scale\\[1\\]", Parameter) ~ "$\\hat{r}_1$", 
      grepl("scale\\[2\\]", Parameter) ~ "$\\hat{r}_2$",
      
      # --- Student Model (eta 代表均值) ---
      grepl("eta\\[1\\]", Parameter) ~ "$\\hat{\\mu}_1$",
      grepl("eta\\[2\\]", Parameter) ~ "$\\hat{\\mu}_2$",
      
      # --- Common Params ---
      grepl("sigma\\[1\\]", Parameter) ~ "$\\hat{\\sigma}_1$",
      grepl("sigma\\[2\\]", Parameter) ~ "$\\hat{\\sigma}_2$",
      grepl("rho\\[1,2\\]", Parameter) ~ "$\\hat{\\rho}$",
      
      # --- 默认保留原名 ---
      TRUE ~ Parameter
    ))
  
  # --- 6. 绘图 (ggplot2 + latex2exp) ---
  p <- ggplot(plot_data, aes(x = Iteration, y = ErgodicMean, color = Chain)) +
    geom_line(size = 0.8) +
    # 关键修改：使用 as_labeller(TeX) 来解析 LaTeX 字符串
    facet_wrap(~Label, scales = "free_y", 
               labeller = as_labeller(latex2exp::TeX, default = label_parsed), 
               ncol = 3) + 
    scale_color_brewer(palette = "Set1") + 
    labs(x = "Iteration", y = "Ergodic mean"#, 
         # title = paste("Posterior Convergence:", model_name)
         ) +
    theme_bw(base_size = 14) + 
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  # --- 7. 保存 ---
  clean_name <- gsub(" ", "_", model_name)
  clean_name <- gsub("[()]", "", clean_name)
  pdf_filename <- paste0("Ergodic_", clean_name, ".pdf")
  full_path <- file.path(output_dir, pdf_filename)
  
  ggsave(full_path, plot = p, width = 9, height = 4, device = cairo_pdf)
  
  message(paste0("   -> Ergodic plot saved to: ", full_path))

}



# ==============================================================================
# 1. 核心仿真函数 (高度封装)
# ==============================================================================

get_fitted_paths <- function(params, time_seq, n_sim = 5000) {
  p <- 2
  # --- A. 参数提取 (利用 unlist 简化) ---
  mu    <- c(as.numeric(params$mu_1), as.numeric(params$mu_2))
  sigma <- c(as.numeric(params$sigma_1), as.numeric(params$sigma_2))
  rho   <- as.numeric(params$rho_12)
  r_vec <- c(as.numeric(params$r_1), as.numeric(params$r_2))
  
  SIG0  <- diag(sigma) %*% matrix(c(1, rho, rho, 1), 2) %*% diag(sigma)
  
  # --- B. 矩阵化仿真 (无显式循环) ---
  # 1. 一次性生成所有随机效应 eta (n_sim x p)
  eta_mat <- abs(mvrnorm(n_sim, mu, SIG0)) 
  
  # 2. 对每个 PC 计算路径统计量
  # lapply 遍历 1 到 p，自动返回合并后的数据框
  res_list <- lapply(1:p, function(j) {
    # 时间项向量 (Length: time_seq)
    lambda_t <- exp(time_seq * r_vec[j]) - 1
    
    # 路径矩阵 (行: n_sim, 列: time_seq)
    # 利用外积: (1/eta) * lambda_t^T
    paths_mat <- (1 / eta_mat[, j]) %*% t(lambda_t)
    
    data.frame(
      Time = time_seq,
      Mean = colMeans(paths_mat),
      Lower = apply(paths_mat, 2, quantile, probs = 0.05),
      Upper = apply(paths_mat, 2, quantile, probs = 0.95),
      PC = paste0("PC", j)
    )
  })
  
  do.call(rbind, res_list)
}
