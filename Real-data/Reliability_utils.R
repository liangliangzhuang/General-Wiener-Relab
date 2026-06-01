# Reliability_utils.R
library(MASS)
library(statmod) # 必须加载，用于 rinvgauss
library(dplyr)

# ==============================================================================
# 1. 经验失效点计算 (保持不变)
# ==============================================================================
# ==============================================================================
# 1. 计算经验失效时间 (Empirical PFT) - 增强鲁棒版
# ==============================================================================
calculate_empirical_PFT <- function(csv_file, thresholds) {
  df <- read.csv(csv_file)
  units <- unique(df$Unit)
  pft_results <- numeric(length(units))
  
  for (i in 1:length(units)) {
    sub_df <- subset(df, Unit == units[i])
    t_vals <- sub_df$Time
    
    # 尝试 NLS 拟合 (首选)
    res_t1 <- NA
    res_t2 <- NA
    
    # --- PC1 ---
    tryCatch({
      fit1 <- nls(sub_df$PC1 ~ c * (exp(a * t_vals) - 1) + sub_df$PC1[1], 
                  start = list(c = 0.1, a = 0.05), 
                  control = nls.control(warnOnly=TRUE, maxiter=100))
      c1 <- coef(fit1)["c"]; a1 <- coef(fit1)["a"]
      res_t1 <- log((thresholds[1] - sub_df$PC1[1]) / c1 + 1) / a1
    }, error = function(e) {
      # 如果 NLS 失败，使用线性插值兜底
      # 找到第一个超过阈值的点，插值计算时间
      idx <- which(sub_df$PC1 >= thresholds[1])[1]
      if (!is.na(idx) && idx > 1) {
        y1 <- sub_df$PC1[idx-1]; y2 <- sub_df$PC1[idx]
        t1 <- t_vals[idx-1];     t2 <- t_vals[idx]
        res_t1 <<- t1 + (t2-t1) * (thresholds[1]-y1) / (y2-y1)
      }
    })
    
    # --- PC2 ---
    tryCatch({
      fit2 <- nls(sub_df$PC2 ~ c * (exp(b * t_vals) - 1) + sub_df$PC2[1], 
                  start = list(c = 0.5, b = 0.05), 
                  control = nls.control(warnOnly=TRUE, maxiter=100))
      c2 <- coef(fit2)["c"]; b2 <- coef(fit2)["b"]
      res_t2 <- log((thresholds[2] - sub_df$PC2[1]) / c2 + 1) / b2
    }, error = function(e) {
      # 线性插值兜底
      idx <- which(sub_df$PC2 >= thresholds[2])[1]
      if (!is.na(idx) && idx > 1) {
        y1 <- sub_df$PC2[idx-1]; y2 <- sub_df$PC2[idx]
        t1 <- t_vals[idx-1];     t2 <- t_vals[idx]
        res_t2 <<- t1 + (t2-t1) * (thresholds[2]-y1) / (y2-y1)
      }
    })
    
    # 取两者最小值，忽略 NA
    times <- c(res_t1, res_t2)
    if (all(is.na(times))) {
      pft_results[i] <- NA
    } else {
      pft_results[i] <- min(times, na.rm = TRUE)
    }
  }
  
  valid_pft <- sort(na.omit(pft_results))
  n <- length(valid_pft)
  
  if (n == 0) {
    warning("无法计算任何单元的失效时间，请检查数据或阈值！")
    return(data.frame(Time = numeric(0), Reliability = numeric(0)))
  }
  
  # 使用中位秩 (Median Rank) 计算经验可靠度，比 (n-1)/n 更准确
  # R(t) = 1 - (i - 0.3) / (n + 0.4)
  rank_probs <- 1 - (1:n - 0.3) / (n + 0.4)
  
  return(data.frame(Time = valid_pft, Reliability = rank_probs))
}
# ==============================================================================
# 2. 方法A: 路径仿真 (对应你代码中的 General 模型部分)
# ==============================================================================
# 适用于: General EE/PP/PE/EP
# 逻辑: for k in time_steps ... cumulative + rnorm(...)
calc_reliability_simulation <- function(params, thresholds, time_seq, drift_type, diff_type, n_sim=5000) {
  
  p <- 2
  # 提取参数
  mu0 <- c(params["mu0[1]"], params["mu0[2]"])
  sigma <- c(params["sigma[1]"], params["sigma[2]"])
  rho <- params["rho[1,2]"]
  csi <- c(params["csi[1]"], params["csi[2]"])
  r_scale <- c(params["true_scale[1]"], params["true_scale[2]"])
  s_scale <- c(params["tau_scale[1]"], params["tau_scale[2]"])
  
  SIG0 <- diag(sigma) %*% matrix(c(1, rho, rho, 1), 2) %*% diag(sigma)
  
  T_sys <- numeric(n_sim)
  
  # 辅助函数：计算时间增量
  get_incr <- function(t_now, t_old, sc, type) {
    if(type == "Exponential") return(exp(t_now * sc) - exp(t_old * sc))
    else return(t_now^sc - t_old^sc)
  }
  
  # --- 仿真循环 (复刻你的代码) ---
  for (q in 1:n_sim) {
    eta <- mvrnorm(1, mu0, SIG0) # eta ~ N(mu0, Sigma)
    cum_deg <- c(0, 0)
    fail_t <- max(time_seq)
    
    for (k in 1:(length(time_seq) - 1)) {
      t_curr <- time_seq[k+1]
      t_prev <- time_seq[k]
      
      for (j in 1:p) {
        dt_true <- get_incr(t_curr, t_prev, r_scale[j], drift_type)
        dt_tau  <- get_incr(t_curr, t_prev, s_scale[j], diff_type)
        
        # [核心公式] General 模型: mean = delta / eta
        mean_inc <- dt_true / eta[j]
        sd_inc   <- csi[j] * sqrt(dt_tau) / abs(eta[j])
        
        cum_deg[j] <- cum_deg[j] + rnorm(1, mean_inc, sd_inc)
      }
      
      if (any(cum_deg >= thresholds)) {
        fail_t <- t_curr
        break
      }
    }
    T_sys[q] <- fail_t
  }
  
  # 计算 R(t)
  R_vals <- sapply(time_seq, function(x) mean(T_sys > x))
  return(data.frame(Time = time_seq, Reliability = R_vals))
}

# ==============================================================================
# 3. 方法B: 逆高斯解析解 (对应你代码中的 R_cal 函数)
# ==============================================================================
# 适用于: Wiener, Student
# 逻辑: rinvgauss() -> inverse_time_function()
calc_reliability_analytical <- function(params, model_type, drift_type, thresholds, time_seq, n_sim=5000) {
  
  p <- 2
  # 提取参数 (注意 Wiener 和 Student 的变量名映射)
  if (model_type == "Wiener") {
    mu_vec <- c(params["mu0[1]"], params["mu0[2]"])
    sigma_vec <- c(params["sigma[1]"], params["sigma[2]"])
    rho <- params["rho[1,2]"]
    # Wiener 代码中 sigma2_err 对应波动率系数，相当于 R_cal 中的 r_delta
    delta_vec <- c(params["sigma2_err[1]"], params["sigma2_err[2]"]) 
    gamma_vec <- c(params["scale[1]"], params["scale[2]"])
    v_param <- Inf # Wiener 对应正态，自由度无穷大
    
  } else { # Student
    mu_vec <- c(params["eta[1]"], params["eta[2]"])
    sigma_vec <- c(params["sigma[1]"], params["sigma[2]"])
    rho <- params["rho[1,2]"]
    delta_vec <- c(params["sig.w[1]"], params["sig.w[2]"])
    gamma_vec <- c(params["scale[1]"], params["scale[2]"])
    v_param <- 5 # 假设自由度固定为 5 (或从参数中提取)
  }
  
  SIG0 <- diag(sigma_vec) %*% matrix(c(1, rho, rho, 1), 2) %*% diag(sigma_vec)
  
  # --- 按照 R_cal 逻辑计算 ---
  ft_star <- numeric(n_sim)
  
  for (b in 1:n_sim) {
    # 1. 生成 tau (Student模型会有随机性，Wiener则为1)
    if (is.infinite(v_param)) {
      r_tau <- 1
    } else {
      r_tau <- rgamma(1, shape = v_param/2, scale = 2/v_param)
    }
    
    # 2. 生成 theta (漂移率)
    # [注意] 你的 R_cal 中用 r_SIG0 / r_tau
    r_theta <- mvrnorm(1, mu_vec, SIG0 / r_tau)
    
    # 3. 生成线性尺度下的失效时间 (Inverse Gaussian)
    # mean = D / theta
    # shape = D^2 * tau / delta^2
    lambda_t <- numeric(p)
    for (j in 1:p) {
      m_val <- thresholds[j] / r_theta[j]
      # 防止负漂移导致均值为负 (IG分布要求均值为正)
      m_val <- abs(m_val) 
      s_val <- (thresholds[j]^2 * r_tau) / (delta_vec[j]^2)
      lambda_t[j] <- rinvgauss(1, mean = m_val, shape = s_val)
    }
    
    # 4. 转换回真实时间
    # T = Lambda^(-1)(T_linear)
    if (drift_type == "Exponential") {
      # t = log(lambda + 1) / gamma
      fail_times <- log(lambda_t + 1) / gamma_vec
    } else { # Power
      # t = lambda ^ (1/gamma)
      fail_times <- lambda_t ^ (1/gamma_vec)
    }
    
    ft_star[b] <- min(fail_times)
  }
  
  # 计算 R(t)
  R_vals <- sapply(time_seq, function(x) mean(ft_star > x))
  return(data.frame(Time = time_seq, Reliability = R_vals))
}