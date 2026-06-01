# Student_model.R
library(rjags)
library(coda)

# ==============================================================================
# 1. 构建 Student 模型字符串 (Power 或 Exponential)
# ==============================================================================
get_student_model_string <- function(type = "Exponential") {
  
  # 根据类型选择时间尺度函数
  if (type == "Power") {
    # Power: t^scale
    time_func_calc <- "t_scale[j,k] <- pow(t[k], scale[j])"
  } else {
    # Exponential: exp(t * scale)
    time_func_calc <- "t_scale[j,k] <- exp(t[k] * scale[j])"
  }
  
  model_str <- paste0("
  model {
    # --- Priors ---
    # 均值 eta 的先验
    eta[1:p] ~ dmnorm(mu0[], prec_eta[,])
    
    # 时间尺度参数 scale 和 波动率 sig.w
    for(j in 1:p){
      scale[j] ~ dgamma(0.001, 0.001)
      sig.w[j] ~ dgamma(0.001, 0.001)
    }
    
    # --- Pre-computation of Time Increments ---
    # 注意：JAGS 中 t 的索引通常需要从 1 开始对应 R 中的向量
    for(j in 1:p){
      for(k in 1:(m+1)){
        ", time_func_calc, "
      }
      for(k in 1:m){
        t_diff[j,k] <- t_scale[j,k+1] - t_scale[j,k]
      }
    }
    
    # 辅助零向量
    for(j in 1:p){
      zeros[j] <- 0
    }
    
    # --- Student-t Process Construction (Scale Mixture of Normals) ---
    for(i in 1:n){
      # Gamma 分布用于模拟 t 分布的权重
      tau[i] ~ dgamma(v/2, v/2)
      
      # 随机效应 z (非中心化参数化)
      z[i,1:p] ~ dmnorm(zeros[], SIG0_inv[,])
      
      for(j in 1:p){
        # 构造 theta
        theta[i,j] <- eta[j] + z[i,j] / sqrt(tau[i])
      }
      
      # 观测模型
      for(k in 1:m){
        for(j in 1:p){
          # 均值 = theta * delta_time
          # 精度 = tau / (sigma^2 * delta_time)
          # 添加 1e-10 保护
          mu_y[i,j,k] <- theta[i,j] * t_diff[j,k]
          prec_y[i,j,k] <- tau[i] / (pow(sig.w[j], 2) * t_diff[j,k] + 1.0E-10)
          
          y.diff[j,k,i] ~ dnorm(mu_y[i,j,k], prec_y[i,j,k])
        }
      }
    }
    
    # --- Derived Quantities ---
    # 计算 SIG0 (协方差矩阵)
    for (j in 1:p){
      sigma[j] <- sqrt(SIG0[j,j])
    }
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        rho[j,k] <- SIG0[j,k] / (sigma[j]*sigma[k])
      }
    }
    
    # --- HIW Prior ---
    for (j in 1:p) {
      xii[j] ~ dgamma(0.5, 1.0E-6)
    }
    Delta[1, 1] <- 4 * xii[1]
    Delta[1, 2] <- 0
    Delta[2, 1] <- 0
    Delta[2, 2] <- 4 * xii[2]
    
    SIG0_inv ~ dwish(inverse(Delta), 3)
    SIG0 <- inverse(SIG0_inv[,])
  }
  ")
  return(model_str)
}

# ==============================================================================
# 2. 运行 Student 推断函数
# ==============================================================================
run_student_inference <- function(data_list, type, inits_params, 
                                  n.chains = 2, n.adapt = 1000, n.iter = 5000) {
  
  model_string <- get_student_model_string(type)
  n <- data_list$n
  p <- data_list$p
  
  # --- 清洗数据 ---
  obs_data <- if(!is.null(data_list$y.diff)) data_list$y.diff else data_list$x_diff
  
  student_data <- list(
    y.diff = obs_data,
    n = n, 
    p = p, 
    m = data_list$m,
    t = data_list$t1,      # Student 模型里叫 t
    v = 5,                 # 自由度通常固定为 5 或作为参数，这里根据你脚本固定为 5
    mu0 = rep(0, p),
    prec_eta = diag(0.001, p)
  )
  
  # --- 初始值 ---
  gen_inits <- function() {
    list(
      eta = inits_params$eta,
      SIG0_inv = inits_params$SIG0_inv_scale * diag(1, p),
      tau = rep(1, n),
      z = matrix(0, nrow = n, ncol = p),
      scale = inits_params$scale,
      sig.w = inits_params$sig_w
    )
  }
  
  inits_list <- replicate(n.chains, gen_inits(), simplify = FALSE)
  
  message(paste0(">>> Initializing Student Model: Type=", type))
  
  jags_model <- jags.model(textConnection(model_string),
                           data = student_data,
                           inits = inits_list,
                           n.chains = n.chains,
                           n.adapt = n.adapt)
  
  message(">>> Burning in...")
  update(jags_model, n.adapt)
  
  message(">>> Sampling Posterior...")
  samples <- coda.samples(jags_model, 
                          variable.names = c("eta", "scale", "sig.w", "sigma", "rho"), 
                          n.iter = n.iter)
  
  dic <- dic.samples(jags_model, n.iter = n.iter)
  
  return(list(samples = samples, dic = dic))
}