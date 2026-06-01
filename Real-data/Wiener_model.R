# Wiener_model.R
library(rjags)
library(coda)

# ==============================================================================
# 1. 构建 Wiener 模型字符串 (Power 或 Exponential)
# ==============================================================================
get_wiener_model_string <- function(type = "Exponential") {
  
  # 根据类型选择时间尺度函数
  if (type == "Power") {
    time_func <- "(pow(t1[k+1], scale[j]) - pow(t1[k], scale[j]))"
  } else {
    time_func <- "(exp(t1[k+1]* scale[j]) - exp(t1[k]* scale[j]))"
  }
  
  model_str <- paste0("
  model {
    # --- Priors ---
    # mu0 给定多元正态先验
    mu0[1:p] ~ dmnorm(m0[], Tau0[,])
    
    # 时间尺度参数 scale
    for (j in 1:p) {
      scale[j] ~ dgamma(0.001, 0.001)
    }
    
    # 观测误差扩散参数先验 (Volatility)
    for (j in 1:p){
      sigma2_err[j] ~ dgamma(0.001, 0.001)
      for (k in 1:m){
        # 精度 = 1 / (sigma^2 * delta_time)
        prec[j,k] <- 1/(sigma2_err[j]^2 * (", time_func, " + 1.0E-10))
      }
    }
    
    # --- Likelihood & Random Effects ---
    for (i in 1:n){
      b[i,1:p] ~ dmnorm(mu0[], Tau_b[,]) # Random drift b around mu0
      
      for (j in 1:p){
        for (k in 1:m){
          # 均值 = delta_time * drift
          mu[i,j,k] <- ", time_func, " * b[i,j]
          y.diff[j,k,i] ~ dnorm(mu[i,j,k], prec[j,k])
        }
      }
    }
    
    # --- Derived Quantities (Covariance & HIW) ---
    # 计算随机效应标准差 sigma 和 相关系数 rho
    for (j in 1:p){
      sigma[j] <- sqrt(SIG_b[j,j])
    }
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        rho[j,k] <- SIG_b[j,k] / (sigma[j]*sigma[k])
      }
    }
    
    # HIW Prior Construction (用于 Tau_b)
    for (j in 1:p) {
      xii[j] ~ dgamma(0.5, 1.0E-6) 
    }
    
    Delta[1, 1] <- 4 * xii[1]
    Delta[1, 2] <- 0
    Delta[2, 1] <- 0
    Delta[2, 2] <- 4 * xii[2]
    
    Tau_b ~ dwish(inverse(Delta), 3)
    SIG_b <- inverse(Tau_b[,])
  }
  ")
  return(model_str)
}

# ==============================================================================
# 2. 运行 Wiener 推断函数 (已修复数据和初始值)
# ==============================================================================
# Wiener_model.R 中的 run_wiener_inference 函数

run_wiener_inference <- function(data_list, type, inits_params, 
                                 n.chains = 2, n.adapt = 1000, n.iter = 5000) {
  
  model_string <- get_wiener_model_string(type)
  n <- data_list$n; p <- data_list$p
  
  # 数据清洗
  obs_data <- if(!is.null(data_list$y.diff)) data_list$y.diff else data_list$x_diff
  wiener_data <- list(y.diff = obs_data, n = n, p = p, m = data_list$m, t1 = data_list$t1,
                      m0 = rep(0, p), Tau0 = diag(0.001, p))
  
  # 初始值
  gen_inits <- function() {
    list(mu0 = rep(inits_params$mu0, p), sigma2_err = rep(inits_params$sigma2_err, p),
         Tau_b = inits_params$Tau_b_scale * diag(1, p),
         b = matrix(inits_params$b_val, nrow = n, ncol = p),
         scale = inits_params$scale)
  }
  inits_list <- replicate(n.chains, gen_inits(), simplify = FALSE)
  
  message(paste0(">>> Initializing Wiener Model: Type=", type))
  jags_model <- jags.model(textConnection(model_string), data = wiener_data, inits = inits_list, n.chains = n.chains, n.adapt = n.adapt)
  
  message(">>> Burning in (Update)...")
  update(jags_model, n.adapt) 
  
  message(">>> Sampling Posterior...")
  # [修复点] 添加 sigma2_err 到监控列表
  samples <- coda.samples(jags_model, 
                          variable.names = c("mu0", "sigma", "rho", "scale", "sigma2_err"), 
                          n.iter = n.iter)
  
  dic <- dic.samples(jags_model, n.iter = n.iter)
  return(list(samples = samples, dic = dic))
}