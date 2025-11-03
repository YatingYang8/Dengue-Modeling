# Run SEI-SEIR model with different rainfall functions ------------------------------
rm(list=ls()) #remove previous variable assignments

setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue")

# load packages
library(deSolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(mcmc)

# load data 
source("02_2_SEI-SEIR_model_THR.R")
source("02_1_SEI-SEIR_simulation_setup.R")


# 定义模型，包含待估计参数
seiseir_model_thr <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dM1 <- EFD(temp[t])*pEA(temp[t])*MDR(temp[t])*mu_th(temp[t], hum[t], timestep)^(-1)*(M1+M2+M3)*max((1-((M1+M2+M3)/K_thr(temp[t], rain[t], Rmax, (S+E+I+R)*2, timestep))),0)-(a(temp[t], c1)*pMI(temp[t])*(I/(S+E+I+R))+mu_th(temp[t], hum[t], timestep))*M1
    dM2 <- a(temp[t], c1)*pMI(temp[t])*(I/(S+E+I+R))*M1-(PDR(temp[t])+mu_th(temp[t], hum[t], timestep))*M2
    dM3 <- PDR(temp[t])*M2-mu_th(temp[t], hum[t], timestep)*M3
    dS <- -a(temp[t], c1)*b(temp[t])*(M3/(S+E+I+R))*S + BR*(S/1000)/365 - DR*(S/1000)/365 + ie*(S+E+I+R) - ie*S
    dE <- a(temp[t], c1)*b(temp[t])*(M3/(S+E+I+R))*S-(1.0/5.9)*E - DR*(E/1000)/365 - ie*E
    dI <- (1.0/5.9)*E-(1.0/5.0)*I - DR*(I/1000)/365 - ie*I
    dR <- (1.0/5.0)*I - DR*(R/1000)/365 - ie*R
    list(c(dM1, dM2, dM3, dS, dE, dI, dR))
  })
}




# Define 'a' function globally
a <- function(temp, c1) {
  briere(temp, c1, 13.35, 40.08)
}

# Updated loss function
lossfunction <- function(c1, temp, rain, hum, timestep, startIC, Rmax, observed_I) {
  parameters <- list(
    EFD = EFD, pEA = pEA, MDR = MDR, K_thr = K_thr, 
    pMI = pMI, mu_th = mu_th, PDR = PDR, b = b, 
    a = a, c1 = c1, timestep = timestep
  )
  
  state <- c(M1 = M0-1, M2 = 0, M3 = 1, S = N-5, E = 0, I = 5, R = 0)
  
  out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters, 
             method = "rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
  out2 <- as.data.frame(out)
  predicted_I <- (1.0/5.9)*out2$E - DR*(out2$I/1000)/365 
  
  loss <- sum((predicted_I - observed_I)^2)
  return(loss)
}



# 定义MCMC采样函数
mcmcSampling <- function(init.c1, n.iter) {
  c1.samples <- numeric(n.iter)
  c1.current <- init.c1
  for (i in 2:n.iter) {
    
    # 调整以保证 alpha 和 beta 的合理性
    #alpha <- 2  # 设定较小值
    #beta <- 10000
    
    alpha <- 2  # 设定较小值
    beta <- 10000
    
    # 使用贝塔分布
    new.c1 <- rbeta(1, alpha, beta)
    
    # 计算旧的和新的损失函数
    loss.old <- lossfunction(c1.current, temp, rain, hum, timestep, startIC, Rmax, observed_I)
    loss.new <- lossfunction(new.c1, temp, rain, hum, timestep, startIC, Rmax, observed_I)
    
    # Metropolis算法
    if (loss.new - loss.old < 0 ) {
      c1.current <- new.c1
    }
    
    c1.samples[i] <- c1.current
    
    cat("n.iter:", i,";","c1.current:", c1.samples[i], "\n")
  }
  return(c1.samples)
}


# 初始化参数和数据
set.seed(123)


i=1
climateData$SVPD <- 1
climateData2 <- climateData %>%
  filter(Site == sites[i] & date >= as.Date("2023-04-01"))
temp <- as.numeric(climateData2$Temperature)
rain <- as.numeric(climateData2$Rainfall)
hum <- as.numeric(climateData2$SVPD)
Rmax <- 123
Date <- climateData2$date
N <- population[i]
city <- sites[i]
BR <- BRs[i]
DR <- DRs[i]
times <- seq(1,length(Date), by=1)
K_thr <- K_thr_inverse
M0 <- K_thr(temp[1], mean(rain), Rmax, N, timestep)
parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep=timestep)
state <- c(M1 = M0-1, M2 = 0, M3 = 1, S = N-5, E = 0, I = 5, R = 0)

timestep <- 1/12


# observed data
# observed_data <- read.csv("./output/SEI-SEIR_simulations_THR_a.csv") %>%
#   filter(Rain_function == "Inverse")
# observed_I <- observed_data$I

observed_data <- read.csv("./data/dengue_caseALL.csv") %>%
  filter( date >= as.Date("2023-04-01"))
observed_I <- observed_data$cases 




# 运行MCMC采样
init.c1 <- 0
n.iter <- 500
c1.samples <- mcmcSampling(init.c1=init.c1, n.iter=n.iter)

# 计算点估计和置信区间
point.estimate <- mean(c1.samples)
ci.lower <- quantile(c1.samples, 0.025)
ci.upper <- quantile(c1.samples, 0.975)

# 输出结果
cat("Point Estimate of c1:", point.estimate, "\n")
cat("95% Credible Interval of c1:", ci.lower, "to", ci.upper, "\n")





# #附：最小二乘法代码，可以先用最小二乘法计算出参数的大概范围，然后再用mcmc获得更准确的估计和置信区间。
# rm(list=ls()) #remove previous variable assignments
# 
# setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue")
# 
# # load packages
# library(deSolve)
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# 
# 
# # load data 
# source("02_2_SEI-SEIR_model_THR.R")
# source("02_1_SEI-SEIR_simulation_setup.R")
# 
# # observed data
# observed_data <- read.csv("./output/SEI-SEIR_simulations_THR_a.csv") %>%
#   filter(Rain_function == "Inverse")
# observed_I <- observed_data$I
# 
# 
# #parameters need to estimate
# # biting rate
# a <- function(temp,c1){
#   briere(temp,c1,13.35,40.08)
# }
# 
# seiseir_model_thr <- function(t, state, parameters) {
#   with(as.list(c(state,parameters)), {
#     dM1 <- EFD(temp[t])*pEA(temp[t])*MDR(temp[t])*mu_th(temp[t], hum[t], timestep)^(-1)*(M1+M2+M3)*max((1-((M1+M2+M3)/K_thr(temp[t], rain[t], Rmax, (S+E+I+R)*2, timestep))),0)-(a(temp[t],c1)*pMI(temp[t])*(I/(S+E+I+R))+mu_th(temp[t], hum[t], timestep))*M1
#     dM2 <- a(temp[t],c1)*pMI(temp[t])*(I/(S+E+I+R))*M1-(PDR(temp[t])+mu_th(temp[t], hum[t], timestep))*M2
#     dM3 <- PDR(temp[t])*M2-mu_th(temp[t], hum[t], timestep)*M3
#     dS <- -a(temp[t],c1)*b(temp[t])*(M3/(S+E+I+R))*S + BR*(S/1000)/365 - DR*(S/1000)/365 + ie*(S+E+I+R) - ie*S
#     dE <- a(temp[t],c1)*b(temp[t])*(M3/(S+E+I+R))*S-(1.0/5.9)*E - DR*(E/1000)/365 - ie*E
#     dI <- (1.0/5.9)*E-(1.0/5.0)*I - DR*(I/1000)/365 - ie*I
#     dR <- (1.0/5.0)*I - DR*(R/1000)/365 - ie*R
#     list(c(dM1, dM2, dM3, dS, dE, dI, dR))
#   })
# } 
# 
# 
# 
# i=1
# climateData$SVPD <- 1
# climateData2 <- climateData %>% filter(Site == sites[i] & dengue_date >= as.Date("2022-09-01"))
# temp <- as.numeric(climateData2$Temperature)
# rain <- as.numeric(climateData2$Rainfall)
# hum <- as.numeric(climateData2$SVPD)
# Rmax <- 123
# Date <- climateData2$dengue_date
# N <- population[i]
# city <- sites[i]
# BR <- BRs[i]
# DR <- DRs[i]
# times <- seq(1,length(Date), by=1)
# K_thr <- K_thr_inverse
# M0 <- K_thr(temp[1], mean(rain), Rmax, N, timestep)
# parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep=timestep)
# state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*N, E = startIC$e*N, I = startIC$i*N, R = startIC$r*N)
# 
# fit_sir_model <- function(params.est) {
#   c1 <- params.est[1]
#   
#   
#   # biting rate
#   a <- function(temp,c1){
#     briere(temp,c1,13.35,40.08)
#   }
#   
#   parameters <- c(EFD, pEA, MDR, K_thr, a, c1 = c1, mu_th, PDR, b, timestep = timestep)
#   
#   times <- 1:length(observed_I)
#   fit <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
#   predicted_I <- fit[, "I"]
#   
#   # 计算损失函数
#   loss <- sum((predicted_I - observed_I)^2)
#   return(loss)
# }  
# 
# # 使用优化函数来估计参数
# optim_result <- optim(par = c(0), fn = fit_sir_model,
#                       method = "Nelder-Mead", lower = c(0), upper = c(0.01))
# 
# # 输出最优参数
# cat("Optimized c1:", optim_result$par[1], "\n") #2.32e-04
# 
# 
# 
# 
# 
# #Demo
# # # 定义SIR模型的微分方程
# # sir_equations <- function(time, state, parameters) {
# #   with(as.list(c(state, parameters)), {
# #     dS <- -beta * S * I / N
# #     dI <- beta * S * I / N - gamma * I
# #     dR <- gamma * I
# #     return(list(c(dS, dI, dR)))
# #   })
# # }
# # 
# # # 初始参数和条件
# # N <- 1e8 # 总人口
# # I0 <- 1   # 初始感染者人数
# # R0 <- 0   # 初始康复者人数
# # S0 <- N - I0 - R0 # 初始易感者人数
# # initial_state <- c(S = S0, I = I0, R = R0)
# # parameters <- c(beta = 0.3, gamma = 0.1) # 初始感染率和恢复率
# # 
# # # 时间序列
# # t <- seq(0, 100, by = 1)
# # 
# # # 模拟一些观测数据
# # set.seed(123)
# # observed_I <- rnorm(100, mean = 50, sd = 10) # 模拟100天的感染者数据
# # 
# # # 拟合模型
# # fit_sir_model <- function(params) {
# #   beta_initial <- params[1]
# #   gamma_initial <- params[2]
# #   times <- 1:length(observed_I)
# #   
# #   fit <- ode(y = initial_state, times = times, func = sir_equations,
# #              parms = c(beta = beta_initial, gamma = gamma_initial))
# #   predicted_I <- fit[, "I"]
# #   
# #   # 计算损失函数
# #   loss <- sum((predicted_I - observed_I)^2)
# #   return(loss)
# # }
# # 
# # # 使用优化函数来估计参数
# # optim_result <- optim(par = c(0.3, 0.1), fn = fit_sir_model,
# #                       method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))
# # 
# # # 输出最优参数
# # cat("Optimized beta:", optim_result$par[1], "\n")
# # cat("Optimized gamma:", optim_result$par[2], "\n")
# # 



