
rm(list=ls()) #remove previous variable assignments

setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue")

# load packages
library(deSolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(mcmc)
library(truncnorm)


# 定义SIR模型的微分方程
seir_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta_t <- beta0 * (1 + beta1 * cos(2 * pi * t / T + phi)) 
    dS <- -beta_t * S * I / N + BR * (S / 1000) / 365 - DR * (S / 1000) / 365 + ie * N - ie * S
    dE <- beta_t * S * I / N - sigma * E - DR * (E / 1000) / 365 - ie * E
    dI <- sigma * E - gamma * I - DR * (I / 1000) / 365 - ie * I
    dR <- gamma * I - DR * (R / 1000) / 365 - ie * R
    list(c(dS, dE, dI, dR))
  })
}


# 初始参数和条件
N <- 1333000
sigma = (1/5.9)
gamma = (1/5.0)
ie <- 0 # set immigration and emmigration rate
BR <- 6.89 # birth rates per 1000 persons
DR <- 7.27 # death rates per 1000 persons
timestep = 1/12 # model timestep

observed_data <- read.csv("./data/dengue_caseALL.csv") %>%
  filter( date >= as.Date("2023-01-01"))
sum(observed_data$cases)

Date <- observed_data$date

times <- seq(1,length(Date), by=1)

#state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*N, E = startIC$e*N, I = startIC$i*N, R = startIC$r*N)
state <- c(S = N-5, E = 0, I = 5, R = 0)

parameters <- c( beta0 = 0.2096247, beta1 = 0.5534963, phi = 3.835708, T = 365,sigma = (1/5.9), gamma = (1/5.0), timestep=timestep)

out <- ode(y = state, times = times, func = seir_model, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)

out2 <- as.data.frame(out)

out2$dS <- c(NA, diff(out2$S))
out2$dE <- c(NA, diff(out2$E))
out2$dI <- c(NA, diff(out2$I))
out2$dR <- c(NA, diff(out2$R))

out2$dI_in <- sigma*out2$E - DR*(out2$I/1000)/365 

out2$Date <- Date

sum(out2$dI_in)

#write.csv(out2, "output/SEIR_simulations_result.csv", row.names = F)

sum(predicted_I)
predicted_I <- sigma*out2$E - DR*(out2$I/1000)/365 
observed_data <- read.csv("./data/dengue_caseALL.csv") %>%
  filter( date >= as.Date("2023-01-01"))
observed_I <- observed_data$cases 

mean(observed_I)+sd(observed_I)

mean(predicted_I)+sd(predicted_I)
plotdata <- data.frame(date=Date, predicted_I=predicted_I, observed_I=observed_I)
plotdata$date <- as.Date(plotdata$date)


library(ggplot2)
ggplot(plotdata, aes(x = date, y = predicted_I, color = "predicted")) +
  geom_line() + # 绘制预测值的线
  geom_point(aes(y = observed_I, color = "observed")) + # 绘制观测值的点
  scale_color_manual(values = c("predicted" = "blue", "observed" = "red")) + # 为预测值和观测值设置颜色
  labs(x = "Date", y = "Value", color = "Legend", title = "Predicted vs Observed") + # 设置坐标轴标签和图例标题
  theme_minimal()


#美化图片
library(ggplot2)
library(extrafont)  # 用于加载正式字体
library(scales)     # 用于日期格式化

# 美化后的代码
ggplot(plotdata, aes(x = date, y = predicted_I, color = "Predicted")) +
  geom_line(size = 1.2) +  # 加粗预测值的线
  geom_point(aes(y = observed_I, color = "Observed"), size = 3, alpha = 0.8) +  # 调整观测值点的大小和透明度
  scale_color_manual(values = c("Predicted" = "#1f77b4", "Observed" = "#d62728")) +  # 使用更柔和的颜色
  labs(
    x = "Date", 
    y = "Value", 
    color = "Legend", 
    title = "Predicted vs Observed Dengue Cases Over Time"
  ) +
  theme_minimal(base_family = "Times New Roman") +  # 使用正式字体
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # 标题居中并加粗
    axis.title.x = element_text(size = 12, face = "bold"),  # X轴标题加粗
    axis.title.y = element_text(size = 12, face = "bold"),  # Y轴标题加粗
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # 调整X轴标签角度
    axis.text.y = element_text(size = 10),  # 调整Y轴标签大小
    legend.position = "bottom",  # 图例放在底部
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题加粗
    legend.text = element_text(size = 10),  # 图例文字大小
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),  # 网格线样式
    panel.grid.minor = element_blank()  # 去掉次要网格线
  ) +
  scale_x_date(labels = date_format("%Y-%m"), breaks = date_breaks("months"))  # 格式化日期轴



head(plotdata)


library(ggplot2)
library(ggpubr)  # 用于添加统计信息
library(extrafont)  # 使用正式字体

# 计算相关系数和 p 值
cor_test_result <- cor.test(plotdata$observed_I, plotdata$predicted_I, method = "pearson")
cor_text <- paste0("r = ", round(cor_test_result$estimate, 3), 
                   ", p = ", format.pval(cor_test_result$p.value, digits = 2))

# 绘制散点图
ggplot(plotdata, aes(x = observed_I, y = predicted_I)) +
  geom_point(size = 3, alpha = 0.8, shape = 21, fill = "#4C72B0", color = "#2A4E6C") +  # 带外框的数据点
  geom_abline(slope = 1, intercept = 0, color = "gray80", linetype = "dashed", size = 1) +  # 浅灰色 y = x 虚线
  annotate("text", x = min(plotdata$observed_I), y = 100, 
           label = cor_text, hjust = 0, vjust = 1, size = 4, color = "black", family = "Times New Roman") +  # 左上角显示统计信息
  labs(
    x = "Observed Values", 
    y = "Predicted Values", 
    title = "Observed vs Predicted Values"
  ) +
  theme_minimal(base_family = "Times New Roman") +  # 使用正式字体
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # 标题居中并加粗
    axis.title.x = element_text(size = 12, face = "bold"),  # X轴标题加粗
    axis.title.y = element_text(size = 12, face = "bold"),  # Y轴标题加粗
    axis.text.x = element_text(size = 10),  # X轴标签大小
    axis.text.y = element_text(size = 10),  # Y轴标签大小
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # 网格线样式
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加外边框
    legend.position = "none"  # 去掉图例
  ) +
  scale_x_continuous(limits = c(0, max(plotdata$observed_I, plotdata$predicted_I))) +  # 设置 X 轴范围
  scale_y_continuous(limits = c(0, max(plotdata$observed_I, plotdata$predicted_I)))  # 设置 Y 轴范围


# #################################################################################
# ######################### 最小二乘法和MCMC估计参数 ##############################
# #################################################################################
# 
# rm(list=ls()) #remove previous variable assignments
# 
# setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue")
# 
# # load packages
# library(deSolve)
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(mcmc)
# 
# ############# 用最小二乘法估计参数 ################################################
# 
# # Demo
# # 定义SIR模型的微分方程
# seir_model <- function(t, state, parameters) {
#   with(as.list(c(state, parameters)), {
#     beta_t <- beta0 * (1 + beta1 * cos(2 * pi * t / T + phi)) 
#     dS <- -beta_t * S * I / N + BR * (S / 1000) / 365 - DR * (S / 1000) / 365 + ie * N - ie * S
#     dE <- beta_t * S * I / N - sigma * E - DR * (E / 1000) / 365 - ie * E
#     dI <- sigma * E - gamma * I - DR * (I / 1000) / 365 - ie * I
#     dR <- gamma * I - DR * (R / 1000) / 365 - ie * R
#     list(c(dS, dE, dI, dR))
#   })
# }
# 
# 
# # 初始参数和条件
# N <- 1333000
# sigma = (1/5.9)
# gamma = (1/5.0)
# ie <- 0 # set immigration and emmigration rate
# BR <- 6.89 # birth rates per 1000 persons
# DR <- 7.27 # death rates per 1000 persons
# timestep = 1/12 # model timestep
# 
# load("data/climateData.RData")
# climateData2 <- climateData %>%
#   filter( date >= as.Date("2023-01-01"))
# Date <- climateData2$date
# 
# times <- seq(1,length(Date), by=1)
# 
# #state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*N, E = startIC$e*N, I = startIC$i*N, R = startIC$r*N)
# state <- c(S = N-5, E = 0, I = 5, R = 0)
# 
# parameters <- c( beta0 = 0.2, beta1 = 0.78, phi = -2.5, T = 365,sigma = (1/5.9), gamma = (1/5.0), timestep=timestep)
# 
# 
# out <- ode(y = state, times = times, func = seir_model, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
# 
# out2 <- as.data.frame(out)
# 
# out2$dS <- c(NA, diff(out2$S))
# out2$dE <- c(NA, diff(out2$E))
# out2$dI <- c(NA, diff(out2$I))
# out2$dR <- c(NA, diff(out2$R))
# 
# out2$dI_in <- sigma*out2$E - DR*(out2$I/1000)/365 
# 
# out2$Date <- Date
# 
# sum(out2$dI_in)
# 
# predicted_I <- sigma*out2$E - DR*(out2$I/1000)/365 
# observed_data <- read.csv("./data/dengue_caseALL.csv") %>%
#   filter( date >= as.Date("2023-01-01"))
# observed_I <- observed_data$cases 
# plotdata <- data.frame(date=Date, predicted_I=predicted_I, observed_I=observed_I)
# 
# ggplot(plotdata, aes(x = date, y = predicted_I, color = "predicted")) +
#   geom_line() + # 绘制预测值的线
#   geom_point(aes(y = observed_I, color = "observed")) + # 绘制观测值的点
#   scale_color_manual(values = c("predicted" = "blue", "observed" = "red")) + # 为预测值和观测值设置颜色
#   labs(x = "Date", y = "Value", color = "Legend", title = "Predicted vs Observed") + # 设置坐标轴标签和图例标题
#   theme_minimal()
# 
# 
# 
# # 观测数据
# observed_data <- read.csv("./data/dengue_caseALL.csv") %>%
#   filter( date >= as.Date("2023-01-01"))
# observed_I <- observed_data$cases 
# 
# 
# # 拟合模型
# fit_sir_model <- function(params) {
#   
#   beta0_initial = params[1]
#   beta1_initial = params[2]
#   phi_initial = params[3]
#   
#   times <- 1:length(observed_I)
#   
#   state <- c(S = N-5, E = 0, I = 5, R = 0)
#   
#   fit <- ode(y = state, times = times, func = seir_model,
#              parms = c( beta0 = beta0_initial, beta1 = beta1_initial, phi = phi_initial, T = 365,sigma = (1/5.9), gamma = (1/5.0), timestep=timestep))
#   predicted_I <- sigma*fit[, "E"]- DR*fit[, "I"]/1000/365
#   
#   
#   # 计算损失函数
#   loss <- sum((predicted_I - observed_I)^2)
#   return(loss)
# }
# 
# # 使用优化函数来估计参数
# optim_result <- optim(par = c(0.25, 0.78, 3.78), fn = fit_sir_model,
#                       method = "L-BFGS-B", lower = c(0.2, 0.5,0), upper = c(0.3, 1,6.28))
# 
# # 输出最优参数
# cat("Optimized beta0:", optim_result$par[1], "\n")
# cat("Optimized beta1:", optim_result$par[2], "\n")
# cat("Optimized phi:", optim_result$par[3], "\n")
# 
# 
# parameters <- c( beta0 = 0.2010281, beta1 = 0.7739768, phi = -2.503413, T = 365,sigma = (1/5.9), gamma = (1/5.0), timestep=timestep)
# 
# 
# out <- ode(y = state, times = times, func = seir_model, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
# 
# out2 <- as.data.frame(out)
# 
# out2$dS <- c(NA, diff(out2$S))
# out2$dE <- c(NA, diff(out2$E))
# out2$dI <- c(NA, diff(out2$I))
# out2$dR <- c(NA, diff(out2$R))
# 
# out2$dI_in <- sigma*out2$E - DR*(out2$I/1000)/365 
# 
# out2$Date <- Date
# 
# sum(out2$dI_in)
# 
# predicted_I <- sigma*out2$E - DR*(out2$I/1000)/365 
# observed_data <- read.csv("./data/dengue_caseALL.csv") %>%
#   filter( date >= as.Date("2023-01-01"))
# observed_I <- observed_data$cases 
# plotdata <- data.frame(date=Date, predicted_I=predicted_I, observed_I=observed_I)
# 
# ggplot(plotdata, aes(x = date, y = predicted_I, color = "predicted")) +
#   geom_line() + # 绘制预测值的线
#   geom_point(aes(y = observed_I, color = "observed")) + # 绘制观测值的点
#   scale_color_manual(values = c("predicted" = "blue", "observed" = "red")) + # 为预测值和观测值设置颜色
#   labs(x = "Date", y = "Value", color = "Legend", title = "Predicted vs Observed") + # 设置坐标轴标签和图例标题
#   theme_minimal()
# 
# 
# 
# ######################### 用MCMC估计参数 ########################################
# # Run SEI-SEIR model with different rainfall functions ------------------------------
# rm(list=ls()) #remove previous variable assignments
# 
# setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue")
# 
# # load packages
# library(deSolve)
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(mcmc)
# 
# 
# 
# # 观测数据
# observed_data <- read.csv("./data/dengue_caseALL.csv") %>%
#   filter( date >= as.Date("2023-01-01"))
# observed_I <- observed_data$cases 
# 
# # 定义模型，包含待估计参数
# seir_model <- function(t, state, parameters) {
#   with(as.list(c(state, parameters)), {
#     beta_t <- beta0 * (1 + beta1 * cos(2 * pi * t / T + phi)) 
#     dS <- -beta_t * S * I / N + BR * (S / 1000) / 365 - DR * (S / 1000) / 365 + ie * N - ie * S
#     dE <- beta_t * S * I / N - sigma * E - DR * (E / 1000) / 365 - ie * E
#     dI <- sigma * E - gamma * I - DR * (I / 1000) / 365 - ie * I
#     dR <- gamma * I - DR * (R / 1000) / 365 - ie * R
#     list(c(dS, dE, dI, dR))
#   })
# }  
# 
# # Updated loss function
# lossfunction <- function(beta0,beta1,phi,timestep,observed_I) {
#   parameters <- list(
#     beta0 = beta0, beta1 = beta1, phi = phi, T = 365,sigma = (1/5.9), gamma = (1/5.0), timestep=timestep
#   )
#   
#   state <- c(S = N-5, E = 0, I = 5, R = 0)
#   
#   out <- ode(y = state, times = times, func = seir_model, parms = parameters, 
#              method = "rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
#   
#   out2 <- as.data.frame(out)
#   predicted_I <- (1.0/5.9)*out2$E - DR*(out2$I/1000)/365 
#   
#   
#   loss <- sum((predicted_I - observed_I)^2)
#   return(loss)
# }
# 
# 
# 
# # 定义MCMC采样函数
# mcmcSampling <- function(init.beta0, init.beta1, init.phi, n.iter) {
#   beta0.samples <- numeric(n.iter)
#   beta1.samples <- numeric(n.iter)
#   phi.samples <- numeric(n.iter)
#   
#   #初始化参数
#   beta0.current <- init.beta0
#   beta1.current <- init.beta1
#   phi.current <- init.phi
#   
#   for (i in 2:n.iter) {
#     
#     new.beta0 <- rtruncnorm(1, a = 0, b = Inf, mean = 0.2, sd = 0.2) # 生成截断正态分布
#     new.beta1 <- rbeta(1, 2, 2)
#     new.phi <- runif(1,0,2*pi)
#     
#     # 计算旧的和新的损失函数
#     loss.old <- lossfunction(beta0.current,beta1.current,phi.current,timestep,observed_I)
#     loss.new <- lossfunction(new.beta0,new.beta1,new.phi,timestep,observed_I)
#     
#     # Metropolis算法
#     if (loss.new - loss.old < 0 ) {
#       beta0.current <- new.beta0
#       beta1.current <- new.beta1
#       phi.current <- new.phi
#     }
#     
#     
#     beta0.samples[i] <- beta0.current
#     beta1.samples[i] <- beta1.current
#     phi.samples[i] <- phi.current
#     
#     cat("n.iter:", i,";","beta0.current:", beta0.samples[i],"beta1.current:", beta1.samples[i],
#         "phi.current:", phi.samples[i],"\n")
#   }
#   return(list(beta0.samples = beta0.samples, beta1.samples = beta1.samples,
#               phi.samples = phi.samples))
#   
# }
# 
# 
# # 初始化参数和数据
# set.seed(1)
# 
# 
# 
# # 初始参数和条件
# N <- 1333000
# sigma = (1/5.9)
# gamma = (1/5.0)
# ie <- 0 # set immigration and emmigration rate
# BR <- 6.89 # birth rates per 1000 persons
# DR <- 7.27 # death rates per 1000 persons
# timestep = 1/12 # model timestep
# 
# Date <- observed_data$date
# 
# times <- seq(1,length(Date), by=1)
# 
# #state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*N, E = startIC$e*N, I = startIC$i*N, R = startIC$r*N)
# state <- c(S = N-5, E = 0, I = 5, R = 0)
# 
# #parameters <- c( beta0 = 0.2010281, beta1 = 0.7739768, phi = 3.78, T = 365,sigma = (1/5.9), gamma = (1/5.0), timestep=timestep)
# 
# #out <- ode(y = state, times = times, func = seir_model, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
# 
# # 运行MCMC采样
# init.beta0 <- 0.2010281
# init.beta1 <- 0.7739768
# init.phi <- 3.78
# 
# init.beta0 <- 0.1
# init.beta1 <- 0.1
# init.phi <- 0
# 
# init.beta0 <- 0.2096247
# init.beta1 <- 0.5534963
# init.phi <- 3.835708
# 
# n.iter <- 1000
# #result: n.iter: 1000 ; beta0.current: 0.2037218 beta1.current: 0.7426666 phi.current: 3.499288
# c1.samples <- mcmcSampling(init.beta0=init.beta0, init.beta1=init.beta1, 
#                            init.phi= init.phi, n.iter=n.iter)
# 
# # 计算点估计和置信区间
# point.estimate <- mean(c1.samples)
# ci.lower <- quantile(c1.samples, 0.025)
# ci.upper <- quantile(c1.samples, 0.975)
# 
# # 输出结果
# cat("Point Estimate of c1:", point.estimate, "\n")
# cat("95% Credible Interval of c1:", ci.lower, "to", ci.upper, "\n")



