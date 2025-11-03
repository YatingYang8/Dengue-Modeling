
rm(list=ls()) #remove previous variable assignments

setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue")

# load packages
library(deSolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(mcmc)

# 用最小二乘法估计参数

# Demo
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

load("data/climateData.RData")
climateData2 <- climateData %>%
  filter( date >= as.Date("2023-04-01"))
Date <- climateData2$date

times <- seq(1,length(Date), by=1)

#state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*N, E = startIC$e*N, I = startIC$i*N, R = startIC$r*N)
state <- c(S = N-5, E = 0, I = 5, R = 0)

parameters <- c( beta0 = 0.2, beta1 = 0.6, phi = 0, T = 365,sigma = (1/5.9), gamma = (1/5.0), timestep=timestep)


out <- ode(y = state, times = times, func = seir_model, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)

out2 <- as.data.frame(out)

out2$dS <- c(NA, diff(out2$S))
out2$dE <- c(NA, diff(out2$E))
out2$dI <- c(NA, diff(out2$I))
out2$dR <- c(NA, diff(out2$R))

out2$dI_in <- sigma*out2$E - DR*(out2$I/1000)/365 

out2$Date <- Date

sum(out2$dI_in)

predicted_I <- sigma*out2$E - DR*(out2$I/1000)/365 
observed_data <- read.csv("./data/dengue_caseALL.csv") %>%
  filter( date >= as.Date("2023-04-01"))
observed_I <- observed_data$cases 
plotdata <- data.frame(date=Date, predicted_I=predicted_I, observed_I=observed_I)

ggplot(plotdata, aes(x = date, y = predicted_I, color = "predicted")) +
  geom_line() + # 绘制预测值的线
  geom_point(aes(y = observed_I, color = "observed")) + # 绘制观测值的点
  scale_color_manual(values = c("predicted" = "blue", "observed" = "red")) + # 为预测值和观测值设置颜色
  labs(x = "Date", y = "Value", color = "Legend", title = "Predicted vs Observed") + # 设置坐标轴标签和图例标题
  theme_minimal()



# 观测数据
observed_data <- read.csv("./data/dengue_caseALL.csv") %>%
  filter( date >= as.Date("2023-04-01"))
observed_I <- observed_data$cases 


# 拟合模型
fit_sir_model <- function(params) {
  
  beta0_initial = params[1]
  beta1_initial = params[2]
  phi_initial = params[3]
  
  times <- 1:length(observed_I)
  
  state <- c(S = N-5, E = 0, I = 5, R = 0)

  fit <- ode(y = state, times = times, func = seir_model,
             parms = c( beta0 = beta0_initial, beta1 = beta1_initial, phi = phi_initial, T = 365,sigma = (1/5.9), gamma = (1/5.0), timestep=timestep))
  predicted_I <- sigma*fit[, "E"]- DR*fit[, "I"]/1000/365
  

  # 计算损失函数
  loss <- sum((predicted_I - observed_I)^2)
  return(loss)
}

# 使用优化函数来估计参数
optim_result <- optim(par = c(0.25, 0.55, 0), fn = fit_sir_model,
                      method = "L-BFGS-B", lower = c(0.2, 0.5,-3.14), upper = c(0.3, 1,3.14))

# 输出最优参数
cat("Optimized beta0:", optim_result$par[1], "\n")
cat("Optimized beta1:", optim_result$par[2], "\n")
cat("Optimized phi:", optim_result$par[3], "\n")


parameters <- c( beta0 = 0.2067163, beta1 = 0.5327438, phi = -0.9244718, T = 365,sigma = (1/5.9), gamma = (1/5.0), timestep=timestep)


out <- ode(y = state, times = times, func = seir_model, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)

out2 <- as.data.frame(out)

out2$dS <- c(NA, diff(out2$S))
out2$dE <- c(NA, diff(out2$E))
out2$dI <- c(NA, diff(out2$I))
out2$dR <- c(NA, diff(out2$R))

out2$dI_in <- sigma*out2$E - DR*(out2$I/1000)/365 

out2$Date <- Date

sum(out2$dI_in)

predicted_I <- sigma*out2$E - DR*(out2$I/1000)/365 
observed_data <- read.csv("./data/dengue_caseALL.csv") %>%
  filter( date >= as.Date("2023-04-01"))
observed_I <- observed_data$cases 
plotdata <- data.frame(date=Date, predicted_I=predicted_I, observed_I=observed_I)

ggplot(plotdata, aes(x = date, y = predicted_I, color = "predicted")) +
  geom_line() + # 绘制预测值的线
  geom_point(aes(y = observed_I, color = "observed")) + # 绘制观测值的点
  scale_color_manual(values = c("predicted" = "blue", "observed" = "red")) + # 为预测值和观测值设置颜色
  labs(x = "Date", y = "Value", color = "Legend", title = "Predicted vs Observed") + # 设置坐标轴标签和图例标题
  theme_minimal()



###### plot the beta_t ####
# 定义参数
# 定义时间序列
times <- 1:length(observed_I)
# 计算 beta_t
parameters <- c( beta0 = 0.2067163, beta1 = 0.5327438, phi = -0.9244718, T = 365)
beta_t <- parameters["beta0"] * (1 + parameters["beta1"] * cos(2 * pi * times  / parameters["T"] + parameters["phi"]))

# 创建数据框用于绘图
beta_data <- data.frame(
  Time = times,
  Beta_t = beta_t,
  I = out2$I
) %>%
  mutate( y=Beta_t*I)

# 绘制 beta_t 随时间变化的图像
ggplot(beta_data, aes(x = Time, y = Beta_t)) +
  geom_line(color = "blue", size = 1) +
  labs(
    title = expression(paste(beta[t], " 随时间变化")),
    x = "时间 (天)",
    y = expression(beta[t])
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )


plot_data <- beta_data %>%
  pivot_longer(cols = c(Beta_t, I, y), names_to = "Variable", values_to = "Value")

# 绘图
ggplot(plot_data, aes(x = Time, y = Value, color = Variable, linetype = Variable)) +
  geom_line(size = 1) +
  labs(
    title = "Beta_t, I 和 y 随时间变化",
    x = "时间 (天)",
    y = "值",
    color = "变量",
    linetype = "变量"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# 绘制 Beta_t 随时间变化的图
plot_beta_t <- ggplot(beta_data, aes(x = Time, y = Beta_t)) +
  geom_line(color = "blue", size = 1) +
  labs(
    title = expression(paste(beta[t], " 随时间变化")),
    x = "时间 (天)",
    y = expression(beta[t])
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# 绘制 I 随时间变化的图
plot_I <- ggplot(beta_data, aes(x = Time, y = I)) +
  geom_line(color = "red", size = 1) +
  labs(
    title = "I 随时间变化",
    x = "时间 (天)",
    y = "I"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# 绘制 y 随时间变化的图
plot_y <- ggplot(beta_data, aes(x = Time, y = y)) +
  geom_line(color = "green", size = 1) +
  labs(
    title = expression(paste(beta[t], " *I随时间变化")),
    x = "时间 (天)",
    y = expression(paste(beta[t], " *I"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# 输出各图
print(plot_beta_t)
print(plot_I)
print(plot_y)




#################################################################################
########################### 根据刚刚的参数范围调整SEIR-SEI模型结果 ##############
#################################################################################

# load data 
source("02_2_SEI-SEIR_model_THR.R")
source("02_1_SEI-SEIR_simulation_setup.R")

# run simulations
traitDF <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R","dS", "dE", "dI", "dR", "dI_in", "Date", "Site", "Rain_function")
traitFileName <- "output/SEI-SEIR_simulations_THR.csv"
write.csv(traitDF, traitFileName, row.names = F)

rfunctions_names <- c("Briere", "Quadratic", "Inverse")
rfunctions <- list(K_thr_briere, K_thr_quadratic, K_thr_inverse)

climateData$SVPD <- 1

for (i in 1:length(sites)){
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
  for (k in 1:length(rfunctions_names)){
    K_thr <- rfunctions[[k]]
    M0 <- K_thr(temp[1], mean(rain), Rmax, N, timestep)
    parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep=timestep)
    
    #state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*N, E = startIC$e*N, I = startIC$i*N, R = startIC$r*N)
    state <- c(M1 = M0-1, M2 = 0, M3 = 1, S = N-5, E = 0, I = 5, R = 0)
    
    
    out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
    out2 <- as.data.frame(out)
    
    out2$dS <- c(NA, diff(out2$S))
    out2$dE <- c(NA, diff(out2$E))
    out2$dI <- c(NA, diff(out2$I))
    out2$dR <- c(NA, diff(out2$R))
    
    out2$dI_in <- (1.0/5.9)*out2$E - DR*(out2$I/1000)/365 
    
    out2$Date <- Date
    out2$Site <- sites[i]
    out2$Rain_function <- rfunctions_names[k]
    traitDF <- rbind(traitDF, out2)
    write.csv(traitDF, traitFileName, row.names = F)
    cat("finished running ode for", sites[i], rfunctions_names[k], "\n")
  }
}






# Load your simulation results
traitDF <- read.csv(traitFileName) %>% filter(Rain_function == "Inverse")

# Convert Date to Date type
traitDF$Date <- as.Date(traitDF$Date)

a1 <- unlist(lapply(temp, a))
b1 <- unlist(lapply(temp,b))

EFD1 <- unlist(lapply(temp, EFD))
pEA1 <- unlist(lapply(temp, pEA))
MDR1 <- unlist(lapply(temp, MDR))
mu_th1 <- mapply(mu_th, temp, hum, MoreArgs = list(timestep = 1/12))
pMI1 <- unlist(lapply(temp, pMI))
K_thr1 <- mapply(K_thr,temp, rain, MoreArgs = list(Rmax=123, N=1333000*2, timestep=1/12))

# dM1 <- EFD(temp[t])*pEA(temp[t])*MDR(temp[t])*mu_th(temp[t], hum[t], timestep)^(-1)*(M1+M2+M3)*
#   max((1-((M1+M2+M3)/K_thr(temp[t], rain[t], Rmax, (S+E+I+R)*2, timestep))),0)-(a(temp[t])*pMI(temp[t])*
#                                                                                   (I/(S+E+I+R))+mu_th(temp[t], hum[t], timestep))*M1


traitDF$temp <- temp
traitDF$a <- a1
traitDF$b <- b1
traitDF <- traitDF %>%
  mutate(infection_rate = a1 * b1 * M3)

traitDF$EFD <- EFD1
traitDF$pEA <- pEA1
traitDF$MDR <- MDR1
traitDF$mu_th <- mu_th1
traitDF$pMI <- pMI1
traitDF$K_thr <- K_thr1

traitDF <- traitDF %>%
  mutate(yy =  ifelse((1-((M1+M2+M3)/K_thr))>0, (1-((M1+M2+M3)/K_thr)),0) )

head(traitDF)

