
rm(list=ls())
# data to set up model simulations -------------------------------------------------------------
# load packages
library(deSolve)
library(dplyr)
library(reshape2)
library(lubridate)
library(tidyr)
library(tidyverse)


setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue")
# load climate data
load("data/city_climatedata.RData")

# load simulated imported cases
load("data/citydata_ALL.RData")

load("data/cityrisk_simdata.RData")
#cityrisk_simdata <- left_join(city_climatedata, citydata_ALL[,c("city","cityid","pop_total","total_cases")])

#cityrisk_simdata <- na.omit(cityrisk_simdata)

# load and set initial conditions
load("data/init.cond.RData")
startIC <- subset(init.cond, IC == "54")
#startIC <- subset(init.cond, IC == "52")
#startIC <- subset(init.cond, IC == "18")

# load traits
load("data/Random_sample_of_posterior_traits.RData")

# set immigration and emmigration rate
#ie <- 0.01
ie <- 0

# set up list of sites
sites <- unique(cityrisk_simdata$Site)

# set human population numbers for each site 
population <- unique(cityrisk_simdata$pop_total)

# set birth and death rates
BRs <- c(0) # birth rates per 1000 persons
DRs <- c(0) # death rates per 1000 persons

# model timestep
timestep = 1/12





############################ Simulation ###################################################
# Run SEI-SEIR model ------------------------------

# load packages
library(deSolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# load data 
source("02_2_SEI-SEIR_model_THR (from literature).R")

initial_cases <- c(1,5,10,15,20)
initial_cases <- c(1:20)
initial_cases <- c(40,60,80,100,200)
initial_cases <- c(30,50,70,90,150)

for (j in 1:length(initial_cases)) {

  # 创建保存模拟结果的数据框
  traitDF <- data.frame(matrix(ncol = 16, nrow = 0))
  colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R",
                         "dS", "dE", "dI", "dR", "dI_in", "Date", "Site")
  
  # 自动命名输出文件
  traitFileName <- paste0("output/Subsequent_city_risks_simulations_", initial_cases[j], ".csv")
  write.csv(traitDF, traitFileName, row.names = FALSE)
  
  # 设定SVPD
  city_climatedata$SVPD <- 1
  
  # 循环每一个城市进行模拟
  for (i in 1:length(sites)) {
    climateData2 <- city_climatedata %>%
      dplyr::filter(Site == sites[i] & date >= as.Date("2023-04-01"))
    
    temp <- as.numeric(climateData2$Temperature)
    rain <- as.numeric(climateData2$Rainfall)
    hum <- as.numeric(climateData2$SVPD)
    Rmax <- 123
    Date <- climateData2$date
    N <- population[i]
    city <- sites[i]
    BR <- BRs
    DR <- DRs
    times <- seq(1, length(Date), by = 1)
    K_thr <- K_thr_inverse
    M0 <- K_thr(temp[1], mean(rain), Rmax, N, timestep)
    parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep = timestep)
    
    # 初始状态设定
    state <- c(M1 = M0 - 1,
               M2 = 0,
               M3 = 1,
               S = N - initial_cases[j],
               E = 0,
               I = initial_cases[j],
               R = 0)
    
    # 求解ODE模型
    out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters,
               method = "rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
    
    out2 <- as.data.frame(out)
    out2$dS <- c(NA, diff(out2$S))
    out2$dE <- c(NA, diff(out2$E))
    out2$dI <- c(NA, diff(out2$I))
    out2$dR <- c(NA, diff(out2$R))
    
    out2$dI_in <- (1.0 / 5.9) * out2$E - DR * (out2$I / 1000) / 365
    
    
    a1 <- unlist(lapply(temp, a))
    b1 <- unlist(lapply(temp,b))
    out2 <- out2 %>%
      mutate(temp=temp,
             a1=a1,
             b1=b1) %>%
      mutate(infection_rate = a1 * b1 * M3)
    
    out2$Date <- Date
    out2$Site <- sites[i]
    out2$initial_cases <- initial_cases[j]
    
    # 添加到总数据框中
    traitDF <- rbind(traitDF, out2)
    
    # 写入CSV
    write.csv(traitDF, traitFileName, row.names = FALSE, fileEncoding = "GB18030")
    
    cat("Finished running ODE for", sites[i], "with", initial_cases[j], "initial cases\n")
  }
}






############## Subsequent risks according to Xsbn simulation results ###############

#load simulation results
load("data/cityrisk_simdata.RData")

source("06_02_SEI-SEIR_model_THR with specific imported cases.R")

# 创建保存模拟结果的数据框
traitDF <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R",
                       "dS", "dE", "dI", "dR", "dI_in", "Date", "Site")

# 自动命名输出文件
traitFileName <- paste0("output/Subsequent_city_risks_simulations_", "xsbn", ".csv")
write.csv(traitDF, traitFileName, row.names = FALSE)

# 设定SVPD
city_climatedata$SVPD <- 1
#i=1

# 循环每一个城市进行模拟
for (i in 1:length(sites)) {
  
  cityrisk_simdata2 <- cityrisk_simdata %>%
    dplyr::filter(city == sites[i] & date >= as.Date("2023-04-01"))
  
  climateData2 <- city_climatedata %>%
    dplyr::filter(Site == sites[i] & date >= as.Date("2023-04-01"))
  
  imported_cases <- as.numeric(cityrisk_simdata2$total_cases)
  temp <- as.numeric(climateData2$Temperature)
  rain <- as.numeric(climateData2$Rainfall)
  hum <- as.numeric(climateData2$SVPD)
  Rmax <- 123
  Date <- climateData2$date
  N <- population[i]
  city <- sites[i]
  BR <- BRs
  DR <- DRs
  times <- seq(1, length(Date), by = 1)
  K_thr <- K_thr_inverse
  M0 <- K_thr(temp[1], mean(rain), Rmax, N, timestep)
  parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep = timestep)
  
  # 初始状态设定
  state <- c(M1 = M0 - 1,
             M2 = 0,
             M3 = 1,
             S = N - 1,
             E = 0,
             I = 1,
             R = 0)
  
  # 求解ODE模型
  out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters,
             method = "rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
  
  out2 <- as.data.frame(out)
  out2$dS <- c(NA, diff(out2$S))
  out2$dE <- c(NA, diff(out2$E))
  out2$dI <- c(NA, diff(out2$I))
  out2$dR <- c(NA, diff(out2$R))
  
  out2$dI_in <- (1.0 / 5.9) * out2$E - DR * (out2$I / 1000) / 365
  out2$Date <- Date
  out2$Site <- sites[i]
  out2$initial_cases <- "xsbn"
  
  # 添加到总数据框中
  traitDF <- rbind(traitDF, out2)
  
  # 写入CSV
  write.csv(traitDF, traitFileName, row.names = FALSE, fileEncoding = "GB18030")
  
  cat("Finished running ODE for", sites[i], "with", "xsbn", "cases\n")
}



#################### Subsequent risks by different first arrival time ################

initial_time <- c("2023-04-01","2023-05-01","2023-06-01","2023-07-01")
initial_time <- c("2023-03-01","2023-08-01","2023-09-01","2023-10-01","2023-11-01")

for (j in 1:length(initial_time)) {
  
  # 创建保存模拟结果的数据框
  traitDF <- data.frame(matrix(ncol = 16, nrow = 0))
  colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R",
                         "dS", "dE", "dI", "dR", "dI_in", "Date", "Site")
  
  # 自动命名输出文件
  traitFileName <- paste0("output/Subsequent_city_risks_simulations_", initial_time[j], ".csv")
  write.csv(traitDF, traitFileName, row.names = FALSE)
  
  # 设定SVPD
  city_climatedata$SVPD <- 1
  
  # 循环每一个城市进行模拟
  for (i in 1:length(sites)) {
    climateData2 <- city_climatedata %>%
      dplyr::filter(Site == sites[i] & date >= as.Date(initial_time[j]))
    
    temp <- as.numeric(climateData2$Temperature)
    rain <- as.numeric(climateData2$Rainfall)
    hum <- as.numeric(climateData2$SVPD)
    Rmax <- 123
    Date <- climateData2$date
    N <- population[i]
    city <- sites[i]
    BR <- BRs
    DR <- DRs
    times <- seq(1, length(Date), by = 1)
    K_thr <- K_thr_inverse
    M0 <- K_thr(temp[1], mean(rain), Rmax, N, timestep)
    parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep = timestep)
    
    # 初始状态设定
    state <- c(M1 = M0 - 1,
               M2 = 0,
               M3 = 1,
               S = N - 1,
               E = 0,
               I = 1,
               R = 0)
    
    # 求解ODE模型
    out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters,
               method = "rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
    
    out2 <- as.data.frame(out)
    out2$dS <- c(NA, diff(out2$S))
    out2$dE <- c(NA, diff(out2$E))
    out2$dI <- c(NA, diff(out2$I))
    out2$dR <- c(NA, diff(out2$R))
    
    out2$dI_in <- (1.0 / 5.9) * out2$E - DR * (out2$I / 1000) / 365
    
    
    a1 <- unlist(lapply(temp, a))
    b1 <- unlist(lapply(temp,b))
    out2 <- out2 %>%
      mutate(temp=temp,
             a1=a1,
             b1=b1) %>%
      mutate(infection_rate = a1 * b1 * M3)
    
    out2$Date <- Date
    out2$Site <- sites[i]
    out2$initial_time <- initial_time[j]
    
    # 添加到总数据框中
    traitDF <- rbind(traitDF, out2)
    
    # 写入CSV
    write.csv(traitDF, traitFileName, row.names = FALSE, fileEncoding = "GB18030")
    
    cat("Finished running ODE for", sites[i], "with", initial_time[j], "initial time\n")
  }
}







############## Subsequent risks according to travel frequency ###############

#load simulation results
load("data/cityrisk_simdata.RData")

cityrisk_simdata <- cityrisk_simdata %>%
  mutate(date = as.Date(date)) %>%
  mutate(month = month(date),
         day = day(date))

# 生成 imported_cases1：每年5~10月的每月1号输入1例
cityrisk_simdata <- cityrisk_simdata %>%
  mutate(imported_cases1 = ifelse(month %in% 5:10 & day == 1, 1, 0))

# 生成 imported_cases2：每年5~10月的每月1号和15号各输入1例
cityrisk_simdata <- cityrisk_simdata %>%
  mutate(imported_cases2 = ifelse(month %in% 5:10 & day %in% c(1, 15), 1, 0))

# 生成 imported_cases3：每年5~10月的每月1、8、15、22号各输入1例
cityrisk_simdata <- cityrisk_simdata %>%
  mutate(imported_cases3 = ifelse(month %in% 5:10 & day %in% c(1, 8, 15, 22), 1, 0))


# 生成 imported_cases1：每年5~10月的每月1号输入2例
cityrisk_simdata <- cityrisk_simdata %>%
  mutate(imported_cases4 = ifelse(month %in% 5:10 & day == 1, 2, 0))

# 生成 imported_cases2：每年5~10月的每月1号和15号各输入2例
cityrisk_simdata <- cityrisk_simdata %>%
  mutate(imported_cases5 = ifelse(month %in% 5:10 & day %in% c(1, 15), 2, 0))

# 生成 imported_cases3：每年5~10月的每月1、8、15、22号各输入2例
cityrisk_simdata <- cityrisk_simdata %>%
  mutate(imported_cases6 = ifelse(month %in% 5:10 & day %in% c(1, 8, 15, 22), 2, 0))

source("06_02_SEI-SEIR_model_THR with specific imported cases.R")



case_types <- c("imported_cases1", "imported_cases2", "imported_cases3")
case_types <- c("imported_cases4", "imported_cases5", "imported_cases6")

# 遍历不同输入频率
for (case_type in case_types) {
  
  # 创建保存模拟结果的数据框
  traitDF <- data.frame(matrix(ncol = 16, nrow = 0))
  colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R",
                         "dS", "dE", "dI", "dR", "dI_in", "Date", "Site")
  
  # 自动命名输出文件
  traitFileName <- paste0("output/Subsequent_city_risks_simulations_", case_type, ".csv")
  write.csv(traitDF, traitFileName, row.names = FALSE)
  
  # 设定SVPD
  city_climatedata$SVPD <- 1
  
  # 循环每一个城市进行模拟
  for (i in 1:length(sites)) {
    
    cityrisk_simdata2 <- cityrisk_simdata %>%
      dplyr::filter(city == sites[i] & date >= as.Date("2023-04-01"))
    
    climateData2 <- city_climatedata %>%
      dplyr::filter(Site == sites[i] & date >= as.Date("2023-04-01"))
    
    # 取不同的imported_cases定义
    imported_cases <- as.numeric(cityrisk_simdata2[[case_type]])
    
    temp <- as.numeric(climateData2$Temperature)
    rain <- as.numeric(climateData2$Rainfall)
    hum <- as.numeric(climateData2$SVPD)
    Rmax <- 123
    Date <- climateData2$date
    N <- population[i]
    city <- sites[i]
    BR <- BRs
    DR <- DRs
    times <- seq(1, length(Date), by = 1)
    K_thr <- K_thr_inverse
    M0 <- K_thr(temp[1], mean(rain), Rmax, N, timestep)
    parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep = timestep)
    
    # 初始状态设定
    state <- c(M1 = M0 - 1,
               M2 = 0,
               M3 = 1,
               S = N - 1,
               E = 0,
               I = 1,
               R = 0)
    
    # 求解ODE模型
    out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters,
               method = "rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
    
    out2 <- as.data.frame(out)
    out2$dS <- c(NA, diff(out2$S))
    out2$dE <- c(NA, diff(out2$E))
    out2$dI <- c(NA, diff(out2$I))
    out2$dR <- c(NA, diff(out2$R))
    
    out2$dI_in <- (1.0 / 5.9) * out2$E - DR * (out2$I / 1000) / 365
    out2$Date <- Date
    out2$Site <- sites[i]
    out2$initial_cases <- case_type
    
    # 添加到总数据框中
    traitDF <- rbind(traitDF, out2)
    
    # 写入CSV
    write.csv(traitDF, traitFileName, row.names = FALSE, fileEncoding = "GB18030")
    
    cat("Finished running ODE for", sites[i], "with", case_type, "cases\n")
  }
}
















# 创建保存模拟结果的数据框
traitDF <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R",
                       "dS", "dE", "dI", "dR", "dI_in", "Date", "Site")

# 自动命名输出文件
traitFileName <- paste0("output/Subsequent_city_risks_simulations_", "xsbn", ".csv")
write.csv(traitDF, traitFileName, row.names = FALSE)

# 设定SVPD
city_climatedata$SVPD <- 1
#i=1

# 循环每一个城市进行模拟
for (i in 1:length(sites)) {
  
  cityrisk_simdata2 <- cityrisk_simdata %>%
    dplyr::filter(city == sites[i] & date >= as.Date("2023-04-01"))
  
  climateData2 <- city_climatedata %>%
    dplyr::filter(Site == sites[i] & date >= as.Date("2023-04-01"))
  
  imported_cases <- as.numeric(cityrisk_simdata2$total_cases)
  temp <- as.numeric(climateData2$Temperature)
  rain <- as.numeric(climateData2$Rainfall)
  hum <- as.numeric(climateData2$SVPD)
  Rmax <- 123
  Date <- climateData2$date
  N <- population[i]
  city <- sites[i]
  BR <- BRs
  DR <- DRs
  times <- seq(1, length(Date), by = 1)
  K_thr <- K_thr_inverse
  M0 <- K_thr(temp[1], mean(rain), Rmax, N, timestep)
  parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep = timestep)
  
  # 初始状态设定
  state <- c(M1 = M0 - 1,
             M2 = 0,
             M3 = 1,
             S = N - 1,
             E = 0,
             I = 1,
             R = 0)
  
  # 求解ODE模型
  out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters,
             method = "rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
  
  out2 <- as.data.frame(out)
  out2$dS <- c(NA, diff(out2$S))
  out2$dE <- c(NA, diff(out2$E))
  out2$dI <- c(NA, diff(out2$I))
  out2$dR <- c(NA, diff(out2$R))
  
  out2$dI_in <- (1.0 / 5.9) * out2$E - DR * (out2$I / 1000) / 365
  out2$Date <- Date
  out2$Site <- sites[i]
  out2$initial_cases <- "xsbn"
  
  # 添加到总数据框中
  traitDF <- rbind(traitDF, out2)
  
  # 写入CSV
  write.csv(traitDF, traitFileName, row.names = FALSE, fileEncoding = "GB18030")
  
  cat("Finished running ODE for", sites[i], "with", "xsbn", "cases\n")
}


















head(traitDF)
traitFileName
# Load your simulation results
traitDF <- read.csv("output/Subsequent_city_risks_simulations_15.csv", fileEncoding = "GB18030") %>% 
  filter(Site == c("广州市"))

traitDF <- read.csv("output/Subsequent_city_risks_simulations_20.csv", fileEncoding = "GB18030") %>% 
  filter(Site == c("上海市"))

traitDF <- read.csv("output/Subsequent_city_risks_simulations_20.csv", fileEncoding = "GB18030") %>% 
  filter(Site == c("北京市"))

# Convert Date to Date type
traitDF$Date <- as.Date(traitDF$Date)

a1 <- unlist(lapply(temp, a))
b1 <- unlist(lapply(temp,b))

traitDF$temp <- temp
traitDF$a <- a1
traitDF$b <- b1
traitDF <- traitDF %>%
  mutate(infection_rate = a1 * b1 * M3)

# Reshape data for plotting
traitDF_long <- traitDF %>%
  pivot_longer(cols = c(S, E, I, R,dI_in), names_to = "State", values_to = "Population")

# Plotting the results with different colors for rainfall functions
ggplot(traitDF_long, aes(x = Date, y = Population)) +
  geom_line() +
  facet_wrap(~ State, scales = "free_y") +
  labs(title = "SEI-SEIR Model Results by Rainfall Function",
       y = "Population",
       x = "Date"
      ) +
  theme_minimal() +
  theme(legend.position = "bottom")

case <- traitDF_long %>%
  filter(
         State == "dI_in")
sum(case$Population)
# #6238

predicted_I <- traitDF$dI_in


plotdata <- data.frame(date=Date, predicted_I=predicted_I)

ggplot(plotdata, aes(x = date, y = predicted_I, color = "predicted")) +
  geom_line() + # 绘制预测值的线
  scale_color_manual(values = c("predicted" = "blue")) + # 为预测值和观测值设置颜色
  labs(x = "Date", y = "Value", color = "Legend", title = "Predicted vs Observed") + # 设置坐标轴标签和图例标题
  theme_minimal()

