# Run SEI-SEIR model with different rainfall functions ------------------------------
rm(list=ls()) #remove previous variable assignments

# load packages
library(deSolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# load data 
#source("02_2_SEI-SEIR_model_THR.R")
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

traitDF$temp <- temp
traitDF$a <- a1
traitDF$b <- b1
traitDF <- traitDF %>%
  mutate(infection_rate = a1 * b1 * M3)

# Reshape data for plotting
traitDF_long <- traitDF %>%
  pivot_longer(cols = c(S, E, I, R,dI_in), names_to = "State", values_to = "Population")

# Plotting the results with different colors for rainfall functions
ggplot(traitDF_long, aes(x = Date, y = Population, color = Rain_function)) +
  geom_line() +
  facet_wrap(~ State, scales = "free_y") +
  labs(title = "SEI-SEIR Model Results by Rainfall Function",
       y = "Population",
       x = "Date",
       color = "Rain Function") +
  theme_minimal() +
  theme(legend.position = "bottom")

case <- traitDF_long %>%
  filter(Rain_function == "Inverse",
         State == "dI_in")
sum(case$Population)
# #6238

predicted_I <- traitDF$dI_in

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
