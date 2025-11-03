library(dplyr)
library(reshape2)
library(lubridate)

setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue")

load()

####======================== 1.Climate dataset ================================#####

########### mean Temperature data from NOAA ############################################################
Tmean <- read.csv("D:/BaiduNetdiskDownload/【立方数据学社】1981-2023年逐日平均气温/【立方数据学社】1981-2023年逐日平均气温/【立方数据学社】1981-2023年逐日平均气温.csv",fileEncoding = "GB18030")  

Tmean=melt(Tmean,
           id=c(names(Tmean)[c(1:5)]),
           measure.vars=c(names(Tmean)[-c(1:5)]))

Tmean <- Tmean %>%
  mutate(variable = sub("^X(\\d{6})", "\\1", variable))

Tmean = Tmean %>%
  subset(select=-c(1))

names(Tmean) <- c("province","provinceid","city","cityid","date","Temperature")

# 使用dplyr的filter函数筛选日期范围内的数据
Tmean_filtered <- Tmean %>%
  filter(date >= "20230101" & date <= "20231231")

Tmean_filtered$date <-parse_date_time(Tmean_filtered$date,orders = "ymd")

colSums(is.na(Tmean_filtered))
Tmean <- na.omit(Tmean_filtered)

Tmean$date <- as.Date(Tmean$date)



################ Rainfall data from NOAA #######################################
Rainfall <- read.csv("D:/BaiduNetdiskDownload/【立方数据学社】1981-2023年逐日降水量/【立方数据学社】1981-2023年逐日降水量.csv",fileEncoding = "GB18030") 

Rainfall=melt(Rainfall,
              id=c(names(Rainfall)[c(1:5)]),
              measure.vars=c(names(Rainfall)[-c(1:5)]))

Rainfall <- Rainfall %>%
  mutate(variable = sub("^X(\\d{6})", "\\1", variable))

Rainfall = Rainfall %>%
  subset(select=-c(1))

names(Rainfall) <- c("province","provinceid","city","cityid","date","Rainfall")

# 使用dplyr的filter函数筛选日期范围内的数据
Rainfall_filtered <- Rainfall %>%
  filter(date >= "20230101" & date <= "20231231")

Rainfall_filtered$date <-parse_date_time(Rainfall_filtered$date,orders = "ymd")

colSums(is.na(Rainfall_filtered))
Rainfall <- na.omit(Rainfall_filtered)

Rainfall$date <- as.Date(Rainfall$date)

################################ combined ######################################
climateData <- left_join(Tmean,Rainfall)

climateData$Site <- climateData$city

climateData <- climateData %>%
  mutate(Temperature = as.numeric(Temperature),
         Rainfall = as.numeric(Rainfall))

climateData$date <- as.Date(climateData$date)

city_climatedata <-  climateData %>%
  select(date,Temperature,Rainfall,Site,cityid)


head(city_climatedata)


save(city_climatedata, file = "./data/city_climatedata.RData")



##################################################################################


#load simulation results
load(file = "./output/travel_simulation_results.RData")
#str(travel_simulation_results)
library(purrr)

travel_month_all <- map(travel_simulation_results, function(sim) {
  sim$transferred_cases %>%
    mutate(month=month(Date)) %>%
    group_by(city,month) %>%
    summarise(total_cases = sum(cases, na.rm = TRUE), .groups = 'drop')
})

travel_month_mean <- bind_rows(travel_month_all) %>%
  group_by(city,month) %>%
  summarise(travel_cases = mean(total_cases, na.rm = TRUE), .groups = 'drop')



load(file = "./output/local_simulation_results.RData")

# 计算每次模拟中每个城市的总迁移病例数
local_month_all <- map(local_simulation_results, function(sim) {
  sim$transferred_cases %>%
    mutate(month=month(Date)) %>%
    group_by(city,month) %>%
    summarise(total_cases = sum(cases, na.rm = TRUE), .groups = 'drop')
})

local_month_mean <- bind_rows(local_month_all) %>%
  group_by(city,month) %>%
  summarise(local_cases = mean(total_cases, na.rm = TRUE), .groups = 'drop')

Total_month_mean <- left_join(travel_month_mean,local_month_mean) %>%
  mutate(total_cases = local_cases + travel_cases)

# Total_case_mean <- Total_month_mean %>%
#   group_by(city) %>%
#   summarise(local_cases = sum(local_cases, na.rm = TRUE),
#             travel_cases = sum(travel_cases, na.rm = TRUE),
#             total_cases = sum(total_cases, na.rm = TRUE))



# 提取1次模拟中的转移病例数据
travel_sim1 <- travel_simulation_results[[3]][[1]]

travel_sim1 <- travel_sim1 %>%
  group_by(city,Date) %>%
  summarise(travel_cases = sum(cases))


local_sim1 <- local_simulation_results[[3]][[1]]

local_sim1 <- local_sim1 %>%
  group_by(city,Date) %>%
  summarise(local_cases = sum(cases))

# 所有可能的城市和月份组合
all_cities <- unique(unlist(map(travel_simulation_results, ~ unique(.x$transferred_cases$city))))

all_dates <- seq.Date(from = as.Date("2023-01-01"), to = as.Date("2023-12-31"), by = "day")

full_grid <- expand.grid(city = all_cities, date = all_dates) %>%
  mutate(month = month(date))


Totalcase_city_day <- full_grid %>%
  left_join(Total_month_mean) %>%
  mutate(day = day(date)) %>%
  mutate(local_cases = replace_na(local_cases, 0)) %>%
  mutate(travel_cases = replace_na(travel_cases, 0)) %>%
  mutate(total_cases = local_cases+travel_cases) %>%
  mutate(
    local_cases = if_else(day == 1, local_cases, 0),
    travel_cases = if_else(day == 1, travel_cases, 0),
    total_cases = if_else(day == 1, total_cases, 0)
  ) %>%
  mutate(time = as.numeric(date - as.Date("2023-01-01")) + 1)


names(Totalcase_city_day)

#########################################################################################

# load climate data
load("data/city_climatedata.RData")

# load simulated imported cases
load("data/citydata_ALL.RData")

cityrisk_simdata <- left_join(city_climatedata, citydata_ALL[,c("city","cityid","pop_total")]) %>%
  left_join(Totalcase_city_day[,c("city","date","total_cases")])

cityrisk_simdata <- na.omit(cityrisk_simdata)


save(cityrisk_simdata, file = "./data/cityrisk_simdata.RData")

