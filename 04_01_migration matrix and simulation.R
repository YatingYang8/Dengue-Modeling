## Mobility matrix
library(readxl) 
library(dplyr)
library(tidyr)
library(reshape2)
library(lubridate)
library(ggplot2)

# Move-out Movement Index
file_path <- "./part2/Raw data/Baidu migration/2023迁徙规模指数.xlsx"

data <- read_excel(file_path, sheet = 1)

head(data)

Moveout_index <- data %>%
  filter(城市 == "西双版纳傣族自治州") %>%
  filter(类型 == "迁出") %>%
  dplyr::select("城市", "日期", "迁徙指数")

colnames(Moveout_index) <- c("city","Date","out_index")

Moveout_index$Date <- as.Date(Moveout_index$Date, format = "%Y%m%d")

head(Moveout_index)


##提取上海的迁徙指数备用
shanghai <- data %>%
  filter(城市 == "上海") %>%
  filter(类型 == "迁出") %>%
  dplyr::select("城市", "日期", "迁徙指数")

colnames(shanghai) <- c("city","Date","out_index")

shanghai$Date <- as.Date(shanghai$Date, format = "%Y%m%d")

head(shanghai)



# Move-out proportion
# 导入数据
file_path <- "E:/fdu/PhD project/Infectious disease/spatial/part2/Raw data/Baidu migration/西双版纳傣族自治州.csv"
data <- read.csv(file_path, stringsAsFactors = FALSE)
head(data)

names(data)[names(data) == "地区A"] <- "city"
names(data)[names(data) == "地区A代码"] <- "cityid"
names(data)[names(data) == "日期"] <- "Date"

data$Date <- as.Date(data$Date, format = "%Y-%m-%d")
data$cityid <- as.character(data$cityid)


# 进行left_join,使用Date和city作为连接键
Moveout_matrix <- left_join(data, Moveout_index)

# 查看合并后的数据
head(Moveout_matrix)

#标准化A迁出B占A迁出总人口的比值使其总和为100 
Moveout_matrix <- Moveout_matrix %>%
  group_by(Date) %>%
  mutate(scaled_proportion = (`A迁出B占A迁出总人口的比值` / sum(`A迁出B占A迁出总人口的比值`)) * 100) %>%
  ungroup()

Moveout_matrix <- Moveout_matrix %>%
  mutate(out_total = scaled_proportion * out_index) %>%
  filter(year(Date) == 2023)

head(Moveout_matrix)

# 加载dplyr包
library(dplyr)


Moveout_matrix <- Moveout_matrix %>%
  mutate(B所属省份 = case_when(
    地区B代码 == 110000 ~ "北京市",  # 如果地区B代码是310000，B所属省份为北京市
    地区B代码 == 120000 ~ "天津市",  # 如果地区B代码是120000，B所属省份为天津市
    地区B代码 == 310000 ~ "上海市", 
    地区B代码 == 500000 ~ "重庆市", 
    TRUE ~ B所属省份               # 其他情况保持原来的B所属省份
  ))

# 查看修改后的数据
head(Moveout_matrix)


Moveout_day <- Moveout_matrix %>%
  group_by(Date) %>%
  summarise(out_total = sum(out_total))

Moveout_area <- Moveout_matrix %>%
  group_by(地区B,地区B代码,B所属省份) %>%
  summarise(out_total = sum(out_total))



##### 重新整理数据，填补缺失的日期，获得完整数据集
full_dates <- seq(as.Date("2023-01-01"), as.Date("2023-12-31"), by = "day")  # 生成完整的日期序列

# 获取所有地区的唯一组合
unique_areas <- Moveout_matrix %>%
  dplyr::select(地区B, 地区B代码, B所属省份) %>%
  distinct()

area_mean <- Moveout_matrix %>%
  group_by(地区B,地区B代码,B所属省份) %>%
  summarise(A迁出B占A迁出总人口的比值mean = mean(A迁出B占A迁出总人口的比值))


# 生成所有日期和地区的组合
complete_data <- expand.grid(Date = full_dates, 地区B = unique_areas$地区B ) %>%
  left_join(Moveout_index[,c("Date","out_index")], by = c("Date")) %>%
  left_join(Moveout_matrix[,c("Date","地区B","A迁出B占A迁出总人口的比值")]) %>%
  left_join(area_mean) %>%
  mutate( proportion = ifelse(is.na(A迁出B占A迁出总人口的比值), A迁出B占A迁出总人口的比值mean ,A迁出B占A迁出总人口的比值)) 

Moveout_matrix_complete <- complete_data %>%
  group_by(Date) %>%
  mutate(scaled_proportion = (proportion / sum(proportion)) * 100) %>%
  mutate(out_total = out_index * scaled_proportion * 1000) %>%
  ungroup()
  
# 按日期排序
Moveout_matrix_complete <- Moveout_matrix_complete %>%
  arrange(Date, 地区B代码)

library(openxlsx)
#write.xlsx(Moveout_matrix_complete, "./data/Moveout_matrix.xlsx")





##############################计算连接概率矩阵###################################
# 定义Sigmoid函数
max_val <- max(Moveout_index$out_index)

sigmoid <- function(x) {
  return(1 / (1 + exp(-x + max_val / 2)))
}

Moveout_index$s_index <- sigmoid(Moveout_index$out_index)
shanghai$s_sh <- sigmoid(shanghai$out_index)

Moveout_index$link_prob <- (Moveout_index$s_index/shanghai$s_sh) * 0.01








################### 本地病例的跨城市传播风险模拟###################################
library(dplyr)

set.seed(123) 

##### load modelling local cases
traitDF <- read.csv("output/SEIR_simulations_result.csv") 

traitDF$Date <- as.Date(traitDF$Date)

local_cases <- traitDF %>%
  dplyr::select("Date","I") %>%
  mutate(local_cases = ifelse(I < 1, 0, round(I)))

sum(local_cases$local_cases)


sim_local_data <- Moveout_index %>%
  left_join(local_cases) 

head(sim_local_data)

# 设定模拟次数
n_sim <- 1000

# 存储每次模拟的结果
local_simulation_results <- list()

for (sim in 1:n_sim) {
  
  # 复制数据，避免修改原数据
  sim_data <- sim_local_data
  
  # 存储转移病例数
  transferred_cases <- data.frame()
  
  # 记录每个城市的首达时间
  first_arrival <- data.frame()
  
  # 遍历每一天的病例数据
  for (i in 1:nrow(sim_data)) {
    
    day_data <- sim_data[i, ]  # 获取当前行数据
    
    # 提取日期、城市、病例数、迁移概率
    date <- day_data$Date
    city <- day_data$city
    local_cases <- day_data$local_cases
    link_prob <- day_data$link_prob
    
    # 确定迁移的病例数（按照 link_prob 计算）
    moving_cases <- rbinom(1, local_cases, link_prob)  
    
    if (moving_cases > 0) {
      # 获取该日期对应的迁移概率数据
      move_probs <- Moveout_matrix_complete %>%
        filter(Date == date, city == city) %>%
        select(地区B, scaled_proportion)
      
      # 归一化概率，确保总和为1
      move_probs$scaled_proportion <- move_probs$scaled_proportion / sum(move_probs$scaled_proportion)
      
      # 按照 scaled_proportion 进行目标城市分配
      target_cities <- sample(move_probs$地区B, moving_cases, replace = TRUE, prob = move_probs$scaled_proportion)
      
      # 统计各城市接收到的病例数
      city_counts <- as.data.frame(table(target_cities))
      names(city_counts) <- c("city", "cases")
      city_counts$Date <- date
      
      # 记录转移病例数
      transferred_cases <- rbind(transferred_cases, city_counts)
      
      # 更新首达时间
      for (j in 1:nrow(city_counts)) {
        target_city <- city_counts$city[j]
        if (!(target_city %in% first_arrival$city)) {
          first_arrival <- rbind(first_arrival, data.frame(city = target_city, first_date = date))
        }
      }
    }
  }
  
  # 存储单次模拟的结果
  local_simulation_results[[sim]] <- list(transferred_cases = transferred_cases, first_arrival = first_arrival)
}

save(local_simulation_results, file = "./output/local_simulation_results.RData")
load(file = "./output/local_simulation_results.RData")

# 计算1000次模拟的平均结果
# 计算每次模拟的总迁移病例数
library(purrr)
simulation_sums <- map_dbl(local_simulation_results, ~ sum(.x$transferred_cases$cases, na.rm = TRUE))

# 计算 1000 次模拟的均值
mean_transferred_cases <- mean(simulation_sums, na.rm = TRUE)

# 计算 1000 次首达时的均值
final_first_arrival <- bind_rows(lapply(local_simulation_results, `[[`, "first_arrival")) %>%
  group_by(city) %>%
  summarise(first_date = mean(first_date), .groups = 'drop')

# 计算每次模拟中每个城市的总迁移病例数
city_simulation_sums <- map(local_simulation_results, function(sim) {
  sim$transferred_cases %>%
    group_by(city) %>%
    summarise(total_cases = sum(cases, na.rm = TRUE), .groups = 'drop')
})

mean_city_transfers <- bind_rows(city_simulation_sums) %>%
  group_by(city) %>%
  summarise(mean_cases = mean(total_cases, na.rm = TRUE), .groups = 'drop')









################### 游客病例的跨城市传播风险模拟###################################
set.seed(123)  # 保持可复现性

##### load modelling traveller cases
travel_cases <- read.csv("output/Travellercases_simulations_result.csv") %>%
  dplyr::select(Date,traveller_cases) %>%
  mutate(Date = as.Date(Date))

head(travel_cases)
sum(travel_cases$traveller_cases)

# 模拟每个游客的停留天数
simulate_departures <- function(travel_cases) {
  
  travel_cases_expanded <- travel_cases %>%
    rowwise() %>%
    mutate(traveller_cases = round(traveller_cases)) %>%  # 四舍五入得到整数
    ungroup() %>%
    uncount(traveller_cases)  # 展开每个病例
  
  # 模拟每个人的停留天数（weibull分布 λ=3）
  travel_cases_expanded <- travel_cases_expanded %>%
    mutate(stay_days = ceiling(rweibull(n(), shape = 2, scale = 3)), 
           departure_date = Date + days(stay_days))  # 计算离开日期
  
  # 统计每天的迁出病例数
  departures_by_date <- travel_cases_expanded %>%
    group_by(departure_date) %>%
    summarise(departed_cases = n(), .groups = "drop")
  
  return(departures_by_date)
}

# 设定模拟次数
n_simulations <- 1000

# 运行 1000 次模拟
traveldays_results <- map(1:n_simulations, ~ simulate_departures(travel_cases))

# 合并所有模拟结果
all_departures <- bind_rows(traveldays_results, .id = "sim_id")


head(all_departures)

#跨城市传播风险模拟

n_sim <- 1000

# 存储每次模拟的结果
travel_simulation_results <- list()

for (sim in 1:n_sim) {
  
  # 复制数据，避免修改原数据
  sim_data <- all_departures %>%
    filter(sim_id == sim)
  
  # 存储转移病例数
  transferred_cases <- data.frame()
  
  # 记录每个城市的首达时间
  first_arrival <- data.frame()
  
  # 遍历每一天的病例数据
  for (i in 1:nrow(sim_data)) {
    
    day_data <- sim_data[i, ]  # 获取当前行数据
  
    # 提取日期、城市、病例数、迁移概率
    date <- day_data$departure_date
    departed_cases <- day_data$departed_cases
    
    # 确定迁移的病例数（全部都迁移）
    moving_cases <- departed_cases  
    
    if (moving_cases > 0) {
      # 获取该日期对应的迁移概率数据
      move_probs <- Moveout_matrix_complete %>%
        filter(Date == date) %>%
        dplyr::select(地区B, scaled_proportion)
      
      # 归一化概率，确保总和为1
      move_probs$scaled_proportion <- move_probs$scaled_proportion / sum(move_probs$scaled_proportion)
      
      # 按照 scaled_proportion 进行目标城市分配
      target_cities <- sample(move_probs$地区B, moving_cases, replace = TRUE, prob = move_probs$scaled_proportion)
      
      # 统计各城市接收到的病例数
      city_counts <- as.data.frame(table(target_cities))
      names(city_counts) <- c("city", "cases")
      city_counts$Date <- date
      
      # 记录转移病例数
      transferred_cases <- rbind(transferred_cases, city_counts)
      
      # 更新首达时间
      for (j in 1:nrow(city_counts)) {
        target_city <- city_counts$city[j]
        if (!(target_city %in% first_arrival$city)) {
          first_arrival <- rbind(first_arrival, data.frame(city = target_city, first_date = date))
        }
      }
    }
  }
  
  # 存储单次模拟的结果
  travel_simulation_results[[sim]] <- list(transferred_cases = transferred_cases, first_arrival = first_arrival)
}

save(travel_simulation_results, file = "./output/travel_simulation_results.RData")
load(file = "./output/travel_simulation_results.RData")

# 计算1000次模拟的平均结果
# 计算每次模拟的总迁移病例数
simulation_sums <- map_dbl(travel_simulation_results, ~ sum(.x$transferred_cases$cases, na.rm = TRUE))

# 计算 1000 次模拟的均值
mean_transferred_cases <- mean(simulation_sums, na.rm = TRUE)

# 计算 1000 次首达时的均值
final_first_arrival <- bind_rows(lapply(travel_simulation_results, `[[`, "first_arrival")) %>%
  group_by(city) %>%
  summarise(first_date = mean(first_date), .groups = 'drop')

# 计算每次模拟中每个城市的总迁移病例数
city_simulation_sums <- map(travel_simulation_results, function(sim) {
  sim$transferred_cases %>%
    group_by(city) %>%
    summarise(total_cases = sum(cases, na.rm = TRUE), .groups = 'drop')
})

mean_city_transfers <- bind_rows(city_simulation_sums) %>%
  group_by(city) %>%
  summarise(mean_cases = mean(total_cases, na.rm = TRUE), .groups = 'drop')








