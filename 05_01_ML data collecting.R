
################# Machine learning data collecting ##############################
library(dplyr); library(raster); library(rgdal); library(Matrix)
library(stringr); library(ggplot2); library(lubridate);library(cowplot);library(colorspace)
library(magrittr);  library(spdep);library(rgeos);library(sf);library(sp);library(ggspatial)
library(readxl)
################### 1.outcome ###################################################
######### Traveller cases migration simulation ##################
load(file = "./output/travel_simulation_results.RData")

# 计算1000次模拟的平均结果
# 计算每次模拟的总迁移病例数
library(purrr)

# 计算每次模拟中每个城市的总迁移病例数
city_simulation_sums <- map(travel_simulation_results, function(sim) {
  sim$transferred_cases %>%
    group_by(city) %>%
    summarise(total_cases = sum(cases, na.rm = TRUE), .groups = 'drop')
})

travel_city_transfers <- bind_rows(city_simulation_sums) %>%
  group_by(city) %>%
  summarise(travel_cases = mean(total_cases, na.rm = TRUE), .groups = 'drop')



load(file = "./output/local_simulation_results.RData")

# 计算每次模拟中每个城市的总迁移病例数
city_simulation_sums <- map(local_simulation_results, function(sim) {
  sim$transferred_cases %>%
    group_by(city) %>%
    summarise(total_cases = sum(cases, na.rm = TRUE), .groups = 'drop')
})

local_city_transfers <- bind_rows(city_simulation_sums) %>%
  group_by(city) %>%
  summarise(local_cases = mean(total_cases, na.rm = TRUE), .groups = 'drop')

Total_city_transfers <- left_join(travel_city_transfers,local_city_transfers) %>%
  mutate(total_cases = local_cases + travel_cases)


############first arrival time
# 计算 1000 次首达时的均值
travel_first_arrival <- bind_rows(lapply(travel_simulation_results, `[[`, "first_arrival")) %>%
  group_by(city) %>%
  summarise(travel_first_date = median(first_date), .groups = 'drop') #不服从正态分布

local_first_arrival <- bind_rows(lapply(local_simulation_results, `[[`, "first_arrival")) %>%
  group_by(city) %>%
  summarise(local_first_date = median(first_date), .groups = 'drop')


Total_first_arrival <- left_join(travel_first_arrival,local_first_arrival) %>%
  mutate(
    first_date = pmin(travel_first_date, local_first_date, na.rm = TRUE),  # 取最早日期
    time = as.integer(first_date - as.Date("2023-01-01")) + 1  # 计算天数，2023-01-01 对应 time = 1
  )


Outcome <- left_join(Total_city_transfers,Total_first_arrival) %>%
  dplyr::select("city","total_cases","time")



################### 2.Features ###################################################
setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/Raw data/city features")

# # 社会经济维度 Socio-Economic Dimensions --------------------------------
# ## 人口基础
# pop_total         # 年末常住人口（万人）    
# pop_density       # 人口密度（人/平方公里）  
# gender_ratio      # 性别比                  
# urbanization      # 城镇化率（%）            

pop_total <- read.csv("./pop_total.csv",fileEncoding ="GB18030" )
urbanization <- read_excel("./urbanization.xlsx")

pop_density <- read_excel("./pop_density.xlsx") %>%
  dplyr::select("地区名称", "市代码","2022") %>%
  dplyr::rename("city" = "地区名称",
                "cityid" = "市代码",
                "pop_density" = "2022")


# ## 经济规模与结构
# gdp_total         # 地区生产总值（亿元）      
# gdp_pct_primary   # 第一产业占GDP比重（%）     
# gdp_pct_secondary # 第二产业占GDP比重（%）     
# gdp_pct_tertiary  # 第三产业占GDP比重（%）     
# gdp_per_capita    # 人均GDP（元）            

gdp_total <- read_excel("./gdp_total.xlsx")
gdp_pct_primary <- read_excel("./gdp_pct_primary.xlsx")

gdp_pct_secondary <- read_excel("./gdp_pct_secondary.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "gdp_pct_secondary" = "2022")

gdp_pct_tertiary <- read_excel("./gdp_pct_tertiary.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "gdp_pct_tertiary" = "2022")

gdp_per_capita <- read_excel("./gdp_per_capita.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "gdp_per_capita" = "2022")



# ## 就业与收入
# employee_avg      # 在岗职工平均人数（万人）    
# wage_total        # 在岗职工工资总额（万元）    
# wage_avg          # 在岗职工平均工资（元）      

employee_avg <- read_excel("./employee_avg.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2019") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "employee_avg" = "2019")

wage_total <- read_excel("./wage_total.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2019") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "wage_total" = "2019")

wage_avg <- read_excel("./wage_avg.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2019") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "wage_avg" = "2019")



# # 空间交互维度 Spatial Interaction --------------------------------
# distance     # 与西双版纳测地线距离（公里）
# passenger_highway # 公路客运量（万人）          
# passenger_water   # 水运客运量（万人）          
# passenger_air     # 民用航空客运量（万人） 
# passenger_railway     # 铁路客运量（万人）
# passenger_total   # 客运总量（万人）            

library(geosphere)

passenger_highway <- read_excel("./passenger_highway.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "passenger_highway" = "2022")


passenger_water <- read_excel("./passenger_water.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2019") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "passenger_water" = "2019")

passenger_air <- read_excel("./passenger_air.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2019") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "passenger_air" = "2019")

passenger_railway <- read_excel("./passenger_railway.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2014") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "passenger_railway" = "2014")

passenger_total <- read_excel("./passenger_total.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2014") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "passenger_total" = "2014")


# # 城市设施与资源维度 Urban Facilities --------------------------------
# ## 市政基础设施
# water_supply      # 供水普及率（%）            
# gas_supply        # 燃气普及率（%）            
# road_area         # 人均道路面积（平方米）      
# sewage_treatment  # 污水处理率（%）            
# waste_treatment   # 生活垃圾处理率（%）        
# green_space       # 人均公园绿地面积（平方米）  

water_supply<- read_excel("./water_supply.xlsx") %>%
  dplyr::select("地区名称", "市代码","2022") %>%
  dplyr::rename("city" = "地区名称",
                "cityid" = "市代码",
                "water_supply" = "2022")

gas_supply <- read_excel("./gas_supply.xlsx") %>%
  dplyr::select("地区名称", "市代码","2022") %>%
  dplyr::rename("city" = "地区名称",
                "cityid" = "市代码",
                "gas_supply" = "2022")

road_area <- read_excel("./road_area.xlsx") %>%
  dplyr::select("地区名称", "市代码","2022") %>%
  dplyr::rename("city" = "地区名称",
                "cityid" = "市代码",
                "road_area" = "2022")

sewage_treatment <- read_excel("./sewage_treatment.xlsx") %>%
  dplyr::select("地区名称", "市代码","2022") %>%
  dplyr::rename("city" = "地区名称",
                "cityid" = "市代码",
                "sewage_treatment" = "2022")

waste_treatment <- read_excel("./waste_treatment.xlsx") %>%
  dplyr::select("地区名称", "市代码","2022") %>%
  dplyr::rename("city" = "地区名称",
                "cityid" = "市代码",
                "waste_treatment" = "2022")

green_space <- read_excel("./green_space.xlsx") %>%
  dplyr::select("地区名称", "市代码","2022") %>%
  dplyr::rename("city" = "地区名称",
                "cityid" = "市代码",
                "green_space" = "2022")


# ## 医疗教育资源
# hospitals         # 医院数（个）              
# hospital_beds     # 医院床位数（张）          
# doctors           # 执业或助理医师数（人）      
# school_uni        # 普通高等学校数（所）        
# school_high       # 普通中学数（所）            
# school_primary    # 普通小学数（所）            
# school_vocational # 中等职教学校数（所）        

hospitals <- read_excel("./hospitals.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "hospitals" = "2022")

hospital_beds <- read_excel("./hospital_beds.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "hospital_beds" = "2022")

doctors <- read_excel("./doctors.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "doctors" = "2022")

school_uni <- read_excel("./school_uni.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "school_uni" = "2022")

school_high <- read_excel("./school_high.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "school_high" = "2022")

school_primary <- read_excel("./school_primary.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "school_primary" = "2022")

school_vocational <- read_excel("./school_vocational.xlsx") %>%
  dplyr::select("CITY" ,"CITY_ID","2022") %>%
  dplyr::rename("city" = "CITY",
                "cityid" = "CITY_ID",
                "school_vocational" = "2022")


All <- left_join(Outcome,pop_total) %>%
  left_join(urbanization) %>%
  left_join(pop_density) %>%
  left_join(gdp_total) %>%
  left_join(gdp_pct_primary) %>%
  left_join(gdp_pct_secondary) %>%
  left_join(gdp_pct_tertiary) %>%
  left_join(gdp_per_capita) %>%
  left_join(employee_avg) %>%
  left_join(wage_avg) %>%
  left_join(wage_total) %>%
  left_join(passenger_highway) %>%
  left_join(passenger_water) %>%
  left_join(passenger_air) %>%
  left_join(passenger_railway) %>%
  left_join(passenger_total) %>%
  left_join(water_supply) %>%
  left_join(gas_supply) %>%
  left_join(sewage_treatment) %>%
  left_join(waste_treatment) %>%
  left_join(green_space) %>%
  left_join(road_area) %>%
  left_join(hospitals) %>%
  left_join(hospital_beds) %>%
  left_join(doctors) %>%
  left_join(school_uni) %>%
  left_join(school_high) %>%
  left_join(school_primary) %>%
  left_join(school_vocational) 
  

colSums(is.na(All))

# distance_xsbn     # 与西双版纳测地线距离（公里）

# 示例：使用tidygeocoder获取坐标
if (!require("tidygeocoder")) install.packages("tidygeocoder")
library(tidygeocoder)

cities <- data.frame(city = All$city)
cities <- geocode(cities, city, method = "arcgis")  # 使用ArcGIS服务

# 西双版纳中心点坐标
xishuangbanna <- c(100.8030, 22.0094)  

# 提取城市坐标矩阵
city_coords <- as.matrix(cities[, c("long", "lat")])


dist_matrix <- distm(
  x = city_coords,
  y = xishuangbanna,  # 如果计算多个点间的距离，这里也可以是矩阵
  fun = distGeo        # 使用最精确的WGS84椭球方法
)

All$distance <- dist_matrix[,1]

ML_origindata <- All

save(ML_origindata, file = "E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue/data/ML_origindata.RData")
write.csv(ML_origindata, file = "./data/ML_origindata.csv",fileEncoding = "GB18030")



####### add attractiveness ranking of city #####################################
load(file = "./data/ML_origindata.RData")
head(ML_origindata)
#####city class data #######
# Select required fields
library(jsonlite)
library(dplyr)
library(readr)

# Read the JSON file
json_data <- fromJSON("E:/fdu/PhD project/Infectious disease/spatial/part2/Raw data/city features/城市.txt")

city_rank <- json_data %>%
  dplyr::select(城市, city_id, 省份, 城市级别, 全国排名) %>%
  dplyr::rename(city = "城市",
         cityid = city_id,
         province = "省份",
         classification = "城市级别",
         rank = "全国排名"
         ) %>%
  dplyr::mutate(cityid = as.numeric(cityid))

citydata_ALL <- left_join(ML_origindata,city_rank[,c("cityid","classification","rank")], by="cityid")

colSums(is.na(citydata_ALL))


####### add south and north area #####################################
head(citydata_ALL)

unique(citydata_ALL$province)

province_area <- tibble(
  province = c(
    # 南方
    "江苏省","浙江省","安徽省","福建省","江西省","上海市","广东省","广西壮族自治区","海南省",
    "四川省","贵州省","云南省","重庆市","湖北省","湖南省",
    # 北方
    "北京市","天津市","河北省","山西省","内蒙古自治区","辽宁省","吉林省","黑龙江省",
    "陕西省","甘肃省","青海省","宁夏回族自治区","新疆维吾尔自治区","河南省","山东省",
    # 其他（西藏通常划为西南，这里算南方）
    "西藏自治区"
  ),
  Area = c(
    rep("South", 15),
    rep("North", 15),
    "South"  # 西藏
  )
)

# 合并到你的数据
citydata_ALL <- citydata_ALL %>%
  left_join(province_area, by = "province")

# 检查结果
table(citydata_ALL$Area, useNA = "ifany")




##very sb area2
# 定义区域分类
east <- c("北京市", "天津市", "河北省", "上海市", "江苏省", "浙江省", "福建省", "山东省", "广东省", "海南省")
central <- c("山西省", "安徽省", "江西省", "河南省", "湖北省", "湖南省")
west <- c("内蒙古自治区", "广西壮族自治区", "重庆市", "四川省", "贵州省", 
          "云南省", "西藏自治区", "陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区")
northeast <- c("辽宁省", "吉林省", "黑龙江省")

# 添加 Area2 变量
citydata_ALL$Area2 <- dplyr::case_when(
  citydata_ALL$province %in% east ~ "East",
  citydata_ALL$province %in% central ~ "Central",
  citydata_ALL$province %in% west ~ "West",
  citydata_ALL$province %in% northeast ~ "North east",
  TRUE ~ NA_character_  # 万一有不在分类里的
)

# 检查结果
table(citydata_ALL$Area2, useNA = "ifany")


save(citydata_ALL, file = "E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue/data/citydata_ALL.RData")
write.csv(citydata_ALL, file = "./data/citydata_ALL.csv",fileEncoding = "GB18030")

load(file = "E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue/data/citydata_ALL.RData")
