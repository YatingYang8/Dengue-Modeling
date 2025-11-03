
################################ Calculate IR(t) ################################
library(dplyr)
library(readxl)
library(lubridate)

##### load modelling cases
traitDF <- read.csv("output/SEIR_simulations_result.csv") 

# Convert Date to Date type
traitDF$Date <- as.Date(traitDF$Date)

# Reshape data for plotting
traitDF_long <- traitDF %>%
  pivot_longer(cols = c(S, E, I, R,dI_in), names_to = "State", values_to = "Population")

# Plotting the results with different colors for rainfall functions
ggplot(traitDF_long, aes(x = Date, y = Population)) +
  geom_line() +
  facet_wrap(~ State, scales = "free_y") +
  labs(title = "SEIR Model Results",
       y = "Population",
       x = "Date") +
  theme_minimal() +
  theme(legend.position = "bottom")

case <- traitDF_long %>%
  filter( State == "dI_in")
sum(case$Population)
# #6238


#### load number of travellers (单位：万人次)
pop <- read_excel("E:/fdu/PhD project/Infectious disease/spatial/part2/Raw data/population.xlsx")
pop$Num_traveller_daily <- pop$Number_traveller/pop$monthdays

head(pop)

library(ggplot2)
library(dplyr)
library(scales)  # 用于日期格式化
library(extrafont)  # 使用正式字体

# 确保 Date 列是日期格式
pop$Date <- as.Date(pop$Date)

# 绘制柱状图
ggplot(pop, aes(x = Date, y = Number_traveller)) +
  geom_bar(stat = "identity", fill = "#4C72B0", color = "white", alpha = 0.8) +  # 柱状图样式
  labs(
    x = "Date (Year-Month)", 
    y = "Number of Travellers", 
    title = "Monthly Number of Travellers to Xishuangbanna",
    caption = "Data source: China Statistical Yearbook"  # 添加数据来源
  ) +
  theme_minimal(base_family = "Times New Roman") +  # 使用正式字体
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # 标题居中并加粗
    axis.title.x = element_text(size = 12, face = "bold"),  # X轴标题加粗
    axis.title.y = element_text(size = 12, face = "bold"),  # Y轴标题加粗
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # 调整X轴标签角度
    axis.text.y = element_text(size = 10),  # Y轴标签大小
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # 网格线样式
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加外边框
    plot.caption = element_text(size = 9, hjust = 1)  # 数据来源注释
  ) +
  scale_x_date(labels = date_format("%Y-%m"), breaks = date_breaks("months")) +  # 格式化日期轴
  scale_y_continuous(expand = c(0, 0))  # 去掉Y轴下方空白
library(ggplot2)
library(dplyr)
library(scales)  # 用于日期格式化
library(extrafont)  # 使用正式字体

# 确保 Date 列是日期格式
pop$Date <- as.Date(pop$Date)

# 绘制柱状图
ggplot(pop, aes(x = Date, y = Number_traveller)) +
  geom_bar(stat = "identity", fill = "#4C72B0", color = "white", alpha = 0.8) +  # 柱状图样式
  geom_text(aes(label = round(Number_traveller)), vjust = -0.5, size = 3.5, color = "black", family = "Times New Roman") +  # 在柱子顶端显示数字
  labs(
    x = "Date (Year-Month)", 
    y = expression(paste("Number of Travellers (", 10^4, ")")),
    title = "Monthly Number of Travellers",
    caption = "Data source: China Statistical Yearbook"  # 添加数据来源
  ) +
  theme_minimal(base_family = "Times New Roman") +  # 使用正式字体
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # 标题居中并加粗
    axis.title.x = element_text(size = 12, face = "bold"),  # X轴标题加粗
    axis.title.y = element_text(size = 12, face = "bold"),  # Y轴标题加粗
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # 调整X轴标签角度
    axis.text.y = element_text(size = 10),  # Y轴标签大小
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # 网格线样式
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加外边框
    plot.caption = element_text(size = 9, hjust = 1)  # 数据来源注释
  ) +
  scale_x_date(labels = date_format("%Y-%m"), breaks = date_breaks("months")) +  # 格式化日期轴
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # 调整Y轴范围，为数字留出空间


## IR
IRt <- traitDF %>%
  dplyr::select(time, Date, S, dI_in) %>%
  dplyr::mutate(IR = dI_in/S ,
                Year = year(Date),
                Month = month(Date)) 

IRt_data <- left_join(IRt,pop[,c("Year","Month","Num_traveller_daily")],by=c("Year","Month"))  
IRt_data$IR <- ifelse(IRt_data$IR > 0, IRt_data$IR, 0)

IRt_data <- IRt_data %>%
  mutate(traveller_cases = IR * Num_traveller_daily *10000) 


#plot
head(IRt_data)
ggplot(IRt_data, aes(x = as.factor(Month), y = traveller_cases)) +
  geom_point() +
  labs(title = "Cases from travellers",
       y = "Number of cases",
       x = "Date"
       ) +
  theme_minimal() +
  theme(legend.position = "bottom")



ggplot(IRt_data, aes(x = as.Date(Date), y = traveller_cases)) +
  geom_point() +
  scale_x_date(date_labels = "%m",date_breaks = "1 month") +
  labs(title = "Cases from travellers",
       y = "Number of cases",
       x = "Month"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# 确保 Month 是因子类型
IRt_data$Month <- as.factor(IRt_data$Month)

IRt_data_month <- IRt_data %>%
  group_by(Month) %>%
  summarise(traveller_cases = sum(traveller_cases)) 

ggplot(IRt_data_month, aes(x = as.factor(Month), y = traveller_cases)) +
  geom_bar(stat = "identity", fill = "#D62728", color = "white", alpha = 0.8) +  # 柱状图样式
  geom_text(aes(label = round(traveller_cases)), vjust = -0.5, size = 3.5, color = "black", family = "Times New Roman") +  # 在柱子顶端显示数字
  labs(
    x = "Month", 
    y = "Cases from travellers",
    title = "Monthly Number of Cases from Travellers"
  ) +
  theme_minimal(base_family = "Times New Roman") +  # 使用正式字体
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # 标题居中并加粗
    axis.title.x = element_text(size = 12, face = "bold"),  # X轴标题加粗
    axis.title.y = element_text(size = 12, face = "bold"),  # Y轴标题加粗
    axis.text.x = element_text(size = 10,  hjust = 1),  # 调整X轴标签角度
    axis.text.y = element_text(size = 10),  # Y轴标签大小
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # 网格线样式
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加外边框
    plot.caption = element_text(size = 9, hjust = 1)  # 数据来源注释
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # 调整Y轴范围，为数字留出空间

# 绘制图表

sum(IRt_data$traveller_cases)

#write.csv(IRt_data, "output/Travellercases_simulations_result.csv", row.names = F)
