library(dplyr)


load(file = "E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue/data/citydata_ALL.RData")
head(citydata_ALL)

province_data <- citydata_ALL %>%
  group_by(province) %>%
  summarise(total_cases_sum = sum(total_cases, na.rm = TRUE)) %>%
  mutate(perc = total_cases_sum / sum(total_cases_sum) ,
         scale_case = scale(total_cases_sum))

head(province_data)


val_data <- read.csv("E:/fdu/PhD project/Infectious disease/spatial/part2/Raw data/validation data.csv",fileEncoding = "GB18030")

plot_data <- left_join(province_data, val_data) %>%
  mutate(
    perc_val = val_case / sum(val_case, na.rm = TRUE),
    scale_val_case = scale(val_case)
  )

head(plot_data)
colnames(plot_data)

plot_data_clean <- na.omit(plot_data)

library(ggplot2)
library(patchwork)

# 计算相关性和p值的函数
get_cor_label <- function(x, y) {
  ct <- cor.test(x, y, use = "complete.obs")
  r <- round(ct$estimate, 3)
  p <- signif(ct$p.value, 3)
  paste0("r = ", r, ", p = ", p)
}

# ---- 1. 实际值 vs 模拟值 ----
label1 <- get_cor_label(plot_data_clean$total_cases_sum, plot_data_clean$val_case)
p1 <- ggplot(plot_data_clean, aes(x = total_cases_sum, y = val_case)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = 2) +
  geom_text(aes(
    x = max(total_cases_sum, na.rm = TRUE) * 0.7,
    y = max(val_case, na.rm = TRUE) * 0.9,
    label = label1
  ), inherit.aes = FALSE) +
  labs(x = "Simulated Cases (total cases)",
       y = "Observed Cases (val_case)",
       title = "Correlation: Observed vs Simulated (Absolute Values)") +
  theme_classic()

p1

# ---- 2. 百分比 vs 百分比 ----
label2 <- get_cor_label(plot_data_clean$perc, plot_data_clean$perc_val)
p2 <- ggplot(plot_data_clean, aes(x = perc, y = perc_val)) +
  geom_point(color = "darkgreen", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = 2) +
  geom_text(aes(
    x = max(perc, na.rm = TRUE) * 0.7,
    y = max(perc_val, na.rm = TRUE) * 0.9,
    label = label2
  ), inherit.aes = FALSE) +
  labs(x = "Simulated Percentage (perc)",
       y = "Observed Percentage (perc_val)",
       title = "Correlation: Observed vs Simulated (Percentages)") +
  theme_minimal()

p2

# ---- 3. 标准化值 vs 标准化值 ----
label3 <- get_cor_label(plot_data_clean$scale_case, plot_data_clean$scale_val_case)
p3 <- ggplot(plot_data_clean, aes(x = scale_case, y = scale_val_case)) +
  geom_point(color = "purple", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = 2) +
  geom_text(aes(
    x = max(scale_case, na.rm = TRUE) * 0.7,
    y = max(scale_val_case, na.rm = TRUE) * 0.9,
    label = label3
  ), inherit.aes = FALSE) +
  labs(x = "Simulated (Standardized)",
       y = "Observed (Standardized)",
       title = "Correlation: Observed vs Simulated (Standardized)") +
  theme_minimal()

p3
# ---- 拼图 ----
p1
p2


library(dplyr)

plot_data_clean <- plot_data_clean %>%
  mutate(province = recode(province,
                           "上海市" = "Shanghai",
                           "北京市" = "Beijing",
                           "吉林省" = "Jilin",
                           "广东省" = "Guangdong",
                           "江苏省" = "Jiangsu",
                           "江西省" = "Jiangxi",
                           "浙江省" = "Zhejiang",
                           "海南省" = "Hainan",
                           "湖北省" = "Hubei",
                           "湖南省" = "Hunan",
                           "福建省" = "Fujian",
                           "重庆市" = "Chongqing",
                           "陕西省" = "Shanxi"
  ))


library(ggplot2)

ggplot(plot_data_clean, aes(x = total_cases_sum, y = val_case, color = province)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  geom_text(aes(label = province), vjust = -0.8, size = 3, show.legend = FALSE) +  # 显示省份名称
  labs(
    title = "Validation of the Model Using 2023 Surveillance \nReported Dengue Importation from Xishuangbanna",
    x = "Simulated Imported Cases",
    y = "Reported Cases",
    color = "Province"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),   # legend标题缩小
    legend.text = element_text(size = 9),                   # legend文字缩小
    legend.key.size = unit(0.4, "cm"),                      # legend小点缩小
    legend.box.margin = margin(0, 0, 0, 0),                 # 减小legend边距
    legend.spacing.y = unit(0.2, "cm")                      # legend项之间距离缩小
  )


  