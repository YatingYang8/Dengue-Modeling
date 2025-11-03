library(dplyr); library(raster); library(rgdal); library(Matrix)
library(stringr); library(ggplot2); library(lubridate);library(cowplot);library(colorspace)
library(magrittr);  library(spdep);library(rgeos);library(sf);library(sp);library(ggspatial)
library(readxl);library(reshape2)
library(tidyverse)
library(caret)     # 用于预处理标准化
library(corrplot)  # 可视化相关系数矩阵
library(car)       # 用于 VIF 计算


setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue")

load(file = "./data/ML_origindata.RData")

summary(ML_origindata)

###查看变量的分布，看是否需要做data transforming

############### Outcome
ggplot(ML_origindata, aes(x = total_cases)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  #geom_density(color = "red", size = 1, alpha = 0.3) +
  labs(title = "Distribution of total_cases", x = "Total Cases", y = "Count") +
  theme_minimal()

ggplot(ML_origindata, aes(x = log(total_cases))) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_density(color = "red",alpha = 0.1, size = 1) +
  labs(title = "Distribution of log total_cases", x = "Total Cases", y = "Count") +
  theme_minimal()

ggplot(ML_origindata, aes(x = log(log(total_cases)))) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_density(color = "red",alpha = 0.1, size = 1) +
  labs(title = "Distribution of total_cases", x = "Total Cases", y = "Count") +
  theme_minimal()

qqnorm(ML_origindata$total_cases, main = "Q-Q Plot of total_cases")
qqline(ML_origindata$total_cases, col = "red")

ggplot(ML_origindata, aes(x = time)) +
  geom_histogram(bins = 30, fill = "orange", color = "black", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  labs(title = "Distribution of time", x = "Time", y = "Count") +
  theme_minimal()

ggplot(ML_origindata, aes(x = log(time))) +
  geom_histogram(bins = 30, fill = "orange", color = "black", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  labs(title = "Distribution of log time", x = "Time", y = "Count") +
  theme_minimal()


###数据高度偏态，且存在明显的离散等级，应该转换为分类变量
summary(ML_origindata$total_cases)


ml_data <- ML_origindata %>%
  mutate(risk_level = case_when(
    total_cases <= 2.44 ~ 0,
    total_cases > 2.44 ~ 1
  ))


table(ml_data$risk_level)

head(ml_data)


###################### data cleaning and EDA ###################################

colSums(is.na(ml_data))

data <- ml_data %>%
  dplyr::select(-passenger_water, -passenger_air, -passenger_railway)

colSums(is.na(data))

# 2. 数据清洗
# 删除缺失值
data <- na.omit(data)
table(data$risk_level)

head(data)

# 选择数值型变量进行建模
num_vars <- data %>% 
  dplyr::select(where(is.numeric))%>% 
  dplyr::select(-risk_level, -total_cases, -time)  # 不包含因变量

# 2. 标准化处理
scaled_data <- scale(num_vars)

cor_matrix <- cor(scaled_data, use = "complete.obs")

corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.6)

# 5. 找出高度相关变量对（|r| > 0.85）
high_cor_pairs <- findCorrelation(cor_matrix, cutoff = 0.7, names = TRUE, verbose = TRUE)
high_cor_pairs

# 6. 删除冗余变量（初步）
data_reduced <- data %>%
  dplyr::select(-all_of(high_cor_pairs))

# 7. 再做 VIF 检查
# 建一个临时线性模型（只用于诊断 VIF）
model_formula <- as.formula(paste("risk_level ~", paste(colnames(data %>% dplyr::select(where(is.numeric), -risk_level,-total_cases,-time,-cityid)), collapse = "+")))
lm_model <- lm(model_formula, data = data)

vif_values <- vif(lm_model)
print(vif_values)

high_vif <- vif_values[vif_values > 10]
print(high_vif)

high_cor_pairs <- c("doctors","hospital_beds","pop_total","gdp_total","employee_avg",  
                     "school_high","urbanization" )
data_reduced <- data %>% dplyr::select(-all_of(high_cor_pairs))

#检查删除后的vif和correlation
model_formula <- as.formula(paste("risk_level ~", paste(colnames(data_reduced %>% dplyr::select(where(is.numeric), -risk_level,-total_cases,-time,-cityid)), collapse = "+")))
lm_model <- lm(model_formula, data = data_reduced)

vif_values <- vif(lm_model)
print(vif_values)



############## Correlation with outcome #########################################
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(corrplot)
# 
# # 保留数值型变量（去除字符型，如 city, province 等）
# numeric_data <- data %>%
#   dplyr::select(where(is.numeric))
# 
# # 计算与 risk_level 的 Pearson 相关性
# cor_matrix <- cor(numeric_data, use = "complete.obs", method = "pearson")
# 
# # 提取与 risk_level 的相关性（排除自身）
# cor_target <- tibble(
#   variable = rownames(cor_matrix),
#   correlation = cor_matrix[, "risk_level"]
# ) %>%
#   filter(variable != "risk_level") %>%
#   arrange(desc(abs(correlation)))
# 
# # 查看前 10 个相关性最大变量
# print(dplyr::slice(cor_target, 1:15))
# 
# var_exclude <- cor_target %>%
#   filter(abs(correlation) < 0.2)
# 
# 保留变量列表
exclude_vars <- c("waste_treatment","gdp_pct_secondary")
exclude_vars
data_reduced <- data_reduced %>% dplyr::select(-any_of(exclude_vars))

head(data_reduced)

# 4. 数据预处理
library(caret)
set.seed(123)

# 将risk_level转为因子（如果是分类任务）
data_reduced <- data_reduced %>%
  dplyr::select(-any_of(c("total_cases","time","cityid", "city", "province")))
data_reduced$risk_level <- as.factor(data_reduced$risk_level)


# 把 high 编码为 1，low 编码为 0（或反之，根据你的实际含义）
table(data_reduced$risk_level)

# data_reduced <- data_reduced %>%
#   mutate(risk_level = ifelse(risk_level == "3", 1, 0)) 

data_reduced$risk_level <- as.factor(data_reduced$risk_level)

# 划分训练/测试集
trainIndex <- createDataPartition(data_reduced$risk_level, p = 0.7, list = FALSE)
train_data <- data_reduced[trainIndex, ]
test_data  <- data_reduced[-trainIndex, ]


table(train_data$risk_level)



# 加载包
library(randomForest)
library(xgboost)
library(e1071)
library(nnet)
library(fastshap)
library(pROC)
library(dplyr)
library(ggplot2)
library(tidyr)
library(caret)
library(fastshap)
library(shapviz)

# 确保分类变量为 factor
train_data$risk_level <- as.factor(train_data$risk_level)
test_data$risk_level <- as.factor(test_data$risk_level)
names(train_data)
names(test_data)

# 特征矩阵
X_train <- train_data %>% dplyr::select(-risk_level)
X_test <- test_data %>% dplyr::select(-risk_level)
y_train <- train_data$risk_level
y_test <- test_data$risk_level

# RF
rf_model <- randomForest(risk_level ~ ., data = train_data, ntree = 500)
rf_prob <- predict(rf_model, test_data, type = "prob")[, 2]
rf_auc <- auc(roc(y_test, rf_prob))

shap_rf <- explain(rf_model, X = X_train, pred_wrapper = function(object, newdata) {
  predict(object, newdata, type = "prob")[, 2]
})

sv_rf <- shapviz(shap_rf, X=X_train)

sv_importance(sv_rf, kind = "bar")       # 平均SHAP值重要性条形图
sv_importance(sv_rf, kind = "beeswarm")  # beeswarm图，显示变量影响分布


# XGBoost
# 转换为 xgboost 所需格式
xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = as.numeric(y_train) - 1)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = as.numeric(y_test) - 1)

xgb_model <- xgboost(data = xgb_train, objective = "binary:logistic", nrounds = 100, verbose = 0)
xgb_prob <- predict(xgb_model, xgb_test)
xgb_auc <- auc(roc(y_test, xgb_prob))

shap_xgb <- explain(xgb_model, X = X_train, pred_wrapper = function(object, newdata) {
  predict(object, as.matrix(newdata))
})

sv_xgb <- shapviz(shap_xgb, X=X_train)

sv_importance(sv_xgb, kind = "bar")       # 平均SHAP值重要性条形图
sv_importance(sv_xgb, kind = "beeswarm")  # beeswarm图，显示变量影响分布



# SVM
svm_model <- svm(risk_level ~ ., data = train_data, probability = TRUE)
svm_prob <- attr(predict(svm_model, test_data, probability = TRUE), "probabilities")[, 2]
svm_auc <- auc(roc(y_test, svm_prob))

shap_svm <- explain(svm_model, X = X_train, pred_wrapper = function(object, newdata) {
  attr(predict(object, newdata, probability = TRUE), "probabilities")[, 2]
})

sv_svm <- shapviz(shap_svm, X=X_train)

sv_importance(sv_svm, kind = "bar")       # 平均SHAP值重要性条形图
sv_importance(sv_svm, kind = "beeswarm")  # beeswarm图，显示变量影响分布


# NN
# 需要把响应变量转成 dummy
y_train_num <- class.ind(y_train)
nn_model <- nnet(x = as.matrix(X_train), y = y_train_num, size = 5, softmax = TRUE, maxit = 500, trace = FALSE)
nn_prob <- predict(nn_model, newdata = as.matrix(X_test))[, 2]
nn_auc <- auc(roc(y_test, nn_prob))

shap_nn <- explain(nn_model, X = X_train, pred_wrapper = function(object, newdata) {
  predict(object, as.matrix(newdata))[, 2]
})

sv_nn <- shapviz(shap_nn, X=X_train)

sv_importance(sv_nn, kind = "bar")       # 平均SHAP值重要性条形图
sv_importance(sv_nn, kind = "beeswarm")  # beeswarm图，显示变量影响分布



#计算AUC
auc_df <- data.frame(
  Model = c("Random Forest", "XGBoost", "SVM", "Neural Network"),
  AUC = c(rf_auc, xgb_auc, svm_auc, nn_auc)
)

ggplot(auc_df, aes(x = Model, y = AUC, fill = Model)) +
  geom_bar(stat = "identity") +
  ylim(0, 1) +
  theme_minimal() +
  ggtitle("AUC Comparison Across Models") +
  theme(legend.position = "none")



# 创建 ROC 曲线对象
roc_rf <- roc(y_test, rf_prob, quiet = TRUE)
roc_xgb <- roc(y_test, xgb_prob, quiet = TRUE)
roc_svm <- roc(y_test, svm_prob, quiet = TRUE)
roc_nn  <- roc(y_test, nn_prob, quiet = TRUE)

# 使用 ggplot2 绘图：整合数据
library(tibble)
library(ggplot2)

roc_df <- bind_rows(
  tibble(
    FPR = 1 - roc_rf$specificities,
    TPR = roc_rf$sensitivities,
    Model = "Random Forest",
    AUC = as.numeric(auc(roc_rf))
  ),
  tibble(
    FPR = 1 - roc_xgb$specificities,
    TPR = roc_xgb$sensitivities,
    Model = "XGBoost",
    AUC = as.numeric(auc(roc_xgb))
  ),
  tibble(
    FPR = 1 - roc_svm$specificities,
    TPR = roc_svm$sensitivities,
    Model = "SVM",
    AUC = as.numeric(auc(roc_svm))
  ),
  tibble(
    FPR = 1 - roc_nn$specificities,
    TPR = roc_nn$sensitivities,
    Model = "Neural Network",
    AUC = as.numeric(auc(roc_nn))
  )
)


# 用 AUC 值作为图例标签
roc_df <- roc_df %>%
  group_by(Model) %>%
  mutate(AUC_Label = paste0(Model, " (AUC = ", round(first(AUC), 3), ")"))

# 绘图
ggplot(roc_df, aes(x = FPR, y = TPR, color = AUC_Label)) +
  geom_line(size = 0.8) +
  geom_abline(linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curve Comparison Across Models",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")





color_palette <- c(
  "Random Forest (AUC = 0.993)" = "#1b9e77",
  "XGBoost (AUC = 0.979)"       = "#d95f02",
  "SVM (AUC = 0.995)"           = "#7570b3",
  "Neural Network (AUC = 0.825)"= "#e7298a"
)

linetype_palette <- c(
  "Random Forest (AUC = 0.993)" = "solid",
  "XGBoost (AUC = 0.979)"       = "dashed",
  "SVM (AUC = 0.995)"           = "dotdash",
  "Neural Network (AUC = 0.825)"= "twodash"
)



# 如果需要动态提取 AUC 标签，请保留这部分逻辑
roc_df <- roc_df %>%
  group_by(Model) %>%
  mutate(AUC_Label = paste0(Model, " (AUC = ", round(first(AUC), 3), ")")) %>%
  ungroup()
unique(roc_df$AUC_Label)
# 强制设定图例顺序和一致性（避免绘图时颜色/线型混乱）
roc_df$AUC_Label <- factor(roc_df$AUC_Label,
                           levels = c("Random Forest (AUC = 0.993)",
                                      "XGBoost (AUC = 0.979)",
                                      "SVM (AUC = 0.995)",
                                      "Neural Network (AUC = 0.825)"))

# 绘图
p3 <- ggplot(roc_df, aes(x = FPR, y = TPR, color = AUC_Label, linetype = AUC_Label)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dotted", color = "gray40", size = 0.8) +
  labs(
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = NULL,
    linetype = NULL
  ) +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(values = linetype_palette) +
  theme(
    panel.background = element_rect(fill = "transparent", color = "black", size = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  coord_equal()  # 保持坐标轴1:1比例

p3


p1 <- sv_importance(sv_rf, kind = "bar")       # 平均SHAP值重要性条形图
p1

p2 <- sv_importance(sv_rf, kind = "beeswarm")  # beeswarm图，显示变量影响分布
p2

p = gridExtra::grid.arrange(p1,p2,p3, ncol = 3, nrow = 1, widths=c(1,1,1))

