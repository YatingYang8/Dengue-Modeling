library(dplyr)
library(reshape2)
library(lubridate)

#read Dengue data
dd <- read.csv("./data/Dengue case data 2023.csv",fileEncoding = "GB18030")

names(dd)

dd <- dd %>%
  dplyr::select("发病日期","现住地址国标")

colnames(dd) <- c("date","areaid")

dd$date <- as.Date(dd$date)

dd <- dd %>%
  mutate(
    year = format(date, "%Y"),
    month = format(date, "%m")
  )

dd$cases <- 1

dd_sum <- dd %>%
  group_by(date) %>%
  summarise(cases=sum(cases,na.rm = TRUE))


date_range <- seq(as.Date("2023-01-01"), as.Date("2023-12-31"), by="day")

total_index <- data.frame(date = date_range)

dd1 <- left_join(total_index,dd_sum, by = c("date") )

dd1$cases[is.na(dd1$cases)] <- 0

dd1 <- dd1 %>%
  mutate(
    date = as.Date(date, format="%Y-%m-%d"),  
    year = format(date, "%Y"),
    month = format(date, "%m"),
    yearmonth = format(date, "%Y-%m"),
  )

dd1$year <- as.integer(dd1$year)

#save(dd1, file = "./data/dengue_caseALL.RData")
#write.csv(dd1, file = "./data/dengue_caseALL.csv",fileEncoding = "GB18030")
