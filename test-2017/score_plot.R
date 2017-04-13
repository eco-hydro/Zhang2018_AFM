rm(list = ls())
library(openxlsx)
source('R/mainfunc/main.R', encoding = 'UTF-8')

df <- read.xlsx("data/成绩.xlsx") %>% subset(totl >= 320 & profession >= 180 & policy >= 55 & english >= 55)
df_opt <- subset(df, code == 105100)



df$direction %<>% gsub("\\(|全日制|\\)|内科学|外科学|（|）", "", .)
df$studentId %<>% as.character()
df$majorType %<>% as.numeric()

df_list <- split(df, df$direction)
writelist_ToXlsx(df_list, "scores.xlsx")

df %>% {
  . %$% which(totl > 320)
}

unique(df[, c("code", "major")])

substr(as.character(df$code), 4, 6) %>% unique()


Id <- 1:400
values <- numeric(length(Id))
for (i in seq_along(Id)){
  x <- df[df$code ==105100 & df$totl > Id[i], ]
  # x <- subset(df, (code == 105100) & (totl > Id[i])) #%>% {split(., .$direction)}
  values[i] <- nrow(x)
}

plot(Id, values)
abline(v = 339, col = "red")



x <- subset(df, code ==105100 & totl >= 320 & policy >= 50 & english >= 50 & profession >= 180)

x <- df[df$code ==105100 & df$totl > 300, ]

x$direction %<>% as.factor()

ggplot(x, aes(x = direction, y = totl, fill = direction)) + 
  geom_boxplot() + geom_jitter(width = 0.1) + 
  ylab("scores") + 
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  guides(fill = "none")