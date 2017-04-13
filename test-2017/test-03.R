library(psych)
library(magrittr)

windowsFonts(Times = windowsFont("Times New Roman"), ST = windowsFont("宋体"), 
             YH = windowsFont("Microsoft Yahei"), Arial = windowsFont("Arial"))


dd <- read.table("124.txt", header = T)
rownames(dd) <- c(1:12, "year")
dd <- dd[, -1]

ww <- cor.plot(dd, numbers = TRUE, show.legend = TRUE, 
               zlim = c(-0.4, 0.4), labels = NULL, family = "Arial", 
               cex = 1.2, cex.axis = 1.3, 
               font.axis = 2)
# , labels = "month", 
# n.legend = 10, main = "NAO")