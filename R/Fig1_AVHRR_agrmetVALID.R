rm(list = ls())
source('R/mainfunc/main.R', encoding = 'UTF-8')

load("data/phenology_list20170105.rda")
load("../TP_phenoAgr/agrmet25_TP.rda")
load("../data/00basemap.rda")

## FUNCTIONS
get_validInfo <- function(pheno, agr.list, Id_clip){
  ## 计算RMSE and MAE
  EOS <- lapply(pheno, `[`, i = Id_clip, j = )

  RE <- list()
  # 14 methods loop
  for (i in seq_along(EOS)){
    x <- EOS[[i]] %>% t %>% set_colnames(NULL) %>% set_rownames(NULL) %>% data.frame(.)
    Id_year <- lapply(agr.list, function(x) x$year - 1981)
    REj <- list()
    for (j in seq_along(agr.list)){
      agri <- agr.list[[j]]
      # Id_year[[j]]
      REj[[j]] <- (x[Id_year[[j]], j] - agri$doy) %>% {.[!is.na(.)]}
    }
    RE[[i]] <- setNames(REj, names(agr.list))
  }
  RE %<>% set_names(names(EOS))

  # 计算RMSE, mean bias, ABE(mean absolute error)
  # NA values in x have been removed, stations order is correct
  info <- lapply(RE, function(xlist)
    ldply(xlist, function(x)
      data.frame(
        bias = mean(x),
        MAE = mean(abs(x)),
        RMSE = sqrt(sum(x ^ 2) / length(x))
      ), .id = NULL))
  infoValid <- lapply(1:3, function(i) do.call(cbind.data.frame, lapply(info, `[[`, i))) %>% 
    set_names(colnames(info[[1]]))
  return(infoValid)
}

show_valid <- function(validInfo, trim = FALSE){
  if (trim) {
    #delete Maturity, SD, Dormancy and RD
    xinterp <- c(4.5, 7.5)
    validInfo %<>% lapply(`[`, , j = c(1:4, seq(5, 16, 2)))
  } else {
    xinterp <- c(4.5, 10.5)
  } 
  suppressMessages(df <- melt(validInfo) %>% set_names(c("var", "value", "indice")))
  # linedf <- data.frame(yintercept = c(15, -15, 15, 15), variable = )
  line_df <- data.frame(indice = c("bias", "bias","MAE", "RMSE"), val = c(15, -15, 15, 15))
  p <- 
    ggplot(df, aes(x = var, y=value, color = var)) +
    geom_hline(aes(yintercept = val), line_df, color = "red", linetype=2, size = 0.7) + 
    geom_hline(yintercept = 30, color = "red", linetype = 1, size = 0.7) +
    # geom_hline(yintercept = 0, color = "black", linetype = 1, size = 0.7) +
    geom_boxplot(outlier.size=2) + geom_jitter(width = 0.15, size = 1.7, alpha = 1) +
    geom_vline(xintercept = xinterp, color = "grey40", linetype=1, size = 0.7) + 
    facet_grid(indice~., scales = "free") +
    # geom_hline(data=data.frame(indice = "bias", stringsAsFactors = F), yintercept = -15, color = "red", linetype=2) +
    theme(axis.text = element_text(size = 14, family = "Times"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          # panel.spacing = unit(0.4, "lines"),
          plot.margin = margin(4, 3, 4, 3),
          strip.text = element_text(face="bold", size = 14, family = "Times")) +
    guides(color = FALSE)
  return(p)
}
## -------------------------------------------------

validInfo.sos <- get_validInfo(SOS_list, agr_SOS$agr, Id_clip)
validInfo.eos <- get_validInfo(EOS_list, agr_EOS$agr, Id_clip)

# df <- suppressMessages(melt(list(SOS = validInfo.sos, EOS = validInfo.eos))) %>% 
#   set_names(c("var", "value", "indice", "pheno"))

p1 <- show_valid(validInfo.sos, trim = T)
p2 <- show_valid(validInfo.eos, trim = T)
p <- gridExtra::grid.arrange(p1, p2, ncol = 2, padding = unit(0.2, "line"))

ggsave("Fig2_agrmet_Valid.pdf", p, width = 12, height = 6.86)
# ggsave("c.png", p, width = 12, height = 6.86)
# used ggsave for temporary
# CairoPNG(file = "b.png", width = 12, height = 6.86, units = "in", dpi = 250, pointsize = 24)
cairo_pdf(file = "Fig2_agrmet_Valid.pdf", width = 12, height = 6.86)
grid::grid.draw(p)
grid::grid.text("(a)", x=unit(0.066, "npc"), y=unit(0.97, "npc"), 
                gp=grid::gpar(fontfamily = "Times", fontsize=16, fontface = "bold"))
grid::grid.text("(b)", x=unit(0.57, "npc"), y=unit(0.97, "npc"), 
                gp=grid::gpar(fontfamily = "Times", fontsize=16, fontface = "bold"))
dev.off()


## ----------------------------------------
## fenxiban
validInfo_trim.eos <- llply(validInfo.eos, `[`, i=, j = c(1:4, seq(5, 16, 2)))
validInfo_trim.sos <- llply(validInfo.sos, `[`, i=, j = c(1:4, seq(5, 16, 2)))
p1 <- show_valid(validInfo_trim.sos, T)
p2 <- show_valid(validInfo_trim.eos, T)
p <- gridExtra::grid.arrange(p1, p2, ncol = 2, padding = unit(0.2, "line"))





stations.sp@data <- EOS.mean[Id_clip, ] - agr_mean

stations.sp@data <- info[[1]]

brks <- c(-50, -30,  -15, 0, 15, 30, 120)

cex = seq(1.3, 0.8, length.out = 3) %>% c(., rev(.))
pch = rep(c(25, 24), c(3,3))
cols <- cols.fun(brks)
cols <- c("firebrick1" ,"orange3", "darkgoldenrod2", 
          brewer.pal(9, "YlGnBu")[c(6, 7)], "green4")

# stations.sp$agr <-  agr[, mean(doys, na.rm = T), by = stationId][, V1]
# x_mean <- subset(df, var == "mean")
# x_sd <- subset(df, var == "sd")


brks <- c(175, seq(210, 260, 15), 300)
x <- subset(df, var == "mean")
x$val %<>% cut(brks)

# shape = as.factor(作物名称)
ggplot(x, aes(long, lat, color = val, size = val)) + 
  geom_point() + 
  geom_path(data = shp, aes(x = long, y = lat, group = group, color=NULL, size=NULL, shape=NULL)) + 
  facet_grid(var~.)

# data$stationId %<>% paste0()
data$stationId %<>% as.numeric()
data <- left_join(data, stations[, c("stationId", "VId")], by = "stationId")
ggplot(data, aes(x = Year, y = doys, color = 作物名称, shape = 作物名称)) + 
  geom_line(size = 1) + geom_point(size = 3) + 
  facet_wrap(~VId, ncol = 5, nrow = 5, scales = "free_y") + 
  scale_x_continuous(breaks = c(1990, 2000, 2010)) + 
  theme(legend.position = c(0.5, 0.04)) + 
  guides(col = guide_legend(nrow = 1, byrow = TRUE)) + 
  theme(legend.position="bottom", 
        # legend.position = c(0.22, 0.02), 
        legend.justification = c(0, 0))

a <- melt(infoValid)
ggplot(a, aes(x = variable, y=value, color = variable)) + 
  geom_boxplot() + facet_grid(.~L1) 