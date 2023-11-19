library(dplyr)
library(ggplot2)
library(ggrepel)
library(geomtextpath)

set_dir = dirname(rstudioapi::getSourceEditorContext()$path)
load(paste(set_dir , "/deco_data/results.rda", sep = ''))

unique(Data$airport) # "lhr" "zrh" "bru" "fra"
city <- "fra"


sel_data <- function(x, city){
  SubData_MSC <- Data %>% filter(airport == city) %>% filter(horizon == x) %>% filter(category == "mcb")
  SubData_DSC <- Data %>% filter(airport == city) %>% filter(horizon == x) %>% filter(category == "dsc")
  SubData_UNC <- Data %>% filter(airport == city) %>% filter(horizon == x) %>% filter(category == "unc")
  SubData_CRPS <- Data %>% filter(airport == city) %>% filter(horizon == x) %>% filter(category == "crps")
  return(list(SubData_UNC$value, SubData_MSC$value, SubData_DSC$value, SubData_CRPS$value))
}

model_names <- c("HCLR" ,  "EMOS"  , "BMA"  ,  "IDR_cw" ,"IDR_st" ,"ENS")

if (city == "fra"){
  city_name = 'Frankfurt'
} else if (city == "bru"){
  city_name = 'Brussels'
}else if (city == "lhr"){
  city_name = 'London'
}else if (city == "zrh"){
  city_name = 'Zurich'
}

deco1 <- sel_data(1, city)
deco2 <-sel_data(2, city)
deco3 <-sel_data(3, city)
deco4 <-sel_data(4, city)
deco5 <-sel_data(5, city)
df_scores <-  data.frame(models = rep(model_names, 5), msc  = c(deco1[[2]] / deco1[[1]][1], deco2[[2]] / deco2[[1]][1], deco3[[2]] / deco3[[1]][1], deco4[[2]] / deco4[[1]][1], deco5[[2]] / deco5[[1]][1]) ,
                         dsc = c(deco1[[3]] / deco1[[1]][1], deco2[[3]] / deco2[[1]][1], deco3[[3]] / deco3[[1]][1], deco4[[3]] / deco4[[1]][1], deco5[[3]] / deco5[[1]][1]))



unc_str <- round(deco1[[1]][1], 2)
df_scores$Lag <- as.character(sort(rep(c(1, 2, 3, 4, 5), 6)))
df_scores$Row <- as.character(c(1:30))

if (city == "bru") {
  crps_isolines <- c(0.86, 0.96, 1.06, 1.16, 1.26, 1.36, 1.46, 1.56, 1.66)
  hjust <- 0.55
} else if (city == "fra") {
  crps_isolines <- c(0.59, 0.65, 0.71, 0.77, 0.83, 0.89, 0.95, 1.01, 1.07)
  hjust <- 0.45
} else if (city == "lhr") {
  crps_isolines <- c(0.57,0.65,0.73,  0.81, 0.89, 0.97, 1.05, 1.13, 1.20)
  hjust <- 0.4
} else if (city == "zrh") {
  crps_isolines <- c(0.91 , 1.02, 1.13, 1.24, 1.35, 1.46, 1.57, 1.68, 1.79)
  hjust <- 0.41
}

intercepts <-(unc_str - crps_isolines) /unc_str
len_1 <- length(intercepts)
diff_small <- intercepts[2] - intercepts[1]
diff_big <- intercepts[(length(intercepts))] - intercepts[(length(intercepts)-1)]
iso = data.frame(intercept = intercepts , slope = rep(1, length(intercepts)))
len_1 <- length(iso$intercept)

crps_iso <- round(crps_isolines, 2)
iso2 <- data.frame(intercept = iso$intercept , slope = iso$slope, crps_iso = crps_iso)

diff_crps <- crps_iso[2] - crps_iso[1]
diff_crps2 <- crps_iso[len_1] - crps_iso[(len_1 - 1)]

p <- ggplot(df_scores) +
  geom_labelabline(
    data = iso2, aes(intercept = intercept, slope = slope, label = crps_iso), color = "gray50",
    hjust = hjust, size = 7 * 0.36, text_only = TRUE, boxcolour = NA, straight = TRUE
  ) + geom_abline(
    data = iso, aes(intercept = intercept, slope = slope), color = "lightgray", alpha = 0.5,
    size = 0.5 ) +geom_point(aes(x = msc, y = dsc, color = Lag), size = 1.5)+
  geom_text_repel(aes(x = msc, y = dsc, label = models, color = Lag),show.legend=FALSE,
                  max.overlaps = NA, size = 8 * 0.36, nudge_x = 0,seed = 4) +
  xlab("MCB") +
  ylab("DSC") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin=grid::unit(c(1,0,1,0), "mm")
  )

p <-  p +  geom_label(data = NULL,x = -Inf, y = Inf, vjust = 1,hjust = 0, label =paste("UNC =" , unc_str)) +
  ggtitle(city_name) + scale_fill_discrete(name = "Prediction horizon")

p

ggsave(paste(set_dir , "/figures/",city,"_deco.pdf", sep=""), width = 120, height = 120, unit = "mm", device = "pdf")
