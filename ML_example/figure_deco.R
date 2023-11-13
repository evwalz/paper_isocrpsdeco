# load msc, dsc, crps and unc for all methods
library(ggplot2)
library(ggrepel)
library(geomtextpath)

set_dir = dirname(rstudioapi::getSourceEditorContext()$path)
data_directory = 'yacht' # bostonHousing, concrete, energy, yacht

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


hhjust <- 0.36
dim_out <- 6
if (data_directory == 'bostonHousing'){
  data_title <- 'Boston'
  dim_out <- 6
  hhjust <- 0.43
} else if (data_directory == 'yacht'){
  data_title <- 'Yacht'
  dim_out <- 7
  hhjust <- 0.53
} else if (data_directory == 'energy'){
  data_title <- 'Energy'
  dim_out <- 7
  hhjust <- 0.6
} else if (data_directory == 'concrete'){
  data_title <- 'Concrete'
  dim_out <- 6
  hhjust <- 0.55
}else {
  data_title <- data_directory
}


methods = c('laplace','mc_dropout', 'single_gaussian', 'easyuq','easyuq_smooth', 'cp', 'cp_smooth')
methods_name <- c( 'Laplace',  'MC Dropout', 'Single Gaussian',  'EasyUQ', 'Smooth EasyUQ','CP', 'Smooth CP')


cols <- gg_color_hue(length(methods))

methods_msc <- rep(0, length(methods))
methods_dsc <- rep(0, length(methods))
methods_Nmsc <- rep(0, length(methods))
methods_Ndsc <- rep(0, length(methods))
methods_crps <-rep(0, length(methods))
methods_unc <-rep(0, length(methods))

k <- 1
for (method in methods) {
  if (method == 'cp_smooth') {
    deco_path <- paste(set_dir, '/deco_data/cp/UCI_datasets/',  data_directory , '/deco_results_smooth/', sep = '')
  } else if (method == 'easyuq_smooth') {
    deco_path <- paste(set_dir, '/deco_data/easyuq/UCI_datasets/',  data_directory , '/deco_results_smooth/', sep = '')
  }  else {
    deco_path <- paste(set_dir, '/deco_data/', method, '/UCI_datasets/',  data_directory , '/deco_results/', sep = '')
  }
  methods_msc[k] <- mean(read.table(paste(deco_path, "msc.txt", sep=""))$V1)
  methods_dsc[k]  <- mean(read.table(paste(deco_path, "dsc.txt", sep=""))$V1)
  methods_unc[k]  <- mean(read.table(paste(deco_path, "unc.txt", sep=""))$V1)
  methods_crps[k]  <- mean(read.table(paste(deco_path, "crps.txt", sep=""))$V1)
  k <- k+1
}

unc_str <- methods_unc[1]
df_scores <- data.frame(models = methods_name, msc  = methods_msc/unc_str, dsc = methods_dsc/unc_str)

unc_val <- round(methods_unc[1], 2)


if (data_directory == 'concrete'){
  crps_isolines <- seq(2, 3.5, 0.1)
} else if (data_directory == 'bostonHousing'){
  crps_isolines <- c(1.0, 1.1,1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3)
} else if (data_directory == 'energy'){
  crps_isolines <- seq(0.19, 0.43, 0.03)
}else if (data_directory == 'yacht'){
  crps_isolines <- seq(0.35, 0.59, 0.03)
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

intercept <- c(intercepts[1]-3*diff_small,intercepts[1]-2*diff_small, intercepts[1]-diff_small ,intercepts , intercepts[len_1]+diff_big, intercepts[len_1]+2*diff_big, intercepts[len_1]+3*diff_big, intercepts[len_1]+4*diff_big)
iso = data.frame(intercept = intercept , slope = rep(1, length(intercept)))


p <- ggplot(df_scores) +
  geom_labelabline(
    data = iso2, aes(intercept = intercept, slope = slope, label = crps_iso), color = "gray50",
    hjust = hhjust, size = 7 * 0.36, text_only = TRUE, boxcolour = NA, straight = TRUE
  ) + geom_abline(
    data = iso, aes(intercept = intercept, slope = slope), color = "lightgray", alpha = 0.5,
    size = 0.5 ) +geom_point(aes(x = msc, y = dsc, color = models), size = 1.5)+
  geom_text_repel(aes(x = msc, y = dsc, label = models, color = models),
                  max.overlaps = NA, size = 11 * 0.36, nudge_x = 0,seed = 4) +
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
  )+
  scale_color_manual(values = cols)


p + geom_label(data = NULL, x = -Inf, y = Inf, vjust = 1,hjust = 0,label =paste("UNC =" , unc_val)) + ggtitle(data_title)

ggsave(paste(set_dir, "/figures/",data_directory,"_deco.pdf", sep=""), width = 120, height = 120, unit = "mm", device = "pdf")
