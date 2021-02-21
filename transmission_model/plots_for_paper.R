library(here)
library(tidyverse)
library(sf)
library(viridis)
library(ggplot2)
library(stringr)
library(data.table)
library(ggforce)
library(cowplot)
source(here("transmission_model/geom-stepribbon.r"))
final_data <- readRDS(here("data/transmission_output_joint.rds"))

### ratio figure main
ggplot(data = final_data) + 
  theme(legend.position = "right") +
  geom_sf(aes(fill = Ratio)) +
  scale_fill_gradientn(colors=leaflet::colorBin(grDevices::colorRamp(c(hsv(0.06, 0.0, 0.9),
                                                                       hsv(0.03, 0.0, 0.8),
                                                                       hsv(0.01, 0.4, 0.8),
                                                                       hsv(0.99, 1.0, 0.7)),
                                                                     bias=1.2),
                                                domain = c(0,1))(seq(0, 1, 0.2))) +
  coord_sf(crs = sf::st_crs("EPSG:3035")) +
  facet_wrap(. ~ epiweek, ncol = 4) +
  theme_void() + 
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.margin=grid::unit(c(0,0,0,0), "mm"), 
    strip.text.x = element_text(size = 20),
    legend.text = element_text(size=20),
    legend.title = element_text(size=20),
    legend.key.size = unit(1, 'cm'))
  
ggsave(here("figures/rt_estimates_ratio_map.pdf"),width=12,height=10, dpi=600)

# Scatter plot for RS(-) and RS(+)
sgss <- readRDS(here("data/sgss_stp_new_43_56_weeks.rds"))
sgss.to.areas <- fread(here("data/ltla-nhser-stp.csv")) %>% select("stp_name","nhser_name")
sgss.to.areas <- sgss.to.areas[, list(nhser_name = nhser_name[1]),by=stp_name]
sgss <- left_join(sgss,sgss.to.areas,by=c("area"="stp_name"))
stp_rt <- left_join(final_data %>%
                     mutate(epiweek = as.double(str_remove_all(final_data$epiweek, "Week " ))) %>%
                     rename(area = name)
                   ,sgss,by=c("epiweek","area"))

lockdown_start <- (45 + 4/7)
lockdown_end <- (49 + 2/7)

stp_rt$freq <- 100*stp_rt$sgss_s_negative_corrected_adj1 / (stp_rt$sgss_s_negative_corrected_adj1+stp_rt$sgss_s_positive_corrected_adj1)
stp_rt <- filter(stp_rt,epiweek >= 45 & epiweek <= 55)
ggplot(stp_rt,aes(x=`R(S+)`,y=`R(S-)`)) +
  geom_point(aes(size=freq,colour=nhser_name,shape=as.factor(epiweek)),alpha=.9) +
  geom_abline(slope=1) +
  geom_vline(xintercept=1,col="gray") +
  scale_x_continuous(limits=c(0,4),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,4),expand=c(0,0)) +
  geom_hline(yintercept=1,col="gray") +
  theme_bw() +
  coord_fixed() +
labs(size="% novel variant", colour="",shape="epiweek") +
scale_shape_manual(values = c(0:2,5:7,9:14)) +
theme(axis.text.x = element_text( angle = 45,
                                  hjust = 1),
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24)
    ) +
geom_smooth(method='lm', se=FALSE,color = 'darkgrey', fullrange = TRUE)
ggsave(here("figures/figure_2e.pdf"),width=12,height=10, dpi=600)
# figure for plotting time varying adavantage
time_varying <- readRDS(here("data/time_varying_advantage.rds"))
ggplot(data=time_varying, aes(x=week, y=`50%`, ymin=`2.5%`, ymax=`97.5%`)) +
  theme_bw() +
  theme(axis.text.x = element_text( angle = 45,
                                    hjust = 1),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 24)
      ) +
  theme(legend.position = "right") +
  geom_step(size=2) + 
  geom_stepribbon(alpha=0.3) +
  xlab(as.expression("Week")) + 
  scale_x_continuous(breaks = (45:55)) + 
  ylab("Multiplicative advantage \n for VOC over non-VOC")
ggsave(here("figures/figure_2d.pdf"),width=12,height=10, dpi=600)

# observation vs prediction plot
facet_plot_data <-
  readRDS(here("data/obs_predicted_data.rds")) %>%
  mutate(facet_name = str_c(Area," (",SGTF,")"))
base <- ggplot(data = facet_plot_data, aes(x = Cases_week, y = Prediction, ymin = Pli, ymax = Pui)) + 
  geom_blank() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
      axis.text = element_text(size = 24),
      axis.title = element_text(size = 24)
      ) 
p <- list()
for(i in seq(1,7)){
  p[[i]] <- base + 
    geom_abline() +
    geom_pointrange(size=0.1, color="black", fill="white", shape=22) + 
    theme_light()+
    theme(strip.text.x = element_text(
      color = "black"),
      strip.background = element_blank()
    ) +  # guides(color = guide_legend(override.aes = list(size = c(2,1,1)))) +
    ylab("Prediction") +
    xlab("Cases per week") +
    # facet_grid_paginate(Area ~ SGTF, scales = "free", ncol = 2, nrow = 6)
    facet_wrap_paginate(. ~ facet_name, scales = "free", ncol = 2, nrow = 6, page = i)
}

pdf(here("figures/obs_pred_each_stp.pdf"))
invisible(lapply(p, print))
dev.off()

ggplot(data = facet_plot_data, aes(x = Cases_week, y = Prediction, ymin = Pli, ymax = Pui)) + 
  geom_abline() +
  geom_pointrange(size=0.1, color="black", fill="white", shape=22) + 
  ylab("Prediction") +
  xlab("Cases per week") +
  facet_wrap(. ~ SGTF, scales = "free", ncol = 2) +
  theme_light()+
  theme(axis.text.x = element_text( angle = 45,hjust = 1)) + 
  theme(strip.text.x = element_text(color = "black", size = 22)) +
  theme(strip.background = element_blank()) + 
  theme(axis.text = element_text(size = 22)) +
  theme(axis.title = element_text(size = 22))
ggsave(here("figures/ED6.pdf"),width=12,height=10, dpi=600)

# comparison plots for original and different tg
df_rt <- readRDS(here("data/df_rt.rds")) %>% mutate(epiweek = as.double(str_remove_all(epiweek,"Week ")))
stp_rt = left_join(df_rt ,sgss,by=c("epiweek","area"))
stp_rt$freq <- 100*stp_rt$sgss_s_negative_corrected_adj1 / (stp_rt$sgss_s_negative_corrected_adj1+stp_rt$sgss_s_positive_corrected_adj1)
stp_rt <- filter(stp_rt,epiweek >= 45 & epiweek <= 55)
stp_rt$epiweek <- as.factor(stp_rt$epiweek)
time_varying <- readRDS(here("data/time_varying_advantage_all.rds"))
plot_ratio_scatter <- function(data,title = 'Original' ,range=3, need_legend = FALSE){
  g <- 
    ggplot(data,aes(x=`R(S+)`,y=`R(S-)`/`R(S+)`)) +
    geom_point(aes(colour=nhser_name,shape=epiweek),size = 2, alpha=.9) +
    scale_shape_manual(values = c(0:2,5:7,9:14)) +
    # geom_abline(slope=1) +
    # geom_vline(xintercept=1,col="gray") +
    # geom_hline(yintercept=1,col="gray") +
    scale_x_continuous(limits=c(0,range),expand=c(0,0)) +
    scale_y_continuous(limits=c(0,range),expand=c(0,0)) +
    theme_bw() +
    coord_fixed()+
    theme(axis.text.x =element_text(angle = 45,hjust = 1)) +
    theme(axis.text = element_text(size = 16)) + 
    theme(axis.title = element_text(size = 16)) +
    theme(legend.text=element_text(size=12)) +
    theme(legend.title = element_text(size=12)) +
    theme(plot.title = element_text(size=20)) +
    theme(legend.position="bottom",legend.direction = "horizontal", legend.box = "horizontal") +
    geom_smooth(method='lm', se=FALSE,color = 'darkgrey', fullrange = TRUE) +
    ggtitle(title) +
    ylab("Ratio") 
  if(need_legend){
    return(get_legend(g))
  }
  return( g + theme(legend.position='none'))
}
g1 <- plot_ratio_scatter(stp_rt %>% filter(`Generation time` == 'original')) 
g2 <- plot_ratio_scatter(stp_rt %>% filter(`Generation time` == '10% reduce'), title = '10% reduce')
g3 <- plot_ratio_scatter(stp_rt %>% filter(`Generation time` == '25% reduce'), title = '25% reduce')
g4 <- plot_ratio_scatter(stp_rt %>% filter(`Generation time` == 'original'), need_legend = TRUE)

plot_time_varying_advantage <-  function(data,title = 'Original'){
  ggplot(data=data %>% filter(`Generation time`==title), aes(x=Epiweek, y=mean, ymin=li, ymax=ui)) +
    theme_bw() +
    theme(axis.text.x =element_text(angle = 45,hjust = 1)) +
    theme(axis.text = element_text(size = 16)) + 
    theme(axis.title = element_text(size = 16)) +
    theme(legend.text=element_text(size=12)) +
    theme(legend.title = element_text(size=12)) +
    theme(plot.title = element_text(size=18)) + 
    theme(legend.position = "right") +
    geom_step(size=2) + 
    coord_fixed(ratio = 6)+
    geom_stepribbon(alpha=0.3) +
    xlab(as.expression("Week")) + 
    scale_x_continuous(breaks = (45:55)) + 
    scale_y_continuous(limits=c(1,2.75),expand=c(0,0)) +
    ylab("Multiplicative advantage \n for VOC over non-VOC")
}

g5 <- plot_time_varying_advantage(time_varying)
g6 <- plot_time_varying_advantage(time_varying, title = '10% reduce')
g7 <- plot_time_varying_advantage(time_varying, title = '25% reduce')

plot_grid( plot_grid(g1,g2,g3, ncol = 3), plot_grid(g4), plot_grid(g5,g6,g7,ncol = 3),
           rel_heights = c(4, 1, 4), nrow = 3)
ggsave(here("figures/ED7.pdf"),width=12,height=10, dpi=600)
