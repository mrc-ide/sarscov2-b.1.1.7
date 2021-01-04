library(here)
library(tidyverse)
library(sf)
library(viridis)
library(stringr)

transmission_data <- readRDS(here("data/transmission_data.rds"))
base <- ggplot(data = transmission_data) + 
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text = ggplot2::element_text(size = 22),
    axis.title = ggplot2::element_text(size = 22)
  ) +
  ggplot2::theme(legend.position = "right") 

base +
  geom_sf(aes(fill = Difference)) +
  scale_fill_gradientn(colors = c("#2ec4b6", "#e9ecef", "#ff9f1c"),
                       limits=c(-1.2, 1.2),
                       breaks=c(-1, -0.5, 0, 0.5, 1)) +
  # scale_fill_manual(values = myColors) +
  coord_sf(crs = sf::st_crs("EPSG:3035")) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
  theme_void() + 
  facet_wrap(. ~ epiweek, ncol = 3) +
  theme(strip.text.x = element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.key.size = unit(1, 'cm')) 

ggsave(here("figures/rt_estimates_differences.svg"),width=12,height=10, dpi=600)
ggsave(here("figures/rt_estimates_differences.png"),width=12,height=10, dpi=600)
ggsave(here("figures/rt_estimates_differences.pdf"),width=12,height=10, dpi=600)

base + 
  geom_sf(aes(fill = Ratio)) +
  # scale_fill_gradientn(values=c(0,0.1,0.8,1),
  #                      colors = c("#e9ecef", "#ff9f1c")) +
  scale_fill_gradientn(colors=leaflet::colorBin(grDevices::colorRamp(c(hsv(0.06, 0.0, 0.9),
                                                                       hsv(0.03, 0.0, 0.8),
                                                                       hsv(0.01, 0.4, 0.8),
                                                                       hsv(0.99, 1.0, 0.7)),
                                                                     bias=1.2),
                                                domain = c(0,1))(seq(0, 1, 0.1))) +
  # scale_fill_manual(values = myColors) +
  coord_sf(crs = sf::st_crs("EPSG:3035")) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
  theme_void() +
  facet_wrap(. ~ epiweek, ncol = 3) +
  theme(strip.text.x = element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.key.size = unit(1, 'cm')) 
ggsave(here("figures/rt_estimates_ratio.svg"),width=12,height=10, dpi=600)
ggsave(here("figures/rt_estimates_ratio.png"),width=12,height=10, dpi=600)
ggsave(here("figures/rt_estimates_ratio.pdf"),width=12,height=10, dpi=600)