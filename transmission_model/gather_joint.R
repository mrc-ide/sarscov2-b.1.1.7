library(here)
library(tidyverse)
library(sf)
library(viridis)
library(stringr)
library(epidemia)
library(matrixStats)
library(lubridate)
library(ggforce)
resdir <- here("transmission_model/fits/")
allf <- list.files(resdir,pattern=".*joint.rds")
##consensus
neg_45 <- NULL
neg_46 <- NULL
neg_47 <- NULL
neg_48 <- NULL
neg_49 <- NULL
neg_50 <- NULL
neg_51 <- NULL
neg_52 <- NULL
neg_53 <- NULL
neg_54 <- NULL
neg_55 <- NULL
neg <- sapply(allf,function(f){
    load(file.path(resdir,f))
    neg <- as.matrix(res_new$fit, pars = "R|negTRUE")
    neg_time_45 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[45,all]") + neg 
    neg_time_46 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[46,all]") + neg_time_45
    neg_time_47 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[47,all]") + neg_time_46
    neg_time_48 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[48,all]") + neg_time_47
    neg_time_49 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[49,all]") + neg_time_48
    neg_time_50 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[50,all]") + neg_time_49
    neg_time_51 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[51,all]") + neg_time_50
    neg_time_52 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[52,all]") + neg_time_51
    neg_time_53 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[53,all]") + neg_time_52
    neg_time_54 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[54,all]") + neg_time_53
    neg_time_55 <- as.matrix(res_new$fit, pars = "R|rw(time = negweek)[55,all]") + neg_time_54
    neg_45 <<- rbind(neg_45, neg_time_45)
    neg_46 <<- rbind(neg_46, neg_time_46)
    neg_47 <<- rbind(neg_47, neg_time_47)
    neg_48 <<- rbind(neg_48, neg_time_48)
    neg_49 <<- rbind(neg_49, neg_time_49)
    neg_50 <<- rbind(neg_50, neg_time_50)
    neg_51 <<- rbind(neg_51, neg_time_51)
    neg_52 <<- rbind(neg_52, neg_time_52)
    neg_53 <<- rbind(neg_53, neg_time_53)
    neg_54 <<- rbind(neg_54, neg_time_54)
    neg_55 <<- rbind(neg_55, neg_time_55)
    
    return(NULL)
})

x <- bind_rows(
    round(quantile(exp(rowMeans(neg_45)),c(0.5,0.025,0.975)),2),
    round(quantile(exp(rowMeans(neg_46)),c(0.5,0.025,0.975)),2),
    round(quantile(exp(rowMeans(neg_47)),c(0.5,0.025,0.975)),2),
    round(quantile(exp(rowMeans(neg_48)),c(0.5,0.025,0.975)),2),
    round(quantile(exp(rowMeans(neg_49)),c(0.5,0.025,0.975)),2), 
    round(quantile(exp(rowMeans(neg_50)),c(0.5,0.025,0.975)),2),
    round(quantile(exp(rowMeans(neg_51)),c(0.5,0.025,0.975)),2),
    round(quantile(exp(rowMeans(neg_52)),c(0.5,0.025,0.975)),2),
    round(quantile(exp(rowMeans(neg_53)),c(0.5,0.025,0.975)),2),
    round(quantile(exp(rowMeans(neg_54)),c(0.5,0.025,0.975)),2),
    round(quantile(exp(rowMeans(neg_55)),c(0.5,0.025,0.975)),2))

x$week <- seq(45,55)
x <- relocate(x, week)
print(x)

facet_plot_data <- tibble()
map_data <- lapply(allf,function(f){
    load(file.path(resdir,f))
    prediction <- posterior_predict(res_new$fit, seed= 12345)
    facet_plot_data <<- 
        bind_rows (
            facet_plot_data,
            left_join(tibble(group = prediction$group,
                             date = prediction$time, 
                             Prediction = matrixStats::colMedians(prediction$draws),
                             Pli = matrixStats::colQuantiles(prediction$draws, probs = 0.025),
                             Pui = matrixStats::colQuantiles(prediction$draws, probs = 0.975)), 
                      res_new$fit$data) %>%
                mutate(SGTF = case_when(neg ~ 'S-', !neg ~ 'S+'),
                       Area = res_new$stp) %>%
                select(Area, SGTF, Prediction, Cases_week, Pli, Pui, week)
        )
    
    rt <- posterior_rt(res_new$fit, seed = 12345)    
    
    tibble(area=res_new$stp,
           Ratio_45=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==45)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==45)],0.5),
           Ratio_46=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==46)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==46)],0.5),
           Ratio_47=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==47)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==47)],0.5),
           Ratio_48=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==48)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==48)],0.5),
           Ratio_49=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==49)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==49)],0.5),
           Ratio_50=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==50)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==50)],0.5),
           Ratio_51=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==51)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==51)],0.5),
           Ratio_52=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==52)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==52)],0.5),
           Ratio_53=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==53)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==53)],0.5),
           Ratio_54=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==1)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==1)],0.5),
           Ratio_55=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==2)] / rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==2)],0.5),
           Difference_45=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==45)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==45)],0.5),
           Difference_46=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==46)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==46)],0.5),
           Difference_47=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==47)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==47)],0.5),
           Difference_48=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==48)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==48)],0.5),
           Difference_49=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==49)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==49)],0.5),
           Difference_50=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==50)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==50)],0.5),
           Difference_51=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==51)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==51)],0.5),
           Difference_52=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==52)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==52)],0.5),
           Difference_53=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==53)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==53)],0.5),
           Difference_54=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==1)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==1)],0.5),
           Difference_55=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==2)] - rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==2)],0.5),
           'R(S-)_45'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==45)],0.5),
           'R(S-)_46'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==46)],0.5),
           'R(S-)_47'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==47)],0.5),
           'R(S-)_48'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==48)],0.5),
           'R(S-)_49'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==49)],0.5),
           'R(S-)_50'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==50)],0.5),
           'R(S-)_51'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==51)],0.5),
           'R(S-)_52'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==52)],0.5),
           'R(S-)_53'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==53)],0.5),
           'R(S-)_54'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==1)],0.5),
           'R(S-)_55'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_negative") & epiweek(rt$time)==2)],0.5),
           'R(S+)_45'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==45)],0.5),
           'R(S+)_46'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==46)],0.5),
           'R(S+)_47'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==47)],0.5),
           'R(S+)_48'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==48)],0.5),
           'R(S+)_49'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==49)],0.5),
           'R(S+)_50'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==50)],0.5),
           'R(S+)_51'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==51)],0.5),
           'R(S+)_52'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==52)],0.5),
           'R(S+)_53'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==53)],0.5),
           'R(S+)_54'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==1)],0.5),
           'R(S+)_55'=quantile(rt$draws[,which(rt$group==str_c(res_new$stp,"_positive") & epiweek(rt$time)==2)],0.5),
    )
})
map_data <- do.call("bind_rows",map_data)
map_data <- pivot_longer(map_data, cols = -area,
                         names_to = c(".value", "epiweek"),
                         names_pattern = "(.+)_(.+)")

load_uk_stp <- function() {
    shape <- st_read(here("transmission_model/map_data/Sustainability_and_Transformation_Partnerships__April_2020__Boundaries_EN_BUC.shp"),
                     quiet = TRUE) %>%
        st_transform("+proj=longlat +datum=WGS84 +init=epsg:3035") %>%
        select(stp20nm) %>%
        rename(name=stp20nm)
    
    shape <- shape %>% 
        mutate(name=gsub("Frimley Health and Care ICS", "Frimley Health", name)) %>% 
        mutate(name=gsub("Cornwall and the Isles of Scilly Health and Social Care Partnership", "Cornwall and the Isles of Scilly", name)) %>% 
        mutate(name=gsub("Surrey Heartlands Health and Care Partnership", "Surrey Heartlands", name)) %>% 
        mutate(name=gsub("Sussex and East Surrey Health and Care Partnership", "Sussex and East Surrey", name))
    
    shape
}

final_data <- inner_join(load_uk_stp(), map_data, by = c("name"= "area")) %>%
    # filter(epiweek == 48 | epiweek == 50) %>%
    filter(epiweek < 56) %>%
    mutate(epiweek = str_c("Week ", epiweek))

ggplot(data = final_data) + 
    theme_bw() +
    theme(axis.text.x =element_text(angle = 45,hjust = 1),
            axis.text =element_text(size = 22),
            axis.title =element_text(size = 22)
          ) +
    theme(legend.position = "right") +
    geom_sf(aes(fill = Ratio)) +
    scale_fill_gradientn(colors=leaflet::colorBin(grDevices::colorRamp(c(hsv(0.06, 0.0, 0.9),
                                                                         hsv(0.03, 0.0, 0.8),
                                                                         hsv(0.01, 0.4, 0.8),
                                                                         hsv(0.99, 1.0, 0.7)),
                                                                       bias=1.2),
                                                  domain = c(0,1))(seq(0, 1, 0.2))) +
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

facet_plot_data <-
    facet_plot_data %>%
    mutate(facet_name = str_c(Area," (",SGTF,")"))
base <- ggplot(data = facet_plot_data, aes(x = Cases_week, y = Prediction, ymin = Pli, ymax = Pui)) + 
    geom_blank() + 
    theme_bw() +
    theme(axis.text.x =element_text(angle = 45,hjust = 1),
          axis.text =element_text(size = 22),
          axis.title =element_text(size = 22)
    )

for(i in seq(1,7)){
    print(
        base + 
        geom_abline() +
        geom_pointrange(size=0.1, color="black", fill="white", shape=22) + 
        theme_light()+
        theme(strip.text.x = element_text(
            color = "black"),
            strip.background = element_blank()
        ) +
        ylab("Prediction") +
        xlab("Cases per week") +
        facet_wrap_paginate(. ~ facet_name, scales = "free", ncol = 2, nrow = 6, page = i)
    )
}

base + 
    geom_abline() +
    geom_pointrange(size=0.1, color="black", fill="white", shape=22) + 
    theme_light()+
    theme(strip.text.x = element_text(
        color = "black"),
        strip.background = element_blank()
    ) +  # guides(color = guide_legend(override.aes = list(size = c(2,1,1)))) +
    ylab("Prediction") +
    xlab("Cases per week") +
    facet_wrap(. ~ SGTF, scales = "free", ncol = 2)
