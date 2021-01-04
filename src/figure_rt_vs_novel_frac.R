library(tidyverse)
library(here)
library(data.table)

stp_rt <- 
  readRDS(here("data/combine_rt_stp_all_neg_pos.rds"))  #%>%
 # select(epiweek, area, value)
sgss = readRDS(here("data/sgss_stp.rds"))

sgss.to.areas = fread(here("data/ltla-nhser-stp.csv")) %>% select("stp_name","nhser_name")
sgss.to.areas = sgss.to.areas[, list(nhser_name = nhser_name[1]),by=stp_name]
sgss = left_join(sgss,sgss.to.areas,by=c("area"="stp_name"))
stp_rt = left_join(stp_rt,sgss,by=c("epiweek","area"))

stp_rt <- filter(stp_rt, (sgss_s_positive > 5 & sgss_s_negative > 5))
lockdown_start = (45 + 4/7)
lockdown_end = (49 + 2/7)

stp_rt$freq = 100*stp_rt$sgss_s_negative_corrected / (stp_rt$sgss_s_negative_corrected+stp_rt$sgss_s_positive_corrected)
#lockdown = filter(stp_rt,epiweek >= lockdown_start & epiweek <= lockdown_end)
stp_rt = filter(stp_rt,epiweek >= 45 & epiweek <= 50)
# stp_rt$epiweek = factor(stp_rt$epiweek)
g = ggplot(stp_rt,aes(x=`R(S+)`,y=`R(S-)`))
g = g + geom_point(aes(size=freq,colour=nhser_name,shape=as.factor(epiweek)),alpha=.9)
#g = g + scale_shape_discrete(values = c(0:3,8,9))
g = g + geom_abline(slope=1)
g = g + geom_vline(xintercept=1,col="gray")
g = g + scale_x_continuous(limits=c(0,2.5),expand=c(0,0))
g = g + scale_y_continuous(limits=c(0,2.5),expand=c(0,0))
g = g + geom_hline(yintercept=1,col="gray")
g = g + theme_bw()
g = g + coord_fixed()
g = g + labs(size="% novel variant", colour="",shape="epiweek") 
#g

g = g + scale_shape_manual(values = c(0:3,5,6)) #,solid=T)
g
ggsave(here("figures/estimates-by-week.pdf"),g,width=10,height=8)
ggsave(here("figures/estimates-by-week.png"),g,width=10,height=8)

stp_rt <- 
  readRDS(here("data/combine_rt_stp_all_neg_pos.rds"))  %>%
  select(epiweek, area, Rt)
stp_geneome <- 
  readRDS(here("data/stp_data_genome.rds")) %>%
  rename(area = stp_name, novel = "B.1.1.7") %>%
  mutate(novel_frac = novel/(novel + other)) %>%
  select(area, novel_frac, epiweek)
  
stp_combined <-
  inner_join(stp_rt, stp_geneome) %>%
  filter(epiweek > 41, epiweek < 51) %>%
  mutate(epiweek = str_c("Week ", epiweek))

base <- ggplot(data = stp_combined) + 
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

g1 <- base + 
  geom_point(mapping = aes(Rt, novel_frac), alpha=0.8, size = 3) +
  facet_wrap(. ~ epiweek, ncol = 4) +
  xlab("Time Varying Reproduction Number") +
  ylab("Percentage of new variant among all genomes") +
  theme(strip.text.x = element_text(size = 18))
g1
ggsave(here('figures/rt_vs_novel_fraction.pdf'), width=11,height=8,dpi=600)
ggsave(here('figures/rt_vs_novel_fraction.png'), width=11,height=8,dpi=600)
ggsave(here('figures/rt_vs_novel_fraction.svg'),width=11,height=8,dpi=600)
