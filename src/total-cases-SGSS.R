library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(data.table)
library(here)

pops = fread(here("data/stp_population_long.csv"))[subgroup == "All"]
setnames(pops,"Y2018","pop")
sgss = readRDS(here("data/sgss_stp_new_43_56_weeks.rds"))

sgss <-
  sgss %>%
  select(-c(sgss_s_positive, sgss_s_negative, total_cases, sgss_s_negative_corrected, sgss_s_positive_corrected)) %>%
  rename(sgss_s_positive = sgss_s_positive_adj1 ,
         sgss_s_negative = sgss_s_negative_adj1,
         total_cases = total_cases_adj1,
         sgss_s_positive_corrected = sgss_s_positive_corrected_adj1,
         sgss_s_negative_corrected  = sgss_s_negative_corrected_adj1)

lockdown_start = (45 + 4/7)
lockdown_end = (49 + 2/7)
third_lockdown = (54 + 1/7)
sgss = left_join(sgss,pops,by=c("area" = "AREA"))

df = data.table(sgss)[, list(neg_corrected = total_cases * sgss_s_negative_corrected / 
                               (sgss_s_negative_corrected + sgss_s_positive_corrected) / pop * 100000 ,
               pos_scaledup = total_cases * sgss_s_positive_corrected / 
                 (sgss_s_negative_corrected + sgss_s_positive_corrected),# / pop * 100000 ,
               neg_scaledup = total_cases * sgss_s_negative / 
                 (sgss_s_negative + sgss_s_positive), # / pop * 100000 ,
               pos_uncorrected = total_cases * sgss_s_positive / 
                 (sgss_s_negative + sgss_s_positive) / pop * 100000 ,
               total = total_cases, # /pop * 100000,
               neg = sgss_s_negative, # / pop * 100000 ,
               pos = sgss_s_positive, # / pop * 100000 ,
               n501y = sum(sgss_s_negative_corrected) / sum(sgss_s_negative_corrected + sgss_s_positive_corrected)),
        by=list(epiweek, area)]

# for disclosure purposes according to PHE, need to not release counts below 5
df$neg[df$neg < 5] = NA
df$pos[df$pos < 5] = NA

n501y.order = df[epiweek == 50,]
n501y.order  = n501y.order$area[order(n501y.order$n501y,decreasing=TRUE)] # order from lowest n501y % to highest

df$area_sorted = factor(df$area,levels=n501y.order)
STPs = c(1,17,33) #2,5,17,25,29,33,38,41)
a <- ggplot(df[area %in% n501y.order[STPs] & epiweek>=45 & epiweek <= 56,], aes(y=total, x=epiweek)) + 
  geom_vline(xintercept = lockdown_start,col="darkred") +
  geom_vline(xintercept = lockdown_end,col="darkred") +
  geom_vline(xintercept = third_lockdown,col="darkred") +
  
  geom_rect(aes(xmin = lockdown_start, xmax = lockdown_end, ymin = -Inf, ymax = Inf),
            fill = "pink", alpha = 0.03) +
  geom_rect(aes(xmin = third_lockdown, xmax = 56, ymin = -Inf, ymax = Inf),
            fill = "pink", alpha = 0.03) +
  
  geom_line(lwd=2, aes(colour=n501y)) +
  #theme(legend.position = "none") +
  scale_colour_gradient(low = "#51c7b7", high = "#f99d38", 
                        name='% S-',
                        limits=c(0,1), labels = scales::percent(0.25*0:4) )  +
#  theme(legend.position="top",
#        legend.justification="right",
#        legend.margin=margin(0,0,0,0),
#        legend.box.margin=margin(-10,-10,-10,-10)) +
  geom_line(aes(y=neg),colour="#f99d38",lwd=1) +
  geom_line(aes(y=pos),colour="#51c7b7",lwd=1) +
  geom_text(data=df[area %in% n501y.order[STPs] & epiweek == 56], aes(x=epiweek, y=neg),label="S-",colour="#f99d38",
            hjust=0,nudge_x=.15) +
  geom_text(data=df[area %in% n501y.order[STPs] & epiweek == 56], aes(x=epiweek, y=total),label="total",colour="#555555",
            hjust=0,nudge_x=.15) +
  geom_text(data=df[area %in% n501y.order[STPs] & epiweek == 56], aes(x=epiweek, y=pos),label="S+",colour="#51c7b7",
            hjust=0,nudge_x=.15) +
  #ggtitle(unique(csv_stp$area))+
  scale_y_continuous(expand=c(0,0)) + #,limits=c(0,800)) +
  scale_x_continuous(expand=c(0,0),labels=c("Week",paste(c(46,48,50,52,1,3))),breaks=seq(44,56,2)) +
  coord_cartesian(xlim=c(44,56),clip='off') +
  theme_light()+ 
  theme(strip.text.x = element_text(
    color = "black"), 
    strip.background = element_blank()
  ) +  # guides(color = guide_legend(override.aes = list(size = c(2,1,1)))) +
  xlab("") +
  ylab("Number of cases") +
  facet_wrap(~ area_sorted, scales="free",labeller = label_wrap_gen()) + 
  theme(panel.spacing = unit(2, "lines"), legend.position="top",
        plot.margin=margin(r=1.5,unit='cm')) #legend.margin=margin(l = 1.5, unit='cm'))
a
ggsave(here("figures/cases-SGSS-43-56-3panel.png"),width=9,height=3)
ggsave(here("figures/cases-SGSS-43-56-3panel.pdf"),width=9,height=3)

df$abbrev = gsub("North East","NE",df$area) 
df$abbrev = gsub("South East","SE",df$abbrev)
df$abbrev = gsub("North West","NW",df$abbrev)
df$abbrev = gsub("South West","SW",df$abbrev)
n501y.order = df[epiweek == 56,]
n501y.order  = n501y.order$abbrev[order(n501y.order$n501y,decreasing=TRUE)] # order from lowest n501y % to highest

df$abbrev_sorted = factor(df$abbrev,levels=n501y.order)

a <- ggplot(df[epiweek>=45 & epiweek <= 56,], aes(y=total, x=epiweek)) + 
  geom_vline(xintercept = lockdown_start,col="darkred") +
  geom_vline(xintercept = lockdown_end,col="darkred") +
  geom_vline(xintercept = third_lockdown,col="darkred") +
  geom_line(lwd=2, aes(colour=n501y)) +
  #theme(legend.position = "none") +
  
  geom_rect(aes(xmin = lockdown_start, xmax = lockdown_end, ymin = -Inf, ymax = Inf),
            fill = "pink", alpha = 0.03) +
  geom_rect(aes(xmin = third_lockdown, xmax = 56, ymin = -Inf, ymax = Inf),
            fill = "pink", alpha = 0.03) +
  
  scale_colour_gradient(low = "#51c7b7", high = "#f99d38", 
                        name='% S-',
                        limits=c(0,1), labels = scales::percent(0.25*0:4) )  +
  theme(legend.position="top",
        legend.justification="right",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) +
  geom_line(aes(y=neg),colour="#f99d38",lwd=1) +
  geom_line(aes(y=pos),colour="#51c7b7",lwd=1) +
  geom_text(data=df[area %in% n501y.order[STPs] & epiweek == 56], aes(x=epiweek, y=neg),label="S-",colour="#f99d38",
            hjust=0,nudge_x=.15) +
  geom_text(data=df[area %in% n501y.order[STPs] & epiweek == 56], aes(x=epiweek, y=total),label="total",colour="#555555",
            hjust=0,nudge_x=.15) +
  geom_text(data=df[area %in% n501y.order[STPs] & epiweek == 56], aes(x=epiweek, y=pos),label="S+",colour="#51c7b7",
            hjust=0,nudge_x=.15) +
  #ggtitle(unique(csv_stp$area))+
  scale_y_continuous(expand=c(0,0)) + #,limits=c(0,800)) +
  scale_x_continuous(expand=c(0,0),labels=c("Week",paste(c(46,48,50,52,1,3))),breaks=seq(44,56,2)) +
  coord_cartesian(xlim=c(44,56),clip='off') +
  theme_light()+ 
  theme(strip.text.x = element_text(
    color = "black"), 
    strip.background = element_blank()
  ) +  # guides(color = guide_legend(override.aes = list(size = c(2,1,1)))) +
  xlab("") +
  ylab("Number of cases") +
  facet_wrap(~ abbrev_sorted, scales="free",labeller = label_wrap_gen()) + 
  theme(panel.spacing.x = unit(.2, "lines")) ## This needs to be fixed

ggsave(here("figures/cases-SGSS-43-56.pdf"),width=16,height=12)
ggsave(here("figures/cases-SGSS-43-56.png"),width=16,height=12)
