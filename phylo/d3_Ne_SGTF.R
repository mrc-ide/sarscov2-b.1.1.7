# Compares the Ne of B.1.1.7 estimated from mlesky against TPP-adjusted SGTF
# Required: 
# - RDS files of mlesky outputs from d2_mlesky.R

require( hrbrthemes )
require( scales )
require( ggplot2 )
require( ggrepel )


# Read in mlesky output
tN = readRDS( "Sample_England_sampler1_B.1.1.7_2021-02-13_n=3000_n_tree_dating_10_mlesky.rds" )

# Calculates the median and 95% HPD of Ne(t)
q_ne = as.data.frame(t(apply( tN$ne, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) )))


# Converts decimal date to epiweek; the first few weeks of 2021 are interpreted as weeks 54-56 of 2020
q_ne$epiweek = lubridate::epiweek(lubridate::date_decimal(tN$time))
q_ne$epiweek = ifelse(q_ne$epiweek < 4, q_ne$epiweek+53, q_ne$epiweek)

# Read SGSS data
sgss_stp_new_43_56_weeks <- readRDS("data/sgss_stp_new_43_56_weeks.rds")


# Creates data frame extracting the total number of TPP-adjusted SGTF across all STPs for each epiweek; Ne (and HPD) at the end of each epi week estimated from mlesky.
pldf <- as.data.frame(do.call(rbind, lapply(unique(sgss_stp_new_43_56_weeks$epiweek), function(week) {
  if(nrow(sgss_stp_new_43_56_weeks[sgss_stp_new_43_56_weeks$epiweek == week, ]) == 42) {# checking there are no duplicates; there are 42 STPs
    total_S_neg = sum(sgss_stp_new_43_56_weeks[sgss_stp_new_43_56_weeks$epiweek == week, "sgss_s_negative_corrected"])
    ne = tail( q_ne[which(q_ne$epiweek == week), ][,"50%"], 1)
    neub = tail( q_ne[which(q_ne$epiweek == week), ][,"97.5%"], 1)
    nelb = tail( q_ne[which(q_ne$epiweek == week), ][,"2.5%"], 1)
    return(c(week = week, total_S_neg = total_S_neg, ne = ne, neub = neub, nelb = nelb))
  }
}
)
)
)

# Remove any rows with NAs
pldf = pldf[!is.na(pldf$ne),] 

# Plot
pl = ggplot(pldf, aes(x = ne, y = total_S_neg)) + geom_point( shape = 15) + 
  geom_errorbarh(aes(xmin = nelb, xmax = neub, y = total_S_neg)) +
  scale_y_continuous(
    trans = "log10",
    breaks = function(x) {
      brks <- extended_breaks(Q = c(1, 5))(log10(x))
      10^(brks[brks %% 1 == 0])
    },
    labels = math_format(format = log10)
  ) +
  scale_x_continuous(
    trans = "log10",
    breaks = function(x) {
      brks <- extended_breaks(Q = c(1, 5))(log10(x))
      10^(brks[brks %% 1 == 0])
    },
    labels = math_format(format = log10)
  )+
  theme_bw() + labs(x = "Effective population size B.1.1.7", y = "TPP-adjusted SGTF case numbers") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=14,face="bold"), panel.grid.minor = element_blank())+ annotation_logticks() +
  geom_label_repel(aes(x = ne, y = total_S_neg, label = week),alpha = 0.8)
pl

ggsave( plot = pl, file = "TPP-adjusted_SGTF_vs_Ne_mlesky_to_week_56.pdf", width = 8, height = 8 )
ggsave( plot = pl, file = "TPP-adjusted_SGTF_vs_Ne_mlesky_to_week_56.png", width = 8, height = 8 )

