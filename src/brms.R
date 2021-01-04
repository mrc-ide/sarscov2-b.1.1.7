
library(dplyr)
library(brms)
library(pracma)
library(here)
## read in #######################################################

x9 <- as.data.frame(readRDS(here('data/sgss_stp.rds')))
x10 <- as.data.frame(readRDS(here('data/stp_data_genome.rds')))
x11 <- as.data.frame(readRDS(here('data/stp_rt.rds')))
area_names <- unique(x9$area)
x11 <- subset(x11,area%in%area_names)
x11 <- left_join(x11,x9,by=c('area','epiweek'))
colnames(x11)[colnames(x11)=='area'] <- 'stp_name'
x10 <- left_join(x11,x10,by=c('stp_name','epiweek'))
colnames(x10)[colnames(x10)=='B.1.1.7'] <- 'B117'

# Rt columns
xtemp <- x10
xtemp$epiweek <- xtemp$epiweek - 1
xtemp$Rtnextweek <- xtemp$value
xtemp$lownextweek <- xtemp$CIlow
xtemp$upnextweek <- xtemp$CIup
x10 <- left_join(x10,select(xtemp,c('stp_name','epiweek','Rtnextweek','lownextweek','upnextweek')),
                 by=c('stp_name','epiweek'))

## stan variables
iters <- 8000
cores <- 4
chains <- 4
thin <- 10

## prepare dataset ################################################

areas <- unique(x10$regcd)
w1 <- 44; w2 <- 50
wks <- w1:w2
lagflags <- c('lag','nolag')
for(lagflag in lagflags){
  
  # shift Rt by one week, so genomes now correspond to Rt one week hence
  # lag or no lag
  d1 <- x10
  if(lagflag=='nolag')
    colnames(d1)[colnames(d1)=='value'] <- 'Rt'
  if(lagflag=='lag'){
    d1$Rt <- d1$Rtnextweek
    d1$CIlow <- d1$lownextweek
    d1$CIup <- d1$upnextweek
    # to set weeks onto Rt scale (not genome data scale)
    d1$epiweek <- d1$epiweek+1
  }
  d1 <- subset(d1,!is.na(Rt))
  
  ## R values
  byweek <- expand.grid(regcd=areas,epiweek=wks,stringsAsFactors = F)
  byweek <- left_join(byweek,d1[,colnames(d1)%in%c("epiweek","Rt","CIlow","CIup",'other','B117',"regcd","sgss_s_negative_corrected","sgss_s_positive_corrected")],by=c('epiweek','regcd'))
  byweek$Rstd <- apply(byweek,1,function(x){
    df <- subset(byweek,epiweek==x[2]&regcd==as.character(x[1])); 
    (df$CIup-df$CIlow)/(qnorm(0.975)*2)
  })
  #var <- byweek$Rstd^2
  #byweek$logRstd <- log(1+var^2/byweek$Rt^2)
  week_ind <- which(colnames(byweek)=='epiweek')
  nm_ind <- which(colnames(byweek)=='regcd')
  
  ## time of min Rt
  minweek <- sapply(areas,function(x){sub=subset(byweek,regcd==x);sub$epiweek[which.min(sub$Rt)]})
  byweek$relweek <- byweek$epiweek - minweek[match(byweek$regcd,names(minweek))]
  #byweek$epiweek <- byweek$relweek
  
  ## Genomic: wt / variant
  byweek$sample_total <- byweek$other + byweek$B117
  
  ## S gene: positive / negative reads
  byweek$negative <- round(byweek$sgss_s_negative_corrected)
  byweek$positive <- round(byweek$sgss_s_positive_corrected)
  byweek$sgene_total <- byweek$positive + byweek$negative
  d4 <- subset(byweek,sgene_total > 0 | sample_total > 0)
  
  d4$sgene_freq <- d4$negative/d4$sgene_total
  d4$sample_freq <- d4$B117/d4$sample_total
  
  
  ## genome data #########################################################
  
  # 6 brm with latent variable for frequency
  d4$fr <- as.numeric(NA)        
  
  ## R -- additive, S gene, assuming a linear effect of frequency on Rt, using sequence samples
  # fixed effect of STP
  bf_B <- bf(B117 | trials(sample_total) ~ 0+  mi(fr) ,  family = binomial(link = logit))
  bf_R <- bf( Rt | mi(Rstd) ~  I(exp( mi(fr)   ) / ( 1 + exp( mi(fr)))) + s(scale(epiweek),k=4)  + (regcd))
  
  d5 <- subset(d4,sample_total>0)
  area_coef <- unique(d4$regcd)
  area_coef <- area_coef[- length(area_coef)]
  
  bf_f <- bf(fr | mi() ~ regcd*scale(epiweek) ,family=gaussian)
  #bf_R <- bf( Rt | mi(Rstd) ~  mi(fr) + s(scale(epiweek),k=4)  + (regcd))
  mod3.4 <- brm(bf_B + bf_R + bf_f, data = d5,              
                chains = chains, cores = cores, iter=iters,thin=thin,
                prior = c(prior(constant(1), coef = mifr, resp = B117),
                          prior(normal(0, 1), class='Intercept',resp='Rt'),
                          prior(normal(0, 1), class='Intercept',resp='fr'),
                          prior(normal(0, 1), coef='scaleepiweek',resp='fr'),
                          prior_string('normal(0, 1)',resp='fr',coef = paste("regcd", area_coef,':scaleepiweek', sep="")),
                          prior_string('normal(0, 1)',resp='Rt',coef = paste("regcd", area_coef, sep="")),
                          prior_string('normal(0, 1)',resp='fr',coef = paste("regcd", area_coef, sep=""))),
                control = list(adapt_delta = 0.99, max_treedepth = 15))
  # lag :  0.67      (0.42     0.93)   /    0.37  0.15  0.62  relweek
  # no lag : 0.54       (0.32     0.77) / 0.20  0.02  0.41 relweek
  saveRDS(mod3.4,paste0(here(),'/results/','mod3.4-',lagflag,w1,w2,'.rds'))
  print(round(quantile(posterior_samples(mod3.4,pars='Rt_IexpmifrD1Pexpmifr')[,1],c(0.5,0.025,0.975)),2))
  
  ## S gene ####################################################
  
  d5 <- subset(d4,sgene_total>0)
  area_coef <- unique(d4$regcd)
  area_coef <- area_coef[-length(area_coef)]
  
  bf_B <- bf(negative | trials(sgene_total) ~ 0+  mi(fr) , family = binomial(link = logit))
  bf_R <- bf( Rt | mi(Rstd) ~  I(exp( mi(fr)   ) / ( 1 + exp( mi(fr)))) + s(scale(epiweek),k=4)  + (regcd))
  #bf_R <- bf( Rt | mi(Rstd) ~  mi(fr) + s(scale(epiweek),k=4)  + (regcd))
  mod4.2 <- brm(bf_B + bf_R + bf_f, data = d5, 
                chains = chains, cores = cores, iter=iters,thin=thin,
                prior = c(prior(constant(1), coef = mifr, resp = negative),
                          prior(normal(0, 1), class='Intercept',resp='Rt'),
                          prior(normal(0, 1), class='Intercept',resp='fr'),
                          prior(normal(0, 1), coef='scaleepiweek',resp='fr'),
                          prior_string('normal(0, 1)',resp='fr',coef = paste("regcd", area_coef,':scaleepiweek', sep="")),
                          prior_string('normal(0, 1)',resp='fr',coef = paste("regcd", area_coef, sep="")),
                          prior_string('normal(0, 1)',resp='Rt',coef = paste("regcd", area_coef, sep=""))),
                control = list(adapt_delta = 0.99, max_treedepth = 15))
  saveRDS(mod4.2,paste0(here(),'/results/','mod4.2-',lagflag,w1,w2,'.Rds'))
  print(round(quantile(posterior_samples(mod4.2,pars='Rt_IexpmifrD1Pexpmifr')[,1],c(0.5,0.025,0.975)),2))
  # lag : 0.48   (0.33     0.64)  /  0.30  0.13  0.46 relweek
  # no lag : 0.42      (0.29     0.54) /  0.26  0.14  0.39 relweek
  
}

for(lagflag in lagflags){
  cat(paste0('\nFor ',lagflag,':\n'))
  
  cat('Binomial genome model:\n')
  x <- readRDS(paste0(here(),'/results/','mod3.4-',lagflag,w1,w2,'.rds'))
  print(round(quantile(posterior_samples(x,pars='Rt_IexpmifrD1Pexpmifr')[,1],c(0.5,0.025,0.975)),2))
  
  cat('Binomial S gene model:\n')
  x <- readRDS(paste0(here(),'/results/','mod4.2-',lagflag,w1,w2,'.rds'))
  print(round(quantile(posterior_samples(x,pars='Rt_IexpmifrD1Pexpmifr')[,1],c(0.5,0.025,0.975)),2))
  
}
