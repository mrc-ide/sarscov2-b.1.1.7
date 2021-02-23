
library(dplyr)
library(brms)
library(pracma)
library(lubridate)
library(here)

## read in #######################################################

x9 <- as.data.frame(readRDS(here('data/sgss_stp_new_43_56_weeks.rds')))
x11 <- readRDS(here('data/stp_rt_all.rds'))
area_names <- unique(x9$area)
x11 <- subset(x11,area%in%area_names)
x11 <- full_join(x11,x9,by=c('area','epiweek'))
ei <- which(colnames(x11)=='epiweek')
si <- which(colnames(x11)=='regcd')

x10 <- x11

# Rt columns
xtemp <- x10
xtemp$epiweek <- xtemp$epiweek - 1
xtemp$Rtnextweek <- xtemp$`R(all)`
xtemp$lownextweek <- xtemp$RallLI
xtemp$upnextweek <- xtemp$RallUI
x10 <- full_join(x10,select(xtemp,c('regcd','epiweek','Rtnextweek','lownextweek','upnextweek')),
                 by=c('regcd','epiweek'))

## stan variables #############################################################
library(rstan)
options(mc.cores = parallel::detectCores())
iters <- 4000
cores <- 2
chains <- 2
thin <- 1
stanvars <- stanvar(scode="  real sigmoid(real x){
  return 1/(1+exp(-x));
  }  
",block='functions') 

fileName <- here('src/stan_models//addvsmult.stan')
stan_code <- readChar(fileName, file.info(fileName)$size)

bf_f <- bf(fr | mi() ~ regcd*scale(epiweek) ,family=gaussian)
bf_B <- bf(negative | trials(sgene_total) ~ 0+  mi(fr) , family = binomial(link = logit))
bf_R <- bf( Rt | mi(Rstd) ~  sigmoid( mi(fr)) + s(scale(epiweek),k=4)  + (regcd))

dts <- c('genome','sgss')
models <- c('multiplicative','additive')
lagflags <- c('nolag')#,'lag')
lagflag='nolag' #for(lagflag in lagflags){

## prepare dataset ################################################

areas <- unique(x10$regcd)
w1s <- c(45)
w2s <- c(50)
w2 <- max(w2s)
w1 <- min(w1s)

# shift Rt by one week, so genomes now correspond to Rt one week hence
# lag or no lag
d1 <- x10
if(lagflag=='nolag'){
  colnames(d1)[colnames(d1)=='R(all)'] <- 'Rt'
  d1$CIlow <- d1$RallLI
  d1$CIup <- d1$RallUI
}
if(lagflag=='lag'){
  d1$Rt <- d1$Rtnextweek
  d1$CIlow <- d1$lownextweek
  d1$CIup <- d1$upnextweek
  # to set weeks onto Rt scale (not genome data scale)
  d1$epiweek <- d1$epiweek+1
}
#d1 <- subset(d1,!is.na(Rt))

## R values
wks <- w1:w2
byweek <- expand.grid(regcd=areas,epiweek=wks,stringsAsFactors = F)
byweek <- left_join(byweek,d1[,colnames(d1)%in%c("Rv","Rnv", "sgss_s_negative" ,"probs","epiweek","Rt","CIlow","CIup",'Other','B117',"regcd","sgss_s_negative_corrected","sgss_s_positive_corrected")],by=c('epiweek','regcd'))
byweek$Rstd <- apply(byweek,1,function(x){
  df <- subset(byweek,epiweek==x[2]&regcd==as.character(x[1])); 
  (df$CIup-df$CIlow)/(qnorm(0.975)*2)
})

## Genomic: wt / variant
#byweek$sample_total <- byweek$Other + byweek$B117

## S gene: positive / negative reads
byweek$negative <- round(byweek$sgss_s_negative_corrected)
byweek$positive <- round(byweek$sgss_s_positive_corrected)
byweek$sgene_total <- byweek$positive + byweek$negative
d4 <- subset(byweek,sgene_total > 0 ) # | sample_total > 0)

d4$sgene_freq <- d4$negative/d4$sgene_total
#d4$sample_freq <- d4$B117/d4$sample_total

d4$fr <- as.numeric(NA)  


## start ########################################################
w1 <- 45; 

results <- list()

#for(w2 in w2s){
wks <- w1:w2

absRerror <- Median <- l95 <- u95 <- nR <- covR <- c()
#for(i in 3)  {
i <- 3 
dt <- dts[ceiling(i/2)]
mod <- models[i%%2 + 1]

## select data ##################################
if(dt=='genome'){
  d5 <- subset(d4,sample_total>0&epiweek%in%wks)
  ## match names to formula bf_B: B117 -> negative, sample_total -> sgene_total
  d5$negative <- d5$B117
  d5$sgene_total <- d5$sample_total
}else if(dt=='sgss'){
  d5 <- subset(d4,sgene_total>0&epiweek%in%wks)
}

## run or load saved ####################################

resultFileName <- here(paste0('results/',dt,mod,lagflag,w1,w2,'2021.Rds'))

if(file.exists(resultFileName)){
  
  resStan <- readRDS(resultFileName)
}else{
  stan_data <- make_standata(bf_B + bf_R + bf_f, data = d5)
  
  stan_data$additive <- as.numeric(mod=='additive')
  
  resStan <- stan(model_code = stan_code, data = stan_data,
                  chains = chains,cores=cores, iter = iters, thin = thin,
                  control = list(adapt_delta = 0.99, max_treedepth = 15)
  )
  
  saveRDS(resStan,resultFileName)
}
print(resultFileName)
print(quantile(extract(resStan, pars = c('bsp_Rt'), permuted = TRUE, inc_warmup = FALSE,include = TRUE)$bsp_Rt))

