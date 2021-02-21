### for epidemia with the new observation format
#### JOINT MODEL FOR NEW VARIANT

library(optparse)
library(here)
library(epidemia) ##new version of epidemia #install.packages("devtools") devtools::install_github("ImperialCollegeLondon/epidemia")
library(tidyverse)
library(readr)
options(mc.cores=parallel::detectCores())

option_list <- list(
  make_option(c("--areaname"),action="store", default="Cumbria_and_North_East",help="Which area to run this for [default \"%default\"]"),
  make_option(c("--input"),action="store", default=here("data/sgtf_transmission_data.rds"),help="Input file [default \"%default\"]"),
  make_option(c("--nchains"),action="store", type="integer", default=4,help="Number of Chains [default \"%default\"]"),
  make_option(c("--iter"),action="store", type="integer", default=2000,help="Number of iterations [default \"%default\"]"),
  make_option(c("--thin"),action="store", type="integer", default=1,help="Amount of thinning of results [default \"%default\"]"),
  make_option(c("--variant"),action="store", default="joint",help="Which varinat to use [default \"%default\"]"),
  make_option(c("--fit_params"),action="store", default=here("data/fit_params.csv"),  help="population file [default \"%default\"]"),
  make_option(c("--output"),action="store", default=here("transmission_model/fits"), help="Output directory [default \"%default\"]"),
  make_option(c("--i2o_rates"),action="store",default=here('data/i2o_rates.rds'),help="i2o rates file [default \"%default\"]"),
  make_option(c("--subgroup"),action="store", default="All",help="Which subgroup to use [default \"%default\"]")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (!is.element(opt$variant,c("joint"))){
  stop("variant is not one of the allowed options")
}
# parsing data
dataorig <- readRDS(opt$input)
testmult <- function(){
  opt$areaname <- unique(dataorig$area)[10:13]
}

data <- 
  dataorig %>%
  ungroup() %>%
  group_by(area) %>%
  arrange(date, .by_group=T)

Caseend <- max(data$date)
cat("Case End=",as.character(Caseend),"\n")

a <- opt$areaname
a <- gsub("_"," ",a)
cat("Running area =",a,"\n")

data <- 
  data %>%
  ungroup() %>%
  filter(is.element(area,a), date <= Caseend) %>%
  mutate(week = pmax(epiweek,42))

datapos <- data%>%mutate(Cases_week=corrected_positive)%>%
  mutate(stp=factor(gsub("\\(|\\)|,","",area)))%>%
  mutate(area=paste(area,"_positive",sep=""))%>%mutate(negweek=NA,neg=FALSE)%>% 
  select(area, date, week,Cases_week,negweek,neg,stp) %>%
  ungroup()

dataneg <- data%>%mutate(Cases_week=corrected_negative)%>%
  mutate(stp=factor(gsub("\\(|\\)|,","",area)))%>%
  mutate(area=paste(area,"_negative",sep=""))%>%mutate(negweek=pmax(epiweek,45),neg=TRUE)%>%
  select(area, date, week,Cases_week,negweek,neg,stp) %>%
  ungroup()
data <- bind_rows(datapos,dataneg)

# her efor historical reasons that data should span at least 3 weeks
if(nrow(data) >= 21){
  
  fileName = paste(opt$output,"/fm-", gsub(" ","_",a),"-",opt$variant,".rds",sep="")
  if(!file.exists(fileName)){
    i2o_rates <- readRDS(opt$i2o_rates) %>% filter(type=='IAR') %>% select(date, value) %>% rename(iar=value)
    data <- left_join(data, i2o_rates)%>%
      fill(iar,.direction="down")
    
    IAR_sd <- readRDS(opt$i2o_rates) %>% filter(type=='IAR_sd') %>% select(value)
    
    obs <- list()
    # weekly distribution
    i2o2week <- function(i2o)
      rowSums(sapply(0:6, function (k) c(rep(0,k),i2o,rep(0,6-k))))
    
    args <- list()
    pops <- read_csv(here("data/stp_population.csv")) %>% filter(subgroup == opt$subgroup, AREA == gsub("_"," ",opt$areaname))
    data$pop <- pops$Y2018
    args$data <- data
    args$obs <- epiobs(formula=Cases_week ~ 0 + iar,
                       link="identity",
                       # family="neg_binom", # if you beleive data is more overdispressed than a multiple of mean use neg
                       family="quasi_poisson",
                       prior_aux = rstanarm::normal(location=3,2),
                       prior=rstanarm::normal(1,IAR_sd$value,autoscale=FALSE),
                       i2o=i2o2week(c(0,0,0,rep(1/10,10))))
    if (opt$variant=="joint"){
      args$rt <- epirt(formula=R(area,date) ~  rw(time=week,prior_scale=.15)+neg+ rw(time=negweek), 
                       prior=rstanarm::normal(0,0.25),
                       prior_intercept=rstanarm::normal(log(1.25),.1),
                       link='log'
      )
    }
    seed_days <- 14
    args$inf=epiinf(gen=EuropeCovid$si,
                    prior_tau=rstanarm::exponential(rate = 1/5000), # uninformative seeding
                    seed_days=seed_days,
                    pop_adjust = TRUE,
                    susceptibles = pop)
    
    args$group_subset <- c(outer(a,c("_positive","_negative"),FUN=function(x,y)paste(x,y,sep="")))
    args$algorithm <- "sampling"
    adapt_delta <- 0.92
    args$data$area <- factor(args$data$area)
    args$iter=opt$iter
    args$chains=opt$nchains
    args$thin=opt$thin
    args$control=list(adapt_delta=adapt_delta,max_treedepth=18)
    cat("iter=",args$iter, "thin= ", args$thin, " adapt_delta=",adapt_delta,"\n")
    time <- system.time({fit <- do.call("epim", args)})
    
    sampler_params <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)
    Rhat <- max(rstan::summary(fit$stanfit)$summary[,"Rhat"])
    divergent <- sum(sapply(sampler_params, function(x) x[,"divergent__"]))
    cat("Rhat=",Rhat," divergent steps=", divergent,"\n")
    
    
    res_new <- list(fit=fit,
                    model=fit$formula,
                    last_obs_date=fit$data$date[max(which(!is.na(fit$data$Cases)))],
                    today=max(fit$data$date),
                    stp=a,
                    time=time)
    
    
    if (Rhat>=1.2||divergent>=opt$nchains*opt$iter/2.*.02){
      if (divergent>=opt$nchains*opt$iter/2.*.02) {
        adapt_delta <- (adapt_delta+1.)/2.
      } else if (Rhat>=1.2 && opt$iter<6000){
        opt$iter <- opt$iter*2
        opt$thin <- opt$thin*2
      }
      write.table(data.frame(Area_name=a,
                             iter=opt$iter,
                             thin=opt$thin,
                             adapt_delta=adapt_delta),
                  sep="&",append=TRUE,col.names=FALSE,row.names=FALSE,
                  file=opt$fit_params)
      dir.create(file.path(opt$output,"failedruns"), showWarnings = FALSE, recursive = TRUE)
      
      # save result to file
      save(res_new, file=paste(opt$output,"/failedruns/fm-", gsub(" ","_",a),Sys.time(), ".rds",sep=""))
      
      stop("Sampling not successful; parameters adjusted")
    }
    warnings()
    # save result to file
    save(res_new, file=paste(opt$output,"/fm-", gsub(" ","_",a),"-",opt$variant,".rds",sep=""))
    
  } else{
    cat("File",fileName," already exists for ",a,"\n")
  }
} else{
  cat("Can't process area ",a," as we do not have at least three weeks of data\n")
  res_new <- NA
  save(res_new, file=paste(opt$output,"/fm-", gsub(" ","_",a),"-",opt$variant,".rds",sep=""))
}
