library(loo)
library(dplyr)
library(rstan)
library(data.table)
library(lubridate)
library(RColorBrewer)
library(here)

## general variables #############################################################

# file name parts
outdir <- here('figures/')
tbit <- c('','Tminus')
databit <- c('sgss','both','genome')

## for stan 
iters <- 4000
cores <- 2
chains <- 2
thin <- 1

# parameters
w1 <- 45
w2 <- 55
wks <- w1:w2
fourages <- F # whether to run with four ages
## define age groups
agebands <- agebands0 <- list(c(0,10),c(11,18),c(19,49),c(50,120))

## read in #######################################################

# ## process data offline
# x15 <- read.csv('data/pillar2_sgtf_deaths_20210130-191401-ab9d0152.csv.xz',stringsAsFactors=F)
# x15 <- subset(x15,!is.na(stp_name))
# stp_to_reg=read.csv("data/stp2region.csv")
# x15=left_join(x15,stp_to_reg,by="stp_name")
# x15$nhser_name=as.factor(x15$nhser_name.y)
# x15$epiweek <- epiweek(x15$specimen_date)
# x15$epiweek[x15$epiweek<30] <- 53 + x15$epiweek[x15$epiweek<30]
# x15$stp_name <- as.factor(x15$stp_name)
# x17 <- read.csv('data/a3.2-prNellyDelOrOtherDel-STP-cr.csv',stringsAs=F)
# colnames(x17)[colnames(x17)=='fac_regcd'] <- 'stp_code'
# colnames(x17)[colnames(x17)=='sample_date'] <- 'specimen_date'
# x17 <- x17 %>% distinct(stp_code, specimen_date, .keep_all = TRUE)
# #x15 <- left_join(x15,x17[,colnames(x17)%in%c('specimen_date','stp_code','PrNellyDel_nellyDelOrOtherdel')],by=c('specimen_date','stp_code'))
# #x15$PrNellyDel_nellyDelOrOtherdel[x15$epiweek>=50] <- 1
# #x15$s_positive_adj1 <- x15$s_positive_adj1 + x15$s_negative_adj1  * (1-x15$PrNellyDel_nellyDelOrOtherdel)
# #x15$s_negative_adj1 <- x15$s_negative_adj1  * x15$PrNellyDel_nellyDelOrOtherdel
# agebands <- agebands0 <- list(c(0,10),c(11,18),c(19,49),c(50,120))
# x15$ageband <- paste(agebands[[1]],collapse='-')
# for(i in 2:length(agebands)) x15$ageband[x15$age>=agebands[[i]][1]] <- paste(agebands[[i]],collapse='-')
# x16 <- data.table(x15)[,lapply(.SD, sum, na.rm=TRUE),by=c('nhser_name','stp_name','epiweek','ageband'),
#                        .SDcols=c('s_positive_adj1','s_negative_adj1','s_na_adj1')]
# saveRDS(x16,'data/stp_age_sgtf.Rds')
# 
x16 <- readRDS(here('data/stp_age_sgtf.Rds'))

if(!fourages){
  agebands <- agebands0 <- list(c(11,18),c(19,120))
  
  x16 <- subset(x16,ageband!='0-10')
  x16$ageband[x16$ageband!='11-18'] <- '19-120'
  x16 <- data.table(x16)[,lapply(.SD, sum, na.rm=TRUE),by=c('nhser_name','stp_name','epiweek','ageband'),
                         .SDcols=c('s_positive_adj1','s_negative_adj1','s_na_adj1')]
  
}

x16 <- x16[x16$epiweek>=w1 & x16$epiweek<=w2,]

stp2reg <- x16[,lapply(.SD, sum, na.rm=TRUE),by=c('nhser_name','stp_name'),.SDcols=c('s_positive_adj1')]
regions <- stp2reg$nhser_name

#genome data

genomes <- read.csv(here("data/B.1.1.7-weightedcounts_nhsregion_date-2021-02-08.csv"),stringsAsFactors=F)
genomes$epiweek <- epiweek(genomes$date)
genomes$epiweek[genomes$epiweek<30] <- 53 + genomes$epiweek[genomes$epiweek<30]
genomes$nhser_name=as.factor(genomes$region)
genomes <- genomes[genomes$epiweek>=w1 & genomes$epiweek<=w2,]


## prepare stan data #############################
stan_data <- list()
## spatial dimensions and indices
N_reg <- length(unique(regions))
# region parameter
stan_data$N_reg <- N_reg
# genome data
stan_data$N_reg_g <- N_reg
# STPs
N_i <- length(regions)
stan_data$N_i <- N_i
# STP map to region - model
stan_data$reg_ind <- as.integer(regions)
# STP map to region - genome
stan_data$reg_ind_g <- as.integer(regions)
## overwrite for data d=3
stan_data$reg_ind_g_d3 <- 1:length(unique(regions))
stan_data$reg_ind_d3 <- 1:length(unique(regions))
stan_data$N_i_d3 <- length(unique(regions))
## time dimensions and indices
N_t <- length(wks)
stan_data$N_t <- N_t
# time parameters
N_t_par <- N_t-1
stan_data$N_t_par <- N_t-1
stan_data$t_ind <- 1:N_t
## age dimensions and indices
N_age_par <- N_age <- length(agebands0)
stan_data$N_age <- stan_data$N_age_par <- N_age
stan_data$a_ind <- 1:(N_age+1)
## overwrite for data d=3
stan_data$N_age_d3 <- stan_data$N_age_par_d3 <- 1
stan_data$a_ind_d3 <- rep(1,1+1)

## aggregate SGTF data
stan_data$Y_Snegative <- round(xtabs(s_negative_adj1 ~ epiweek + stp_name + ageband, x16))
stan_data$Y_Spositive <- round(xtabs(s_positive_adj1 ~ epiweek + stp_name + ageband, x16))
stan_data$Y_Sna <- xtabs(s_na_adj1 ~ epiweek + stp_name + ageband, x16)

## overwrite for data d=3
stan_data$Y_Snegative_d3 <- array(round(xtabs(s_negative_adj1 ~ epiweek + nhser_name, x16)),dim=c(N_t,N_reg,1),dimnames=list(wks,levels(regions),list('all')))
stan_data$Y_Spositive_d3 <- array(round(xtabs(s_positive_adj1 ~ epiweek + nhser_name, x16)),dim=c(N_t,N_reg,1),dimnames=list(wks,levels(regions),list('all')))
stan_data$Y_Sna_d3 <- array(xtabs(s_na_adj1 ~ epiweek + nhser_name, x16),dim=c(N_t,N_reg,1),dimnames=list(wks,levels(regions),list('all')))

## aggregate genome data
stan_data$genomes_voc <- round(xtabs(B.1.1.7 ~ epiweek + nhser_name , genomes))
stan_data$genomes_nonvoc <- round(xtabs(Other ~ epiweek + nhser_name , genomes))

## stan code flags
stan_data$prior_only <- F
stan_data$Tonly <- 0
stan_data$infTs <- 1 # infer generation time Tg and CV s
stan_data$Tminus <- 0 # set Tg for S- equal to that for S+?
stan_data$NB <- 1 # fix NB var ~ mean?
stan_data$Gonly <- 1 # use only genome data?

## parameter values
stan_data$mean_generation_time_T <- 6.4
stan_data$log_Tg_sd <- 0.5 # gives SD in Tg of 1
stan_data$sq_CV_generation_time_s <- 0.4
stan_data$log_s_sd <- 0.04  # gives SD in CV of Tg of 0.04


## initial values
init_fun <- function(...) list(b_t=array(-0.0,dim=c(N_t_par,N_reg,N_age_par)),
                               b_std=0.01,b_t_std=0.01, 
                               s_var=0.0,Tg_var=0.0,log_T_ratio=0.0,
                               gamma=0.01, reciprocal_phi=100,
                               firstobspos=array(-0.0,dim=c(N_i,N_age)),
                               firstobsneg=array(-0.0,dim=c(N_i,N_age)),
                               rti=array(-0.0,dim=c(N_t-1,N_i,N_age)))


## omit covariates ###################################################################
# construct models by omitting covariates

stan_data0 <- stan_data
fileName <- here('growth_factor_model/growthfactor_model.stan')

modnames <- c('noageortimeorreg','notimeorreg','noage','all')

stan_data_list <- list()
for(md in modnames) stan_data_list[[md]] <- stan_data0

## no age models
for(md in modnames[c(1,3)]){
  stan_data_list[[md]]$a_ind <- rep(1,N_age+1)
  stan_data_list[[md]]$N_age_par <- 1
}

## no time or reg models
for(md in modnames[c(1,2)]){
  stan_data_list[[md]]$N_reg <- 1
  stan_data_list[[md]]$reg_ind <- rep(1,stan_data0$N_i)
  stan_data_list[[md]]$reg_ind_d3 <- rep(1,stan_data0$N_i_d3)
  stan_data_list$notimeorreg$N_t_par <- 1
  stan_data_list$notimeorreg$t_ind <- rep(1,stan_data0$N_t)
}

## specify models #############################################################
# 1. no time, age, region or y
# 2. no time, age, region but with y
# 3. as 1 but with age
# 4. as 2 but with age
# 5. time and region, but no age or y
# 6. time and region, no age but with y
# 7. time, region and age, no y
# 8. time, region and age, with y

model_list <- list()
for(md in 1:length(modnames)){
  for(tt in 0:1){ # whether generation time for VOC will be different or not
    index <- (md-1)*2 + tt + 1
    ab <- agebands0
    if(md%in%c(1,3)) ab <- '11-120'
    model_list[[index]] <- list(modname=modnames[md],tt=tt,agebands=ab) 
  }
}

## start ######################################################################
stan_model <- rstan::stan_model(
  file = fileName,
  model_name= 'stan_model'
)

mod2run=1:8

for(d in 1:3){
  for(mi in mod2run){
    # read in model-specific variables
    for(mm in 1:length(model_list[[mi]])) assign(names(model_list[[mi]])[mm],model_list[[mi]][[mm]])
    # update stan data 
    stan_data <- stan_data_list[[modname]]
    stan_data$Gonly <- d-2
    stan_data$Tminus <- tt
    # replace if genome only - does not use age or STP geography
    if(stan_data$Gonly==1) 
      for(v in c('a_ind','N_age','N_age_par','N_i','reg_ind','reg_ind_g','Y_Spositive','Y_Snegative','Y_Sna')) 
        stan_data[[v]] <- stan_data[[paste0(v,'_d3')]]
    # update dimensions for init
    N_t_par <- stan_data$N_t_par; N_age_par <- stan_data$N_age_par; N_reg <- stan_data$N_reg; N_i <- stan_data$N_i; N_age <- stan_data$N_age
    # create file name from variables
    rname <- paste0(tbit[tt+1],modname,databit[d])
    fn <- paste0(here("results/"),rname,paste(range(wks),collapse=''),'.Rds')
    # run stan, save
    resStan <- sampling(stan_model, data = stan_data,init=init_fun,
                        chains = chains,cores=cores, iter = iters, thin = thin,
                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                        save_warmup = FALSE)
    saveRDS(list(resStan=resStan,stan_data=stan_data),fn)
    # compute LOO-CV, save
    loo_1 <- loo(resStan, cores = 2, moment_match=T)
    saveRDS(loo_1,paste0('loo',fn))
    
    ## results ################################################################
    ## set up variables to use
    outfile.base <- rname
    ##
    ## build labels
    ##
    lbls <- expand.grid(week=wks,
                        loc=levels(regions),
                        age_band=unlist(sapply(agebands0,paste,collapse='-')),stringsAsFactors=F)#,
    lbls$loc_idx <- match(lbls$loc,unique(lbls$loc))
    lbls$week_idx <- match(lbls$week,unique(lbls$week))
    lbls$age_idx <- match(lbls$age_band,unique(lbls$age_band))

    ## posterior predictive ################################################
    ## posterior predictive check on S- S+ 
    po <- rstan::extract(resStan)
    pp <- as.data.table(reshape2::melt( po$S_negative_pr) )
    setnames(pp, 2:4, c('week_idx','loc_idx','age_idx'))
    pp[, variable:= 'Snegative']
    tmp <- as.data.table(reshape2::melt( po$S_all_pr) )
    setnames(tmp, 2:4, c('week_idx','loc_idx','age_idx'))
    tmp[, variable:= 'Total cases']
    pp <- rbind(pp,tmp)
    pp <- pp[, list(value=quantile(value, p=c(.5,.025,.975)), stat=c('PP_M','PP_CL','PP_CU')), by=c('week_idx','loc_idx','age_idx','variable')]
    pp <- dcast.data.table(pp, variable+age_idx+loc_idx+week_idx~stat, value.var='value')
    tmp <- as.data.table(reshape2::melt( stan_data$Y_Snegative) )
    setnames(tmp, 1:4, c('week_idx','loc_idx','age_idx','observed'))
    tmp[, variable:= 'Snegative']
    tmp2 <- as.data.table(reshape2::melt( stan_data$Y_Snegative+stan_data$Y_Spositive+stan_data$Y_Sna) )
    setnames(tmp2, 1:4, c('week_idx','loc_idx','age_idx','observed'))
    tmp2[, variable:= 'Total cases']
    tmp <- rbind(tmp, tmp2)
    
    pp <- cbind(tmp,pp[,5:7])

    print(subset(pp, observed<PP_CL | observed>PP_CU))
    
    print(pp[, list(PP_OBS_OUT_CRI = mean(observed<PP_CL | observed>PP_CU)), by='week_idx'])
    
    ## time share #########################################################
    if(d<3&length(agebands0)==4){
      {
        # load age population counts
        dp <- as.data.table(read.csv('data/ukmidyearestimates20192020ladcodes.csv', stringsAsFactors = FALSE))
        dp <- subset(dp, select=-All.ages)
        setnames(dp, 4:94, paste0('age-',0:90))
        setnames(dp, 'Geography1', 'Geography')
        dp <- melt(dp, id.vars=c('Code','Name','Geography'), variable.name = 'age', value.name = 'pop_count')
        set(dp, NULL, 'age', dp[,as.integer(gsub('age-','',age))])
        set(dp, NULL, 'pop_count', dp[,as.integer(gsub(',','',pop_count))])
        
        # make age bands
        tmp <- data.table(age_from= 19L, age_to=49L)
        tmp2 <- data.table(age_from= 0L, age_to=10L)
        tmp <- rbind(tmp2, tmp)
        tmp2 <- data.table(age_from= 11L, age_to=18L)
        tmp <- rbind(tmp2, tmp)
        tmp2 <- data.table(age_from= 50L, age_to=99L)
        tmp <- rbind(tmp2, tmp)
        tmp <- tmp[, list(age=seq.int(age_from, age_to)), by=c('age_from','age_to')]
        tmp[, age_band:= paste0(age_from,'-',age_to)]
        
        tmp <- tmp[, list(age=seq.int(age_from, age_to)), by=c('age_from','age_to')]
        tmp[, age_band:= paste0(age_from,'-',age_to)]
        dp <- merge(dp, tmp, by='age')
        dp <- dp[, list(pop_count=sum(pop_count)), by=c('Code','Name','Geography','age_from','age_to','age_band')]
        
        # aggregate to STP areas
        tmp <- as.data.table(read.csv(file.path('data/sgtf_wes.csv')))
        tmp <- unique(subset(tmp, select=c(LTLA_code, LTLA_name, STP_code, STP_name )))
        tmp$LTLA_code %in% dp$Code #  cannot map all LTLAs to pop count data.table
        setnames(tmp, 'LTLA_code', 'Code')
        tmp <- subset(tmp, select = c(Code, STP_code, STP_name))
        dp2 <- merge(tmp, dp, by='Code')
        dp2 <- dp2[, list(pop_count=sum(pop_count)), by=c('STP_code','STP_name','age_from','age_to','age_band')]
        
        # make population fractions
        tmp <- dp2[, list(pop_total=sum(pop_count)), by='STP_code']
        dp2 <- merge(dp2, tmp, by='STP_code')
        dp2[, pop_prop := pop_count / pop_total]
        set(dp2, NULL, 'STP_name', dp2[, as.character(STP_name)])
        
        # make reference population
        dpref <- subset(dp, Name=='UNITED KINGDOM')
        tmp <- dpref[, list(pop_total=sum(pop_count)), by='Code']
        dpref <- merge(dpref, tmp, by='Code')
        dpref[, pop_prop := pop_count / pop_total]
        setnames(dpref, c('pop_count', 'pop_total', 'pop_prop'), c('popref_count', 'popref_total', 'popref_prop'))
        
      }
      ##
      ## make share of age groups over time
      ##
      # dp2 <- readRDS('data/populationbyage.Rds')
      dp2$loc <- stp_to_reg$nhser_name[match(dp2$STP_name,stp_to_reg$stp_name)]
      dp3 <- dp2[,lapply(.SD,sum,na.rm=T),by=.(loc,age_band),.SDcols=('pop_count')]
      dp3[,pop_total:=sum(pop_count),by=loc]
      # dpref <- readRDS('data/referencepopulationbyage.Rds')
      dp4 <- dp3[dpref,on='age_band']
      dp4$age_band[dp4$age_band=='50-99'] <- '50-120'
      dp4$age_to[dp4$age_to==99] <- 120
      dp4[,reg_prop:=pop_count/pop_total]
      
      dsa <- as.data.table(reshape2::melt(po$expected_Snegative)) # pi_negative)) # 
      setnames(dsa, 2:4, c('week_idx','loc_idx','age_idx'))
      dsa[, variable := 'S_negative']
      tmp <- as.data.table(reshape2::melt(po$expected_Spositive)) # pi_positive)) # 
      setnames(tmp, 2:4, c('week_idx','loc_idx','age_idx'))
      tmp[, variable := 'S_positive']
      dsa <- rbind(dsa, tmp)
      dsa <- merge(dsa, lbls, by=c('week_idx','loc_idx','age_idx'))
      dsa <- dsa[dp4,on=c('age_band','loc')]
      dsa[,foi:=value/pop_count]
      dsa[,cases := pop_total*popref_prop*foi]
      dsa[,total_reg_cases:=sum(cases),by=c('loc','week','iterations','variable')]
      dsa[,value:=cases/total_reg_cases]
      dsa <- dsa[!is.na(dsa$value),]
      # dsa[,value:=value/reg_prop*popref_prop] # shortcut using pi
      
      dsa <- dsa[, list(value=quantile(value, p=c(.5,.025,.975)), stat=c('M','CL','CU')), by=c('week_idx','loc_idx','age_idx','variable')]
      dsa <- dcast.data.table(dsa, week_idx+loc_idx+age_idx+variable~stat, value.var='value')
      
      dsa <- merge(dsa, lbls, by=c('week_idx','loc_idx','age_idx'))
      dp <- dsa#subset(dsa, loc%in%c('London',"East of England","South East"))
      p <- ggplot(dp, aes(x=lubridate::ymd( "2020-01-01" ) + lubridate::weeks(week) - 1 )) +
        geom_ribbon(data=subset(dp, variable=='S_positive'), aes(ymin=CL, ymax=CU, fill=age_band, group=interaction(variable,age_band)), alpha=0.3) +
        geom_ribbon(data=subset(dp, variable=='S_negative'), aes(ymin=CL, ymax=CU, fill=age_band, group=interaction(variable,age_band)), alpha=0.3) +
        geom_line(aes(y=M, linetype=factor(variable, levels=c('S_positive','S_negative'), labels=c('S+','S-')), group=interaction(variable,age_band)), show.legend= FALSE) +
        geom_point(aes(y=M, fill=age_band, pch=factor(variable, levels=c('S_positive','S_negative'), labels=c('S+','S-'))), size=2, colour='black') +
        ggsci::scale_colour_npg() +
        ggsci::scale_fill_npg() +
        scale_x_date(breaks='2 weeks',expand=c(0,0)) +
        scale_y_continuous(labels= scales::percent, expand=c(0,0), lim=c(0,.8), breaks=seq(0,1,.2)) +
        coord_cartesian(ylim=c(0,0.8)) +
        scale_linetype_manual(values=c('S-'='11','S+'='solid')) +
        scale_shape_manual(values=c('S-'=23, 'S+'=21)) +
        theme_bw() +
        labs(x='', y='Proportion of all cases within age band', colour='Age band', fill='Age band', pch='') +
        facet_grid(~loc) +
        guides(colour=FALSE) +
        theme(text=element_text(family="sans"),
              legend.position='bottom',
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              panel.spacing = unit(2, "lines"),
              strip.background = element_blank()
        )
      p
      ggsave(file=file.path(outdir, paste0(outfile.base,'std_shareage_time.pdf')), p, limitsize = FALSE, w=15, h=6)
      
      
      
    }
  }
  
  
  ## plot two advantages ##############################
  ## once all d=1,2 models have run, put two from d=1 and d=2 side by side
  if(d==2){
    ##
    ## make transmission advantage over time
    ##
    # select models
    rlistsgss <- readRDS('./noagesgss4555.Rds')
    rlistboth <- readRDS('./noageboth4555.Rds')
    ## get first model advantage
    resStan <- rlistsgss[[1]]
    po <- rstan::extract(resStan)
    dsa <- as.data.table(reshape2::melt(po$b_t))
    setnames(dsa, 2:4, c('week_idx','loc_idx','age_idx'))
    set(dsa, NULL, 'value', dsa[,exp(value)])
    dsa <- dsa[, list(value=quantile(value, p=c(.5,.025,.975)), stat=c('M','CL','CU')), by=c('week_idx','age_idx')]
    dsa <- dcast.data.table(dsa, week_idx+age_idx~stat, value.var='value')
    dsa <- merge(dsa, lbls, by=c('week_idx','age_idx'))
    dsa$dataset <- 'SGTF'
    ## get second model advantage
    resStan <- rlistboth[[1]]
    po <- rstan::extract(resStan)
    dsa2 <- as.data.table(reshape2::melt(po$b_t))
    setnames(dsa2, 2:4, c('week_idx','loc_idx','age_idx'))
    set(dsa2, NULL, 'value', dsa2[,exp(value)])
    dsa2 <- dsa2[, list(value=quantile(value, p=c(.5,.025,.975)), stat=c('M','CL','CU')), by=c('week_idx','age_idx')]
    dsa2 <- dcast.data.table(dsa2, week_idx+age_idx~stat, value.var='value')
    dsa2 <- merge(dsa2, lbls, by=c('week_idx','age_idx'))
    dsa2$dataset <- 'SGTF+genome'
    dsa <- rbind(dsa,dsa2)
    ## plot
    p <- ggplot(dsa, aes(x=week)) +
      geom_ribbon(aes(ymin=CL, ymax=CU), fill='navyblue', alpha=.4) +
      geom_line(aes(y=M), colour='navyblue') +
      geom_point(aes(y=M), colour='navyblue') +
      ggsci::scale_colour_npg() +
      ggsci::scale_fill_npg() +
      scale_x_continuous(expand=c(0,0),breaks=wks[seq(1,length(wks),by=2)]) +
      labs(x='Week',y='VOC transmission advantage') +
      theme_bw(base_size = 15) + facet_wrap(~dataset) +
      theme(legend.position = '')
    p
    ggsave(file=file.path(outdir, paste0('sgssandboth_VOCtrmadvtime_overall.pdf')), p, limitsize = FALSE, w=8, h=4)
    
    
    
    
  }

}




## loo ################################################################
# http://mc-stan.org/loo/articles/loo2-with-rstan.html
if(0){
  resdir <- here("results/")
  regex_pars = "b_t|reciprocal_phi|rti|log_growth_Spos_t_sd|b_std|b_t_std|log_T_ratio|s_var|Tg_var"
  for(d in 1:3){
  mod2run <- 1:8
  if(d==3) mod2run <- c(1,2,5,6)
  tt <- 1 # 0=no T minus, 1=T minus = y T plus
  loos <- loo_list <- waics <- Ts <- bs <- list()
  for(i in 1:length(mod2run)){
    # extract model values
    mi <- mod2run[i]
    for(mm in 1:length(model_list[[mi]])) assign(names(model_list[[mi]])[mm],model_list[[mi]][[mm]])
    # make file name
    rname <- paste0(tbit[tt+1],modname,databit[d])
    fn <- paste0(resdir,rname,paste(range(wks),collapse=''),'.Rds')
    # read in
    resStanlist <- readRDS(fn)
    resStan <- resStanlist$resStan
    # print diagnostics
    su <- summary(resStan)$summary
    su <- su[grepl(regex_pars,rownames(su)), ]
    print(c(mi,min(su[,'n_eff']),max(su[,'Rhat'])))
    log_lik_1 <- extract_log_lik(resStan, merge_chains = FALSE)
    # r_eff <- relative_eff(exp(log_lik_1), cores = 2)
    # read in LOO
    fn <- paste0(resdir,'loo',rname,paste(range(wks),collapse=''),'.Rds')
    loo_list[[i]] <- readRDS(fn)
    loos[[i]] <- loo_list[[i]]$estimates[1,1]
    # compute WAIC
    waics[[i]] <- waic(log_lik_1)
    # extract posteriors
    Tqs <- round(quantile(extract(resStan, pars = c('T_ratio'), permuted = TRUE, inc_warmup = FALSE,include = TRUE)$T_ratio,c(0.5,0.025,0.975)),2)
    Ts[[i]] <-  paste0(Tqs[1],' (',Tqs[2],',',Tqs[3],')')
    bqs <- round(quantile(exp(extract(resStan, pars = c('b_t'), permuted = TRUE, inc_warmup = FALSE,include = TRUE)$b_t),c(0.5,0.025,0.975)),2)
    bs[[i]] <-  paste0(bqs[1],' (',bqs[2],',',bqs[3],')')
  }
  ## print
  print(loo_list)
  
  comp <- loo_compare(loo_list)
  print(round(apply(comp,1,function(x)pnorm(0,x[1],x[2],lower=F)),3))
  comp2 <- loo_compare(waics)
  print(round(apply(comp2,1,function(x)pnorm(0,x[1],x[2],lower=F)),3))
  
  mdnames <- paste0('model',mod2run)
  print(cbind(T_ratio=unlist(Ts),
        Effect_size=unlist(bs),
        loo=round(unlist(loos)),
        waic=round(sapply(waics,function(x)x$estimates[1,1]))))#,
  }
  
}



