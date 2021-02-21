library(dplyr)
library(pracma)
library(rstan)
library(lubridate)

maxweek = 56

outdir = '.'
outfile.base='d4'


## mcmc parms #################################################

fileName <- 'd4.stan'
iters <- 8000
cores <- 6
chains <- 2
thin <- 2

## read #######################################################


x10 <- read.csv('../data/B.1.1.7-weightedcounts_nhsregion_date-2021-02-08.csv',stringsAs=F)
colnames(x10)[colnames(x10)=='B.1.1.7'] <- 'B117'

x10$date <- as.Date( x10$date )
x10$epiweek <- epiweek( x10$date )
x10$epiweek <- ifelse( year(x10$date) == 2021 & x10$epiweek<53, 53 + x10$epiweek, x10$epiweek)
x10$day <- as.numeric(round(difftime(x10$date,min(x10$date),units='days')+1))

x10 <- subset(x10,epiweek>42)
x10 <- subset(x10,epiweek<=maxweek)

cat( 'Sample size \n' )
print( with( x10, c( sum( B117 ) , sum( Other )) ) ) 


## dataset ################################################

areas <- unique(x10$region)
w1 <- 43
w2 <- max( x10$epiweek )
wks <- w1:w2

regcds <- unique(x10$region)
x10$total <- x10$B117 + x10$Other

byweek <- expand.grid(regcd=regcds,time=wks,stringsAsFactors = F)
byweek$B117  <- apply(byweek,1,function(x)sum(subset(x10,epiweek==x[2]&region==x[1])$B117))
byweek$total <- apply(byweek,1,function(x)sum(subset(x10,epiweek==x[2]&region==x[1])$total))
byweek$epiweek <- byweek$time

d8 <- byweek 


#################################################################

stan_data <- list()
stan_data$N_w <- length(unique(d8$epiweek))
stan_data$N_i <- length(unique(d8$regcd))
stan_data$N_t <- length(unique(d8$time))
time_to_week <- sapply(1:stan_data$N_t,function(x)which(sort(unique(d8$epiweek))==subset(d8,time==unique(d8$time)[x])$epiweek[1]))
stan_data$t_to_w <- time_to_week
stan_data$i_ind <- match(d8$regcd,unique(d8$regcd))
stan_data$t_ind <- match(d8$time,unique(d8$time))
stan_data$N_B117 <- nrow(d8)
stan_data$N_total <- nrow(d8)
stan_data$Y_B117 <- sapply(unique(d8$regcd),function(x)
  sapply(unique(d8$time),function(y){
    sub=subset(d8,regcd==x&time==y)
    ifelse(nrow(sub)==0,0,round(sub$B117))
    }))
stan_data$Y_total <- sapply(unique(d8$regcd),function(x)
  sapply(unique(d8$time),function(y){
    sub=subset(d8,regcd==x&time==y)
    ifelse(nrow(sub)==0,0,round(sub$total))
    }))

stan_data$first_offset <- -7 
stan_data$prior_only <- F

init_fun <- function(...) list(firstobs=array(-0.0,dim=c(stan_data$N_i)),
                               rho_0=array(-0.0,dim=c(stan_data$N_w-1,stan_data$N_i)),
                               overall_rho_0=array(-0.0,dim=c(stan_data$N_w-1)))


stan_model_freq <- rstan::stan_model(
  file = fileName,
  model_name= 'stan_model_freq'
)


resStan <- sampling(stan_model_freq, data = stan_data,init=init_fun,
                chains = chains,cores=cores, iter = iters, thin = thin,
                control = list(adapt_delta = 0.99, max_treedepth = 15)
)




ps <- t(apply((extract(resStan, pars = c('overall_rho'), permuted = TRUE, inc_warmup = FALSE,include = TRUE)$overall_rho),2,quantile,c(0.5,0.025,0.975)))
print( 'rho:' )
print( cbind(week=sort(unique(d8$epiweek))[-1],round(ps,2)) )

freq <- c(sapply(1:stan_data$N_w,function(x)
  apply(extract(resStan, pars = c('freq'), permuted = TRUE, inc_warmup = FALSE,include = TRUE)$freq[,x,],2,mean)))
freqdat <- data.frame(pred=freq,
                      epiweek=rep(sort(unique(d8$epiweek)),each=stan_data$N_i),
                      regcd=rep(unique(d8$regcd),length((unique(d8$epiweek)))),
                      stringsAsFactors=F)


dat <- as.data.frame(cbind(time=sort(unique(d8$time))[-1],ps))
dat$median <- as.numeric(dat$`50%`)
dat$upper <- as.numeric(dat$`97.5%`)
dat$lower <- as.numeric(dat$`2.5%`)
p0 <- ggplot(dat,aes(x=time)) + 
  geom_line(aes(y=median),size=1) + theme_minimal(base_size = 15) +
  ylab('Growth difference (1/week)') + geom_ribbon(aes(ymin=lower, ymax=upper), 
              alpha=0.2, size=0) + ylim(0,max(dat$upper)+.01)

ggsave(p0,file=file.path(outdir, paste0(outfile.base,'rho.png')),width=4,height=3)
ggsave(p0,file=file.path(outdir, paste0(outfile.base,'rho.pdf')),width=4,height=3)



## checks ######################################################
su <- summary(resStan)$summary
trms <- "rho|overall_rho|firstobs|freq|B_count"
su <- su[grepl(trms,rownames(su)), ]
min(su[,'n_eff'])
max(su[,'Rhat'])
which(su[,'n_eff']<1000)

su <- su[grepl(trms,rownames(su)), ]
write.csv(su, file=file.path(outdir, paste0(outfile.base,'_margposu.csv')))

##
## build labels
##
lbls <- unique(subset(d8, select=c( epiweek, regcd)))
lbls$region <- match(lbls$regcd,colnames(stan_data$Y_B117))
lbls$rgn <- colnames(stan_data$Y_B117)
lbls$rgn[lbls$rgn=='Yorkshire_and_The_Humber'] <- 'Yorkshire'
lbls$rgn <- gsub('_',' ',lbls$rgn)
lbls$week <- lbls$epiweek + 1 - min(lbls$epiweek)


##
## make growth diff over time
##
library(data.table)

po <- rstan::extract(resStan)
dsa <- as.data.table(reshape2::melt(po$rho))
setnames(dsa, 2:3, c('week','region'))
set(dsa, NULL, 'value', dsa[,(value)])
dsa <- dsa[, list(value=quantile(value, p=c(.5,.025,.975)), stat=c('M','CL','CU')), by=c('week','region')]
dsa <- dcast.data.table(dsa, week+region~stat, value.var='value')
lbls1 <- lbls; lbls1$week <- lbls1$week-1
dsa <- merge(dsa, lbls1, by=c('week','region'))
wksub <- wks[c(2,5,8,11)]
p1 <- ggplot(dsa, aes(x=epiweek)) +
  geom_ribbon(aes(ymin=CL, ymax=CU), alpha=.5) +
  geom_line(aes(y=M)) +
  geom_point(aes(y=M)) +
  ggsci::scale_colour_npg() +
  ggsci::scale_fill_npg() +
  facet_wrap(~rgn) +
  scale_x_continuous(expand=c(0,0)) +
  labs(x='Week',y='Growth difference (1/week)') +
  scale_x_continuous('Week', wksub, wksub, range(wksub)) +
  theme_bw(base_size = 14) +
  theme(legend.position = 'bottom')# +
p1
ggsave(file=file.path(outdir, paste0(outfile.base,'_VOCtrmadvtime.png')), p1, limitsize = FALSE, w=5, h=5)
ggsave(file=file.path(outdir, paste0(outfile.base,'_VOCtrmadvtime.pdf')), p1, limitsize = FALSE, w=5, h=5)
#facet_grid(~loc)
#ggsave(file=file.path(outdir, paste0(outfile.base,'_VOCtrmadvtime.pdf')), p, limitsize = FALSE, w=12, h=5)
#ggsave(file=file.path(outdir, paste0(outfile.base,'_VOCtrmadvtime.png')), p, limitsize = FALSE, w=12, h=5)


## posterior predictive check
po <- rstan::extract(resStan)
pp <- as.data.table(reshape2::melt( po$B_count) )
setnames(pp, 2:3, c('time','region'))
pp <- pp[, list(value=quantile(value, p=c(.5,.025,.975)), stat=c('PP_M','PP_CL','PP_CU')), by=c('time','region')]
pp <- dcast.data.table(pp, region+time~stat, value.var='value')
tmp <- as.data.table(reshape2::melt( unname(stan_data$Y_B117)) )
setnames(tmp, 1:3, c('time','region','observed'))
pp <- merge(pp,tmp,by=c('time','region'))

print( 'posterior predictive check' )
print(
	round(nrow(subset(pp, !observed<PP_CL & !observed>PP_CU))/nrow(pp)*100)
)




##
##  posterior predictive scatter plot
##
p <- ggplot(pp, aes(x=observed, y=PP_M)) +
  geom_errorbar(aes(ymin=PP_CL, ymax=PP_CU), colour='grey80') +
  geom_point(aes(colour=factor(region))) +
  geom_abline(slope=1, intercept=0) +
  ggsci::scale_colour_npg() +
  scale_x_log10(expand=c(0,0)) +
  scale_y_log10(expand=c(0,0)) +
  #facet_grid(age_band~variable) +
  theme_bw() +
  theme(legend.position='') +
  labs(x='Observed',y='Predicted',colour='region')
ggsave(file=file.path(outdir, paste0(outfile.base,'rhoscatter.png')),p,width=5,height=5)




##  frequency plot ###############################################

dpred = d8 
freq <- c(sapply(1:stan_data$N_w,function(x)
  apply(extract(resStan, pars = c('freq'), permuted = TRUE, inc_warmup = FALSE,include = TRUE)$freq[,x,],2,mean)))
freqdat <- data.frame(pred=freq,epiweek=rep(wks,each=stan_data$N_i),regcd=rep(unique(d8$regcd),length(wks)),stringsAsFactors=F)
dpred <- left_join(dpred,freqdat,by=c('epiweek','regcd'))
dpred$raw <- dpred$B117/dpred$total

dpred$region = gsub( d8$regcd, patt = '_', rep = ' ' )
dpred$region [ grepl( dpred$region, patt = 'Yorkshire') ] <- 'Yorkshire' 
dpred$date <- as.Date( '2020-01-01' ) + (dpred$epiweek-1)*7 

dpred$flb = as.vector(sapply(1:stan_data$N_w,function(x)
  apply(extract(resStan, pars = c('freq'), permuted = TRUE, inc_warmup = FALSE,include = TRUE)$freq[,x,],2
  ,function(x) as.numeric( quantile(na.omit(x), prob=.025 ) )
  )))
dpred$fub = as.vector(sapply(1:stan_data$N_w,function(x)
  apply(extract(resStan, pars = c('freq'), permuted = TRUE, inc_warmup = FALSE,include = TRUE)$freq[,x,],2
  ,function(x) as.numeric( quantile(na.omit(x), prob=.975 ) )
  )))
dpred <- na.omit( dpred )
rawci = with( dpred, {
	t(sapply( 1:nrow(dpred), function(k) {
		bt = binom.test( round(B117[k]), round(total[k]), 1 ) 
		as.numeric( bt$conf.int ) 
	}))
})
colnames( rawci ) <- c( 'rawmin', 'rawmax' )
dpred <- cbind( dpred, rawci ) 

saveRDS( list(
  drho = dat 
  , dpred = dpred
  , drho2 = dsa 
), file = paste0( outfile.base, '-results.rds' ) 
)

dg = 2
p = ggplot( data=dpred, aes( x=date, y = pred, ymin = flb, ymax = fub , colour = region, fill=region, pch = region  )  )  + geom_path() + 
geom_point( aes( date , raw, colour = region), position=position_dodge(dg)  )  + theme_classic()  + geom_ribbon( alpha = .25, colour = NA )
p = p + geom_errorbar( aes(x=date, ymin = rawmin, ymax = rawmax ),  width=.2, position=position_dodge(dg))
p = p + xlab('') + ylab( 'Estimated frequency VOC' ) + theme(legend.title = element_blank()) + scale_shape_manual(values=seq(0,length( unique( dpred$region))))
ggsave( p , file=file.path(outdir, paste0(outfile.base,'estfreq.png')),width=4.75,height=3.5)
ggsave( p , file=file.path(outdir, paste0(outfile.base,'estfreq.pdf')),width=4.75,height=3.5)

pldf <- merge( dpred, dsa, by = c('epiweek', 'regcd' ))
pldf$frequency = pldf$pred 
pldf$rho =pldf$M
p3 = ggplot( aes(x = frequency,y = rho, colour = regcd), data = pldf ) + geom_point() + geom_path() + theme_minimal() + theme( legend.title=element_blank() ) + xlab('Sample frequency') + ylab('Growth rate difference (1/week)')
ggsave( p3 , file=file.path(outdir, paste0(outfile.base,'freqXrho.png')),width=4.75,height=3.5)
ggsave( p3 , file=file.path(outdir, paste0(outfile.base,'freqXrho.pdf')),width=4.75,height=3.5)

# corresponding plot by time 
p4 = ggplot( aes(x = epiweek,y = rho, colour = regcd), data = pldf ) + geom_point() + geom_path() + theme_minimal() + theme( legend.title=element_blank() ) + xlab('Week') + ylab('Growth rate difference (1/week)')
ggsave( p4 , file=file.path(outdir, paste0(outfile.base,'weekXrho.png')),width=4.75,height=3.5)
ggsave( p4 , file=file.path(outdir, paste0(outfile.base,'weekXrho.pdf')),width=4.75,height=3.5)

# check for frequency dependence 
print( summary( lm( rho ~ epiweek  , data = pldf )  )  )
print( summary( lm( rho ~ frequency  , data = pldf )  )  )
print( summary( lm( rho ~ frequency + epiweek  , data = pldf )  )  )

