library( lubridate ) 
library( ggplot2 )
library( mgcv ) 


drho  = readRDS( 'd4-results.rds' )$drho 
dpred = readRDS( 'd4-results.rds' )$dpred

# compute growthrate (1/day) for voc/non-voc 
## growth based on ML phylodyn
skyfit0 = readRDS('B.1.1.7_2021-02-13_n_tree_dating_10_dated_trees_mlesky.rds_mlesky.rds')
ne = skyfit0$ne
gr = apply( ne, MAR=2, FUN=function(x)  c( NA, diff( log(x) ) / (365*diff(skyfit0$time)[1]) )  ) # 1/days
## same for lineage 177 
skyfit = readRDS('B.1.177_2021-02-13_n_tree_dating_10_dated_trees_mlesky.rds_mlesky.rds')
mane = skyfit$ne 
magr = apply( mane, MAR=2, FUN=function(x)  c( NA, diff( log(x)) / (365*diff(skyfit$time)[1]))  ) # 1/days


## time limits for generation times
alen = 43 ## max days of truncated inf period  (days)
ataxis <- 1:alen 
si_tlen = alen/365
## gamma scale parameter for generation times 
gamscale <- (6.4* 0.444)^2 / 6.4 # 

## find dates corresponding to weeks, put all variables on same time axis 
drho$week <- drho$time 
.d <- as.Date( '2020-01-01' ) + 1:500 
.w <- epiweek( .d ) 
.w [  year(.d ) > 2020  &  .w<53 ] <- .w [  year(.d ) > 2020  &  .w<53 ] + 53
.ww = unique( .w )
.dd = lapply( .ww , function(w) max(.d[ .w == w]) ) # last day of each epiweek 
drho$date <-  as.Date( unlist( .dd[ match( drho$week ,  .ww )  ]  )  , origin = as.Date('1970-01-01' ) ) 
drho$time <- decimal_date ( drho$date )
drho[ , c('median', 'lower', 'upper') ] <- drho[ , c('median', 'lower', 'upper') ]/7 # change to daily growth rate 
drho$sd = with ( drho , (upper - lower ) / (2*1.96) )
rho_interp <- approxfun(  drho$time, drho$median, rule=2 )
rho_lb_interp <- approxfun(  drho$time, drho$lower, rule=2 )
rho_ub_interp <- approxfun(  drho$time, drho$upper, rule=2 )


#' compute selcoef (R ratio - 1) at given time
#' 
#' @param tt time 
#' @param rhoquantile quantile 0-1 of posterior growth rate difference 
#' @param Tg mean of generation time distribution 
time2selcoef <- function(tt, rhoquantile, Tg = NULL )
{
	# compute discretised gamma 
	gamshape = Tg / gamscale 
	si1 <- diff(pgamma( 0:alen, shape=gamshape, scale = gamscale ) ) 
	si1 = si1 / sum( si1 )
	
	magr_interpfun <- approxfun( skyfit$time, magr[, sample( ncol(magr), size=1) ] , rule=2)
	gr_interpfun <- approxfun( skyfit0$time, gr[, sample( ncol(gr), size=1) ] , rule=2)
	
	qrho = qnorm( rhoquantile, mean = drho$median, sd = drho$sd )
	rhofun = approxfun( drho$time, qrho, rule = 2)
	
	rataxis <- rev ( seq(  tt-si_tlen+1/365, tt, by = 1/365 ) ) 
	
	Rnv0= 1/sum( si1 * exp( -cumsum(magr_interpfun(rataxis)*1) ) * 1 )
	Rnv1= 1/sum( si1 * exp( -cumsum(gr_interpfun(rataxis)*1 - rhofun(rataxis)*1) ) * 1 )
	
	Rvoc0 =  1/sum( si1 * exp( -(cumsum(gr_interpfun(rataxis))*1) ) * 1 ) 
	Rvoc1 =  1/sum( si1 * exp( -(cumsum(magr_interpfun(rataxis) + rhofun(rataxis))*1) ) * 1 ) 
	
	s0 =  (Rvoc1/Rnv0 + Rvoc0/Rnv1)/2 - 1
	
	c( s0, Rvoc0, Rnv0, Rvoc1, Rnv1 )
}

# simulation replicates 
nrep = 3000
montecarlo_inputs = data.frame( rhoquantile = runif( nrep , 0, 1 )
 , Tg = rnorm( nrep, 6.4, sd = (4/(2*1.96) ))
)
time_axis = drho$time 
res0 =   lapply( 1:nrow( montecarlo_inputs ), function(k) sapply( time_axis, function(tt) time2selcoef( tt, montecarlo_inputs$rhoquantile[k], montecarlo_inputs$Tg[k]  ) ) ) 
a0 = array( do.call( c, res0 ) , dim = c( dim(res0[[1]]), length( res0 ) ) )

O = a0[1,,]

## estimated selcoef (ratio R - 1 ), computed over entire time axis 
print( quantile( O , c( .5, .025, .975 )) )

'
      50%      2.5%     97.5% 
0.8838412 0.4889253 1.4499576 

      50%      2.5%     97.5% 
0.8818959 0.4217651 1.5587428 

'





#~ --------------
# recompute to to test sensitivity to difference in generation times between voc and non-voc 
time2selcoef2 <- function(tt, rhoquantile=.5, Tg1 = 6.5, Tg2 = 6.5 )
{
	gamshape1 = Tg1 / gamscale 
	si1 = diff(pgamma( 0:alen, shape=gamshape1, scale = gamscale ) ) 
	si1 = si1 / sum( si1 )
	
	gamshape2 = Tg2 / gamscale 
	si2 = diff(pgamma( 0:alen, shape=gamshape2, scale = gamscale ) ) 
	si2 = si2 / sum( si2 )
	
	magr_interpfun <- approxfun( skyfit$time, magr[, sample( ncol(magr), size=1) ] , rule=2)
	gr_interpfun <- approxfun( skyfit0$time, gr[, sample( ncol(gr), size=1) ] , rule=2)
	
	qrho = qnorm( rhoquantile, mean = drho$median, sd = drho$sd )
	rhofun = approxfun( drho$time, qrho, rule = 2)
	
	rataxis <- rev ( seq(  tt-si_tlen+1/365, tt, by = 1/365 ) ) 
	
	s0 = sum( si1 * exp( -cumsum(magr_interpfun(rataxis)*1) ) * 1 ) /  sum( si2 * exp( -(cumsum(magr_interpfun(rataxis) + rhofun(rataxis))*1) ) * 1 )  - 1
	s1 = sum( si1 * exp( -cumsum(gr_interpfun(rataxis)*1 - rhofun(rataxis)*1 ) * 1 ) ) /  sum( si2 * exp( -(cumsum(gr_interpfun(rataxis) )*1) ) * 1 )  - 1
	(s0 + s1 ) / 2
}


fTg = 1-.25 # prop reduction in Tg in voc
O25 = sapply( 1:nrow( montecarlo_inputs ), function(k) sapply( time_axis, function(tt) time2selcoef2( tt, montecarlo_inputs$rhoquantile[k], Tg1=montecarlo_inputs$Tg[k], Tg2=fTg*montecarlo_inputs$Tg[k]  )[1] ) )
fTg50 = 1-.50 # prop reduction in Tg in voc
O50 = sapply( 1:nrow( montecarlo_inputs ), function(k) sapply( time_axis, function(tt) time2selcoef2( tt, montecarlo_inputs$rhoquantile[k], Tg1=montecarlo_inputs$Tg[k], Tg2=fTg50*montecarlo_inputs$Tg[k]  )[1] ) )


print( quantile( O25 , c( .5, .025, .975 )) )
print( quantile( O50 , c( .5, .025, .975 )) )


#~ ---------------------------------------------------------------------
# Figures 

library( glue ) 
cols = rev( c( '#cc993373' , '#3366cc64' )) # background
cols1 = rev( c( '#f8766dff' , '#00bfc4ff' ) ) #lines

###
# R ratio 
library( pammtools ) # for step ribbon 

spldf0  = as.data.frame( 1+ t( apply( O, 1, function(x) quantile(x, c(.5, .025, .975))  ) ) )
spldf = spldf0 
colnames( spldf ) <- c( 'median', 'lower', 'upper' )
spldf$Scenario = 'Baseline'
spldf$date <- as.Date( date_decimal( time_axis ) )

spldf1 = as.data.frame( 1+  t( apply( O25, 1, function(x) quantile(x, c(.5, .5, .5))  ) ) )
colnames( spldf1 ) <- c( 'median', 'lower', 'upper' )
spldf1$Scenario = glue('-{round((1-fTg)*100)}% Tg')
spldf1$date <- as.Date( date_decimal( time_axis ) )

spldf2 = as.data.frame( 1+  t( apply( O50 , 1, function(x) quantile(x, c(.5, .5, .5))  ) ) )
colnames( spldf2 ) <- c( 'median', 'lower', 'upper' )
spldf2$Scenario = glue('-{round((1-fTg50)*100)}% Tg')
spldf2$date <- as.Date( date_decimal( time_axis ) )

spldf <- rbind( spldf, spldf1,  spldf2 )

prr = ggplot( aes( x = date, y = median, ymin = lower, ymax = upper, colour = Scenario, fill=Scenario) ,  data = spldf ) + 
#~ geom_path(aes(colour=Scenario), lwd=1) + 
#~ geom_ribbon( lwd = 0, alpha = .25) + 
geom_step( aes(colour=Scenario), lwd=1) + 
geom_stepribbon(lwd = 0, alpha = .25) + 
theme_minimal() + 
ylab('Multiplicative advantage for VOC over non-VOC') + xlab('' ) + 
theme_classic() +
theme(legend.pos='top')

ggsave( prr, file = 'd5-ReffRatio.png', width = 12, height = 8 , units = 'cm' )
ggsave( prr, file = 'd5-ReffRatio.pdf', width = 12, height = 8 , units = 'cm' )


####
# R from phylodyn method 
Rvoc = a0[2,,]
Rnv = a0[3,,] 

Rvocdf = as.data.frame( t( apply( Rvoc, 1, function(x) quantile(x, c(.5, .025, .975))  ) ) )
colnames( Rvocdf ) <- c( 'median', 'lower', 'upper' )
Rvocdf$Lineage = 'B.1.1.7'
Rvocdf$date <- as.Date( date_decimal( time_axis ) )

Rnvdf = as.data.frame( t( apply( Rnv, 1, function(x) quantile(x, c(.5, .025, .975))  ) ) )
colnames( Rnvdf ) <- c( 'median', 'lower', 'upper' )
Rnvdf$Lineage = 'B.1.177'
Rnvdf$date <- as.Date( date_decimal( time_axis ) )

Rdf = rbind( Rvocdf , Rnvdf ) 

pr = ggplot( aes( x = date, y = median, ymin = lower, ymax = upper, colour = NULL, fill=Lineage, group=Lineage) ,  data = Rdf )  + 
geom_path(aes(colour=Lineage), lwd=0.5) + 
geom_ribbon( lwd = 0, alpha = 0.35) + theme_minimal() + ylab('Rt') + xlab('' ) + 
scale_color_manual( name=NULL, values = cols1 ) + 
scale_fill_manual( name='Lineage', values=cols ) + 
geom_hline( yintercept = 1, colour = 'red', lwd = .35) + 
guides(color = FALSE)  + 
theme(axis.text.x = element_text( color="#000000", size=11), axis.text.y = element_text( color="#000000",  size=11), legend.pos='none') + 
scale_y_log10() +
theme(legend.position='none')
ggsave( pr, file = 'd5-R.png', width = 6, height = 6 , units = 'cm' )
ggsave( pr, file = 'd5-R.pdf', width = 6, height = 6 , units = 'cm' )



###
# Ne(t) 
q_ne = t(apply( ne, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ))
q_mane  = t(apply( mane, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ))
colnames( q_ne ) = colnames( q_mane ) = c( 'y', 'ylb', 'yub' )
pldf0 = as.data.frame( q_ne ) ; pldf0$Lineage = 'B.1.1.7'; pldf0$time = skyfit0$time 
pldf1 = as.data.frame( q_mane ); pldf1$Lineage = 'B.1.177'; pldf1$time = skyfit$time 
pldf = rbind( pldf0, pldf1 )

pne = ggplot( aes(x = as.Date( date_decimal( time)), y = y, colour = Lineage, fill = Lineage , ymin = ylb, ymax = yub ) , data = na.omit( pldf ) ) + 
geom_path(lwd = .5) + 
geom_ribbon( alpha = .35 , lwd = 0) + 
guides(color = FALSE)  + 
scale_color_manual( name=NULL, values = cols1 ) + 
scale_fill_manual( name='', values=cols ) +
#~ scale_fill_manual( name='Lineage', values=cols ) +
xlab('') + ylab('Effective population size' ) + 
theme_classic()  + 
theme(axis.text.x = element_text( color="#000000", size=11), axis.text.y = element_text( color="#000000",  size=11), legend.pos='top') 

ggsave( pne, file = 'd5-Ne.png', width = 6, height = 6 , units = 'cm' )
ggsave( pne, file = 'd5-Ne.pdf', width = 6, height = 6 , units = 'cm')


