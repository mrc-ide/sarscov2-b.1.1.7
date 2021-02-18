# Maximum likeihood estimation of Ne(t) from trees dated with treedater
# Required: 
# - RDS files containing dated trees from d1_date_trees.R

# devtools::install_github('emvolz-phylodynamics/mlesky')


library( ape )
library( lubridate )
library( glue )
library( mlesky )
library( treedater ) 
library( sarscov2 ) 
library( ggplot2 )
library( grid )
library( gridExtra )
library( ggtree )
library( alakazam )
library( stringi )


#' Runs mlesky on dated trees from treedater
#' @param tds_list A list of a list of dated trees. Output from d1_date_trees.R
#' @param ofn Name of mlesky output RDS file to save
#' @param taxis Time period over which to plot mlesky
#' @return Output from mlesky: Effective population size over time (Ne(t))
run_mlesky <- function(tds_list, ofn, taxis = taxis) {
  
  if(!inherits(tds_list, what = c("list")))
    tds_list = readRDS(tds_list)
  
  
  res_mlesky_list = lapply(tds_list, function(tds) {
    
    weeks = round(as.numeric((date_decimal(max(tds[[1]]$sts))-date_decimal(min(tds[[1]]$sts)))/7))
    res <- weeks * 2
    class( tds ) <- 'multiPhylo' 
    
    tds = lapply( tds , 
                  function(x) {x$tip.label = unlist(lapply(strsplit(x$tip.label, '[|]'), function(y) paste0(y[1])))
                  return(x)}
    )
    
    NeStartTimeBeforePresent =  max( tds[[1]]$sts ) - decimal_date(as.Date('2020-10-25'))
    
    # run mlesky. WARNING: number of cores demanded is ncpu * mc.cores; in this case it is 3*10 = 30.
    sgs = parallel::mclapply( tds, function(td) {
      mlskygrid(td, tau = NULL, tau_lower=.001, tau_upper = 10 , sampleTimes = td$sts , res = res, ncpu = 3, NeStartTimeBeforePresent = NeStartTimeBeforePresent)
    }, mc.cores = 10 )
    
    
    out = lapply(sgs, function(sg) {
      with( sg, approx( time, ne, rule = 1, xout = taxis )$y )
    })
    
    out
    
  })
  
  res_mlesky <-  do.call( cbind, lapply(res_mlesky_list, function(x) do.call( cbind, x ) ))
  
  saveRDS( list( time = taxis, ne = res_mlesky ) , file=paste0(ofn, "_mlesky", '.rds' ))
  
  res_mlesky
}



# time period over which to estimate Ne
taxis = decimal_date( seq( as.Date( '2020-10-15') , as.Date('2021-01-24'), by = 1) )


res_mlesky_B.1.1.7 = run_mlesky(tds_list = "Sample_England_sampler1_B.1.1.7_2021-02-13_n=3000_n_tree_dating_10_dated_trees.rds",
                        ofn = "Sample_England_sampler1_B.1.1.7_2021-02-13_n=3000_n_tree_dating_10", taxis = taxis)

res_mlesky_B.1.177 = run_mlesky(tds_list = "Sample_England_sampler1_B.1.177_2021-02-13_n=3000_n_tree_dating_10_dated_trees.rds",
                                ofn = "Sample_England_sampler1_B.1.177_2021-02-13_n=3000_n_tree_dating_10", taxis = taxis)

res_mlesky_control = run_mlesky(tds_list = "Sample_England_matchSample_control_2021-02-13_n_tree_dating_10_dated_trees.rds",
                                ofn = "Sample_England_matchSample_control_2021-02-13_n_tree_dating_10", taxis = taxis)
