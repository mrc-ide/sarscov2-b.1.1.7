# Date lists of ML trees using treedater.
# Required: 
# - nwk files for lineages: B.1.1.7, B.1.177 and control 
# - metadata for each sequence: metadata.csv
#   - column names: central_sample_id (name of sequence which matches the tip.label of nwk trees), sample_date (date of sampling in format YYYY-MM-DD)


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




#' Dates a single ML tree using treedater and a dataframe containing sampling times and clock rate
#' @param mltr A single ML tree
#' @param metadata A dataframe containing sequence name and sampling time
#' @param meanrate Clock rate
#' @return A dated tree
datetree <- function(mltr, metadata, meanrate)
{
  sts <- setNames( metadata$sample_time[  match( mltr$tip.label , metadata$central_sample_id ) ], mltr$tip.label )
  tr <- di2multi(mltr, tol = 1e-05)
  tr = unroot(multi2di(tr))
  tr$edge.length <- pmax(1/29000/5, tr$edge.length)
  dater(unroot(tr), sts[tr$tip.label], s = 29000, omega0 = meanrate, numStartConditions = 0, meanRateLimits = c(meanrate, meanrate + 1e-6), ncpu = 6)
}



#' Dates a list of ML trees using treedater and a dataframe containing sampling times and clock rate
#' @param mltr_fn Path to nwk file containing a list of ML trees
#' @param ofn Name of treedater output RDS file to save
#' @param n_tree_dating For each ML tree, this is the number of times to sample clock rate and re-date tree
#' @param metadata A dataframe containing sequence name and sampling time
#' @param meanrate Mean clock rate
#' @param meanratesd Standard deviation of mean clock rate
#' @param ncpu Number of cores demanded for parallelisation
#' @return A list of a list of dated trees. First level = number of ML trees (in nwk file). Second level = number of times each ML is dated with treedater (n_tree_dating).
date_trees <- function(mltr_fn, ofn, n_tree_dating = 10, metadata, meanrate, meanratesd, ncpu = 4, ...)
{
  mltr = read.tree(mltr_fn)
  
  # checking all samples have metadata attached... removing tips that aren't able to be matched
  mltr = lapply(mltr, function(tr) ape::drop.tip(tr, tr$tip.label[!tr$tip.label %in% metadata$central_sample_id]))
  
  
  tds = parallel::mclapply(mltr, function(tr) {
    tmp = lapply(1:n_tree_dating, function(x) {
      td = datetree( tr , metadata = metadata, meanrate =  max( 0.0001, rnorm( 1, meanrate, sd = meanratesd ) ) )
      td$tip.label =  paste0(td$tip.label, '|', as.Date(date_decimal(td$sts)), '|', td$sts ) # Rendering the tip.label in format for phylodynamic analyses
      td
    })
    tmp
  }, mc.cores = ncpu)
  
  saveRDS( tds , file=paste0(ofn, "_dated_trees", '.rds' ))
  
  tds
}



# Sample times 
metadata = read.csv( "metadata.csv" , stringsAs = FALSE , header=TRUE )
metadata$sample_date <- as.Date( metadata$sample_date )
metadata$sample_time <- decimal_date( metadata$sample_date ) # converting to decimal date for use in treedater




# Distribution of clock rate to sample
#95% HPD interval	[5.2065E-4, 6.7144E-4]
mr = 5.9158E-4
mrci = 	c( 5.2065E-4, 6.7144E-4)
mrsd = diff( mrci ) / 4 / 1.96


# make treedater trees for each lineage
tds_list_B.1.1.7 = date_trees(mltr_fn = "sampler1_B.1.1.7_2021-02-13_n=3000.nwk",
                      ofn = paste0('Sample_England_sampler1_B.1.1.7_2021-02-13_n=3000_n_tree_dating_10'), 
                      metadata = metadata, 
                      meanrate = mr,
                      n_tree_dating = 10,
                      meanratesd = mrsd, 
                      ncpu = 4)

# make treedater trees for each lineage
tds_list_B.1.177 = date_trees(mltr_fn = "matchSample_notB.1.1.7_lineageB.1.177_2021-02-13.nwk",
                              ofn = paste0('Sample_England_sampler1_B.1.177_2021-02-13_n=3000_n_tree_dating_10'), 
                              metadata = metadata, 
                              meanrate = mr,
                              n_tree_dating = 10,
                              meanratesd = mrsd, 
                              ncpu = 4)


# make treedater trees for each lineage
tds_list_control = date_trees(mltr_fn = "matchSample_notB.1.1.7_2021-02-13.nwk",
                              ofn = paste0('Sample_England_matchSample_control_2021-02-13_n_tree_dating_10'), 
                              metadata = metadata, 
                              meanrate = mr,
                              n_tree_dating = 10,
                              meanratesd = mrsd, 
                              ncpu = 4)






