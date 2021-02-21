library(optparse)
library(here)
library(epidemia)
library(dplyr)

option_list <- list(
    make_option(c("--input"),action="store", default=here("data/targets.csv"),help="Input file [default \"%default\"]"),
    make_option(c("--dInput"),action="store", default=here("data/sgtf_transmission_data.rds"),help="Data Input file [default \"%default\"]"),
    make_option(c("--fitscript"),action="store", default=here("transmission_model/models/joint_model.R"),help="Name of file used for fitting regions [default \"%default\"]"),
    make_option(c("--variant"),action="store", default="joint",help="Which varinat to use [default \"%default\"]"),
    make_option(c("--output"),action="store", default="fits/",help="Output directory [default \"%default\"]")
)
opt <- parse_args(OptionParser(option_list=option_list))


targets <- read.csv(opt$input,stringsAsFactors = FALSE)
targets$filename <- paste("fm-",gsub(" ","_",targets$AREA),"-",targets$subgroup,".rds",sep="")
dir.create(opt$output)
dir.create(paste(opt$output,"log/",sep=""))


sink(paste(opt$output,"Makefile",sep=""))

cat("SHELL=/bin/bash\n")

cat("all:  ",      sample(targets$filename))
cat("\n\n")

#cat("\n\ninclude regions.makefile\n\n")

cat("skip= test -e \"$@\" -a ! -s\"$@\"\n") ## to skip starting a simulation if it already exists as a file
cat("lock= > \"$@\" && for x in $^; do while test ! -s $$x; do sleep 5;done; done\n\n") ## create target as empty file and then wait for all dependencies to be non-empty
##  usage $(skip)||(($lock)&& what needs doing)

cat("clean:\n")
cat("\trm *.rds\n\n")

cat("clean0:\n")
cat("\tfind . -size  0 -print -delete\n\n")



for (j in 1:(dim(targets)[1])){
    i <- targets[j,]
    cat(i$filename,": ",opt$fitscript,
        "\n",sep="")
    cat("\techo \"starting ",i$filename,"\"\n",sep="")
    cat("\tfor i in 1 2 3 4 5; do \\\n")
    cat("\t  echo Sampling Attempt $$i ; \\\n")
    cat("\t  (timeout 150m nice Rscript $< --areaname=\"",i$AREA,"\"",
        " --subgroup=",i$subgroup,
        " --input=",opt$dInput,
        " --iter=3000",
        " --variant=",opt$variant,
        " --nchains=3",
        " --thin=4",
        " --output=.",
        " 1>>\"log/",i$AREA,".log\"",
        " 2>&1",
        ")&&break;\\\n",
        sep="")
    cat("\tdone \n")
    cat("\techo \"completed sampling for ",i$AREA,"\"\n",sep="")
    cat("\n")
}
sink()

