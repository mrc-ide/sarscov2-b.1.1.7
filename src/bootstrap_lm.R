#################
## main code for the frequentis regression results
####################

library(tidyverse)
library(parallel)
library(here)
##full data range for bootstrapping
weekstart <-43
weekend <- 50 

##subsets for specific models
week_same <- list(start=44,end=50) ##ranges for standard model
week_next <- list(start=43,end=49) ##ranges for lagged model  

source(here("src/bootstrap_functions.R"))

datadir <- here("data")

models <- list(lmadd=function(dat){
    dat <- dat %>% filter(epiweek>=week_same$start,epiweek<=week_same$end)
    fit <- lm(Rt~factor(epiweek)+area+novel_frac,data=dat)
    fit$coeff[["novel_frac"]]
},
lmeradd=function(dat){
    dat <- dat %>% filter(epiweek>=week_same$start,epiweek<=week_same$end)
    fit <- lme4::lmer(Rt~factor(epiweek)+(1|area)+novel_frac,data=dat)
    lme4::fixef(fit)["novel_frac"]
},
lmaddnext=function(dat){
    dat <- dat %>% filter(epiweek>=week_next$start,epiweek<=week_next$end)
    fit <- lm(Rtnext~factor(epiweek)+area+novel_frac,data=dat)
    fit$coeff[["novel_frac"]]
},
lmeraddnext=function(dat){
    dat <- dat %>% filter(epiweek>=week_next$start,epiweek<=week_next$end)
    fit <- lme4::lmer(Rtnext~factor(epiweek)+(1|area)+novel_frac,data=dat)
    lme4::fixef(fit)["novel_frac"]
}
)

standardCI <- list(lmadd=function(dat){
    dat <- dat %>% filter(epiweek>=week_same$start,epiweek<=week_same$end)
    fit <- lm(Rt~factor(epiweek)+area+novel_frac,data=dat)
    confint(fit, "novel_frac")
},
lmeradd=function(dat){
    dat <- dat %>% filter(epiweek>=week_same$start,epiweek<=week_same$end)
    fit <- lme4::lmer(Rt~factor(epiweek)+(1|area)+novel_frac,data=dat)
    confint(fit, "novel_frac")
},lmaddnext=function(dat){
    dat <- dat %>% filter(epiweek>=week_next$start,epiweek<=week_next$end)
    fit <- lm(Rtnext~factor(epiweek)+area+novel_frac,data=dat)
    confint(fit, "novel_frac")
},
lmeraddnext=function(dat){
    dat <- dat %>% filter(epiweek>=week_next$start,epiweek<=week_next$end)
    fit <- lme4::lmer(Rtnext~factor(epiweek)+(1|area)+novel_frac,data=dat)
    confint(fit, "novel_frac")
})


### prepLTLA 
la_rt<- readRDS(file.path(datadir,"la_rt.rds"))%>%
    group_by(area)%>%mutate(Rtnext=lead(Rt,order_by=epiweek),
                            CIupnext=lead(CIup,order_by=epiweek),
                            CIlownext=lead(CIlow,order_by=epiweek)
                            )%>%ungroup

LTLA_sgss <- readRDS(file.path(datadir,"sgss_la.rds"))%>%
    mutate(novel_frac=sgss_s_negative_corrected/(sgss_s_negative_corrected+sgss_s_positive_corrected))%>%
    inner_join(la_rt,by=c("area","epiweek"))%>%
    rename(novel=sgss_s_negative,other=sgss_s_positive)

### prepare Rt
stp_rt <- readRDS(file.path(datadir,"stp_rt.rds"))%>%rename(Rt=value)%>%
    group_by(area)%>%mutate(Rtnext=lead(Rt,order_by=epiweek),
                            CIupnext=lead(CIup,order_by=epiweek),
                            CIlownext=lead(CIlow,order_by=epiweek)
                            )%>%ungroup
stp_genome_raw <- readRDS(file.path(datadir,"stp_data_genome.rds"))%>%
    rename(area=stp_name)%>%
    rename(novel="B.1.1.7")%>%
    select(area,epiweek,novel,other)%>%
    mutate(novel_frac=novel/(novel+other))%>%   
    inner_join(stp_rt,by=c("area","epiweek"))

stp_genome <- stp_genome_raw%>%filter(novel+other>=5)  ## exclude small counts


stp_sgss <- readRDS(file.path(datadir,"sgss_stp.rds"))%>%
    mutate(novel_frac=sgss_s_negative_corrected/(sgss_s_negative_corrected+sgss_s_positive_corrected))%>%
    inner_join(stp_rt,by=c("area","epiweek"))%>%
    rename(novel=sgss_s_negative,other=sgss_s_positive)
    

alldata <- list("LA SGSS"=list(dat=LTLA_sgss),
                "STP gen"=list(dat=stp_genome),
                "STP SGSS"=list(dat=stp_sgss))
    
res <- tibble()
for (d in names(alldata))
    res <- res%>%
        bind_rows(
            runstandardCI(alldata[[d]]$dat,weekstart=weekstart,weekend=weekend)%>%
            mutate(data=d,bootstrap="standard CI")
        )

nboot <- 1000
RNGkind("L'Ecuyer-CMRG")
set.seed(120310 )
mc.reset.stream()
res <- res%>%
    bind_rows(
        runbootstrap(LTLA_sgss,weekstart=weekstart,weekend=weekend,nboot=nboot,bdat=bdat_areanovelfrac)%>%
        mutate(data= "LTLA SGSS",bootstrap = "area,Rt,novel_frac")
       ,
        runbootstrap(stp_sgss,weekstart=weekstart,weekend=weekend,nboot=nboot,bdat=bdat_areanovelfrac)%>%
        mutate(data= "stp SGSS",bootstrap = "area,Rt,novel_frac")
       ,
        runbootstrap(stp_genome,weekstart=weekstart,weekend=weekend,nboot=nboot,bdat=bdat_areanovelfracgenome)%>%
        mutate(data= "STP genome",bootstrap = "area,Rt,novel_frac")
    )

res$startweek[is.element(res$model,c("lmadd","lmeradd"))] <- week_same$start
res$endweek[is.element(res$model,c("lmadd","lmeradd"))] <- week_same$end
res$startweek[is.element(res$model,c("lmaddnext","lmeraddnext"))] <- week_next$start
res$endweek[is.element(res$model,c("lmaddnext","lmeraddnext"))] <- week_next$end

## for publ
TableR <- res %>% filter(bootstrap=="area,Rt,novel_frac")%>%
    filter(is.element(model,c("lmadd","lmeradd")))
o <- c(5,6,1,2,3,4)
TableR <- TableR[o,]
TableR <- TableR%>%
    mutate(res=paste(round(effect,2)," [",round(CIlow,4),", ",round(CIup,4),"]",sep=""))%>%
    relocate(res)%>%select(res,model,data)



TableR_lag <- res %>% filter(bootstrap=="area,Rt,novel_frac")%>%
    filter(is.element(model,c("lmaddnext","lmeraddnext")))
o <- c(5,6,1,2,3,4)
TableR_lag<- TableR_lag[o,]
TableR_lag <- TableR_lag%>%
    mutate(res=paste(round(effect,2)," [",round(CIlow,4),", ",round(CIup,4),"]",sep=""))%>%
    relocate(res)%>%select(res,model,data)


publ <- list(TableR=TableR,
             TableR_lag=TableR_lag,
             dimgenome=list(all=dim(stp_genome_raw%>%filter(epiweek>=week_next$star,epiweek<=week_next$end))[1],
                            exclsmall=dim(stp_genome%>%filter(epiweek>=week_next$star,epiweek<=week_next$end))[1]),
             dimstp_sgss=dim(stp_sgss%>%filter(epiweek>=week_next$star,epiweek<=week_next$end))[1]
             )

countdat <- stp_genome%>%filter(epiweek>=week_next$start,epiweek<=week_next$end)
publ$countsGenome <-
    list(novel=sum(countdat$novel),  
         total=sum(countdat$novel+countdat$other), 
         meannovel=mean(countdat$novel), 
         meantoal=mean(countdat$novel+countdat$other)) 


countdat <- stp_sgss%>%filter(epiweek>=week_next$start,epiweek<=week_next$end)

publ$countsSGSS <- list(totalcorrected=sum(countdat$sgss_s_positive_corrected+countdat$sgss_s_negative_corrected),
                        total=sum(countdat$novel+countdat$other),
                        novelcorrected=sum(countdat$sgss_s_negative_corrected))


### classic CI - not in paper - for comparison purposes
TableR_asympt <- res %>% filter(bootstrap=="standard CI")%>%
    filter(is.element(model,c("lmadd","lmeradd")))%>%
    filter(is.element(data,c("LA SGSS", "STP SGSS","STP gen")))
o <- c(3,4,1,2,5,6)
TableR_asympt<- TableR_asympt[o,]
TableR_asympt <- TableR_asympt%>%
    mutate(res=paste(round(effect,2)," [",round(CIlow,4),", ",round(CIup,4),"]",sep=""))%>%
    relocate(res)%>%select(res,model,data)

TableRlag_asympt <- res %>% filter(bootstrap=="standard CI")%>%
    filter(is.element(model,c("lmaddnext","lmeraddnext")))%>%
    filter(is.element(data,c("LA SGSS", "STP SGSS","STP gen")))
o <- c(3,4,1,2,5,6)
TableRlag_asympt<- TableRlag_asympt[o,]
TableRlag_asympt <- TableRlag_asympt%>%
    mutate(res=paste(round(effect,2)," [",round(CIlow,4),", ",round(CIup,4),"]",sep=""))%>%
    relocate(res)%>%select(res,model,data)

TableR_asympt
TableRlag_asympt

print(publ)

saveRDS(file=here("results/respubl.rds"),publ)
