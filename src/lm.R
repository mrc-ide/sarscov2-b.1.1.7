library(tidyverse)
library(here)
library(gridExtra)

##full data range for bootstrapping
weekstart <-45
weekend <- 50

datadir <- here("data")

formulas <- list(lmadd=Rt~factor(epiweek)+area+novel_frac,
                 lmeradd=Rt~factor(epiweek)+(1|area)+novel_frac
)

models <- list(
  lmadd=function(dat){
    dat <- dat %>% filter(epiweek>=weekstart,epiweek<=weekend)
    fit <- lm(formulas[["lmadd"]],data=dat)
    fit$coeff[["novel_frac"]]
  },
  lmeradd=function(dat){
    dat <- dat %>% filter(epiweek>=weekstart,epiweek<=weekend)
    fit <- lme4::lmer(formulas[["lmeradd"]],data=dat)
    lme4::fixef(fit)["novel_frac"]
  }
)

standardCI <- list(lmadd=function(dat){
  dat <- dat %>% filter(epiweek>=weekstart,epiweek<=weekend)
  fit <- lm(formulas[["lmadd"]],data=dat)
  confint(fit, "novel_frac")
},
lmeradd=function(dat){
  dat <- dat %>% filter(epiweek>=weekstart,epiweek<=weekend)
  fit <- lme4::lmer(formulas[["lmeradd"]],data=dat)
  confint(fit, "novel_frac")
})

la_rt<- readRDS(file.path(datadir,"la_rt_new_43_56_weeks.rds"))

LTLA_sgss <- readRDS(file.path(datadir,"sgss_la_new_43_56_weeks.rds"))%>%
  mutate(novel_frac=sgss_s_negative_corrected/(sgss_s_negative_corrected+sgss_s_positive_corrected))%>%
  inner_join(la_rt,by=c("area","epiweek"))%>%
  rename(novel=sgss_s_negative,other=sgss_s_positive)

stp_rt <- readRDS(file.path(datadir,"stp_rt_all.rds"))%>%
  rename(Rt="R(all)",CIlow=RallLI,CIup=RallUI)%>%
  select(area,epiweek,Rt,CIlow,CIup)

stp_sgss <- readRDS(file.path(datadir,"sgss_stp_new_43_56_weeks.rds"))%>%
  mutate(novel_frac=sgss_s_negative_corrected/(sgss_s_negative_corrected+sgss_s_positive_corrected))%>%
  inner_join(stp_rt,by=c("area","epiweek"))%>%
  rename(novel=sgss_s_negative,other=sgss_s_positive)

alldata <- list("LA SGSS"=list(dat=LTLA_sgss),
                "STP SGSS"=list(dat=stp_sgss))


##### Code for running bootstraps

######### independent sampling 
bdat_indep<- function(dat){
  dat[sample(1:dim(dat)[1],replace=TRUE),]
}

## area based sampling
bdat_area<- function(dat){
  a <- unique(dat$area)
  l <- length(a)
  s <- sample(a,size=l,replace=TRUE)
  bind_rows(lapply(1:l,function(x) dat%>%filter(area==s[x])%>%mutate(area=x)))%>%
    mutate(area=paste(area))    
}

bdat_areanovelfrac <- function(dat){
  res <- dat%>%
    bdat_area()
  total <- res$novel+res$other
  positive <- rep(0,length(total))
  positive[total>0] <- rbinom(sum(total>0),
                              size=total[total>0],
                              prob=(res$novel/total)[total>0])
  res$novel_frac <- positive*res$probs/total
  res%>%mutate(area=paste(area))
}




runbootstrap <- function(dat,weekstart,weekend,nboot=1e3,bdat=bdat_area){
  dat <- dat%>%filter(epiweek>=weekstart )%>%filter(epiweek<=weekend)
  samples <- mclapply(1:nboot,function(i) bdat(dat),mc.cores=50)
  res <- tibble()
  for (m in names(models)){
    
    T <- models[[m]]
    
    Tboot <- mclapply(samples, T,mc.cores=50)
    Tboot <- unlist(Tboot)
    Tobs <- T(dat)
    
    q <- quantile(Tboot-Tobs,c(0.975, .025))
    CI <- Tobs - q
    add <- tibble(startweek=weekstart,endweek=weekend,model=m,effect=Tobs,CIlow=CI[1],CIup=CI[2])
    res <- bind_rows(res,add)
  }
  res
}


runstandardCI <- function(dat,weekstart,weekend){
  dat <- dat%>%filter(epiweek>=weekstart )%>%filter(epiweek<=weekend)
  res <- tibble()
  for (m in names(models)){
    Tobs <- models[[m]](dat)
    CI <- standardCI[[m]](dat)
    add <- tibble(startweek=weekstart,endweek=weekend,model=m,effect=Tobs,CIlow=CI[1,1],CIup=CI[1,2])
    res <- bind_rows(res,add)
  }
  res
}

###############
res <- tibble()
for (d in names(alldata))
  res <- res%>%
  bind_rows(
    runstandardCI(alldata[[d]]$dat,weekstart=weekstart,weekend=weekend)%>%
      mutate(data=d,bootstrap="standard CI")
  )

nboot <- 1000
library(parallel)
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
  )


## for publ

TableR <- res %>% filter(bootstrap=="area,Rt,novel_frac")%>%
  filter(is.element(model,c("lmadd","lmeradd")))
TableR <- TableR%>%
  mutate(res=paste(round(effect,2)," [",round(CIlow,4),", ",round(CIup,4),"]",sep=""))%>%
  relocate(res)%>%select(res,model,data)

TableR
pdf(here("figures/ED_table2.pdf"), height=11, width=8.5)
grid.table(TableR)
dev.off()
