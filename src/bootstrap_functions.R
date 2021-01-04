## independent sampling 
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



bdat_novelfracgenome <- function(dat){
    dat %>% 
        mutate(total=novel+other)%>%
        mutate(novel= rbinom(length(total),
                             size=total,
                             prob=(novel_frac)))%>%
        mutate(novel_frac = (novel)/(total))
}


bdat_areanovelfracgenome <- function(dat){
    dat %>%
    bdat_area    %>%
        bdat_novelfracgenome()
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
        Tbiascor <- Tobs-mean(Tboot-Tobs)
        Mediancor <- Tobs-median(Tboot-Tobs)
        CI <- Tobs - q
        if (length(grep("mult",m))==1){
            Tobs <- exp(Tobs)
            Tbiascor <- exp(Mediancor)
            CI <- exp(CI)
        }
        add <- tibble(startweek=weekstart,endweek=weekend,model=m,effect=Tobs,CIlow=CI[1],CIup=CI[2],Tbiascor=Tbiascor,Mediancor=Mediancor)
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
