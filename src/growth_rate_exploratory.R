# Growth Rate Ratio 
library(here)
library(ggplot2)
library(ggpubr)
library(viridis)
library(matrixStats)
library(mgcv)
d = readRDS(here("data/sgss_stp.rds"))
d=as.data.frame(d)
# parameters
RESOLUTION="area"
choice_weeks = 46:50
CORRECT=TRUE

un=sort(unique(d[,RESOLUTION])) # unique LTLA names ~300
weeks=sort(unique(d$epiweek))

S_neg = matrix(NA, nrow=length(un),ncol=length(weeks))
S_pos = matrix(NA, nrow=length(un),ncol=length(weeks))

#d$sgss_s_negative[d$sgss_s_negative<5]=NA
#d$sgss_s_negative[d$sgss_s_positive<5]=NA



rownames(S_neg)=rownames(S_pos)=un
colnames(S_neg)=colnames(S_pos)=weeks

if(CORRECT){
	d$S_negative = (d$total_cases * d$sgss_s_negative) / (d$sgss_s_negative+ d$sgss_s_positive)
	d$S_positive = (d$total_cases * d$sgss_s_positive) / (d$sgss_s_negative+ d$sgss_s_positive)
}else {
	d$S_negative =  d$sgss_s_negative 
	d$S_positive =  d$sgss_s_positive
}

d$S_negative[d$S_negative==0]=NA
d$S_positive[d$S_positive==0]=NA

# pivot table explicit
for(i in 1:length(un)){
	for(j in 1:length(weeks)){
		S_neg[i,j] = sum(d$S_negative[d$epiweek%in%weeks[j] & d[,RESOLUTION]%in%un[i]])
		S_pos[i,j] = sum(d$S_positive[d$epiweek%in%weeks[j] & d[,RESOLUTION]%in%un[i]])
	}
}
# subset to the number of weeks of interest
S_neg=S_neg[,weeks%in%choice_weeks]
S_pos=S_pos[,weeks%in%choice_weeks]

S_neg_rate = matrix(NA, nrow=nrow(S_neg),ncol=ncol(S_neg)-1)
S_pos_rate = matrix(NA, nrow=nrow(S_neg),ncol=ncol(S_neg)-1)

for(i in 1:nrow(S_neg)){
	for(j in 1:(ncol(S_neg)-1)){
		S_neg_rate[i,j] = (S_neg[i,j+1]/S_neg[i,j])^(6.5/7)
		S_pos_rate[i,j] = (S_pos[i,j+1]/S_pos[i,j])^(6.5/7)
	}
}


S_rate_ratio = S_neg_rate/S_pos_rate
# frequencies corresponding to numerator week
S_proportion_ratio = S_neg[,1:(ncol(S_neg)-1)]/(S_neg[,1:(ncol(S_neg)-1)]+S_pos[,1:(ncol(S_pos)-1)])

df=list(x=as.vector(t(S_proportion_ratio)),
y=as.vector(t(S_rate_ratio)),
week=rep(choice_weeks[1:(length(choice_weeks)-1)],length(un)) # numerator week
)

knots=df$x[!is.na(df$y)]
l <- gam( (y)~ s(x), data=data.frame(x=df$x,y=df$y),method = "REML")
p <- predict(l,se.fit = TRUE)
upr <- (p$fit + (2 * p$se.fit))
lwr <- (p$fit - (2 * p$se.fit))
mn = (p$fit)


g=ggplot()
g=g+geom_point(data=data.frame(x=df$x,y=df$y,Week=as.factor(df$week)),aes(x=x,y=y,shape=Week,colour=Week),alpha=0.9,size=2.4)
g=g+geom_line(data=data.frame(x=knots,y=mn),aes(x=x,y=y),col='skyblue',size=1.5)
g=g+geom_ribbon(data=data.frame(x=knots,ymin=lwr,ymax=upr,y=mn),aes(x=x,y=y,ymin=ymin,ymax=ymax),fill='grey',alpha=0.2)

g=g+geom_hline(yintercept=1, linetype="dashed", color = "purple", size=1)
g=g + theme_bw() +
	ylim(0,4)+
	xlab("Proportion S gene negative")+
	ylab("Multiplicative advatage of variant") +
  theme(legend.position="none")

ggsave(here('figures/funnelplot_ratio.png'),g,width=7,height=4)

# Growth rate difference
d = readRDS(here("data/sgss_stp.rds"))
d=as.data.frame(d)

un=sort(unique(d[,RESOLUTION])) # unique LTLA names ~300
weeks=sort(unique(d$epiweek))

S_neg = matrix(NA, nrow=length(un),ncol=length(weeks))
S_pos = matrix(NA, nrow=length(un),ncol=length(weeks))

#d$sgss_s_negative[d$sgss_s_negative<5]=NA
#d$sgss_s_positive[d$sgss_s_positive<5]=NA


rownames(S_neg)=rownames(S_pos)=un
colnames(S_neg)=colnames(S_pos)=weeks

if(CORRECT){
  d$S_negative = (d$total_cases * d$sgss_s_negative) / (d$sgss_s_negative+ d$sgss_s_positive)
  d$S_positive = (d$total_cases * d$sgss_s_positive) / (d$sgss_s_negative+ d$sgss_s_positive)
}else {
  d$S_negative =  d$sgss_s_negative 
  d$S_positive =  d$sgss_s_positive
}

d$S_negative[d$S_negative==0]=NA
d$S_positive[d$S_positive==0]=NA

# pivot table explicit
for(i in 1:length(un)){
  for(j in 1:length(weeks)){
    S_neg[i,j] = sum(d$S_negative[d$epiweek%in%weeks[j] & d[,RESOLUTION]%in%un[i]])
    S_pos[i,j] = sum(d$S_positive[d$epiweek%in%weeks[j] & d[,RESOLUTION]%in%un[i]])
  }
}
# subset to the number of weeks of interest
S_neg=S_neg[,weeks%in%choice_weeks]
S_pos=S_pos[,weeks%in%choice_weeks]

S_neg_rate = matrix(NA, nrow=nrow(S_neg),ncol=ncol(S_neg)-1)
S_pos_rate = matrix(NA, nrow=nrow(S_neg),ncol=ncol(S_neg)-1)

for(i in 1:nrow(S_neg)){
  for(j in 1:(ncol(S_neg)-1)){
    S_neg_rate[i,j] = (S_neg[i,j+1]/S_neg[i,j])^(6.5/7)
    S_pos_rate[i,j] = (S_pos[i,j+1]/S_pos[i,j])^(6.5/7)
  }
}

S_rate_ratio = (S_neg_rate-S_pos_rate) 
# frequencies corresponding to numerator week
S_proportion_ratio = S_neg[,1:(ncol(S_neg)-1)]/(S_neg[,1:(ncol(S_neg)-1)]+S_pos[,1:(ncol(S_pos)-1)])

df=list(x=as.vector(t(S_proportion_ratio)),
        y=as.vector(t(S_rate_ratio)),
        week=rep(choice_weeks[1:(length(choice_weeks)-1)],length(un)) # numerator week
)


knots=df$x[!is.na(df$y)]
l <- gam( (y)~ s(x), data=data.frame(x=df$x,y=df$y),method = "REML")
p <- predict(l,se.fit = TRUE)
upr <- (p$fit + (2 * p$se.fit))
lwr <- (p$fit - (2 * p$se.fit))
mn = (p$fit)



g2=ggplot()
g2=g2+geom_point(data=data.frame(x=df$x,y=df$y,Week=as.factor(df$week)),aes(x=x,y=y,shape=Week,colour=Week),alpha=0.9,size=2.4)
g2=g2+geom_line(data=data.frame(x=knots,y=mn),aes(x=x,y=y),col='skyblue',size=1.5)
g2=g2+geom_ribbon(data=data.frame(x=knots,ymin=lwr,ymax=upr,y=mn),aes(x=x,y=y,ymin=ymin,ymax=ymax),fill='grey',alpha=0.2)

g2=g2+geom_hline(yintercept=0, linetype="dashed", color = "purple", size=1)
g2=g2 + theme_bw() +
  ylim(0,3)+
  xlab("Proportion S gene negative")+
  ylab("Additive advatage of variant")+
  theme(legend.position="none")
# scale_color_viridis(discrete = TRUE, option = "D")
g2
ggsave(here('figures/funnelplot_difference.png'),g2,width=7,height=4)
