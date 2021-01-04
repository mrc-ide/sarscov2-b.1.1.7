library(here)
library(kernlab)
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

rownames(S_neg)=rownames(S_pos)=un
colnames(S_neg)=colnames(S_pos)=weeks

if(CORRECT){
	d$S_negative = (d$total_cases * d$sgss_s_negative) / (d$sgss_s_negative+ d$sgss_s_positive)
	d$S_positive = (d$total_cases * d$sgss_s_positive) / (d$sgss_s_negative+ d$sgss_s_positive)
}else {
	d$S_negative =  d$sgss_s_negative 
	d$S_positive =  d$sgss_s_positive
}

d$S_negative[d$S_negative==0]=1
d$S_positive[d$S_positive==0]=1

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

growth_rates_p = as.vector(S_pos_rate)
growth_rates_n = as.vector(S_neg_rate)

ks.test(matrix(growth_rates_p),matrix(growth_rates_n))

kmmd(matrix(growth_rates_p),matrix(growth_rates_n))

mat =  matrix(c(sum(growth_rates_p>=1), sum(growth_rates_n>=1), sum(growth_rates_p<1), sum(growth_rates_n<1)), nrow = 2,
	              dimnames =
	       list(c("S_positive", "S_negative"),
		    c("r>=1", "r<1")))
t(mat)		    
fisher.test(t(mat), alternative = "less")

mat2 =  matrix(c(sum(growth_rates_p>1 & growth_rates_n>1), sum(growth_rates_p>1 & growth_rates_n<=1)	, sum(growth_rates_p<=1 & growth_rates_n>1)	, sum(growth_rates_p<=1 & growth_rates_n<=1)), nrow = 2,
	              dimnames =
	       list(c("Sneg>1", "Sneg<=1"),
		    c("Spos>1", "Spos<=1")))
		    
print(t(mat2))
print(mcnemar.test(mat2))
