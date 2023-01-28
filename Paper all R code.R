### R function for random effect model meta-analysis

random_effect_meta<-function(B,V){
  W<-1/V                                #Inverse variance weight
  Q<-sum(W*B^2)-(sum(W*B))^2/sum(W)     #total variance or Q statistics
  df<-length(B)-1                       #degree of freedom, number of studies minus 1
  c<-sum(W)-sum(W^2)/sum(W)             #c is the scaling factor

  tau_square<-ifelse(Q>df,(Q-df)/c,0)   #between study variance

  V_star=V+tau_square
  W_star=1/V_star                        #new weight
  Beta_star<-sum(W_star*B)/sum(W_star)   #overall random effect coefficients
  Var_star=1/sum(W_star)                 #ovall variance
  upper_b<-Beta_star+1.96*sqrt(Var_star) #upper 95% CI for overall beta
  lower_b<-Beta_star-1.96*sqrt(Var_star) #lower 95% CI for overall beta

  resultlist<-list('Overal beta'=Beta_star,'Overall variance'=Var_star,'HR'=exp(Beta_star),'HR lower 95% CI'=exp(lower_b),'HR upper 95% CI'=exp(upper_b))#change to hazard ratio
  return(resultlist)
}

# Random meta analyis for beta1
studies<-read.csv('https://raw.githubusercontent.com/enwuliu/meta-analysis/main/random_effect_meta_sim.csv',header=T)
b1<-studies$b1 #beta1 of the 10 cohorts
var_b1<-studies$var_b1 #variances of the coefficients
meta1<-random_effect_meta(b1,var_b1) # random effect meta-analysis for coefficient beta1
meta1

#### using the meta package to conduct random meta for beta1
library(meta)
m.gen<- metagen(TE = b1,
                seTE = sqrt(var_b1),
                studlab = LETTERS[1:10],
                data = ,
                sm = "HR",
                method.tau = "DL",
                fixed = FALSE,
                random = TRUE,
                title = "Use R meta package")

summary(m.gen)

#forest plot for beta1
forest.meta(m.gen, layout = "RevMan5")

# random effect meta-analysis for beta2 using function written by the authors

b2<-studies$b2 #beta2 of the 10 cohorts
var_b2<-studies$var_b2 #variances of the coefficients
meta2<-random_effect_meta(b2,var_b2) # random effect meta-analysis for coefficient beta1
meta2

#### using the meta package to conduct random meta for beta2

m.gen2<- metagen(TE = b2,
                 seTE = sqrt(var_b2),
                 studlab = LETTERS[1:10],
                 data = ,
                 sm = "HR",
                 method.tau = "DL",
                 fixed = FALSE,
                 random = TRUE,
                 title = "Use R meta package")

summary(m.gen2)

#forest plot for beta2

forest.meta(m.gen2, layout = "RevMan5")

#function for condducting random effect meta-analysis for correlation coefficient r

random_effect_meta_r<-function(V11,V22,Cov12,S_size){
  r<-Cov12/sqrt(V11*V22)           #correlation coefficient
  z<-0.5*log((1+r)/(1-r))          #Fisher z transformation
  V<-1/(S_size-3)                  #variance for z
  W<-1/V                            #weight
  Q<-sum(W*z^2)-(sum(W*z))^2/sum(W) #total variance or Q statistics
  df<-length(z)-1                   #degree of freedom, number of studies minus 1
  c<-sum(W)-sum(W^2)/sum(W)         #c is the scaling factor

  tau_square<-ifelse(Q>df,(Q-df)/c,0) #between study variance

  V_star=V+tau_square
  W_star=1/V_star                     #new weight
  z_star<-sum(W_star*z)/sum(W_star)   #overall random effect of z
  r_overall<-(exp(2*z_star)-1)/(exp(2*z_star)+1) #transform overall z* back to overall correlation coefficient
  return(r_overall)
}

V11<-studies$var_b1
V22<-studies$var_b2
V12<-studies$cov_b1b2
S_size<-studies$sample_size

meta3<-random_effect_meta_r(V11,V22,V12,S_size)
meta3  #overall correlation coefficient

## using meta package to perfomr random effect meta-analysis on correlation coefficient r

r<-V12/sqrt(V11*V22) #
m.cor <- metacor(cor = r, 
                 n = S_size,
                 studlab = LETTERS[1:10],
                 data = ,
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "DL",
                 hakn = FALSE,
                 title = "Use R meta package")
summary(m.cor)

#calculate covariance of beta1 and beta2
cov_b1b2<-meta3*sqrt(meta1$`Overall variance`*meta2$`Overall variance`)
cov_b1b2

#plot figure1 
age<-seq(50,90,1)
HR<-exp(meta1$`Overal beta`+meta2$`Overal beta`*age)
HR_low<-exp((meta1$`Overal beta`+meta2$`Overal beta`*age)-1.96*sqrt(meta1$`Overall variance`+age^2*meta2$`Overall variance`+2*age*cov_b1b2))
HR_upper<-exp((meta1$`Overal beta`+meta2$`Overal beta`*age)+1.96*sqrt(meta1$`Overall variance`+age^2*meta2$`Overall variance`+2*age*cov_b1b2))
plot(age,HR,type='l',ylim=c(0,3.5),ylab='HR(fall vs. no fall)')
lines(age,HR_low,type='l',col='blue',lty=2)
lines(age,HR_upper,type='l',col='blue',lty=2)

#overall beta1, beta2 and covaraince by multivariate meta-analysis
library(mvmeta)
library(dplyr)
S<-as.matrix(select(sim_data,var_b1,cov_b1b2, var_b2))
model<-mvmeta(cbind(sim_data$b1,sim_data$b2),S,method="ml")
summary(model)
model$vcov

# figure 2

cov<-model$vcov[1,2]
beta<-model$coefficients
HR_mv<-exp(beta[1]+beta[2]*age)
HR_mv_low<-exp((beta[1]+beta[2]*age)-1.96*sqrt(model$vcov[1,1]+age^2*model$vcov[2,2]+2*age*cov))
HR_mv_upper<-exp((beta[1]+beta[2]*age)+1.96*sqrt(model$vcov[1,1]+age^2*model$vcov[2,2]+2*age*cov))
plot(age,HR_mv,type='l',ylim=c(0,8),ylab='HR(fall vs. no fall)',col='red',lwd=1.0,bty="n")
lines(age,HR_mv_low,type='l',col='red',lty=4)
lines(age,HR_mv_upper,type='l',col='red',lty=4)
lines(age,HR,type='l',ylim=c(0,3.5),ylab='HR(fall vs. no fall)',col='blue',lwd=2.0)
lines(age,HR_low,type='l',col='blue',lty=2)
lines(age,HR_upper,type='l',col='blue',lty=2)
legend(70, 8, c("Multivariate meta-analysis", "Univariate meta-analysis using r", "95% CI multivariate meta-analysis", "95% CI using r"), col=c("red", "blue","red","blue"), lty=c(1,1,4,2),cex=0.8, box.lty=0)


