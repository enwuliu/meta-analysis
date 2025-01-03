studies <- read.csv('https://raw.githubusercontent.com/enwuliu/meta-analysis/main/random_effect_meta_sim.csv', header = TRUE)
b1 <- studies$b1  # beta1 from the 10 cohorts
var_b1 <- studies$var_b1  # variances of the coefficients

fixed_effect_meta <- function(B, V) {
  W <- 1 / V
  Beta <- sum(W * B) / sum(W)
  Var <- 1 / sum(W)
  resultlist <- list('Overall beta' = Beta, 'Overall variance' = Var)
  return(resultlist)
}

random_effect_meta <- function(B, V) {
  W <- 1 / V
  Q <- sum(W * B^2) - (sum(W * B))^2 / sum(W)
  df <- length(B) - 1
  c <- sum(W) - sum(W^2) / sum(W)
  tau_square <- ifelse(Q > df, (Q - df) / c, 0)
  V_star <- V + tau_square
  W_star <- 1 / V_star
  Beta_star <- sum(W_star * B) / sum(W_star)
  Var_star <- 1 / sum(W_star)
  resultlist <- list('Overall beta' = Beta_star, 'Overall variance' = Var_star)
  return(resultlist)
}

fixed_b1 <- fixed_effect_meta(b1, var_b1)
fixed_b1

# $`Overall beta`
# [1] 1.040411
#
# $`Overall variance`
# [1] 0.6841809

random_b1 <- random_effect_meta(b1, var_b1)
random_b1

# $`Overall beta`
# [1] 1.014141
#
# $`Overall variance`
# [1] 0.7308388

library(meta)
b1.meta<- metagen(TE = b1,
                seTE = sqrt(var_b1),
                studlab = LETTERS[1:10],
                data = ,
                sm = "",
                method.tau = "DL",
                fixed = TRUE,
                random = TRUE,
                title = "Use R meta package")

summary(b1.meta)
forest(b1.meta)
#################beta 2

b2 <- studies$b2  # beta2 from the 10 cohorts
var_b2 <- studies$var_b2  # variances of the coefficients

fixed_b2<-fixed_effect_meta(b2, var_b2)
fixed_b2
# $`Overall beta`
# [1] -0.01140614
#
# $`Overall variance`
# [1] 0.0001403961

random_b2<-random_effect_meta(b2, var_b2)
random_b2
# $`Overall beta`
# [1] -0.01140614
#
# $`Overall variance`
# [1] 0.0001403961


b2.meta<- metagen(TE = b2,
                  seTE = sqrt(var_b2),
                  studlab = LETTERS[1:10],
                  data = ,
                  sm = "",
                  method.tau = "DL",
                  fixed = TRUE,
                  random = TRUE,
                  title = "Use R meta package")

summary(b2.meta)
forest(b2.meta)

###meta analysis on r

fixed_effect_meta_r<-function(v1,v2,cov12,sample_size){
  r<-cov12/sqrt(v1*v2)                   #correlation coefficient
  z<-0.5*log((1+r)/(1-r))                #Fisher z transformation
  v<-1/(sample_size-3)                   #variance for z
  W<-1/v                                 #weight
  z_overall<-sum(W*z)/sum(W)             #overall random effect of z
  r_overall<-(exp(2*z_overall)-1)/(exp(2*z_overall)+1) #transform back to r
  return(r_overall)
}
random_effect_meta_r<-function(v1,v2,cov12,sample_size){
  r<-cov12/sqrt(v1*v2)                   #correlation coefficient
  z<-0.5*log((1+r)/(1-r))                #Fisher z transformation
  v<-1/(sample_size-3)                   #variance for z
  W<-1/v                                 #weight
  Q<-sum(W*z^2)-(sum(W*z))^2/sum(W)      #total variance or Q statistics
  c<-sum(W)-sum(W^2)/sum(W)              #c is the scaling factor
  df=length(v1)-1                        #degree of freedom
  tau_square<-ifelse(Q>df,(Q-df)/c,0)    #between study variance

  v_star=v+tau_square
  W_star=1/v_star                        #new weight
  z_star<-sum(W_star*z)/sum(W_star)      #overall random effect of z
  r_overall<-(exp(2*z_star)-1)/(exp(2*z_star)+1) #transform back to r
  return(r_overall)
}

v1<-studies$var_b1
v2<-studies$var_b2
v12<-studies$cov_b1b2
sample_size<-studies$sample_size
r<-v12/sqrt(v1*v2)
fixed_meta_r<-fixed_effect_meta_r(v1,v2,v12,sample_size)
fixed_meta_r
#[1] -0.9600286
random_meta_r<-random_effect_meta_r(v1,v2,v12,sample_size)
random_meta_r
#[1] -0.9493409

r.meta <- metacor(cor = r, 
                 n = sample_size,
                 studlab = LETTERS[1:10],
                 data = ,
                 fixed = TRUE,
                 random = TRUE,
                 method.tau = "DL",
                 hakn = FALSE,
                 title = "Use R meta package")
summary(r.meta)

####plot

library(ggplot2)
library(ggpubr)

plot_with_interaction <- function(age, b1, v1, b2, v2, r) {
  cov_b1b2 <- r * sqrt(v1 * v2)
  RR <- exp(b1 + b2 * age)
  RR_lower <- exp((b1 + b2 * age) - 1.96 * sqrt(v1 + age^2 * v2 + 2 * age * cov_b1b2))
  RR_upper <- exp((b1 + b2 * age) + 1.96 * sqrt(v1 + age^2 * v2 + 2 * age * cov_b1b2))
  ndata <- as.data.frame(cbind(RR, RR_lower, RR_upper, age))
  ggplot(data = ndata, aes(x = age)) +
    geom_line(aes(y = RR)) +
    geom_line(aes(y = RR_lower), color = "steelblue", linetype = "dashed") +
    geom_line(aes(y = RR_upper), color = "steelblue", linetype = "dashed") +
    scale_y_continuous(name = "Rate Ratio (RR): Fallers vs. Non-Fallers",
                       limits = c(0, 4), expand = c(0, 0)) +
    theme_bw() +
    theme(axis.line = element_line(color = 'black'),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())
}

# Fixed-effect model
age <- seq(50, 90, 1)
b1 <- fixed_b1$`Overall beta`
v1 <- fixed_b1$`Overall variance`
b2 <- fixed_b2$`Overall beta`
v2 <- fixed_b2$`Overall variance`
r <- fixed_meta_r

p1 <- plot_with_interaction(age, b1, v1, b2, v2, r)
p1
# Random-effect model
b1_rand <- random_b1$`Overall beta`
v1_rand <- random_b1$`Overall variance`
b2_rand <- random_b2$`Overall beta`
v2_rand <- random_b2$`Overall variance`
r_rand <- random_meta_r

p2 <- plot_with_interaction(age, b1_rand, v1_rand, b2_rand, v2_rand, r_rand)
p2
p3 <- ggarrange(p1, p2,
                labels = c("Fixed Effect", "Random Effect"),
                ncol = 2, nrow = 1)
p3

#ggsave(filename = "E:\\cov_paper\\covmeta.jpeg",p3,units="mm", width = 250, height = 100,device='jpeg', dpi = 300)

#multivariate meta-analysis

library (mvmeta) 
library (dplyr)
S <-as.matrix(select(studies ,var_b1 , cov_b1b2 , var_b2)) 
model <- mvmeta (cbind ( studies$b1,studies $b2),S, method ="ml")
summary(model)

# Call:  mvmeta(formula = cbind(studies$b1, studies$b2) ~ 1, S = S, method = "ml")
# 
# Multivariate random-effects meta-analysis
# Dimension: 2
# Estimation method: ML
# 
# Fixed-effects coefficients
# Estimate  Std. Error        z  Pr(>|z|)  95%ci.lb  95%ci.ub   
# y1    1.1642      0.7770   1.4982    0.1341   -0.3588    2.6871   
# y2   -0.0193      0.0128  -1.5060    0.1321   -0.0445    0.0058   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# Between-study random-effects (co)variance components
# Structure: General positive-definite
# Std. Dev  Corr
# y1    0.8462    y1
# y2    0.0201     1
# 
# Multivariate Cochran Q-test for heterogeneity:
#   Q = 77.9106 (df = 18), p-value = 0.0000
# I-square statistic = 76.9%
# 
# 10 studies, 20 observations, 2 fixed and 3 random-effects parameters
# logLik      AIC      BIC  
# -5.0854  20.1708  25.1495  

model$vcov

# y1.(Intercept) y2.(Intercept)
# y1.(Intercept)    0.603797649  -0.0058182353
# y2.(Intercept)   -0.005818235   0.0001650336

b1_m<-model$coefficients[1,1]
b1_m
b2_m<-model$coefficients[1,2]
b2_m
v1_m<-model$vcov[1,1]
v1_m
v2_m<-model$vcov[2,2]
v2_m
cov_m<-model$vcov[1,2]
cov_m

plot_with_multivariate <- function(age, b1, v1, b2, v2, r, b1_m, v1_m, b2_m,V2_m,cov_m) {
  cov_b1b2 <- r * sqrt(v1 * v2)
  RR <- exp(b1 + b2 * age)
  RR_lower <- exp((b1 + b2 * age) - 1.96 * sqrt(v1 + age^2 * v2 + 2 * age * cov_b1b2))
  RR_upper <- exp((b1 + b2 * age) + 1.96 * sqrt(v1 + age^2 * v2 + 2 * age * cov_b1b2))
  
  RR_m <- exp(b1_m + b2_m * age)
  RR_lower_m <- exp((b1_m + b2_m * age) - 1.96 * sqrt(v1_m + age^2 * v2_m + 2 * age * cov_m))
  RR_upper_m <- exp((b1_m + b2_m * age) + 1.96 * sqrt(v1_m + age^2 * v2_m + 2 * age * cov_m))
  
  
  ndata <- as.data.frame(cbind(RR, RR_lower, RR_upper, RR_m, RR_lower_m, RR_upper_m, age))
  ggplot(data = ndata, aes(x = age)) +
    geom_line(aes(y = RR, color="Correlation coefficient method")) +
    geom_line(aes(y = RR_lower,color = "Correlation coefficient method"), linetype = "dashed") +
    geom_line(aes(y = RR_upper,color = "Correlation coefficient method"),, linetype = "dashed") +
    
    geom_line(aes(y = RR_m,color="Multivariate meta-analysis")) +
    geom_line(aes(y = RR_lower_m, color = "Multivariate meta-analysis"), linetype = "longdash") +
    geom_line(aes(y = RR_upper_m, color = "Multivariate meta-analysis"), linetype = "longdash") +
    
    scale_y_continuous(name = "Rate Ratio (RR): Fallers vs. Non-Fallers",
                       limits = c(0, 5), expand = c(0, 0)) +
    scale_color_manual(values = c("red", "blue"))+
    guides(color = guide_legend(title = "Methods"))+
    
    theme_bw() +
    theme(axis.line = element_line(color = 'black'),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())
  
}

p4<-plot_with_multivariate(age, b1, v1, b2, v2, r,b1_m, v1_m, b2_m, v2_m, cov_m)
p4

#ggsave(filename = "E:\\cov_paper\\multivariate.jpeg",p4,units="mm", width = 250, height = 100,device='jpeg', dpi = 300)



###### use the package
library(devtools)
install_github("enwuliu/covmeta")
library(covmeta)

# Load the dataset
studies <- read.csv('https://raw.githubusercontent.com/enwuliu/meta-analysis/main/random_effect_meta_sim.csv', header = TRUE)

# Function arguments
b1 <- studies$b1          # Beta1 coefficients from the 10 cohorts
v1 <- studies$var_b1      # Variances of Beta1
b2 <- studies$b2          # Beta2 coefficients from the 10 cohorts
v2 <- studies$var_b2      # Variances of Beta2
cov_b1b2 <- studies$cov_b1b2  # Covariance between Beta1 and Beta2
sample_size <- studies$sample_size  # Sample sizes

# Calculate the overall main effect, interaction effect, and covariance
# Fixed-effect meta-analysis
cov_meta(b1, v1, b2, v2, cov_b1b2, sample_size, 'fixed')

# Random-effect meta-analysis
cov_meta(b1, v1, b2, v2, cov_b1b2, sample_size, 'random')


