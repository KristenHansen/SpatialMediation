ammodel = "x ~ m + z"
ymmodel = "y ~ x + z + z2"
ymodel = "y~ m*x + z + z2 "
ysat = "y~x*m + z + z2"
mmodel = "m ~x + z + z2"
library(medflex)
#RMSE's are standradized by the sd of the bootstrap estimates so they are comparable across scenarios and methods 
Mediation <- function(n, a, b, c, xl, lm, ly, d){
  source("iorw.R")
  ab = a*b #indirect effect #.5625
  
  cp= c-ab #cp direct effect
  
  propmed = ab/ (cp + ab) #.6367
  
  z = rbinom(n, size = 1, prob = 0.3)
  
  x = rbinom(n, size = 1,prob = 0.5*z + 0.4)
  if(xl != 0) l = xl*x + (1 - ab^2)*rnorm(n) else l = rep(0, n)
  
  z2 = rnorm(n, 0, 1)
  m = a*x + lm*l + 0.3*z + 0.5*z2 + rnorm(n)
  ey = 1-(cp^2 +b^2 + 2*a*cp*b)
  
  y = cp*x + b*m + d*x*m + ly*l + 0.3*z + 0.5*z2 + ey*rnorm(n)
  data = as.data.frame(cbind(x,m,y, l, z, z2))
  return(data)
}

#Large Total Effect
#3
#Low Total Effect
#0.5
#low mediated = 5% or 1%
#high mediated 60%

#Add a confounder# to data generating function
high_tot = 3
low_tot = 0.5

low_med = 0.05
high_med = 0.60

sample_size = c(rep(100,12), rep(10000,12))
total_effects = c(rep(3,6), rep(0.5, 6), rep(3,6), rep(0.5,6))
proportion_med = c(rep(c(0.05, 0.05, 0.05, 0.6, 0.6, 0.6),4))
Z = c(rep(c(0, 0, 1), 8))
AMInt = c(rep(c(0, 1, 0),8))


dgm_table = as.data.frame(cbind(sample_size,total_effects,proportion_med,Z,AMInt))
dgm_table$n = dgm_table$sample_size
dgm_table$ab = total_effects*proportion_med
dgm_table$a = sqrt(dgm_table$ab)
dgm_table$b = sqrt(dgm_table$ab)
dgm_table$c = dgm_table$total_effects
dgm_table$xl = dgm_table$Z*0.25
dgm_table$lm = dgm_table$Z*0.3
dgm_table$ly = dgm_table$Z*0.3
dgm_table$d = dgm_table$AMInt*0.4
dgm_table$true_nde = dgm_table$c - dgm_table$ab
dgm_table$true_nie = dgm_table$ab

Datasets = list()
for(i in 1:24){
  print(i)
  Datasets[[i]] = Mediation(dgm_table$n[i], dgm_table$a[i], dgm_table$b[i], dgm_table$c[i], dgm_table$xl[i], dgm_table$lm[i], dgm_table$ly[i], dgm_table$d[i])
}


bootstrap_mediation <- function(obsdat, true_nde, true_nie){
  source("MedGFormula.R")
  gform <- GComputationMediation(obsdat, ymodel = ysat, ymmodel = ymmodel, mmodel = mmodel, m_fam = "gaussian")
  
  gform_nde <- gform[[1]]
  gform_nie <- gform[[2]]
  
  gform_boot <- function(x){
    boot_df <- obsdat[sample(nrow(obsdat), replace = T),]
    
    gform <- GComputationMediation(boot_df, ymodel = ysat, ymmodel = ymmodel, mmodel = mmodel, m_fam = "gaussian")
    
    gform_nde <- gform[[1]]
    gform_nie <- gform[[2]]
    
    gform_b <- c(NDE = gform_nde, NIE = gform_nie)
    return(gform_b)
  }
  
  gform_b_list <- lapply(1:250, function(x) gform_boot(x))
  gform_b_df <- data.frame(do.call(rbind, gform_b_list))
  
  gform_nde <- data.frame(nde = gform_nde, 
                         nde_lb = gform_nde - 1.96*sqrt(var(gform_b_df$NDE)), 
                         nde_ub = gform_nde + 1.96*sqrt(var(gform_b_df$NDE)))
  gform_nie <- data.frame(nie = gform_nie, 
                         nie_lb = gform_nie - 1.96*sqrt(var(gform_b_df$NIE)), 
                         nie_ub = gform_nie + 1.96*sqrt(var(gform_b_df$NIE)))
  gform_rmse_nde=sqrt(mean((gform_b_df$NDE - true_nde)^2))
  gform_rmse_nie=sqrt(mean((gform_b_df$NIE - true_nie)^2))
  
  gform <- cbind(gform_nde, gform_nie, gform_rmse_nde, gform_rmse_nie)
  print("gformcomplete")
  
  source("iorw.R")
  iorw_res = iorw(obsdat,expose = obsdat$x, "gaussian", ymmodel, ammodel)
  #Same problem as product method, makes sense this method does not account for this type of confounding
  iorw_nde=iorw_res[[1]]
  iorw_nie = iorw_res[[2]]
  
  iorw_boot<- function(x){
    boot_df <-  obsdat[sample(nrow(obsdat), replace=T),]
    iorw_res = iorw(boot_df,expose = boot_df$x, "gaussian", ymmodel, ammodel)
    #Same problem as product method, makes sense this method does not account for this type of confounding
    iorw_nde=iorw_res[[1]]
    iorw_nie = iorw_res[[2]]
    
    iorw_b <- c(NDE = iorw_nde, NIE = iorw_nie)
    return(iorw_b)
  }
  iorw_b_list <- lapply(1:250, function(x) iorw_boot(x))
  iorw_b_df <- data.frame(do.call(rbind, iorw_b_list))
  #Need to calculate squared error too
  
  iorw_nde <- data.frame(nde = iorw_nde, 
                         nde_lb = iorw_nde - 1.96*sqrt(var(iorw_b_df$NDE)), 
                         nde_ub = iorw_nde + 1.96*sqrt(var(iorw_b_df$NDE)))
  iorw_nie <- data.frame(nie = iorw_nie, 
                         nie_lb = iorw_nie - 1.96*sqrt(var(iorw_b_df$NIE)), 
                         nie_ub = iorw_nie + 1.96*sqrt(var(iorw_b_df$NIE)))
  iorw_rmse_nde=sqrt(mean((iorw_b_df$NDE - true_nde)^2))
  iorw_rmse_nie=sqrt(mean((iorw_b_df$NIE - true_nie)^2))
  
  iorw <- cbind(iorw_nde, iorw_nie, iorw_rmse_nde, iorw_rmse_nie)
  
  print("iorw complete")
  #RMPW

  expData<-neWeight(m ~ x + z + z2, family = "gaussian", data = obsdat)
  #w<-weights(expData)
  #head(w)
  
  neMod1<-neModel(y~x0 + x1,
                  family="gaussian",expData=expData,se = "bootstrap", nBoot = 250)
  # rmpw_de <- summary(neMod1)$coef[2,1]
  # rmpw_ie <- summary(neMod1)$coef[3,1]
  # 
  
  rmpw_rmse_nde=sqrt(mean((neMod1$bootRes$t[,2] - true_nde)^2))
  rmpw_rmse_nie=sqrt(mean((neMod1$bootRes$t[,3] - true_nie)^2))
  
  
  rmpw_nde <- data.frame(nde = summary(neMod1)$coef[2,1], 
                         nde_lb = summary(neMod1)$coef[2,1] - 1.96*summary(neMod1)$coef[2,2], 
                         nde_ub = summary(neMod1)$coef[2,1] + 1.96*summary(neMod1)$coef[2,2])
  rmpw_nie <- data.frame(nie = summary(neMod1)$coef[3,1], 
                         nie_lb = summary(neMod1)$coef[3,1] - 1.96*summary(neMod1)$coef[3,2], 
                         nie_ub = summary(neMod1)$coef[3,1] + 1.96*summary(neMod1)$coef[3,2])
  rmpw <- cbind(rmpw_nde, rmpw_nie, rmpw_rmse_nde, rmpw_rmse_nie)
  
  print("rmpw complete")
  #Product Method 
  m.mod=glm(formula = mmodel ,data=obsdat,family="gaussian")
  
  y.mod=lm(formula = ymodel ,data=obsdat)
  
  
  #summary(lm(y~x,data=data)) 
  #prod_toteff = coef(summary(lm(y~x,data=data)))[2,1] 
  prod_de = coef(summary(y.mod))[3,1]
  prod_ie = coef(summary(y.mod))[2,1]*coef(summary(m.mod))[2,1]
  product_boot<- function(x){
    boot_df <-  obsdat[sample(nrow(obsdat), replace=T),]
    m.mod=glm(formula = mmodel ,data=boot_df,family="gaussian")
    
    y.mod=lm(formula = ymodel ,data=boot_df)
    
    
    #summary(lm(y~x,data=data)) 
    #prod_toteff = coef(summary(lm(y~x,data=data)))[2,1] 
    prod_de = coef(summary(y.mod))[3,1]
    prod_ie = coef(summary(y.mod))[2,1]*coef(summary(m.mod))[2,1]
    
    prod_b <- c(NDE = prod_de, NIE = prod_ie)
    return(prod_b)
  }
  prod_b_list <- lapply(1:250, function(x) product_boot(x))
  prod_b_df <- data.frame(do.call(rbind, prod_b_list))
  #Need to calculate squared error too
  prod_rmse_nde=sqrt(mean((prod_b_df$NDE - true_nde)^2))
  prod_rmse_nie=sqrt(mean((prod_b_df$NIE - true_nie)^2))
  
  
  prod_nde <- data.frame(nde = prod_de, 
                         nde_lb = prod_de - 1.96*sqrt(var(prod_b_df$NDE)), 
                         nde_ub = prod_de + 1.96*sqrt(var(prod_b_df$NDE)))
  prod_nie <- data.frame(nie = prod_ie , 
                         nie_lb = prod_ie - 1.96*sqrt(var(prod_b_df$NIE)), 
                         nie_ub = prod_ie + 1.96*sqrt(var(prod_b_df$NIE)))
  prod <- cbind(prod_nde, prod_nie, prod_rmse_nde, prod_rmse_nie)
  print("product complete")
  #IPTW
  library(mediation)
  library(twang)
  library(survey)
  library(ipw)
  library(geepack)
  
  exposure=ipwpoint(x,family="binomial",link="logit",denominator=~z + z2,data=obsdat)
  obsdat$med = as.numeric(cut(obsdat$m, 
                       quantile(obsdat$m, probs=c(0, 0.2, 0.4, 0.6, 0.8, 1) ) , 
                       include.lowest=TRUE))
  mediator=ipwpoint(med,family="ordinal",numerator=~x,denominator=~x + z + z2 ,data=obsdat, link = "logit")
  
  totalweight=exposure$ipw.weights*mediator$ipw.weights
  
  medmod = geeglm(m ~ x , weights = mediator$ipw.weights, family = "gaussian", id = 1:nrow(obsdat), data = obsdat)
  
  msmIPW=geeglm(y~x*m ,weights=totalweight,family="gaussian",id=1:nrow(obsdat),data=obsdat)
  
  msmexp=geeglm(y ~ x,weights=exposure$ipw.weights,family="gaussian",id=1:nrow(obsdat),data=obsdat)
  
  iptw_nde = coef(summary(msmIPW))[2,1]
  iptw_nie = coef(summary(msmIPW))[3,1]*coef(summary(medmod))[2,1]
  iptw_boot<- function(b){
    boot_df <-  obsdat[sample(nrow(obsdat), replace=T),]
    
    exposure=ipwpoint(x,family="binomial",link="logit",denominator=~z + z2,data=boot_df)
    obsdat$med = as.numeric(cut(boot_df$m, 
                      quantile(boot_df$m, probs=c(0, 0.2, 0.4, 0.6, 0.8, 1) ) , 
                      include.lowest=TRUE))
    mediator=ipwpoint(med,family="ordinal",numerator=~x,denominator=~x + z + z2 ,data=boot_df, link = "logit")
    
    totalweight=exposure$ipw.weights*mediator$ipw.weights
    
    medmod = geeglm(m ~ x , weights = mediator$ipw.weights, family = "gaussian", id = 1:nrow(boot_df), data = boot_df)
    
    msmIPW=geeglm(y~x*m,weights=totalweight,family="gaussian",id=1:nrow(boot_df),data=boot_df)
    
    msmexp=geeglm(y ~ x,weights=exposure$ipw.weights,family="gaussian",id=1:nrow(boot_df),data=boot_df)
    
    iptw_nde = coef(summary(msmIPW))[2,1]
    iptw_nie = coef(summary(msmIPW))[3,1]*coef(summary(medmod))[2,1]
    
    iptw_b <- c(NDE = iptw_nde, NIE = iptw_nie)
    return(prod_b)
  }
  iptw_b_list <- lapply(1:250, function(b) product_boot(b))
  iptw_b_df <- data.frame(do.call(rbind, iptw_b_list))
  #Need to calculate squared error too
  iptw_rmse_nde=sqrt(mean((iptw_b_df$NDE - true_nde)^2))
  iptw_rmse_nie=sqrt(mean((iptw_b_df$NIE - true_nie)^2))
  
  iptw_nde <- data.frame(nde = iptw_nde, 
                         nde_lb = iptw_nde - 1.96*sqrt(var(iptw_b_df$NDE)), 
                         nde_ub = iptw_nde + 1.96*sqrt(var(iptw_b_df$NDE)))
  iptw_nie <- data.frame(nie = iptw_nie , 
                         nie_lb = iptw_nie - 1.96*sqrt(var(iptw_b_df$NIE)), 
                         nie_ub = iptw_nie + 1.96*sqrt(var(iptw_b_df$NIE)))
 iptw <- cbind(iptw_nde, iptw_nie, iptw_rmse_nde, iptw_rmse_nie)
  
  
  return(list(gform = gform, iorw = iorw, prod = prod, iptw = iptw, rmpw = rmpw))
}

results_list = list()
for(i in 1:24){
  print(i)
  results_list[[i]] =  bootstrap_mediation(obsdat = Datasets[[i]], dgm_table$true_nde[i], dgm_table$true_nie[i])
}


results_df <-(data.frame(matrix(unlist(results_list), nrow=length(results_list), byrow = T)))
names(results_df) <- c("gform_nde", "gform_nde_lb", "gform_nde_ub", "gform_nie", "gform_nie_lb", "gform_nie_ub", "gform_rmse_nde", "gform_rmse_nie",
                       "iorw_nde", "iorw_nde_lb", "iorw_nde_ub", "iorw_nie", "iorw_nie_lb", "iorw_nie_ub", "iorw_rmse_nde", "iorw_rmse_nie",
                       "prod_nde", "prod_nde_lb", "prod_nde_ub", "prod_nie", "prod_nie_lb", "prod_nie_ub", "prod_rmse_nde", "prod_rmse_nie",
                       "iptw_nde", "iptw_nde_lb", "iptw_nde_ub", "iptw_nie", "iptw_nie_lb", "iptw_nie_ub", "iptw_rmse_nde", "iptw_rmse_nie",
                       "rmpw_nde", "rmpw_nde_lb", "rmpw_nde_ub", "rmpw_nie", "rmpw_nie_lb", "rmpw_nie_ub", "rmpw_rmse_nde", "rmpw_rmse_nie")
write.csv(results_df, "results_mediation_withw_sat.csv")  
