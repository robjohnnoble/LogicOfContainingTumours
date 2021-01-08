############# HEAT MAP
#http://www.r-graph-gallery.com/27-levelplot-with-lattice/

# growth from N1 to N2:
growth_gompertz <- function(K, N1, N2, rho) {
  top <- log(K / N1)
  bot <- log(K / N2)
  return(as.numeric(1/rho * log(top / bot)))
}
# containment at Nref, with initial resistant population R0:
contain_gompertz <- function(K, R0, N0, Nref, rho) {
  top <- log(N0 / R0)
  bot <- rho * log(K / Nref)
  return(as.numeric(top / bot))
}
#function heat_state_treatment returns the time to state for a given N0 and ratio r0/N0 for the specified treatment
heat_prog_id_cont<-function(N0,ratio,parms=pars_fig){
  contain_gompertz(parms["K"], ratio * N0, N0, N0, parms["rhoG"])
}
# start at N0; contain at Nref; record time until Nacc
heat_fail_id_cont<-function(N0,ratio,parms=pars_fig){
  
  if(N0 > parms["Nref"]) stop(paste("N0 > Nref", N0, parms["Nref"]))
  if(parms["Nref"] > parms["Nacc"]) return(NA) #stop("Nref > Nacc")
  
  #growth from N0 to Nref
  if(N0 == parms["Nref"]) t1 <- 0
  else t1 <- growth_gompertz(parms["K"], N0, parms["Nref"], parms["rhoG"])
  
  #containment at Nref
  t2 <- contain_gompertz(parms["K"], ratio * N0, N0, parms["Nref"], parms["rhoG"])
  
  #growth to Nacc
  t3 <- growth_gompertz(parms["K"], parms["Nref"], parms["Nacc"], parms["rhoG"])
  
  #print(c(t1=t1,t2=t2,t3=t3))
  
  if(min(t1, t2, t3) < 0) return(NA)
  return(t1+t2+t3)
}
heat_prog_id_agg<-function(N0,ratio,parms=pars_fig){
  x<-log(log(parms["K"]/(ratio*N0))/log(parms["K"]/N0))
  if(x==0){
    return(NA)
  }
  return(x/parms["rhoG"])
}
heat_fail_id_agg<-function(N0,ratio,parms=pars_fig){
  #r0 to Nacc
  x<-log(log(parms["K"]/(ratio*N0))/log(parms["K"]/parms["Nacc"]))
  if(x<=0){
    return(NA)
  }
  return(x/parms["rhoG"])
}
heat_prog_cont<-function(N0,ratio,parms=pars_fig){
  first <- heat_prog_id_cont(N0, ratio, parms)
  gN0 <- parms["rhoG"] * log(parms["K"] / N0)
  second <- log(parms["Cmax"]/(parms["Cmax"]-1)) / gN0
  return(max(first - second, 0))
}
heat_fail_cont<-function(N0,ratio,parms=pars_fig){
  first <- heat_fail_id_cont(N0, ratio, parms)
  gNref <- parms["rhoG"] * log(parms["K"] / parms["Nref"])
  second <- log(parms["Cmax"]/(parms["Cmax"]-1)) / gNref
  mini <- heat_fail_id_agg(N0, 1, parms)
  return(max(first - second, mini))
  
  # Nref<-max(parms["Nref"],N0)
  # 
  # t<-log(log(parms["K"]/N0)/log(parms["K"]/Nref))
  # Cmax<-parms["Cmax"]
  # t2<-((1/log(parms["K"]/N0))*(log(1/ratio)-log(Cmax/(Cmax-1))))
  # 
  # t2<-max(0,t2)
  # print(c(t,t2)/parms["rhoG"])
  # return((t+t2)/parms["rhoG"])
}
heat_fail_id_cont_N0 <- function(N0,ratio,parms=pars_fig) {
  t1 <- contain_gompertz(parms["K"], ratio * parms["N0"], N0, N0, parms["rhoG"])
  t2 <- growth_gompertz(parms["K"], N0, parms["Nacc"], parms["rhoG"])
  return(t1 + t2)
}
heat_surv_id_cont_Nacc <- function(N0,ratio,parms=pars_fig) {
  t1 <- growth_gompertz(parms["K"], N0, parms["Nacc"], parms["rhoG"])
  t2 <- contain_gompertz(parms["K"], ratio * parms["N0"], parms["N0"], parms["Nacc"], parms["rhoG"])
  t3 <- growth_gompertz(parms["K"], parms["Nacc"], parms["Ncrit"], parms["rhoG"])
  return(t1 + t2 + t3)
}
heat_prog_agg<-function(N0,ratio,parms=pars_fig){
  parms["r0"]<-ratio*N0
  parms["N0"]<-N0
  yini<-init_yini(parms)
  yini["C"]<-parms["Cmax"]
  #no need to go too high
  t <-seq(from = 0, to=200,length=1000)
  #simulation of mtd
  mtd<-ode(y=yini,parms=parms,times=t,func=MonroGaffney)
  #manually searching for the progression time
  t2<-dic_find(N0,mtd[,2]+mtd[,3],0,length(mtd[,2]))
  return(t[t2])
}
heat_fail_agg<-function(N0,ratio,parms=pars_fig){
  parms["r0"]<-ratio*N0
  parms["N0"]<-N0
  yini<-init_yini(parms)
  yini["C"]<-parms["Cmax"]
  #goes a bit higher than for tprog
  t <-seq(from = 0, to=500,length=1000)
  mtd<-ode(y=yini,parms=parms,times=t,func=MonroGaffney)
  t2<-dic_find(parms["Nref"],mtd[,2]+mtd[,3],0,length(mtd[,2]))
  return(t[t2])
}
heat_fail_agg_extended_model<-function(N0,ratio,parms,beta=NA,gamma=NA){
  parms["r0"]<-ratio*N0
  parms["N0"]<-N0
  if(is.na(beta)) beta <- parms["beta"]
  else parms["beta"] <- beta
  if(is.na(gamma)) gamma <- parms["Ks"] / parms["Kr"]
  else parms["Kr"] <- parms["Ks"] / gamma
  yini<-init_yini(parms)
  yini["C"]<-parms["Cmax"]
  #no need to go too high
  t_max <- 2 * heat_fail_agg_extended_approx(N0,ratio,parms,beta,gamma)
  t <- seq(0, t_max, length = 1e4)
  #simulation of mtd
  mtd<-ode(y=yini,parms=parms,times=t,func=MonroGaffneyExtended)
  # plot(mtd[,2]+mtd[,3] ~ mtd[,1], type = "l", log = "y", ylim = c(1e6, 2e12), ylab = "")
  #manually searching for the progression time
  t2<-dic_find(parms["Nacc"],mtd[,2]+mtd[,3],0,length(mtd[,2]))
  # abline(v = t[t2], lty = 2)
  # abline(h = parms["Nacc"], lty = 2)
  # print(paste("MTD", beta, gamma, t[t2]))
  return(t[t2])
}
heat_fail_cont_extended_model<-function(N0,ratio,parms,beta=NA,gamma=NA,level=NA,endpoint=NA){
  parms["r0"]<-ratio*N0
  parms["N0"]<-N0
  if(is.na(beta)) beta <- parms["beta"]
  else parms["beta"] <- beta
  if(is.na(gamma)) gamma <- parms["Ks"] / parms["Kr"]
  else parms["Kr"] <- parms["Ks"] / gamma
  
  if(is.na(level)) level <- parms["Nacc"]
  if(is.na(endpoint)) endpoint <- parms["Nacc"]
  
  parms["Nref"] <- level
  # print(parms)
  
  yini<-init_yini(parms)
  yini["C"]<-0
  #no need to go too high
  t_max <- 2 * heat_fail_cont_extended_approx(N0,ratio,parms,beta,gamma,level,endpoint)
  t_max <- max(t_max, 1000)
  #print(t_max)
  t <- seq(0, t_max, length = 1e4)
  #print(paste("cont", beta, gamma))
  #print(yini)
  #print(parms)
  #simulation of cont
  untreat<-lsoda(y=yini,parms = parms, times = t, func = MonroGaffneyExtended, maxsteps = 1e4, rtol = 1e-4)
  #print(tail(untreat))
  cont<-continuous(level,t,parms,MonroGaffneyExtendedCont,untreat,method="lsoda")
  # print(tail(cont))
  # plot(cont[,2]+cont[,3] ~ cont[,1], type = "l", col = "black", log = "y",
  #      ylim = c(1e6, 2e12), ylab = "", main = c(beta, gamma))
  # lines(cont[,2] ~ cont[,1], col = "blue")
  # lines(cont[,3] ~ cont[,1], col = "red")
  #manually searching for the progression time
  t2<-dic_find((1+1e-3)*endpoint,cont[,2]+cont[,3],0,length(cont[,2]))
  # abline(h = endpoint, lty = 2)
  # abline(v = t[t2], lty = 2)
  # ww <- cont[,2] + cont[,3]
  # diff1 <- ww - lag(ww)
  # diff2 <- lag(diff1)
  # print(max(diff1, na.rm = TRUE))
  # t2 <- min(which(ww > parms["Nacc"] & diff1 > 20 * diff2))
  # print(paste("cont", beta, gamma, t[t2]))
  # print(c(cont[t2,2] + cont[t2,3], cont[t2+1,2] + cont[t2+1,3], cont[t2+2,2] + cont[t2+2,3]))
  #plot(dose(cont[,2], cont[,3], parms["Kr"], parms["Ks"], parms["alpha"], parms["beta"], parms["Nacc"], parms["Cmax"]) ~ cont[,1])
  if(t2 == length(length(cont[,2]))) return(NA)
  return(t[t2])
}
heat_fail_agg_extended_approx <- function(N0,ratio,parms,beta=NA,gamma=NA){
  parms["r0"]<-ratio*N0
  parms["N0"]<-N0
  if(is.na(gamma)) gamma <- parms["Ks"] / parms["Kr"]
  else parms["Kr"] <- parms["Ks"] / gamma
  
  t <- 1/parms["rhoG"] * log((log(parms["Kr"]/parms["r0"]) / (log(parms["Kr"]/parms["Nacc"]))))
  return(t)
}

# tumour size is contained at "level" and time is measured until tumour size reaches "endpoint"
heat_fail_cont_extended_approx <- function(N0,ratio,parms,beta=NA,gamma=NA,level=NA,endpoint=NA){
  parms["r0"]<-ratio*N0
  parms["N0"]<-N0
  if(is.na(beta)) beta <- parms["beta"]
  else parms["beta"] <- beta
  if(is.na(gamma)) gamma <- parms["Ks"] / parms["Kr"]
  else parms["Kr"] <- parms["Ks"] / gamma
  
  if(is.na(level)) level <- parms["Nacc"]
  if(is.na(endpoint)) endpoint <- parms["Nacc"]
  
  t1 <- 1/parms["rhoG"] * log((log(parms["Ks"]/parms["N0"]) / (log(parms["Ks"]/level))))
  t1 <- as.numeric(t1)
  
  R1 <- parms["r0"] * level / parms["N0"] * exp(-t1 * parms["rhoG"] * log(beta * parms["Ks"]/parms["Kr"]))
  # print(paste("R1 components:", parms["r0"], level, parms["N0"], -t1, parms["rhoG"], log(beta * parms["Ks"]/parms["Kr"])))
  # print(paste(parms["Ks"], parms["Kr"]))
  R1 <- as.numeric(R1)
  
  integrand <- function(R) {1 / (R * log(parms["Kr"] / (R + beta * (level - R))))}
  t2 <- 1/parms["rhoG"] * integrate(integrand, lower = R1, upper = endpoint)$value
  t2 <- as.numeric(t2)
  
  return(t1 + t2)
}
#logistic model:
heat_prog_logistic_id_cont<-function(N0,ratio,parms=pars_fig){
  1/parms["rho"] * log(1/ratio) / (1 - N0 / parms["K"])
}
heat_prog_logistic_id_agg<-function(N0,ratio,parms=pars_fig){
  1/parms["rho"] * log((ratio * N0 - parms["K"]) / (ratio * (N0 - parms["K"])))
}
heat_fail_logistic_id_cont<-function(N0,ratio,parms=pars_fig){
  t1 <- 1/parms["rho"] * log((parms["Nacc"] * (N0 - parms["K"])) / (N0 * (parms["Nacc"] - parms["K"])))
  t2 <- 1/parms["rho"] * log(1/ratio) / (1 - parms["Nacc"] / parms["K"])
  return(t1 + t2)
}
heat_fail_logistic_id_agg<-function(N0,ratio,parms=pars_fig){
  1/parms["rho"] * log(parms["Nacc"] * (ratio * N0 - parms["K"]) / (ratio * N0 * (parms["Nacc"] - parms["K"])))
}
heat_fail_logistic_untreat<-function(N0,ratio,parms=pars_fig){
  1/parms["rho"] * log((parms["Nacc"] * (N0 - parms["K"])) / (N0 * (parms["Nacc"] - parms["K"])))
}
#von Bertalanffy:
heat_fail_vonB_id_cont<-function(N0,ratio,parms,gamma = 1/3){
  t1 <- 1 / (parms["rho"] * parms["K"]^-gamma * gamma) * log((1 - (N0 / parms["K"])^gamma ) / (1 - (parms["Nacc"] / parms["K"])^gamma))
  t2 <- log(1 / ratio) / (parms["rho"] * (parms["Nacc"]^-gamma - parms["K"]^-gamma))
  return(t1 + t2)
}
heat_fail_vonB_id_agg<-function(N0,ratio,parms,gamma = 1/3){
  1 / (parms["rho"] * parms["K"]^-gamma * gamma) * log((1 - (ratio * N0 / parms["K"])^gamma ) / (1 - (parms["Nacc"] / parms["K"])^gamma))
}
heat_fail_vonB_untreated<-function(N0,ratio,parms,gamma = 1/3){
  1 / (parms["rho"] * parms["K"]^-gamma * gamma) * log((1 - (N0 / parms["K"])^gamma ) / (1 - (parms["Nacc"] / parms["K"])^gamma))
}
#power law (including superexponential):
heat_fail_power_id_cont<-function(N0,ratio,parms,gamma = -1/3){
  t1 <- 1 / (parms["rho"] * gamma) * (parms["Nacc"]^gamma - N0^gamma)
  t2 <- log(1 / ratio) / (parms["rho"] * parms["Nacc"]^-gamma)
  return(t1 + t2)
}
heat_fail_power_id_agg<-function(N0,ratio,parms,gamma = -1/3){
  1 / (parms["rho"] * gamma) * (parms["Nacc"]^gamma - (ratio * N0)^gamma)
}
heat_fail_power_untreated<-function(N0,ratio,parms,gamma = -1/3){
  1 / (parms["rho"] * gamma) * (parms["Nacc"]^gamma - N0^gamma)
}
#exponential:
heat_fail_exp_id_cont<-function(N0,ratio,parms,gamma = -1/3){
  t1 <- 1 / parms["rho"] * log(parms["Nacc"] / N0)
  t2 <- log(1 / ratio) / parms["rho"]
  return(t1 + t2)
}
heat_fail_exp_id_agg<-function(N0,ratio,parms,gamma = -1/3){
  1 / parms["rho"] * log(parms["Nacc"] / (ratio * N0))
}
heat_fail_exp_untreated<-function(N0,ratio,parms,gamma = -1/3){
  1 / parms["rho"] * log(parms["Nacc"] / N0)
}

#given two functions of type heat_state_treatment and the length of the heatmap
#returns a data frame with the ratio of time
#increasing the length increases the precision but also greatly increases the computation time
heat_fun1_fun2<-function(length,fun1,fun2,parms=pars_fig){
  #N0<-seq(from=1e6,to=parms["Nacc"],length=length)
  #r0_on_N0<-seq(from=1e-6,to=1,length=length)
  N0<-seq(from=log10(1e6),to=log10(0.99*parms["Nacc"]),length=length)
  r0_on_N0<-seq(from=log10(1e-6),to=log10(0.1),length=length)
  N0<-exp(log(10)*N0)
  r0_on_N0<-exp(log(10)*r0_on_N0)
  res<-NULL
  #create a grid of 2 columbs and length*length row
  data_res<-expand.grid(X=N0,Y=r0_on_N0)
  for(i in 1:length){
    for(j in 1:length){
      #calculate time for each function
      data_res$Z[(i-1)*length+j]<-fun1(N0[j],r0_on_N0[i],parms=parms)
      data_res$V[(i-1)*length+j]<-fun2(N0[j],r0_on_N0[i],parms=parms)
    }
  }
  #ratio
  data_res$Z<-data_res$Z/data_res$V
  data_res$V<-NULL
  
  return(data_res)
}

# variation of beta and gamma:
heat_fun1_fun2_extended<-function(length,fun1,fun2,parms=pars_fig){
  beta<-seq(from=0,to=4,length=length)
  gamma<-seq(from=1e-9,to=10,length=length)
  N0 <- parms["N0"]
  r0_on_N0 <- parms["r0"] / parms["N0"]
  res<-NULL
  #create a grid of 2 columbs and length*length row
  data_res<-expand.grid(X=beta,Y=gamma)
  for(i in 1:length){
    for(j in 1:length){
      #calculate time for each function
      if(beta[j] > 1 && gamma[i] < parms["Ks"]/parms["Nacc"] && beta[j] * gamma[i] > parms["Ks"]/parms["Nacc"]) {
        data_res$Z[(i-1)*length+j]<-NA
        data_res$V[(i-1)*length+j]<-NA
      } else {
      data_res$Z[(i-1)*length+j]<-fun1(N0,r0_on_N0,parms=parms,beta[j],gamma[i])
      data_res$V[(i-1)*length+j]<-fun2(N0,r0_on_N0,parms=parms,beta[j],gamma[i])
      }
    }
  }
  #ratio
  data_res$Z<-data_res$Z/data_res$V
  data_res$V<-NULL
  
  return(data_res)
}

colfunc <- colorRampPalette(c("#f2f2ff", "#0033dd"))

##data_res must be a dataframe with a column X, Y and Z, Z~Y*X
#plot the heatmap, log is for both X and Y
plot_heatmap<-function(data_res,xlabel,ylabel,title,scales=list(),logged=FALSE,min=1,max=3,
                       contour_gap=0.1, myseq = NULL, mylabs = NULL){
  if(logged == TRUE){
    data_res$X<-log10(data_res$X)
    data_res$Y<-log10(data_res$Y)
  } else if(logged == "x") {
    data_res$X<-log10(data_res$X)
  } else if(logged == "y") {
    data_res$Y<-log10(data_res$Y)
  }
  myTheme <- BTCTheme()
  myTheme$panel.background$col = '#ddff00'
  colk <- list(height = 0.7)
  if(!is.null(myseq)) colk <- list(height = 0.7, at = myseq, labels = mylabs)
  pp <- levelplot(Z ~ Y*X,
                  data=data_res,
                  par.settings = myTheme,
                  xlab=xlabel,
                  ylab=ylabel,
                  contour = TRUE,
                  col.regions = colfunc(100),
                  at=seq(min, max, by=contour_gap),
                  scales=scales,
                  colorkey = colk,
                  legend=list(top=list(fun=grid::textGrob("relative\nbenefit", y=-0.9, x=1.2))))
  print(pp)
}
