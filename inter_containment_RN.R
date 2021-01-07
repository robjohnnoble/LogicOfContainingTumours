#time to progression, failure and survival for intermittent treatment
#based on simulations (usually it's better to use the exact formulas time_prog_one etc.)
time_intermittent<-function(target_min,times=times_MG, parms=pars_fig,Cmax=2,func=MonroGaffney, num_steps = 400, return_df = FALSE){
  parms["target_max"] <- parms["N0"]
  parms["target_min"] <- target_min
  parms["Cmax"] <- Cmax
  yini_MG<-init_yini(parms)
  #begin slightly below N0 or there's some weird things happening
  #we give it time to reach nmax before starting the treatment
  yini_MG["y1"]<-yini_MG["y1"]*0.99999
  df_MG <- seq_solve(yini_MG, times, func, parms, dose_update_MG, num_steps = num_steps)
  n_vec <- df_MG[which(df_MG$variable=="n"),"value"]
  t_vec <- df_MG[which(df_MG$variable=="n"),"time"]
  n<- length(n_vec)
  #manually find tprog tfail and tsurv 
  tprog <- dic_find(1.0001 * parms["N0"],n_vec,0,n)
  tfail <- dic_find(parms["Nacc"],n_vec,0,n)
  tsurv <- dic_find(parms["Ncrit"],n_vec,0,n)
  #find tprog, tfail and tsurv for a dose C
  if(return_df) return(list(data = df_MG, times = c(t_vec[tprog], t_vec[tfail], t_vec[tsurv])))
  return(c(t_vec[tprog], t_vec[tfail], t_vec[tsurv]))
}

#size of tumour for intermittent treatment, for cmax =500, it's an approximated idealized intermittent treatment
intermittent_treat_a<-function(parms=pars_fig,log=TRUE,cmax=500,times=times_MG){
  parms["Cmax"] <- cmax
  #idealized agg treatment,
  yini_MG <- init_yini(parms)
  #set target max and min for the intermittent containment
  parms["target_max"] <- parms["Ncrit"]*1.2
  parms["target_min"] <- parms["target_max"]/2
  #we want a idAgg so no sensitive cells at the start
  yini_MG["y1"]<-0
  df_MG <- seq_solve(yini_MG, times, MonroGaffney, parms, dose_update_MG)
  #taking only the variable n = tumour size
  b<-df_MG[which(df_MG$variable=="n"),]
  #just to eliminate other unused variables
  b$variable<-factor(b$variable)
  #renaming the variable to r0 (Nmax=r0)
  b$variable<-revalue(b$variable,c("n"="r0"))
  parms["target_max"] <- parms["N0"]
  parms["target_min"] <- parms["target_max"]/2
  yini_MG<-init_yini(parms)
  #begin slightly below N0 or there's some weird things happening
  #we give it time to reach nmax before starting the treatment
  yini_MG["y1"]<-yini_MG["y1"]*0.98
  df_MG <- seq_solve(yini_MG, times, MonroGaffney, parms, dose_update_MG)
  b<-rbind(b,df_MG[which(df_MG$variable=="n"),])
  b$variable<-factor(b$variable)
  b$variable<-revalue(b$variable,c("n"="N0"))
  #Nacc
  parms["target_max"] <- parms["Nacc"]
  parms["target_min"] <- parms["target_max"]/2
  df_MG <- seq_solve(yini_MG, times, MonroGaffney, parms, dose_update_MG)
  b<-rbind(b,df_MG[which(df_MG$variable=="n"),])
  b$variable<-factor(b$variable)
  b$variable<-revalue(b$variable,c("n"="Nacc"))
  #Ncrit
  parms["target_max"] <- parms["Ncrit"]
  parms["target_min"] <- parms["target_max"]/2
  df_MG <- seq_solve(yini_MG, times, MonroGaffney, parms, dose_update_MG)
  b<-rbind(b,df_MG[which(df_MG$variable=="n"),])
  b$variable<-factor(b$variable)
  b$variable<-revalue(b$variable,c("n"="Ncrit"))
  #untreated
  parms["target_max"] <- parms["Ncrit"]*1.2
  parms["target_min"] <- parms["target_max"]/2
  parms["Cmax"]<-0
  df_MG <- seq_solve(yini_MG, times, MonroGaffney, parms, dose_update_MG)
  b<-rbind(b,df_MG[which(df_MG$variable=="n"),])
  b$variable<-factor(b$variable)
  b$variable<-revalue(b$variable,c("n"="More"))
  #just to reset the row numbering
  rownames(b) <- NULL 
  #fig 2.a
  #setting them to NA so we dont have the horizontal lines
  parms["target_max"]<-NA
  parms["target_min"]<-NA
  #the legend and colors follows the same order as the variable appears
  # my_plot<-plot_mod(df = b, xvar = "time", labels = c("r0","N0", "Nacc", "Ncrit", "More"), 
  #                   my_cols = c("black","green", "red", "blue","yellow"),
  #                   logged = TRUE, parms) + 
  #   xlim(times[1],times[2]) + 
  #   scale_y_log10(limits=c(parms["N0"]/3, 8e11))
  # if(log==FALSE){
  #   my_plot<-my_plot + 
  #     ylim(0,parms["Ncrit"]*1.3)
  # }
  # my_plot
  par(mar = c(3, 4, 1, 1))
  plot_mod_base(df = b, labels = c("r0","N0", "Nacc", "Ncrit", "More"),
                logged = TRUE, parms)
}
#Size of tumour and dose spike for intermittent treatment
intermittent_treat_b<-function(parms=pars_fig,log=TRUE,cmax=500,times=times_MG){
  y<-init_yini(pars_fig)
  ###fig 2.b
  parms["Cmax"] <- cmax
  parms["target_max"] <- parms["Nacc"]
  parms["target_min"] <- parms["target_max"]/2
  df_MG <- seq_solve(y, times, MonroGaffney, parms, dose_update_MG)
  #same as bin intermitten_treat_a execpt that we keep the resistant cells too
  c<-df_MG[which(df_MG$variable=="n"),];
  c<-rbind(c,df_MG[which(df_MG$variable=="y2"),])
  c$variable<-revalue(c$variable,c("y2"="R"))
  c$variable<-factor(c$variable)
  ##put the dose spike in the plot  
  c<-rbind(c,df_MG[which(df_MG$variable=="C"),])
  c$variable<-factor(c$variable)
  #scaling the dose spikes
  c[which(c$variable == "C"), "value"] <- 2e8 * c[which(c$variable == "C"), "value"]
  parms["target_max"] <- NA
  parms["target_min"] <- NA
  # plot3 <- plot_mod(c,"time",c("R","n","C"),logged=log,parms)+ggtitle(sprintf("logged, cmax =%s, idealized intermittent",cmax))
  # if(log==FALSE){
  #   plot3<-plot3+ylim(0,8e11)+ggtitle(sprintf("not logged, cmax =%s, idealized intermittent",cmax))
  # }
  # plot3<-plot3+xlim(0,230)
  # plot3
  par(mar = c(3, 4, 1, 1))
  par(mfrow = c(1, 1))
  plot_mod_base(df = b, labels = c("r0","N0", "Nacc", "Ncrit", "More"),
                logged = TRUE, parms)
  
}
calc_tfailmax<-function(Nmax,Nmin,N0,parms=pars_fig){
  #N0 is the value AFTER IdAgg
  b<-c<-0
  tfailmax <- log(log(parms["K"]/N0)/log(parms["K"]/Nmax))
  q<- floor(log(N0/parms["r0"])/log(Nmax/Nmin))
  rqplusone<-parms["r0"]*(Nmax/N0)*(Nmax/Nmin)^q
  if(rqplusone==0){
    b<-(log(N0/parms["r0"])/log(Nmax/Nmin))
    b<-(b-log(log(parms["K"]/Nmin)/log(parms["K"]/Nmax)))
  }else{
    b<-(q*log(log(parms["K"]/Nmin)/log(parms["K"]/Nmax)))
    c<-(log(log(parms["K"]/rqplusone)/log(parms["K"]/Nmax)))
  }
  tfailmax<-tfailmax+b+c
  tfailmax<-tfailmax/parms["rhoG"]
  tfailmax
}
#time to progression for on value of ymax and ymin (doesnt work with a vector)
time_prog_one<-function(ymax,ymin,parms =pars_fig){
  if(ymax<=parms["r0"]){
    return(time_from_to(parms["r0"],parms["N0"]))
  }
  ymin<-max(ymin,parms["r0"])
  if(ymax>parms["N0"]){
    return(0)
  }
  if(ymax==parms["N0"]){
    tx <-calc_tfailmax(ymax,ymin,parms["N0"],parms)
    return(tx)
  }
  tx <-calc_tfailmax(ymax,ymin,ymax,parms)
  
  #Here we treat at maximum dose until we reach Nmax and then it goes like an untreated tumor
  return(tx+time_from_to(ymax,parms["N0"]))
}
#same for tfail
time_fail_one<-function(ymax,ymin,parms=pars_fig){
  if(ymax<=parms["r0"]){
    return(time_from_to(parms["r0"],parms["Nacc"]))
  }
  if(ymax>parms["Nacc"]){
    return(time_from_to(parms["N0"],parms["Nacc"]))
  }
  if(ymax<=parms["N0"]){
    tx <-calc_tfailmax(ymax,ymin,ymax,parms)
    return(tx+time_from_to(ymax,parms["Nacc"]))
  }
  tx <-calc_tfailmax(ymax,ymin,parms["N0"],parms)
  return(tx+time_from_to(ymax,parms["Nacc"]))
}
#same for tsurv
time_surv_one<-function(ymax,ymin,parms=pars_fig){
  if(ymax<=parms["r0"]){
    return(time_from_to(parms["r0"],parms["Ncrit"]))
  }
  if(ymax>parms["Ncrit"]){
    return(time_from_to(parms["N0"],parms["Ncrit"]))
  }
  if(ymax<=parms["N0"]){
    tx <-calc_tfailmax(ymax,ymin,ymin,parms)
    return(tx+time_from_to(ymax,parms["Ncrit"]))
  }
  tx <-calc_tfailmax(ymax,ymin,parms["N0"],parms)
  return(tx+time_from_to(ymax,parms["Ncrit"]))
}
time_inter<-function(ymax,ymin,parms=pars_fig,fun){
  #apply the function to the whole vector
  res<-rep(x=0,times=length(ymax))
  for(i in 1:length(ymax)){
    res[i]<-fun(ymax[i],ymin[i])
  }
  res
}
#plot time to prog/fail/surv for intermittent treatment
time_treat<-function(nref,parms=pars_fig,logged=FALSE){
  nref_t <- nref
  line_N0<-parms["N0"]
  line_ncrit<-parms["Ncrit"]
  line_nacc<-parms["Nacc"]
  title<-"not logged,"
  log<-""
  if(logged == TRUE){
    title<-"logged,"
    log<-"x"
  }
  plot(nref_t,time_inter(nref,nref/2,fun=time_prog_one),log=log,type="l",col="green",yaxs="i",main = paste(title,"idealized intermittent"),ylim=c(1,500))
  lines(nref_t,time_inter(nref,nref/2,fun=time_fail_one),type="l",col="red")
  lines(nref_t,time_inter(nref,nref/2,fun=time_surv_one),type="l",col="blue")
  abline(v=line_N0,lty=2)
  abline(v=line_ncrit,lty=2)
  abline(v=line_nacc,lty=2)
  legend("topleft",lty=c("solid","solid","solid"), col=c("blue","red","green"),legend = c("Tsurv","Tfail","Tprog"),lwd=2);
  
 
  
}
#time to failure for a varying nmin and Nmax = Nacc
time_fail_var_nmin<-function(parms=pars_fig){
  #Fig 2.d
  ymin<-seq(from=0,to=parms["Nacc"],length=100000)
  ymax<-rep(x=parms["Nacc"],times=length(ymin))
  y<-time_inter(ymax,ymin,fun=time_fail_one)
  plot(x=ymin,y,type="l",xlab="nref",ylab="time",main="Tfail with Nmin from 0 to Nmax = Nacc")
  legend("bottomright",lty=c("solid"), col=c("black"),legend = c("Tfail"),lwd=2);
}
#time to progression for a varying nmin and Nmax = Nacc
time_prog_var_nmin<-function(parms=pars_fig){
  ymin<-seq(from=parms["r0"],to=parms["N0"],length=100000)
  #ymax<-rep(x=parms["N0"],times=length(ymin))
  ymax <- parms["N0"]
  #y<-time_inter(ymax,ymin,fun=time_prog_one)
  y <- sapply(ymin, time_prog_one, ymax = ymax)
  yy <- expression(paste("time to progression (", italic("t"[prog]), ")"))
  l1 <- contain_gompertz(parms["K"], parms["r0"], parms["N0"], parms["N0"], parms["rhoG"])
  l1 <- rep(l1, length(y))
  l2_func <- function(Nref) {
    l2_pars <- parms
    l2_pars["N0"] <- Nref
    l2_pars["Nref"] <- Nref
    l2_pars["Nacc"] <- parms["N0"]
    return(heat_fail_id_cont(l2_pars["N0"], l2_pars["r0"] /l2_pars["N0"], l2_pars))
  }
  l2 <- sapply(ymin, l2_func)
  
  par(mgp=c(2.5, 1, 0), las=1)
  par(mar = c(4,4,1,1))
  plot(x=ymin/parms["N0"],y,type="l",xlab = expression(paste(italic("N"[min]), " / ", italic("N")[0])),ylab=yy,
       ylim = c(0, max(l1) / 0.95), lwd = 2)
  legend("bottomright", bty = "n", lty = c(2, 1, 3), lwd = 2, 
         legend = c(expression(paste("ideal containment at ", italic("N")[0])), 
                    "ideal intermittent", 
                    expression(paste("ideal containment at ", italic("N"[min])))))
  lines(x=ymin/parms["N0"], y = l1, lty = 2, lwd = 2)
  lines(x=ymin/parms["N0"], y = l2, lty = 3, lwd = 2)
}
