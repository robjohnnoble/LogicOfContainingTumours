#continuous containment at target with cmax = parms["Cmax"]
#func is the model used, untreat the untreated tumour
continuous<-function(target,x,parms,func,untreat,method="lsoda",tau1 = 0, tau2 = 0){
  ##create the untreated growth part to containment level (target)
  #find the index at which we reach containment level (target)
  result <- find_time(untreat,target)
  # print(c("result ", result))
  # print(target)
  # print(untreat[result-1, ])
  #take the corresponding untreated part
  new <- untreat[1:(result-1),]
  #treat the cases if new is 1 or 2 dimensional (could probably be done without the if)
  if(length(new)>5){
    #take the first 4 columns (5th one is the dose for resistant cell, not used right now)
    new <- new[,1:4]
    #if target is reached before the end of our time window
    if(length(new)<length(x)){
      #stitch the untreated part with the ideal containment one
      yini_cont <- c(y1=0,y2=0,C = 0, ideal_contain = 0)
      #initialisation of yini_cont(sensible = y1 and resistant = y2)
      #new is : time; sensible; resistant; treatment dose. in terms of columns
      yini_cont["y1"]<-new[result-1,2]
      yini_cont["y2"]<-new[result-1,3]
      #initialisation of dose C
      C<-min(parms["Cmax"],(1+ yini_cont["y2"]/yini_cont["y1"]))
      yini_cont["C"]<-C
      #continuous containment part
      if(identical(func, MonroGaffneyExtendedCont)) out <- lsoda(y=yini_cont, parms = parms, times = x[result:length(x)], func =func,maxsteps = 1e5, rtol = 1e-4, level = target)
      else if(tau1 > 0 | tau2 > 0) out <- lsoda(y=yini_cont, parms = parms, times = x[result:length(x)], func =func,maxsteps = 1e5, rtol = 1e-4, tau1 = tau1, tau2 = tau2)
      else out <- lsoda(y=yini_cont, parms = parms, times = x[result:length(x)], func =func,maxsteps = 1e5, rtol = 1e-4)
      out <- out[,1:4]
      #stitching to the untreated part
      out<-rbind(new,out)
    }else{
      out<-new
    }
  }else{
    #do the same but on 1 dimension
    new <- new[1:4]
    yini_cont <- c(y1=0,y2=0, C = 0, ideal_contain = 0)
    yini_cont["y1"]<-new[2]
    yini_cont["y2"]<-new[3]
    C<-min(parms["Cmax"],1+new[3]/new[2])
    yini_cont["C"]<-C
    
    if(identical(func, MonroGaffneyExtendedCont)) out <- lsoda(y=yini_cont, parms = parms, times = x[result:length(x)], func = func,maxsteps = 1e5, rtol = 1e-4, level = target)
    else if(tau1 > 0 | tau2 > 0) out <- lsoda(y=yini_cont, parms = parms, times = x[result:length(x)], func = func,maxsteps = 1e5, rtol = 1e-4, tau1 = tau1, tau2 = tau2)
    else out <- lsoda(y=yini_cont, parms = parms, times = x[result:length(x)], func = func,maxsteps = 1e5, rtol = 1e-4)
    out <- rbind(new,out[,1:4])
  }
  out
}
#size of tumour for continuous containment
continuous_treat<-function(x,cmax1,cmax2,parms = pars_fig,logged = FALSE,png=FALSE){
  par(mfrow = c(2, 2))
  #initialisation of yini
  yini_untreat<-init_yini(parms)
  ###Idealized continuous containment, cmax1 is usually at 500  
  parms["Cmax"]<-cmax1
  untreat<-ode(y=yini_untreat,parms = parms, times = x, func = MonroGaffney)
  #for the MTD
  yini_untreat["C"]<-cmax2
  mtd<-ode(y=yini_untreat,parms=parms,times=x,func=MonroGaffney)
  #for Idealized containment
  crit<-continuous(parms["Ncrit"],x,parms,MonroGaffneyCont,untreat)
  acc<-continuous(parms["Nacc"],x,parms,MonroGaffneyCont,untreat)
  N0<-continuous(parms["N0"],x,parms,MonroGaffneyCont,untreat)
  yini_untreat["y1"]<-0
  #IdAgg
  zero<-ode(y=yini_untreat,parms = parms, times = x, func = MonroGaffney)
  #no treatment
  max<-continuous(parms["Ncrit"]*10,x,parms,MonroGaffneyCont,untreat)
  #continuous containment
  parms["Cmax"]<-cmax2
  acc_cmax<-continuous(parms["Nacc"],x,parms,MonroGaffneyCont,untreat)
  #adding the sensible and resistant together
  mtd_b<-mtd[,2]+mtd[,3]
  crit_a<-(crit[,2]+crit[,3])
  acc_bis<-(acc[,3])
  acc_4<-(acc[,4])
  acc_a<-(acc[,2]+acc[,3])
  N0_a<-( N0[,2]+ N0[,3])
  zero_a<-(zero[,2]+zero[,3])
  acc_cmax_bis<-(acc_cmax[,3])
  acc_cmax_4<-acc_cmax[,4]
  acc_cmax_a<-(acc_cmax[,2]+acc_cmax[,3])
  max_a<-(max[,2]+max[,3])
  title <- "not logged"
  log<-""
  
  if(logged == TRUE){
    title <- "logged"
    #axis we want to be log
    log<-"y"
  }
  ###Fig 1.a
  
  #if png ==TRUE, create automatically the figure. Useful when you have to make a lot of r0 and cmax variations
  if(png==TRUE){
    text<-sprintf("fig_1_a%s_%s_%s.png",if(logged==TRUE) "_log" else "",parms["r0"],cmax2)
    png(text,width=300,height=220)
  }
  #Ncrit
  par(mar = c(3, 4, 1, 1))
  
  plot(x,crit_a,type="l",col="orange",
       log=log, 
       ylim=c(parms["N0"]/2,2e12), 
       ylab = "number of tumour cells", 
       xlab = "", 
       lwd=2, 
       xaxt = "n", 
       yaxt = "n")
  mtext("time (days)", 1, 2, cex = 0.8)
  mtext("a", adj=-0.25, line=0)
  abline(h = pars_fig["Ncrit"], col = "red", lty = 3)
  axis(2, pars_fig["Ncrit"], labels = "Ncrit", las = 2, col = "red", col.axis = "red")
  axis(2, c(10^10, 10^11, 10^12), labels = parse(text=c("10^10", "10^11", "10^12")), las = 2)
  axis(1, 200*(0:3))
  #Nacc
  lines(x,acc_a,type="l",col="magenta", lwd=2)
  #N0
  lines(x,N0_a,type="l",col="cyan3", lwd=2)
  lines(x,zero_a,type="l",col="black", lty=2, lwd=2)
  #over Ncrit
  lines(x,max_a,type="l",col="black", lwd=2)
  #you can change the position of the legend by changing bottomright to bottomleft, topright or topleft
  legend("bottomright", col=c("black","black","cyan3","magenta","orange"), 
                              lty = c(1, 2, rep(1, 3)),
         legend = c("untreated",
                    "ideal MTD",
                    "contain at N0",
                    "contain at Nacc",
                                 "contain at Ncrit"),
         lwd=2, 
         bty = "n");
  if(png==TRUE){
    dev.off()
  }
  
  #### fig 1_b
  if(png==TRUE){
    text<-sprintf("fig_1_d%s_%s_%s.png",if(logged==TRUE) "_log" else "",parms["r0"],cmax2)
    png(text,width=300,height=220)
  }
  par(mar = c(3, 4, 1, 1))
  plot(x,(acc_cmax_4),type="l",
       xlim=c(0,1000), 
       ylim = c(0, 2.5),
       ylab = "dose (C)", 
       xlab = "", 
       xaxt = "n", 
       yaxt = "n",
       lwd = 2)
  axis(2, 0:2, las = 2)
  axis(1, c(0, 500, 1000))
  mtext("time (days)", 1, 2, cex = 0.8)
  mtext("b", adj=-0.25, line=0)
  lines(x,(acc_4),type="l", lty = 2, lwd = 2)
  legend("bottomright", lty = 1:2,
         legend = c(sprintf("Cmax = %s",cmax2),sprintf("Cmax = %s",cmax1)),
         lwd=2, bty = "n")
  if(png==TRUE){
    dev.off()
  }
  ###Fig 1.c Nacc = 7e10
  if(png==TRUE){
    text<-sprintf("fig_1_b%s_%s_%s.png",if(logged==TRUE) "_log" else "",parms["r0"],cmax2)
    png(text,width=300,height=220)
  }
  par(mar = c(3, 4, 1, 1))
  plot(x,acc_cmax[,2],log=log,type="l",
       lty=2,
       col="grey50",
       xlim=c(0,600),
       ylim=c(parms["N0"]/2,2e12),
       lwd=2, 
       ylab = "number of tumour cells", 
       xlab = "", 
       xaxt = "n", 
       yaxt = "n")
  mtext("time (days)", 1, 2, cex = 0.8)
  mtext("c", adj=-0.25, line=0)
  abline(h = pars_fig["Ncrit"], col = "red", lty = 3)
  axis(2, pars_fig["Ncrit"], labels = "Ncrit", las = 2, col = "red", col.axis = "red")
  axis(2, 10^(10:12), labels = parse(text=c("10^10", "10^11", "10^12")), las = 2)
  axis(1, 200*(0:3))
  lines(x,acc_cmax_bis,lty=4,col="grey75", lwd=2)
  lines(x,acc_a,col="magenta", lwd=2)
  lines(x,acc_cmax_a,col="darkmagenta", lwd=2)
  # lines(x,(mtd_b),lty="solid",col="yellow")
  # lines(x,(zero_a),lty="solid",col="purple")
  # abline(h=(parms["Ncrit"]))
  legend("bottomright",
         col=c("magenta","darkmagenta","grey50","grey75"),
         lty = c(1, 1, 2, 4),
         lwd=2,
         bty = "n",
         legend = c("ideal containment",
                    sprintf("Cmax = %s",cmax2),
                    sprintf("S (Cmax = %s)",cmax2),
                    sprintf("R (Cmax = %s)",cmax2)))
  if(png==TRUE){
    dev.off()
  } 
  
  #### fig 1_d
  if(png==TRUE){
    text<-sprintf("fig_1_c%s_%s.png",if(logged==TRUE) "_log" else "",parms["r0"])
    png(text,width=300,height=220)
  }
  time_cont(nref,logged=TRUE, png = FALSE)
  mtext("d", adj=-0.25, line=0)
  if(png==TRUE){
    dev.off()
  } 
}

# containment for three different values of cmax (c1, c2, c3)
# similar to the "constant_treat" function
continuous_treat_alt<-function(parms=pars_fig,logged=FALSE,c1=0.5,c2=1,c3=2, c2extratext = "", x_cex = 0.7){
  x=seq(from=0,to=2200,length=1000)
  yini<-init_yini(parms)
  yini_untreat<-init_yini(parms)
  parms["Cmax"] <- 500
  yini["C"]<-0
  
  par(mfrow = c(1, 3))
  layout(matrix(c(1,1,2), 1, 3, byrow = TRUE))
  par(mar = c(3, 4, 1, 1))
  
  untreat<-ode(y=yini_untreat,parms = parms, times = x, func = MonroGaffney)
  yini_untreat["C"]<-500
  parms_alt <- parms
  parms_alt["Cmax"] <- 500
  N0<-continuous(parms["N0"],x,parms_alt,MonroGaffneyCont,untreat)
  N0_a<-( N0[,2]+ N0[,3])
  
  #no treatment
  untreated<-ode(y=yini,parms = parms, times = x, func = MonroGaffney)
  #IdAGg
  yini["C"]<-500
  iat <- ode(y=yini,parms = parms, times = x, func = MonroGaffney)
  #containment subject to Cmax
  
  parms_alt["Cmax"] <- c1
  const1<-continuous(parms["N0"],x,parms_alt,MonroGaffneyCont,untreat)
  
  parms_alt["Cmax"] <- c2
  const2<-continuous(parms["N0"],x,parms_alt,MonroGaffneyCont,untreat)
  
  parms_alt["Cmax"] <- c3
  const3<-continuous(parms["N0"],x,parms_alt,MonroGaffneyCont,untreat)
  
  title<-"not logged"
  log<-""
  if(logged == TRUE){
    title<-"logged"
    log<-"y"
  }
  #taking the tumour size for each treatment
  untreated<-untreated[,2]+untreated[,3]
  iat<-iat[,2]+iat[,3]
  N1<-const1[,2]+const1[,3]
  N2<-const2[,2]+const2[,3]
  N3<-const3[,2]+const3[,3]
  
  # first panel:
  plot(x=x,y=untreated, type = "l",
       log = "y", 
       xlim = c(0, 1000), 
       ylim=c(parms["N0"]/2, 2e12),
       lwd = 2, 
       xaxt = "n", 
       yaxt = "n",
       xlab = "",
       ylab = "number of tumour cells",
       col = "black")
  axis(2, pars_fig["Ncrit"], labels = expression(paste(italic("N" [crit]))), las = 2, col = "red", col.axis = "red")
  axis(2, c(10^10, 10^11, 10^12), labels = parse(text=c("10^10", "10^11", "10^12")), las = 2)
  axis(1, 500*(0:3))
  mtext("time (days)", 1, 2, cex = x_cex)
  mtext("a", adj=-0.12, line=-1)
  lines(x=x[-1],y=iat[-1],type="l",col="black",lty=2,lwd = 2)
  lines(x=x, y = N0_a,lty="solid",col="cyan3",lwd = 2)
  lines(x=x,y=N3,lty="solid",col="yellowgreen",lwd = 2)
  lines(x=x,y=N2,lty="solid",col="cyan4",lwd = 2)
  lines(x=x,y=N1,lty="solid",col="blue",lwd = 2)
  abline(h = parms["Ncrit"], col = "red", lty = 3)
  substitute(expr = paste(italic("C"), " = ", c1))
  legend("bottomright",
         bty = "n",
         lty=c("solid","dashed","solid","solid","solid","solid"), 
         col=c("black","black","blue","cyan4","yellowgreen","cyan3"),
         legend = c("untreated",
                    expression("ideal MTD"),
                    substitute(expr = paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = ", c1, ")"), env = base::list(c1 = c1)),
                    substitute(expr = paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = ", c2, ")"), env = base::list(c2 = c2)),
                    substitute(expr = paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = ", c3, ")"), env = base::list(c3 = c3)),
                    substitute(expr = paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = ", infinity, ")"))
         ),
         lwd=2)
  
  #### dose panel:
  par(mar = c(3, 4, 1, 1))
  plot(x[-1],N0[-1,4],type="l",
       xlim=c(0,1000), 
       ylim = c(0, 2.25),
       ylab = "dose (C)", 
       xlab = "", 
       xaxt = "n", 
       yaxt = "n",
       lwd = 2, 
       col = "cyan3")
  axis(2, 0:2, las = 2)
  axis(1, c(0, 500, 1000))
  mtext("time (days)", 1, 2, cex = x_cex)
  mtext("b", adj=-0.3, line=-1)
  lines(x[-1],const3[-1,4],type="l", lwd = 2, col = "yellowgreen")
  lines(x[-1],const2[-1,4],type="l", lwd = 2, col = "cyan4")
  lines(x[-1],const1[-1,4],type="l", lwd = 2, col = "blue")
  legend("bottomright",
         bty = "n",
         lty=c("blank","solid","solid","solid","solid"),
         col=c("white","blue","cyan4","yellowgreen","cyan3"),
         legend = c(expression(""),
                    substitute(expr = paste(italic("C"[max]), " = ", c1), env = base::list(c1 = c1)),
                    substitute(expr = paste(italic("C"[max]), " = ", c2), env = base::list(c2 = c2)),
                    substitute(expr = paste(italic("C"[max]), " = ", c3), env = base::list(c3 = c3)),
                    substitute(expr = paste(italic("C"[max]), " = ", infinity))
         ),
         lwd=2)
}


#In the 3 following functions, y = Nref
#Tsurv for continuous containment
time_surv_cont <-function(y,parms=pars_fig, true_parms = pars_fig){
  a<-b<-c<-0
  #growth from N0 to y = Nref
  a<-log(log(parms["K"]/parms["N0"])/log(parms["K"]/y))
  #there's no growth from N0 to Nref if Nref < N0
  a[which(y<=parms["N0"])] <- 0
  #if Nref <r0 we force it at r0
  y[y<parms["r0"]]<-parms["r0"]
  #if r0<=y<N0, the "new" N0 after initial IdAgg is y
  N0 <- y
  N0[y>=parms["N0"]] <- parms["N0"]
  b<-log(N0/parms["r0"])/log(parms["K"]/y)#stabilization at nref
  c<-log(log(parms["K"]/y)/log(parms["K"]/parms["Ncrit"]))# growth from nref to Ncrit
  res<-(a+b+c)/parms["rhoG"]
  #if we stabilise at >Ncrit, it's the same as a untreated tumour
  res[y>parms["Ncrit"]]<-log(log(parms["K"]/parms["N0"])/log(parms["K"]/parms["Ncrit"]))/parms["rhoG"]
  return(res)
}
#Tfail
time_fail_cont<-function(y,parms=pars_fig){
  #We just replace Ncrit by Nacc
  #I use Ninit and R0 as the one we used before
  parms_prime <- parms
  parms_prime["Ncrit"]<-parms_prime["Nacc"]
  return(time_surv_cont(y,parms_prime, parms))
}
#Tprog
time_prog_cont<-function(y,parms=pars_fig){
  #same as before, replacing Ncrit by N0 = N0
  parms_prime <- parms
  parms_prime["Ncrit"]<-parms_prime["N0"]
  return(time_surv_cont(y,parms_prime))
}
#plotting of Time to prog/fail/surv for continuous containment
time_cont<-function(nref,parms=pars_fig,logged=FALSE, png = FALSE){
  #the distinction was useful before but not anymore
  nref_t<-nref
  title <- "not logged"
  log<-""
  if(logged == TRUE){
    log<-"x"
    title <- "logged"
  }
  par(mar = c(3, 4, 1, 1.8))
  plot(nref_t, log=log, time_surv_cont(nref,parms),
       type="l", col="red", yaxs="i",
       #xlim=c(pars_fig["r0"], pars_fig["Ncrit"]),
       xlim = c(1e8, 2e12),
       ylim=c(-10, 1800), lwd=2, 
       xaxt = "n", 
       yaxt = "n", 
       ylab = "time (days)")
  abline(v = pars_fig["Ncrit"], col = "red", lty = 3)
  axis(1, pars_fig["Ncrit"], labels = expression(paste(italic("N" [crit]), "         ")), col = "red", col.axis = "red")
  axis(1, 10^(2*3:6), labels = parse(text=c("10^6", "10^8", "10^10", "10^12")))
  axis(2, 500*(0:5), las = 2)
  mtext("containment size (number of tumour cells)", 1, 2, cex = 0.7)
  lines((nref_t),time_fail_cont(nref,parms),type="l",col="green3", lwd=2)
  lines((nref_t),time_prog_cont(nref,parms),type="l",col="blue", lwd=2)
  legend("topleft",
         lty=c("solid","solid","solid"), 
         col=c("blue","green3","red"),
         legend = c(expression(paste(italic("t" [prog]))),
                    expression(paste(italic("t" [fail]))),
                    expression(paste(italic("t" [surv])))
                    ),
         lwd=2, 
         bty = "n");
  if(png==TRUE){
    dev.off()
  } 
}
#plot tumour size for IdAgg_N and MTD_NRS
iat_vs_mtd<-function(x,cmax,parms=pars_fig,logged=FALSE,xlim = c(0,230),ylim= NA){
  ###here IAT is IdAgg
  yini<-init_yini(parms)
  #a general definition of ylim so that the figure isnt too hard to decrypt
  if(is.na(ylim)){
    ylim<-c((parms["r0"]/2),(5e11))
  }
  ###mtd => just put C to cmax and use the normal MonroGaffney one
  yini["C"]<-cmax
  mtd <-ode(y=yini ,parms = parms, times = x, func = MonroGaffney)
  #It's not a "real" IdAgg treatment, it's a MTD with a very high dose.
  #Real IdAgg can be made by putting yini["y1"] to 0 and let it grow like an untreated tumour
  yini["C"]<-500
  iat <- ode(y=yini,parms = parms, times = x, func = MonroGaffney)
  #naming our vectors
  iat_n <- iat[,2]+iat[,3]
  mtd_n <- mtd[,2]+mtd[,3]
  mtd_r <- mtd[,3]
  mtd_s <- mtd[,2]
  title<-"not logged"
  log<-""
  if(logged == TRUE){
    title<-"logged"
    log<-"y"  
  }
  plot(x,(iat_n),log=log,lty="dashed",type="l",col="black",main = paste(title, "IAT vs MTD"),xlim=xlim,ylim=ylim)
  lines(x,(mtd_n),lty="dotted",col="red")
  lines(x,(mtd_s),lty="dotted",col="blue")
  lines(x,(mtd_r),lty="dotted",col="green")
  abline(h=(parms["Ncrit"]))
  legend("bottomright",lty=c("dashed","dashed","dotted","dotted"), col=c("black","red","blue","green"),legend = c(sprintf("N_IAT;C=%s",500),sprintf("N_MTD;C=%s",cmax),sprintf("S_MTD;C=%s",cmax),sprintf("R_MTD;C=%s",cmax)),lwd=2);
}
#Comparison of time to prog/fail/surv for IdCont and IdAgg 
time_comp<-function(from,to,length,parms=pars_fig,logged=FALSE){
  #create a vector from 'from' to 'to' with the length 'length' 
  r0<-seq(from = from, to=to,length = length)
  t_prog<-NULL
  t_fail<-NULL
  t_surv<-NULL
  #stabilisation at Nacc for both tfail and tsurv
  for(r in r0){
    parms["r0"]<-r
    #tprog/tfail and tsurv all have two columns, the first one for the continuous containment, the second for the idagg treatment
    t_prog <-rbind(t_prog,c(time_prog_cont(parms["N0"],parms=parms),time_from_to(r,parms["N0"],parms=parms)))
    t_fail <-rbind(t_fail,c(time_fail_cont(parms["Nacc"],parms=parms),time_from_to(r,parms["Nacc"],parms=parms)))
    t_surv <-rbind(t_surv,c(time_surv_cont(parms["Nacc"],parms=parms),time_from_to(r,parms["Ncrit"],parms=parms)))
  }
  title<-"not logged"
  log<-""
  if(logged){
    title<-"logged"
    log<-"x"  
  }
  plot(x=(r0),log=log,y=t_prog[,1],type="l",main=paste(title,"Idealized Containment vs IAT"),col="red")
  lines((r0),y=t_prog[,2],type="l",col="blue")
  legend("topright",lty=c("solid","solid"), col=c("red","blue"),legend = c("Continuous","IAT"),lwd=2);
  
  plot(x=(r0),log=log,y=t_fail[,1],type="l",main=paste(title,"Idealized Containment vs IAT"),col="red")
  lines((r0),y=t_fail[,2],type="l",col="blue")
  legend("topright",lty=c("solid","solid"), col=c("red","blue"),legend = c("Continuous","IAT"),lwd=2);
  
  plot(x=(r0),log=log,y=t_surv[,1],type="l",main=paste(title,"Idealized Containment vs IAT"),col="red")
  lines((r0),y=t_surv[,2],type="l",col="blue")
  legend("topright",lty=c("solid","solid"), col=c("red","blue"),legend = c("Continuous","IAT"),lwd=2);
  
}
#Formula for containment at Cmax
formula_cont<-function(r,N0,cmax,nref,parms=pars_fig){
  #Formula from 1.4 in useful computation
  #growth no Nref
  a <- log(log(parms["K"]/N0)/log(parms["K"]/nref))
  if(cmax<=1){
    #repeat the value x
    return(rep(x=a/parms["rhoG"],times=length(r)))
  }
  #stabilization at Nref
  b<-1/log(parms["K"]/nref)
  b<-b*(log(N0/r)-log(cmax/(cmax-1)))
  b[which(b<0)]<-0
  return((a+b)/parms["rhoG"])
}
#Comparison of tsurv,tfail and tprog for containment at Nacc and MTD for C = Cmax for various r0
fail_surv<-function(from,to,length,cmax,x,parms=pars_fig,logged=FALSE){
  ###TSURV AND TFAIL FOR
  #CONTAINMENT at Nacc AND MTD AT CMAX
  r0 <- seq(from=from,to=to,length=length)
  parms["Cmax"]<-cmax
  mtd<-NULL
  contain<-NULL
  for(r in r0){
    parms["r0"]<-r
    #init_yini uses the parameters of parms to deduce the sensible portion of the tumour from N0 and r0
    yini<-init_yini(parms)
    untreat<-ode(y=yini,parms = parms, times = x, func = MonroGaffney)
    yini["C"]<-cmax
    #mtd
    m<-ode(y=yini,parms=parms,times=x,func=MonroGaffney)
    #containment
    c<-continuous(parms["Nacc"],x,parms,MonroGaffneyCont,untreat,method="adams")
    
    m<-m[,2]+m[,3]
    c<-c[,2]+c[,3]
    #x is the time vector, dic_find find the index just before we cross the target (parms["N0] or parms["nacc"])
    mtd<-rbind(mtd,c(x[dic_find(parms["N0"],m,0,length(m))],x[dic_find(parms["Nacc"],m,0,length(m))],x[dic_find(parms["Ncrit"],m,0,length(m))]))
    contain<-rbind(contain,x[dic_find(parms["Ncrit"],c,0,length(m))])
  } 
  log<-""
  title<-"not logged"
  ####using the formula ?
  if(logged==TRUE){
    title<-"logged"
    log<-"x"
  }
  #mtd[,1] is the time to prog, mtd[,2] to fail, mtd[,3] to surv
  #contain is time to surv since we have a formula for containment
  
  #prog
  a<-formula_cont(r0,parms["N0"],cmax,parms["N0"])
  b<-mtd[,1]
  ylim<-c(min(a,b),max(a,b))
  plot(x=r0,b,type="l",col="red",main=paste(title,"prog"),log=log,ylim=ylim)
  lines(x=r0,a,type="l",col="green")
  legend("bottomleft",lty=c("solid","solid"), col=c("red","green"),legend = c(sprintf("MTD;C=%s",cmax),sprintf("cont;C=%s",cmax)),lwd=2);
  
  #fail
  a<-formula_cont(r0,parms["N0"],cmax,parms["Nacc"])
  b<-mtd[,2]
  ylim<-c(min(a,b),max(a,b))
  plot(x=r0,log=log,b,type="l",col="red",main=paste(title,"fail"),ylim=ylim)
  lines(x=r0,a,type="l",col="green")
  legend("bottomleft",lty=c("solid","solid"), col=c("red","green"),legend = c(sprintf("MTD;C=%s",cmax),sprintf("cont;C=%s",cmax)),lwd=2);
  
  #surv
  a<-contain[,1]
  b<-mtd[,3]
  ylim<-c(min(a,b),max(a,b))
  plot(x=r0,b,type="l",col="red",main=paste(title,"surv"),log=log,ylim=ylim)
  lines(x=(r0),a,type="l",col="green")
  legend("bottomleft",lty=c("solid","solid"), col=c("red","green"),legend = c(sprintf("MTD;C=%s",cmax),sprintf("cont;C=%s",cmax)),lwd=2);
  
}
# #Tsurv, tfail and tprgog for containment, idCont, MTD, idAgg at Nacc for cmax
# #fusion of fig e and fig f
# fig_ef<-function(from,to,length,cmax,x,parms=pars_fig,logged=FALSE){
#   r0 <- seq(from=from,to=to,length=length)
#   parms["Cmax"]<-cmax
#   mtd<-NULL
#   idagg<-NULL
#   idcontain<-NULL
#   contain<-NULL
#   title<-"not logged"
#   log<-""
#   if(logged==TRUE){
#     title<-"logged"
#     log<-"x"
#   }
#   for(r in r0){
#     parms["r0"]<-r
#     yini<-init_yini(parms)
#     untreat<-ode(y=yini,parms = parms, times = x, func = MonroGaffney)
#     yini["C"]<-cmax
#     m<-ode(y=yini,parms=parms,times=x,func=MonroGaffney)
#     c<-continuous(parms["Nacc"],x,parms,MonroGaffneyCont,untreat,method="adams")
#     m<-m[,2]+m[,3]
#     c<-c[,2]+c[,3]
#     idcontain <-rbind(idcontain,c(time_prog_cont(parms["N0"],parms=parms),time_fail_cont(parms["Nacc"],parms=parms),time_surv_cont(parms["Nacc"],parms=parms)))
#     idagg <-rbind(idagg,c(time_from_to(r,parms["N0"],parms=parms),time_from_to(r,parms["Nacc"],parms=parms),time_from_to(r,parms["Ncrit"],parms=parms)))
#     mtd<-rbind(mtd,c(x[dic_find(parms["N0"],m,0,length(m))],x[dic_find(parms["Nacc"],m,0,length(m))],x[dic_find(parms["Ncrit"],m,0,length(m))]))
#     contain<-rbind(contain,c(x[dic_find(parms["N0"],c,0,length(m))],x[dic_find(parms["Nacc"],c,0,length(m))],x[dic_find(parms["Ncrit"],c,0,length(m))]))
#   } 
#   a<-formula_cont(r0,parms["N0"],cmax,parms["N0"])
#   b<-mtd[,1]
#   c<-idagg[,1]
#   d<-idcontain[,1]
#   ylim<-c(min(a,b,c,d),max(a,b,c,d))
#   plot(x=(r0),b,log=log,type="l",col="red",main=paste(title,"prog"),ylim=ylim)
#   lines(x=(r0),a,type="l",col="green")
#   lines(x=(r0),c,type="l",lty="dashed",col="blue")
#   lines(x=(r0),d,type="l",lty="dashed",col="magenta")
#   legend("bottomleft",lty=c("solid","solid","dashed","dashed"), col=c("red","green","blue","magenta"),legend = c(sprintf("MTD;C=%s",cmax),sprintf("cont;C=%s",cmax),"IdAgg","IdContain"),lwd=2);
#   
#   
#   a<-formula_cont(r0,parms["N0"],cmax,parms["Nacc"])
#   b<-mtd[,2]
#   c<-idagg[,2]
#   d<-idcontain[,2]
#   ylim<-c(min(a,b,c,d),max(a,b,c,d))
#   plot(x=log10(r0),b,type="l",col="red",main=paste(title,"fail"),ylim=ylim)
#   lines(x=log10(r0),a,type="l",col="green")
#   lines(x=log10(r0),c,type="l",lty="dashed",col="blue")
#   lines(x=log10(r0),d,type="l",lty="dashed",col="magenta")
#   legend("bottomleft",lty=c("solid","solid","dashed","dashed"), col=c("red","green","blue","magenta"),legend = c(sprintf("MTD;C=%s",cmax),sprintf("cont;C=%s",cmax),"IdAgg","IdContain"),lwd=2);
#   
#   
#   a<-contain[,3]
#   b<-mtd[,3]
#   c<-idagg[,3]
#   d<-idcontain[,3]
#   ylim<-c(min(a,b,c,d),max(a,b,c,d))
#   plot(x=log10(r0),b,type="l",col="red",main=paste(title,"surv"),ylim=ylim)
#   lines(x=log10(r0),a,type="l",col="green")
#   lines(x=log10(r0),c,type="l",lty="dashed",col="blue")
#   lines(x=log10(r0),d,type="l",lty="dashed",col="magenta")
#   legend("bottomleft",lty=c("solid","solid","dashed","dashed"), col=c("red","green","blue","magenta"),legend = c(sprintf("MTD;C=%s",cmax),sprintf("cont;C=%s",cmax),"IdAgg","IdContain"),lwd=2);
#   
# }

# times based on a mix of simulations and formulas:
fig_ef<-function(from,to,length,cmax,x,parms=pars_fig,logged=FALSE, objective = "N0", x_cex = 0.7){
  r0 <- seq(from=from,to=to,length=length)
  parms["Cmax"]<-cmax
  mtd<-NULL
  idagg<-NULL
  idcontain<-NULL
  contain<-NULL
  title<-"not logged"
  log<-""
  if(logged==TRUE){
    title<-"logged"
    log<-"x"
  }
  for(r in r0){
    parms["r0"]<-r
    yini<-init_yini(parms)
    untreat<-ode(y=yini,parms = parms, times = x, func = MonroGaffney)
    yini["C"]<-cmax
    m<-ode(y=yini,parms=parms,times=x,func=MonroGaffney)
    c<-continuous(parms["Ncrit"],x,parms,MonroGaffneyCont,untreat,method="adams")
    m<-m[,2]+m[,3]
    c<-c[,2]+c[,3]
    idcontain <-rbind(idcontain, 
                      c(
                        time_prog_cont(parms["N0"],parms=parms), 
                        time_fail_cont(parms["Nacc"], parms=parms), 
                        time_surv_cont(parms["Ncrit"],parms=parms)
                        )
                      )
    idagg <-rbind(idagg, 
                  c(
                    time_from_to(r,parms["N0"],parms=parms), 
                    time_from_to(r,parms["Nacc"],parms=parms), 
                    time_from_to(r,parms["Ncrit"],parms=parms)
                    )
                  )
    mtd<-rbind(mtd, 
               c(x[dic_find(parms["N0"],m,0,length(m))], 
                 x[dic_find(parms["Nacc"],m,0,length(m))], 
                 x[dic_find(parms["Ncrit"],m,0,length(m))] 
                 )
               )
    contain<-rbind(contain, 
                   c(x[dic_find(parms["N0"],c,0,length(m))], 
                     x[dic_find(parms["Nacc"],c,0,length(m))], 
                     x[dic_find(parms["Ncrit"],c,0,length(m))]
                     )
                   )
  } 
  if(objective == "N0") {
    a<-formula_cont(r0,parms["N0"],cmax,parms["N0"])
    ind <- 1
    yy <- expression(paste("time to progression (", italic("t"[prog]), ")"))
    objective_txt1 <- expression(paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = 2)"))
    objective_txt2 <- expression(paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = ", infinity, ")"))
  }
  else if(objective == "Nacc") {
    a<-formula_cont(r0,parms["N0"],cmax,parms["Nacc"])
    ind <- 2
    yy <- expression(paste("time to treat. failure (", italic("t"[fail]), ")"))
    objective_txt1 <- expression(paste("contain at ", italic("N"[tol]), " (", italic("C"[max]), " = 2)"))
    objective_txt2 <- expression(paste("contain at ", italic("N"[tol]), " (", italic("C"[max]), " = ", infinity, ")"))
  }
  else if(objective == "Ncrit") {
    a<-contain[,3]
    ind <- 3
    yy <- expression(paste("survival time (", italic("t"[surv]), ")"))
    objective_txt1 <- expression(paste("contain at ", italic("N"[crit]), " (", italic("C"[max]), " = 2)"))
    objective_txt2 <- expression(paste("contain at ", italic("N"[crit]), " (", italic("C"[max]), " = ", infinity, ")"))
  } else stop("Unknown objective")
  b<-mtd[,ind]
  c<-idagg[,ind]
  d<-idcontain[,ind]
  ylim<-c(0, 1.25 * max(a,b,c,d))
  plot(x=(r0/parms["N0"]),y=b,
       type="l",col="grey75",
       ylim=ylim, log = "x", lwd=2, # mtd (c=2)
       xaxt = "n", 
       yaxt = "n", 
       xlab = "", 
       ylab = yy)
  axis(1, 10^(-6:0), labels = parse(text=c("10^-6", "10^-5", "10^-4", "10^-3", "10^-2", "10^-1", "1")))
  if(max(ylim) > 800) axis(2, 200*(0:9), las = 2)
  else if(max(ylim) > 1000) axis(2, 400*(0:9), las = 2)
  else axis(2, 100*(0:9), las = 2)
  lines(x=(r0/parms["N0"]),a,type="l",col="black", lwd=2) # cont (c=2)
  lines(x=(r0/parms["N0"]),c,type="l",lty="dashed",col="grey75", lwd=2) # MTD (ideal)
  lines(x=(r0/parms["N0"]),d,type="l",lty="dashed",col="black", lwd=2) # contain (ideal)
  legend("topright", bty = "n",
         lty=c("solid","dashed"), 
         col=c("black","black"),
         legend = c(
           objective_txt1,
           objective_txt2
           ),
         lwd=2)
  legend("bottomleft", bty = "n",
         lty=c("solid","dashed"), 
         col=c("grey75","grey75"),
         legend = c(
           expression(paste("MTD (", italic("C"[max]), " = 2)")),
           expression(paste("MTD (", italic("C"[max]), " = ", infinity, ")"))
         ),
         lwd=2)
  mtext(expression(paste("initial frequency of resistance (", italic("R") [0], " / ", italic("N") [0], ")")), 1, 2, cex = x_cex)
}
# alternative version:
fig_ef_alternative<-function(from,to,length,cmax,x,parms=pars_fig,logged=FALSE, objective = "N0", x_cex = 0.7, level = NA){
  r0 <- seq(from=from,to=to,length=length)/parms["N0"]
  parms["Nref"] <- parms[objective]
  parms["Cmax"] <- cmax
  if(is.na(level)) level <- objective
  
  if(level == "N0" && objective == "N0") {
    a <- sapply(r0, heat_prog_cont, N0 = parms["N0"], parms=parms)
    b <- sapply(r0, heat_prog_agg, N0 = parms["N0"], parms=parms)
    c <- sapply(r0, heat_prog_id_agg, N0 = parms["N0"], parms=parms)
    d <- sapply(r0, heat_prog_id_cont, N0 = parms["N0"], parms=parms)
    ind <- 1
    yy <- expression(paste("time to progression (", italic("t"[prog]), ")"))
    objective_txt1 <- expression(paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = 2)"))
    objective_txt2 <- expression(paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = ", infinity, ")"))
  }
  else if(level == "Nacc" && objective == "Nacc") {
    a <- sapply(r0, heat_fail_cont, N0 = parms["N0"], parms=parms)
    b <- sapply(r0, heat_fail_agg, N0 = parms["N0"], parms=parms)
    c <- sapply(r0, heat_fail_id_agg, N0 = parms["N0"], parms=parms)
    d <- sapply(r0, heat_fail_id_cont, N0 = parms["N0"], parms=parms)
    ind <- 2
    yy <- expression(paste("time to treat. failure (", italic("t"[fail]), ")"))
    objective_txt1 <- expression(paste("contain at ", italic("N"[tol]), " (", italic("C"[max]), " = 2)"))
    objective_txt2 <- expression(paste("contain at ", italic("N"[tol]), " (", italic("C"[max]), " = ", infinity, ")"))
  }
  else if(level == "N0" && objective == "Nacc") {
    a <- sapply(r0, heat_fail_cont_extended_model, N0 = parms["N0"], parms=parms, level = parms["N0"])
    b <- sapply(r0, heat_fail_agg, N0 = parms["N0"], parms=parms)
    c <- sapply(r0, heat_fail_id_agg, N0 = parms["N0"], parms=parms)
    d <- sapply(r0, heat_fail_id_cont_N0, N0 = parms["N0"], parms=parms)
    ind <- 2
    yy <- expression(paste("time to treat. failure (", italic("t"[fail]), ")"))
    objective_txt1 <- expression(paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = 2)"))
    objective_txt2 <- expression(paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = ", infinity, ")"))
  }
  else if(level == "Ncrit" && objective == "Ncrit") {
    parms_a <- parms
    parms_a["Nacc"] <- parms_a["Ncrit"]
    a <- sapply(r0, heat_fail_cont, N0 = parms_a["N0"], parms=parms_a)
    b <- sapply(r0, heat_fail_agg, N0 = parms_a["N0"], parms=parms_a)
    c <- sapply(r0, heat_fail_id_agg, N0 = parms_a["N0"], parms=parms_a)
    d <- sapply(r0, heat_fail_id_cont, N0 = parms_a["N0"], parms=parms_a)
    ind <- 3
    yy <- expression(paste("survival time (", italic("t"[surv]), ")"))
    objective_txt1 <- expression(paste("contain at ", italic("N"[crit]), " (", italic("C"[max]), " = 2)"))
    objective_txt2 <- expression(paste("contain at ", italic("N"[crit]), " (", italic("C"[max]), " = ", infinity, ")"))
  }
  else if(level == "Nacc" && objective == "Ncrit") {
    parms_a <- parms
    parms_a["Nacc"] <- parms_a["Ncrit"]
    a <- sapply(r0, heat_fail_cont_extended_model, N0 = parms_a["N0"], parms=parms_a, level = parms["Nacc"])
    b <- sapply(r0, heat_fail_agg, N0 = parms_a["N0"], parms=parms_a)
    c <- sapply(r0, heat_fail_id_agg, N0 = parms_a["N0"], parms=parms_a)
    d <- sapply(r0, heat_surv_id_cont_Nacc, N0 = parms["N0"], parms=parms)
    ind <- 3
    yy <- expression(paste("survival time (", italic("t"[surv]), ")"))
    objective_txt1 <- expression(paste("contain at ", italic("N"[tol]), " (", italic("C"[max]), " = 2)"))
    objective_txt2 <- expression(paste("contain at ", italic("N"[tol]), " (", italic("C"[max]), " = ", infinity, ")"))
  } 
  else if(level == "N0" && objective == "Ncrit") {
    parms_a <- parms
    parms_a["Nacc"] <- parms_a["Ncrit"]
    a <- sapply(r0, heat_fail_cont_extended_model, N0 = parms["N0"], parms=parms_a, level = parms["N0"])
    b <- sapply(r0, heat_fail_agg, N0 = parms["N0"], parms=parms_a)
    c <- sapply(r0, heat_fail_id_agg, N0 = parms["N0"], parms=parms_a)
    d <- sapply(r0, heat_fail_id_cont_N0, N0 = parms["N0"], parms=parms_a)
    ind <- 2
    yy <- expression(paste("survival time (", italic("t"[surv]), ")"))
    objective_txt1 <- expression(paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = 2)"))
    objective_txt2 <- expression(paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = ", infinity, ")"))
  } else stop("Unknown objective")
  
  ylim<-c(0, min(1.2 * max(a,b,c,d, na.rm = TRUE), 5000))
  plot(x=r0,y=b,
       type="l",col="grey75",
       ylim=ylim, log = "x", lwd=2, # mtd (c=2)
       xaxt = "n", 
       yaxt = "n", 
       xlab = "", 
       ylab = yy)
  axis(1, 10^(-6:0), labels = parse(text=c("10^-6", "10^-5", "10^-4", "10^-3", "10^-2", "10^-1", "1")))
  if(max(ylim) > 800) axis(2, 200*(0:9), las = 2)
  else if(max(ylim) > 1000) axis(2, 400*(0:9), las = 2)
  else axis(2, 100*(0:9), las = 2)
  lines(x=r0,a,type="l",col="black", lwd=2) # cont (c=2)
  lines(x=r0,c,type="l",lty="dashed",col="grey75", lwd=2) # MTD (ideal)
  lines(x=r0,d,type="l",lty="dashed",col="black", lwd=2) # contain (ideal)
  legend("topright", bty = "n",
         lty=c("solid","dashed"), 
         col=c("black","black"),
         legend = c(
           objective_txt1,
           objective_txt2
         ),
         lwd=2)
  legend("bottomleft", bty = "n",
         lty=c("solid","dashed"), 
         col=c("grey75","grey75"),
         legend = c(
           expression(paste("MTD (", italic("C"[max]), " = 2)")),
           expression(paste("MTD (", italic("C"[max]), " = ", infinity, ")"))
         ),
         lwd=2)
  mtext(expression(paste("initial frequency of resistance (", italic("R") [0], " / ", italic("N") [0], ")")), 1, 2, cex = x_cex)
}
#Cont vs MTD, tumour size
cont_vs_mtd<-function(r0,cmax,x,nref,parms=pars_fig,logged=FALSE,png=FALSE){
  parms["r0"]<-r0
  parms["Cmax"]<-cmax
  mtd<-NULL
  contain<-NULL
  yini<-init_yini(parms)
  #untreated part 
  untreat<-ode(y=yini,parms = parms, times = x, func = MonroGaffney)
  yini["C"]<-cmax
  #mtd
  m<-ode(y=yini,parms=parms,times=x,func=MonroGaffney)
  #containment
  c<-continuous(nref,x,parms,MonroGaffneyCont,untreat,method="adams")
  m<-m[,2]+m[,3]
  c<-c[,2]+c[,3]
  title<-"not logged"
  log<-""
  if(logged==TRUE){
    title<-"logged"
    log<-"y"
  }
  if(png==TRUE){
    text<-sprintf("fig_1_g%s_%s_%s.png",if(logged==TRUE) "_log" else "",parms["r0"],cmax)
    png(text,width=800,height=600)
  }
  plot(x,m,type="l",col="red",main=paste(title,sprintf("MTD vs Cont at %s with cmax =%s r0=%s",nref,cmax,r0)),log= log,xlim=c(0,250),ylim=c(min(m,c),parms["K"]))
  lines(x,c,type="l",col="blue")
  abline(h=parms["N0"],lty="dashed")
  abline(h=parms["Nacc"],lty="dashed")
  abline(h=parms["Ncrit"],lty="dashed")
  legend("bottomright",lty=c("solid","solid"), col=c("red","blue"),legend = c(sprintf("MTD;C=%s",cmax),sprintf("cont;C=%s",cmax)),lwd=2);
  if(png==TRUE){
    dev.off()
  }
}