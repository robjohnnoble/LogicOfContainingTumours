#treatment of tumour at constant dose = c1,c2,c3
constant_treat<-function(parms=pars_fig,logged=FALSE,c1=0.5,c2=1,c3=2,func=MonroGaffney, c2extratext = "", x_cex = 0.7, delayed = FALSE){
  x=seq(from=0,to=2200,length=10000)
  yini<-init_yini(parms)
  yini_untreat<-init_yini(parms)
  parms["Cmax"] <- 500
  yini["C"]<-0
  
  untreat<-ode(y=yini_untreat,parms = parms, times = x, func = MonroGaffney)
  yini_untreat["C"]<-500
  parms_alt <- parms
  parms_alt["Cmax"] <- 2
  
  # containment:
  if(delayed) {
    N0<-continuous(parms["Nacc"],x,parms_alt,MonroGaffneyCont,untreat)
  } else {
    N0<-continuous(parms["N0"],x,parms_alt,MonroGaffneyCont,untreat)
  }
  N0_a<-( N0[,2]+ N0[,3])
  
  #no treatment
  untreated<-ode(y=yini,parms = parms, times = x, func = func)
  
  #IdAGg
  yini["C"]<-500
  iat <- ode(y=yini,parms = parms, times = x, func = func)
  
  #constant dose
  yini["C"]<-c1
  if(delayed) {
    const1 <- delayed_treat(c1, parms["Nacc"], times = x, parms = parms, func = func)
  } else {
    const1 <- ode(y=yini,parms = parms, times = x, func = func)
  }
  
  yini["C"]<-c2
  if(delayed) {
    const2 <- delayed_treat(c2, parms["Nacc"], times = x, parms = parms, func = func)
  } else {
    const2 <- ode(y=yini,parms = parms, times = x, func = func)
  }
  
  yini["C"]<-c3
  if(delayed) {
    const3 <- delayed_treat(c3, parms["Nacc"], times = x, parms = parms, func = func)
  } else {
    const3 <- ode(y=yini,parms = parms, times = x, func = func)
  }
  
  title<-"not logged"
  log<-""
  if(logged == TRUE){
    title<-"logged"
    log<-"y"
  }
  #taking the tumour size for each treatment
  untreated<-untreated[,2]+untreated[,3]
  iat<-iat[,2]+iat[,3]
  const1<-const1[,2]+const1[,3]
  const2<-const2[,2]+const2[,3]
  const3<-const3[,2]+const3[,3]
  
  #par(mar = c(3, 4, 1, 1.5))
  plot(x=x,y=untreated, type = "l",
       log = "y", 
       xlim = c(0, 1500), 
       ylim=c(parms["N0"]/2, 2e12),
       lwd = 2, 
       xaxt = "n", 
       yaxt = "n",
       xlab = "",
       ylab = "number of tumour cells",
       col = "black")
  axis(2, pars_fig["Ncrit"], labels = expression(paste(italic("N" [crit]))), las = 2, col = "red", col.axis = "red")
  axis(2, c(10^10, 10^11, 10^12), labels = parse(text=c("10^10", "10^11", "10^12")), las = 2)
  axis(1, 500*(0:4))
  mtext("time (days)", 1, 2, cex = x_cex)
  lines(x=x[-1],y=iat[-1],type="l",col="black",lty=2,lwd = 2)
  lines(x=x, y = N0_a,lty="solid",col="cyan3",lwd = 2)
  lines(x=x,y=const3,lty="solid",col="yellowgreen",lwd = 2)
  lines(x=x,y=const2,lty="solid",col="cyan4",lwd = 2)
  lines(x=x,y=const1,lty="solid",col="blue",lwd = 2)
  abline(h = parms["Ncrit"], col = "red", lty = 3)
  substitute(expr = paste(italic("C"), " = ", c1))
  legend("bottomright",
         bty = "n",
         lty=c("solid","dashed","solid","solid","solid","solid"), 
         col=c("black","black","cyan3","blue","cyan4","yellowgreen"),
         legend = c("untreated",
                    expression("ideal MTD"),
                    expression(paste("contain (", italic("C"[max]), " = 2)")),
                    substitute(expr = paste(italic("C"), " = ", c1), env = base::list(c1 = c1)),
                    substitute(expr = paste(italic("C"), " = ", c2, c2extratext), env = base::list(c2 = c2, c2extratext = c2extratext)),
                    substitute(expr = paste(italic("C"), " = ", c3), env = base::list(c3 = c3))
         ),
         lwd=2);
  
}
#supplementary variant of constant_treat
suppl_constant_treat<-function(parms=pars_fig,logged=FALSE,Cmax,r0,title,func=MonroGaffney, 
                               c2extratext = "", x_cex = 0.7,whichlines,add_legend=FALSE, Cmax_in_legend = FALSE, add_dose_plot = FALSE, 
                               tau1 = 0, tau2 = 0, xlim = c(0, 1000), ylim = c(1e7, 2e12), fine_axes = FALSE){
  x=seq(from=0,to=1.1*max(xlim),length=1000)
  parms["r0"] <- r0
  yini<-init_yini(parms)
  yini_untreat<-init_yini(parms)
  
  cont_func <- MonroGaffneyCont
  
  par(mar = c(3, 4.5, 1.5, 1))
  if(fine_axes) par(mar = c(3, 6.5, 1.5, 1))
  
  #no treatment:
  parms["Cmax"] <- Cmax
  yini["C"]<-0
  untreat<-ode(y=yini_untreat,parms = parms, times = x, func = func)
  #ideal containment at N0:
  parms["Cmax"] <- 500
  yini["C"]<-500
  idcontN0<-continuous(parms["N0"],x,parms,func = cont_func,untreat)
  #containment at N0 with / without mutations:
  parms["Cmax"] <- Cmax
  yini["C"]<-Cmax
  contN0 <- continuous(parms["N0"],x,parms,func = cont_func,untreat)
  contN0WithMutations <- continuous(parms["N0"],x,parms,func = cont_func,untreat, tau1 = tau1, tau2 = tau2)
  #ideal containment at Nacc:
  parms["Cmax"] <- 500
  yini["C"]<-500
  idcontNacc<-continuous(parms["Nacc"],x,parms,func = cont_func,untreat)
  #containment at Nacc:
  parms["Cmax"] <- Cmax
  yini["C"]<-Cmax
  contNacc <- continuous(parms["Nacc"],x,parms,func = cont_func,untreat)
  #idMTD:
  parms["Cmax"] <- 500
  yini["C"]<-500
  idMTD <- ode(y=yini,parms = parms, times = x, func = func)
  #MTD with / without mutations:
  parms["Cmax"] <- Cmax
  yini["C"]<-Cmax
  MTDout <- ode(y=yini,parms = parms, times = x, func = func)
  MTDoutWithMutations <- ode(y=yini,parms = parms, times = x, func = func, tau1 = tau1, tau2 = tau2)
  
  log<-""
  if(logged == TRUE){
    log<-"y"
  }
  #taking the tumour size for each treatment
  idcontN0_N<-(idcontN0[,2]+idcontN0[,3])
  contN0_N<-contN0[,2]+contN0[,3]
  contN0_R<-contN0[,3]
  contN0WithMutations_N<-contN0WithMutations[,2]+contN0WithMutations[,3]
  idcontNacc_N<-(idcontNacc[,2]+idcontNacc[,3])
  contNacc_N<-contNacc[,2]+contNacc[,3]
  idMTD_N<-idMTD[,2]+idMTD[,3]
  MTD_S<-MTDout[,2]
  MTD_R<-MTDout[,3]
  MTD_N<-MTD_S+MTD_R
  MTDWithMutations_N <- MTDoutWithMutations[,2]+MTDoutWithMutations[,3]
  
  print(as.numeric(contN0WithMutations_N[1000] / contN0_N[1000]))
  print(as.numeric(MTDWithMutations_N[1000] / MTD_N[1000]))
  
  #par(mar = c(3, 4, 1, 1.5))
  plot(x=x,y=MTD_N,type="l",
       log = "y", 
       xlim = xlim, 
       ylim=ylim,
       lwd = 2, 
       xaxt = "n", 
       yaxt = "n",
       xlab = "",
       ylab = "",
       col = "red", 
       main = title)
  axis(2, pars_fig["Ncrit"], labels = expression(paste(italic("N" [crit]))), las = 2, col = "red", col.axis = "red")
  
  if(fine_axes) axis(2, 99:101 * 10^8, labels = c(
    expression(paste(9.9, " \u00D7 ", 10^9)),
    expression(paste(10^10)),
    expression(paste(1.01, " \u00D7 ", 10^10))
    ), las = 2)
  else axis(2, 10^(6:12), labels = parse(text=c("10^6", "10^7", "10^8", "10^9", "10^10", "10^11", "10^12")), las = 2)
  
  if(fine_axes) axis(1, 10*c(0:100))
  else axis(1, 500*(0:3))
  
  mtext("time (days)", 1, 2, cex = x_cex)
  
  if(fine_axes) mtext("number of tumour cells", 2, 5, cex = x_cex)
  else mtext("number of tumour cells", 2, 3, cex = x_cex)
  
  lty_vec <- c("solid")
  col_vec <- c("red")
  if("contN0WithMutations" %in% whichlines) leg_vec <- c("without mutation, MTD")
  else leg_vec <- c("MTD")
  if("MTDWithMutations" %in% whichlines) {
    lty_vec <- c("dashed", lty_vec)
    col_vec <- c("black", col_vec)
    leg_vec <- c(leg_vec, "with mutation, MTD")
    lines(x=x,y=MTDWithMutations_N,lty="dashed",col="black",lwd = 2)
  }
  if("idMTD" %in% whichlines) {
    lty_vec <- c("dashed", lty_vec)
    col_vec <- c("black", col_vec)
    leg_vec <- c("ideal MTD", leg_vec)
    lines(x=x[-1],y=idMTD_N[-1],lty="dashed",col="black",lwd = 2)
  }
  if("MTD_S" %in% whichlines) {
    lty_vec <- c(lty_vec, "dashed")
    col_vec <- c(col_vec, "grey50")
    leg_vec <- c(leg_vec, expression(paste(italic("S"), " (MTD)")))
    lines(x=x, y = MTD_S,col="grey50",lwd = 2,lty=2)
  }
  if("MTD_R" %in% whichlines) {
    lty_vec <- c(lty_vec, "dotdash")
    col_vec <- c(col_vec, "grey75")
    leg_vec <- c(leg_vec, expression(paste(italic("R"), " (MTD)")))
    lines(x=x, y = MTD_R,col="grey75",lwd = 2,lty=4)
  }
  if("idcontN0" %in% whichlines) {
    lty_vec <- c(lty_vec, "solid")
    col_vec <- c(col_vec, "cyan3")
    leg_vec <- c(leg_vec, expression(paste("ideal contain at ", italic("N")[0])))
    lines(x=x, y = idcontN0_N,lty="solid",col="cyan3",lwd = 2)
  }
  if("contN0" %in% whichlines) {
    lty_vec <- c(lty_vec, "solid")
    col_vec <- c(col_vec, "darkcyan")
    if(Cmax_in_legend) {
      leg_vec <- c(leg_vec, substitute(expr = paste("contain at ", italic("N")[0], " (", italic("C"[max]), " = ", Cmax, ")"), env = base::list(Cmax = Cmax)))
    } else if("contN0WithMutations" %in% whichlines) {
      leg_vec <- c(leg_vec, expression(paste("without mutation, contain at ", italic("N")[0])))
    } else {
      leg_vec <- c(leg_vec, expression(paste("contain at ", italic("N")[0])))
    }
    lines(x=x, y = contN0_N,lty="solid",col="darkcyan",lwd = 2)
  }
  if("contN0WithMutations" %in% whichlines) {
    lty_vec <- c(lty_vec, "dotdash")
    col_vec <- c(col_vec, "black")
    if(Cmax_in_legend) {
      leg_vec <- c(leg_vec, substitute(expr = paste("with mutation, contain at ", italic("N")[0], " (", italic("C"[max]), " = ", Cmax, ")"), env = base::list(Cmax = Cmax)))
    } else {
      leg_vec <- c(leg_vec, expression(paste("with mutation, contain at ", italic("N")[0])))
    }
    lines(x=x, y = contN0WithMutations_N,lty="dotdash",col="black",lwd = 2)
  }
  if("contN0_R" %in% whichlines) {
    lty_vec <- c(lty_vec, "dotted")
    col_vec <- c(col_vec, "grey75")
    leg_vec <- c(leg_vec, expression(paste(italic("R"), " (contain)")))
    lines(x=x, y = contN0_R,col="grey75",lwd = 2,lty=3)
  }
  if("idcontNacc" %in% whichlines) {
    lty_vec <- c(lty_vec, "solid")
    col_vec <- c(col_vec, "magenta")
    leg_vec <- c(leg_vec, expression(paste("ideal contain at ", italic("N"[tol]))))
    lines(x=x, y = idcontNacc_N,lty="solid",col="magenta",lwd = 2)
  }
  if("contNacc" %in% whichlines) {
    lty_vec <- c(lty_vec, "solid")
    col_vec <- c(col_vec, "darkmagenta")
    if(Cmax_in_legend) {
      leg_vec <- c(leg_vec, substitute(expr = paste("contain at ", italic("N"[tol]), " (", italic("C"[max]), " = ", Cmax, ")"), env = base::list(Cmax = Cmax)))
    } else {
      leg_vec <- c(leg_vec, expression(paste("contain at ", italic("N"[tol]))))
    }
    lines(x=x, y = contNacc_N,lty="solid",col="darkmagenta",lwd = 2)
  }
  
  abline(h = parms["Ncrit"], col = "red", lty = 3)
  substitute(expr = paste(italic("C"), " = ", c1))
  if(add_legend) {
    legend("bottomright",
         bty = "n",
         lty=lty_vec, 
         col=col_vec,
         legend = leg_vec,
         lwd=2)
    
    if(add_dose_plot) {
      par(mar = c(3, 4, 1, 1))
      plot(x[-1],idcontNacc[-1,4],type="l",
           xlim=c(0,1000), 
           ylim = c(0, 2.25),
           ylab = "dose (C)", 
           xlab = "", 
           xaxt = "n", 
           yaxt = "n",
           lwd = 2, 
           col = "magenta")
      axis(2, 0:2, las = 2)
      axis(1, c(0, 500, 1000))
      mtext("time (days)", 1, 2, cex = x_cex)
      lines(x[-1],contNacc[-1,4],type="l", lwd = 2, col = "darkmagenta")
      legend("bottomright",
             bty = "n",
             lty=c("blank", "solid", "solid"),
             col=c("white", "magenta", "darkmagenta"),
             legend = c(expression(""),
                        expression(paste("ideal contain at ", italic("N"[tol]))),
                        expression(paste("contain at ", italic("N"[tol])))
             ),
             lwd=2)
    }
  }
  
}
#time to progression, failure and survival for a constant dose C
time_constant<-function(C,times=times_MG, parms=pars_fig,func=MonroGaffney, delayed = FALSE){
  yini<-init_yini(parms)
  yini["C"]<-C
  times_vec <- seq(from = times[1],to=times[2],length=1e4)
  
  if(delayed) {
    df <- delayed_treat(C, parms["Nacc"], times = times_vec, parms = parms, func = func)
  } else {
    df <-ode(y = yini, times = times_vec, func = func, parms = parms)
  }
  
  df[,2]<-df[,2]+df[,3]
  n<- length(df[,2])
  #manually find tprog tfail and tsurv 
  tprog <- dic_find(parms["N0"],df[,2],0,n)
  tfail <- dic_find(parms["Nacc"],df[,2],0,n)
  tsurv <- dic_find(parms["Ncrit"],df[,2],0,n)
  #find tprog, tfail and tsurv for a dose C
  c(df[tprog,1],df[tfail,1],df[tsurv,1])
}
#time to prog, failure and surv for a constant dose varying from dosemin to dosemax
#and same but for two patients with different sensitivity and the average time to prog/fail/surv
time_dose<-function(dosemin=0,dosemax,step=0.1,parms=pars_fig,func=MonroGaffney, delayed = FALSE){
  time_const<-NULL
  time_patient1<-NULL
  time_patient2<-NULL
  parms["target_max"]<-parms["N0"]
  parms["target_min"]<-0
  c<- seq(from = dosemin, to = 2 * dosemax, by = step)
  
  time_const <- as.data.frame(t(sapply(c,time_constant, func=MonroGaffney, delayed = delayed)))
  colnames(time_const) <- c("prog", "fail", "surv")
  time_const$c <- c
  time_const$c_patient1 <- c*2
  time_const$c_patient2 <- c/2
  
  return(time_const)
}

plot_time_dose <- function(parms=pars_fig, time_const, outcome, delayed = FALSE, x_cex = 0.7) {
  time_const$mean <- NA
  len <- floor(dim(time_const)[1] / 2)
  time_const[1:len, "mean"] <- (time_const[1:len, outcome] + time_const[2*(1:len), outcome])/2
  
  if(outcome == "prog") {
    if(!delayed) {
      top <- time_prog_cont(parms["N0"])
    } else {
      top <- -100
    }
    ylab <- expression(paste("time to progression (", italic("t"[prog]), ")"))
    position <- "bottomright"
    } else if(outcome == "fail") {
      top <- time_fail_cont(parms["Nacc"])
      ylab <- expression(paste("time to treat. failure (", italic("t"[fail]), ")"))
      position <- "bottomright"
    } else if(outcome == "surv") {
      top <- -100 # time_surv_cont(parms["Ncrit"])
      ylab <- expression(paste("survival (", italic("t"[surv]), ")"))
      position <- "bottomright"
    } else stop("Invalid outcome type")
  
  #fig3.c
  ylim<-c(0, max(time_const[, outcome],top))
  #par(mar = c(3, 4.25, 1, 1.5))
  plot(time_const[,"c"], time_const[, outcome], 
       type="l", col="white", xlim = c(0, 1.2 * max(time_const[,"c"]) / 2), ylim=ylim,
       #xlab = expression(paste("dose (", italic("C"), ")")),
       xlab = "",
       ylab = ylab,
       main="", 
       xaxt = "n", 
       yaxt = "n", lwd=2)
  lines(time_const[,"c_patient1"], time_const[, "mean"], 
        type="l", col="gold", lwd=2)
  lines(time_const[,"c"], time_const[, outcome], 
        type="l", col="blue", lwd=2)
  lines(time_const[,"c_patient1"], time_const[, outcome], 
        type="l", col="red", lwd=2)
  # lines(time_const[,"c_patient2"], time_const[, outcome], 
  #       type="l", col="green", lwd=2)
  abline(h = top, lty = 2, lwd=2)
  axis(1, 0:10, labels = 0:10)
  if(max(ylim) > 1000) interval <- 400
  else if(max(ylim) > 600) interval <- 200
  else interval <- 100
  axis(2, interval*(0:9), las = 2)
  mtext(expression(paste("dose (", italic("C"), ")")), 1, 2, cex = x_cex)
  if(outcome == "prog" & !delayed) legend(position, bty = "n",
         lty=c("solid", "solid", "solid", "dashed"), 
         col=c("blue","red","gold","black"),
         legend = c(
           expression(paste(italic(lambda), " = 1")),
           expression(paste(italic(lambda), " = 0.5")),
           "mean",
           expression("ideal cont.")
         ),
         lwd=2)
  else  legend(position, bty = "n",
               lty=c("solid", "solid", "solid"), 
               col=c("blue","red","gold"),
               legend = c(
                 expression(paste(italic(lambda), " = 1")),
                 expression(paste(italic(lambda), " = 0.5")),
                 "mean"
               ),
               lwd=2)
}
