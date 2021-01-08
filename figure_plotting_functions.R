Figure1fgh <- function(x,cmax1,cmax2,parms = pars_fig,logged = FALSE){
  par(mfrow = c(1, 3))
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
  max<-continuous(2*parms["K"],x,parms,MonroGaffneyCont,untreat)
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
  ###Fig 1.f
  #Ncrit
  par(mar = c(3, 4, 1, 1))
  
  plot(x,crit_a,type="l",col="orange",
       log=log, 
       ylim=c(parms["N0"]/2,2e12), 
       ylab = "number of tumor cells", 
       xlab = "", 
       lwd=2, 
       xaxt = "n", 
       yaxt = "n")
  mtext("time (days)", 1, 2, cex = 0.7)
  mtext("f", adj=-0.3, line=-0.5)
  abline(h = pars_fig["Ncrit"], col = "red", lty = 3)
  axis(2, pars_fig["Ncrit"], labels = expression(paste(italic("N" [crit]))), las = 2, col = "red", col.axis = "red")
  axis(2, c(10^10, 10^11, 10^12), labels = parse(text=c("10^10", "10^11", "10^12")), las = 2)
  axis(1, 500*(0:5))
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
                    expression(paste("id. cont. at ", italic("N") [0])),
                    expression(paste("id. cont. at ", italic("N" [tol]))),
                    expression(paste("id. cont. at ", italic("N" [crit])))
         ),
         lwd=2, 
         bty = "n")
  
  #### fig 1_g
  intermittent_treat_a(log=TRUE,cmax=2,times=c(0,2200))
  mtext("g", adj=-0.3, line=-0.5)
  
  #### fig 1_h
  time_cont(nref,logged=TRUE, png = FALSE)
  mtext("h", adj=-0.3, line=-0.5)
}

Figure2heatmaps <- function(parms, func1, func2, length = 150, gap = 0.1, min_z = NA, max_z = NA) {
  parms["Nref"]<-parms["Nacc"]
  
  data_res3<-heat_fun1_fun2(length, func1, func2, parms)
  print(head(data_res3))
  xlab = expression(paste("initial frequency of resistance (", italic("R") [0], " / ", italic("N") [0], ")"))
  ylab = expression(paste("initial number of tumor cells (", italic("N") [0], ")"))
  if(is.na(min_z)) min_z <- floor(min(data_res3$Z) / gap) * gap
  if(is.na(max_z)) max_z <- ceiling(max(data_res3$Z) / gap) * gap
  scales <- list(x=list(at = -6:-1,
                        labels = parse(text=c("10^-6", "10^-5", "10^-4", "10^-3", "10^-2", "10^-1"))), 
                 y=list(at = 6:10,
                        labels = parse(text=c("10^6", "10^7", "10^8", "10^9", "10^10"))))
  # botseq <- floor(min(data_res3$Z) / gap)
  # topseq <- ceiling(max(data_res3$Z) / gap)
  # print(c(botseq, topseq))
  # myseq <- gap*botseq:topseq
  par(mfrow = c(1, 1))
  print(min_z)
  print(max_z)
  plot_heatmap(data_res3,xlab,ylab, 
               NULL, scales=scales, logged=TRUE, min = min_z, max=max_z, contour_gap=gap)
}

FigureFromSupplTable5 <- function(length = 10) {
  afn <- function(pair) 1 + log(1-1/pair["X"]) / log(1/pair["Y"])
  ratio_vec <- 10^seq(-6, -1, length = length)
  lambda_Cmax_vec <- seq(1.5, 3, length = length)
  m <- expand.grid(X = lambda_Cmax_vec, Y = ratio_vec)
  m$Z <- unlist(afn(m))
  my_scales <- list(x=list(at = -6:-1,
                           labels = parse(text=c("10^-6", "10^-5", "10^-4", "10^-3", "10^-2", "10^-1"))), 
                    y=list(at = 0.5*(2:10), labels = 0.5*(2:10))
  )
  print(c(min(m$Z), max(m$Z)))
  xlab <- expression(paste("initial frequency of resistance (", italic("R") [0], " / ", italic("N") [0], ")"))
  ylab <- expression(paste("maximum dose (", italic("C") [max], ")"))
  plot_heatmap(m, xlab, ylab, "", min = 0.5, max = 1, 
               scales = my_scales, logged = "y", contour_gap = 0.05)
}

Figure4heatmap <- function(data_res3, gap = 1, max_z = NA) {
  xlab = expression(paste(italic("K" [s]), " / ", italic("K" [r])))
  ylab = expression(paste(italic(beta)))
  scales <- list(x=list(at = 2*0:5,
                        labels = 2*0:5), 
                 y=list(at = 0:5,
                        labels = 0:5))
  if(is.na(max_z)) {
    max_z <- ceiling(max(data_res3$Z, na.rm = TRUE) / gap) * gap
    print(max_z)
  }
  
  myseq <- seq(0, max_z, by = gap)
  mylabs <- as.character(2^myseq)
  mylabs[length(mylabs)] <- expression(infinity)
  
  par(mfrow = c(1,1))
  
  plot_heatmap(data_res3,xlab,ylab, 
               NULL, scales, logged=FALSE, min = 0, max=max_z, contour_gap = gap, myseq = myseq, mylabs = mylabs)
}

heatmap_data <- function(parms, approx = FALSE, length = 10, logZ = FALSE) {
  parms["Nref"]<-parms["N0"]
  pars_fig_exended <- parms
  pars_fig_exended["Ks"] <- pars_fig_exended["K"]
  pars_fig_exended["Kr"] <- pars_fig_exended["K"]
  pars_fig_exended["alpha"] <- 1
  pars_fig_exended["Cmax"] <- 500
  
  par(mfrow = c(5,5))
  
  if(approx) data_res3<-heat_fun1_fun2_extended(length,heat_fail_cont_extended_approx,heat_fail_agg_extended_approx,pars_fig_exended)
  else data_res3<-heat_fun1_fun2_extended(length,heat_fail_cont_extended_model,heat_fail_agg_extended_model,pars_fig_exended)
  
  if(logZ) data_res3$Z <- log10(data_res3$Z)/log10(logZ)
  
  # print(min(data_res3$Z, na.rm = TRUE))
  # print(max(data_res3$Z, na.rm = TRUE))
  
  return(data_res3)
}

figure_int_versus_contain <- function(df_list, N0) {
  par(mar = c(4, 4.5, 1, 1), las = 1)
  plot(value ~ time, data = filter(df_M, variable == "n"), log = "y", type = "l", 
       ylim = c(5e9, 1.2e10), xlim = c(0, 350), col = "black", 
       xlab = "time (days)", 
       ylab = expression(paste("number of tumor cells", " \u00D7 ", 10^9)), 
       yaxt = "n")
  lines(value ~ time, data = filter(df_M, variable == "y1"), col = "blue")
  lines(value ~ time, data = filter(df_M, variable == "y2"), col = "red")
  abline(v = df_list[["times"]][1], lty = 2)
  t_cont <- formula_cont(pars_fig["r0"], pars_fig["N0"], 2, pars_fig["N0"])
  abline(v = t_cont, lty = 2, col = "grey")
  lines(y1 ~ time, data = N0, col = "skyblue")
  lines(y2 ~ time, data = N0, col = "pink")
  lines(y1 + y2 ~ time, data = N0, col = "grey")
  legend("bottomleft", bty = "n", lty = rep(1, 6),
         col = c("black", "blue", "red", "grey", "skyblue", "pink"), 
         legend = c(expression(paste(italic("N"), " (intermittent containment)")), 
                    expression(paste(italic("S"), " (intermittent containment)")), 
                    expression(paste(italic("R"), " (intermittent containment)")), 
                    expression(paste(italic("N"), " (containment)")), 
                    expression(paste(italic("S"), " (containment)")), 
                    expression(paste(italic("R"), " (containment)"))
         ))
  axis(2, 5:12 * 1e9, 5:12)
}

suppl_dynamics_with_mutations <- function() {
  par(mfrow = c(1, 2))
  add_legend <- TRUE
  par(mfrow = c(2, 2))
  add_legend <- TRUE
  get_r0 <- function(tau1, tau2, n0) tau1 / (tau1 + tau2) * n0 * (1 - n0 ^ (-tau1 - tau2))
  tau_vec <- c(1e-6, 1e-4)
  for(row in 1:2) for(col in 1:2) {
    tau1 <- tau_vec[col]
    tau2 <- 0
    par(mar = c(4,6,1.5,1), las = 0)
    mut_pars <- pars_fig
    r0 <- get_r0(tau1, tau2, mut_pars["N0"])
    print(paste("R0 / N0 =", as.numeric(r0 / mut_pars["N0"])))
    if(row == 1) {
      xlim <- c(0, 1000)
      ylim <- c(1e7, 2e12)
      fine_axes <- FALSE
      adj <- -0.2
    }
    else if(col == 2) {
      xlim <- c(150, 180)
      ylim <- c(9.9e9, 1.01e10)
      fine_axes <- TRUE
      adj <- -0.35
    }
    else {
      xlim <- c(220, 330)
      ylim <- c(9.9e9, 1.01e10)
      fine_axes <- TRUE
      adj <- -0.35
    }
    if(col == 1) title <- substitute(expr = paste(italic("\u03C4"[1]), " = ", " 10"^-6))
    else title <- substitute(expr = paste(italic("\u03C4"[1]), " = ", " 10"^-4))
    suppl_constant_treat(pars_fig, TRUE, 2, r0, title = title, 
                         whichlines = c("contN0", "contN0WithMutations", "MTDWithMutations"), add_legend=add_legend, x_cex = 1, 
                         tau1 = tau1, tau2 = tau2, xlim=xlim, ylim=ylim, fine_axes = fine_axes)
    mtext(letters[(row - 1) * 2 + col], adj=adj, line=-0.5)
    add_legend <- FALSE
  }
}


