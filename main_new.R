library(plyr)
library("lattice")#for heatmap
library(deSolve) # for solving ODEs
library(ggplot2) # for plotting
library(reshape2) # for reshaping dataframes
library(gridExtra) # for plotting
library(viridisLite)
library(dplyr)
library(rasterVis)
source("constant_RN.R")
source("continuous_containment_RN.R")
source("delayed_RN.R")
source("heatmap_new_RN.R")
source("inter_containment_RN.R")
source("mode_func_RN.R")
pars_fig <- init_pars_fig()
yini_MG <- init_yini(pars_fig)
x<-seq(from=0,to=2200,length=1000)
nref<-seq(from=pars_fig["r0"]*0.8,to=3e12,length=10000)

Figure1ab <- function(x,cmax1,cmax2,parms = pars_fig,logged = FALSE){
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
  ###Fig 1.a
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
  mtext("a", adj=-0.3, line=-0.5)
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
  
  #### fig 1_b
  time_cont(nref,logged=TRUE, png = FALSE)
  mtext("b", adj=-0.3, line=-0.5)
  
  plot.new()
  mtext("c", adj=-0.3, line=-0.5)
}

#### fig 1_c
Figure1c <- function(parms, func1, func2, length = 150, gap = 0.1, min_z = NA, max_z = NA) {
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

Figure2 <- function(x,cmax1,cmax2,parms = pars_fig,logged = FALSE,c1=0.67,c2=1,c3=1.5) {
  par(mfrow = c(2, 4))
  layout(matrix(c(1,2,3,4, 5,6,6,7), 2, 4, byrow = TRUE))
  yini_untreat<-init_yini(parms)#initialisation of yini
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
  print("aha1")
  #max<-continuous(parms["Ncrit"]*10,x,parms,MonroGaffneyCont,untreat)
  print("aha2")
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
  #max_a<-(max[,2]+max[,3])
  title <- "not logged"
  log<-""
  
  if(logged == TRUE){
    title <- "logged"
    #axis we want to be log
    log<-"y"
  }
  par(mar = c(3, 4, 1, 1.5))
  
  ###Fig 2_a Nacc = 7e10
  plot(x,acc_cmax[,2],log=log,type="l",
       lty=2,
       col="grey50",
       xlim=c(0,2200),
       ylim=c(parms["N0"]/2,2e12),
       lwd=2, 
       ylab = "number of tumor cells", 
       xlab = "", 
       xaxt = "n", 
       yaxt = "n")
  mtext("time (days)", 1, 2, cex = 0.7)
  mtext("a", adj=-0.3, line=-0.5)
  abline(h = pars_fig["Ncrit"], col = "red", lty = 3)
  axis(2, pars_fig["Ncrit"], labels = expression(paste(italic("N" [crit]))), las = 2, col = "red", col.axis = "red")
  axis(2, 10^(10:12), labels = parse(text=c("10^10", "10^11", "10^12")), las = 2)
  axis(1, 500*(0:5))
  lines(x,acc_cmax_bis,lty=4,col="grey75", lwd=2)
  lines(x,acc_a,col="magenta", lwd=2)
  lines(x,acc_cmax_a,col="darkmagenta", lwd=2)
  # lines(x,(mtd_b),lty="solid",col="yellow")
  # lines(x,(zero_a),lty="solid",col="purple")
  # abline(h=(parms["Ncrit"]))
  legend("bottomright",
         col=c("darkmagenta","grey50","grey75","magenta"),
         lty = c(1, 2, 4, 1),
         lwd=2,
         bty = "n",
         legend = c(expression(paste(italic("N"), " (", italic("C"[max]), " = 2)")),
                    expression(paste(italic("S"), " (", italic("C"[max]), " = 2)")),
                    expression(paste(italic("R"), " (", italic("C"[max]), " = 2)")),
                    expression(paste(italic("N"), " (", italic("C"[max]), " = ", infinity, ")"))
         ))
  
  #### fig 2_b
  plot(x,(acc_cmax_4),type="l",
       xlim=c(0,1000), 
       ylim = c(0, 2.5),
       ylab = expression(paste("dose (", italic(C), ")")), 
       xlab = "", 
       xaxt = "n", 
       yaxt = "n",
       lwd = 2)
  axis(2, 0:2, las = 2)
  axis(1, 500*(0:2))
  mtext("time (days)", 1, 2, cex = 0.7)
  mtext("b", adj=-0.3, line=-0.5)
  lines(x,(acc_4),type="l", lty = 2, lwd = 2)
  legend("bottomright", lty = 1:2,
         legend = c(expression(paste(italic("C"[max]), " = 2")),
                    expression(paste(italic("C"[max]), " = ", infinity))),
         lwd=2, bty = "n")
  
  fig_ef(1e-4*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE, objective = "Nacc")
  mtext("c", adj=-0.3, line=-0.5)
  
  plot.new()
  mtext("d", adj=-0.3, line=-0.5)
  
  intermittent_treat_a(log=TRUE,cmax=2,times=c(0,2200))
  mtext("e", adj=-0.3, line=-0.5)
  
  constant_treat(c1=c1,c2=c2,c3=c3,func=MonroGaffney,logged=TRUE)
  mtext("f", adj=-0.125, line=-0.5)
  
  time_const <- time_dose(dosemin=0,dosemax=4,step=0.01,parms=pars_fig,func=MonroGaffney)
  plot_time_dose(pars_fig, time_const, "prog")
  mtext("g", adj=-0.3, line=-0.5)
}

# # alternative version of Fig 1c, using simulations instead of formulas:
Figure1c_simulations <- function(parms) {
  parms["Nref"]<-parms["N0"]
  length<-5
  pars_fig_exended <- parms
  pars_fig_exended["Ks"] <- pars_fig_exended["K"]
  pars_fig_exended["Kr"] <- pars_fig_exended["K"]
  pars_fig_exended["alpha"] <- 1
  pars_fig_exended["beta"] <- 1
  pars_fig_exended["Cmax"] <- 500
  
  par(mfrow = c(5, 5))

  data_res3<-heat_fun1_fun2(length,heat_fail_cont_extended_model,heat_fail_agg_extended_model,pars_fig_exended)
  xlab = expression(paste("R0 / N0"))
  ylab = expression(paste("N0"))
  max_z <- ceiling(max(data_res3$Z) / 0.1) * 0.1
  print(min(data_res3$Z))
  print(max(data_res3$Z))
  scales <- list(x=list(at = -6:-1,
                        labels = parse(text=c("10^-6", "10^-5", "10^-4", "10^-3", "10^-2", "10^-1"))),
                 y=list(at = 6:10,
                        labels = parse(text=c("10^6", "10^7", "10^8", "10^9", "10^10"))))
  par(mfrow = c(1, 1))
  plot_heatmap(data_res3,xlab,ylab,
               NULL, scales, logged=TRUE, max=max_z)
}

Figure3a_data <- function(parms, approx = FALSE, length = 10, logZ = FALSE) {
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
  
  print(min(data_res3$Z, na.rm = TRUE))
  print(max(data_res3$Z, na.rm = TRUE))
  
  return(data_res3)
}

Figure3a_plot <- function(data_res3, gap = 1, max_z = NA) {
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

logistic_gain <- function(N0, R0, K) {
  top <- log(N0/R0) / (1 - N0 / K)
  bot <- log(N0*(R0-K) / (R0*(N0-K)))
  return(top/ bot)
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

###########

png("FigST5_taller.png", width=750, height=800, res = 200)
FigureFromSupplTable5(150)
dev.off()

logistic_pars <- pars_fig
logistic_pars["K"] <- 6.4e11
logistic_pars["rho"] <- 2.4e-2

png("Fig_logistic_heatmap.png", width=850, height=730, res = 200)
Figure1c(logistic_pars, heat_fail_logistic_id_cont, heat_fail_logistic_id_agg, gap = 0.01)
dev.off()

vonB_pars <- pars_fig
vonB_pars["a"] <- 90
vonB_pars["K"] <- 5e13

png("Fig_vonB_heatmap.png", width=850, height=730, res = 200)
Figure1c(vonB_pars, heat_fail_vonB_id_cont, heat_fail_vonB_id_agg, gap = 0.5)
dev.off()

heat_prog_agg_extended_model(1e10, 1e-6, pars_fig_exended)

pdf("Fig1ab.pdf", width=7, height=2)
Figure1ab(x,500,2,logged=TRUE)
dev.off()

#pdf("Fig1c_blues.pdf", width=7 * 0.6, height=2*3 * 0.6)
png("Fig1c_blues.png", width=850, height=730, res = 200)
Figure1c(pars_fig, heat_fail_id_cont, heat_fail_id_agg, gap = 0.2, min_z = 1, max_z = 2.8)
dev.off()

png("Fig1d_blues.png", width=850, height=730, res = 200)
Figure1c(pars_fig, heat_prog_id_cont, heat_prog_id_agg, gap = 0.2, min_z = 1, max_z = 2.8)
dev.off()

# alternative version of Fig 1c, using simulations instead of formulas:
# Figure1c_simulations(pars_fig)

pdf("Fig2_wider.pdf", width=9.3, height=4)
Figure2(x,500,2,logged=TRUE)
dev.off()

pdf("Fig2_suppl_version_larger.pdf", width=8, height=5)
par(mfrow = c(2, 3))
par(mar = c(4,5,1,1))
layout(matrix(c(1,2,3, 4,5,5), 2, 3, byrow = TRUE))
fig_ef(1e-4*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,objective="N0", x_cex = 0.7)
mtext("a", adj=-0.3, line=-0.5)
fig_ef_alternative(1e-4*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,
                    objective="Nacc", x_cex = 0.7,parms = pars_fig_alt, level = "N0")
mtext("b", adj=-0.3, line=-0.5)
fig_ef_alternative(1e-4*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,
                   objective="Ncrit", x_cex = 0.7,parms = pars_fig_alt, level = "N0")
mtext("c", adj=-0.3, line=-0.5)
fig_ef_alternative(1e-4*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,
                   objective="Ncrit", x_cex = 0.7,parms = pars_fig_alt, level = "Nacc")
mtext("d", adj=-0.3, line=-0.5)
constant_treat(c1=0.74,c2=1.09,c3=2,func=MonroGaffney,logged=TRUE, x_cex = 0.7)
mtext("e", adj=-0.1, line=-0.5)
dev.off()

#data_res3 <- Figure3a_data(pars_fig, approx = FALSE, logZ = 2)
#data_approx <- Figure3a_data(pars_fig, approx = TRUE, logZ = 2)
data_approx_long <- Figure3a_data(pars_fig, approx = TRUE, length = 200, logZ = 2)
data_approx_long$Z <- pmin(data_approx_long$Z, 8)
data_approx_long$Z <- pmax(data_approx_long$Z, 0)

# Figure3a_plot(data_res3, 0.25)
# Figure3a_plot(data_approx, 0.25)

# plot(data_approx$Z ~ data_res3$Z)
# lines(1 * data_res3$Z ~ data_res3$Z, lty = 2, col = "red")

png("Fig3a_blues.png", width=3*75*7 * 0.6, height=3*75*2*3 * 0.6, res = 200)
Figure3a_plot(data_approx_long)
dev.off()

# takes a long time!
time_const <- time_dose(0,4,0.01,func=MonroGaffney)
time_const_delayed <- time_dose(0,4,0.01,func=MonroGaffney, delayed = TRUE)

png("Fig_two_patients_prog.png", width=810, height=3*75*2*3 * 0.6, res = 200)
plot_time_dose(pars_fig, time_const, "prog")
dev.off()
png("Fig_two_patients_fail.png", width=810, height=3*75*2*3 * 0.6, res = 200)
plot_time_dose(pars_fig, time_const, "fail")
dev.off()
png("Fig_two_patients_surv.png", width=810, height=3*75*2*3 * 0.6, res = 200)
plot_time_dose(pars_fig, time_const, "surv")
dev.off()

png("Fig_two_patients_delyaed_prog.png", width=810, height=3*75*2*3 * 0.6, res = 200)
plot_time_dose(pars_fig, time_const_delayed, "prog", delayed = TRUE)
dev.off()
png("Fig_two_patients_delyaed_fail.png", width=810, height=3*75*2*3 * 0.6, res = 200)
plot_time_dose(pars_fig, time_const_delayed, "fail", delayed = TRUE)
dev.off()
png("Fig_two_patients_delyaed_surv.png", width=810, height=3*75*2*3 * 0.6, res = 200)
plot_time_dose(pars_fig, time_const_delayed, "surv", delayed = TRUE)
dev.off()

png("Fig_two_patients_combined.png", width=1500, height=800 * 0.6, res = 200)
par(mfrow = c(1, 3))
plot_time_dose(pars_fig, time_const, "prog")
plot_time_dose(pars_fig, time_const, "fail")
plot_time_dose(pars_fig, time_const_delayed, "fail", delayed = TRUE)
dev.off()

################

png("BestOutcome.png", width=1100, height=600, res = 200)
plot_best_outcome_orientation2(pars_fig)
dev.off()

################

# supplementary dynamical figures:

pdf("suppl_dynamics1.pdf", width=7, height=4.5)
par(mfrow = c(2, 3))
add_legend <- TRUE
for(r0 in c(2.3e5, 1e8)) for(Cmax in c(0.8, 1.5, 5)) {
  par(mar = c(4,4,1.5,1))
  if(r0 == 2.3e5) title <- substitute(expr = paste(italic("R"[0]), " = ", "2.3 ", "\u00D7", " 10"^5, ", ", 
                                                   italic("C"[max]), " = ", Cmax), 
                                      env = base::list(Cmax = Cmax))
  if(r0 == 1e8) title <- substitute(expr = paste(italic("R"[0]), " = ", " 10"^8, ", ", 
                                                 italic("C"[max]), " = ", Cmax), 
                                    env = base::list(Cmax = Cmax))
  suppl_constant_treat(pars_fig, TRUE, Cmax, r0, title = title, 
                       whichlines = c("idMTD", "idcontNacc", "contNacc"), add_legend=add_legend)
  add_legend <- FALSE
}
dev.off()

pdf("suppl_dynamics2.pdf", width=3.5, height=3.5)
par(mfrow = c(1, 1))
par(mar = c(4,4,1.5,1))
suppl_constant_treat(pars_fig, TRUE, 2, 1e6, 
                     whichlines = c("idMTD", "MTD_S", "MTD_R"), title = "", add_legend=TRUE, x_cex = 1)
dev.off()

pdf("suppl_dynamics3.pdf", width=7, height=3.5)
par(mfrow = c(1, 2))
add_legend <- TRUE
for(r0 in c(2.3e5, 1e8)) {
  par(mar = c(4,4,1.5,1))
  if(r0 == 2.3e5) title <- substitute(expr = paste(italic("R"[0]), " = ", "2.3 ", "\u00D7", " 10"^5))
  if(r0 == 1e8) title <- substitute(expr = paste(italic("R"[0]), " = ", " 10"^8))
  suppl_constant_treat(pars_fig, TRUE, 2, r0, title = title, 
                       whichlines = c("contN0"), add_legend=add_legend, x_cex = 1)
  add_legend <- FALSE
}
dev.off()

# containment with three different values of Cmax:
pdf("containment_three_Cmax_values.pdf", width=6.75, height=2)
continuous_treat_alt(c1=0.9,c2=1.1,c3=2,logged=TRUE, x_cex = 0.7)
dev.off()

# for new Figure 1c, 1d, 1e:
pdf("const_Cmax2.pdf", width=7.25, height=2)
par(mfrow = c(1, 3))
suppl_constant_treat(pars_fig, TRUE, 2, 2.3e5, 
                     whichlines = c("idMTD", "MTD_R", "contN0", "contN0_R"), title = "", add_legend=TRUE, x_cex = 0.7)
suppl_constant_treat(pars_fig, TRUE, 2, 2.3e5, title = "", 
                     whichlines = c("idMTD", "idcontNacc", "contNacc"), add_legend=TRUE, x_cex = 0.7, 
                     Cmax_in_legend = FALSE, add_dose_plot = TRUE)
dev.off()

pars_fig_alt <- pars_fig
pars_fig_alt["Kr"] <- pars_fig_alt["K"]
pars_fig_alt["Ks"] <- pars_fig_alt["K"]
pars_fig_alt["alpha"] <- 1
pars_fig_alt["beta"] <- 1
pdf("New_Fig2_panels.pdf", width=8, height=5)
par(mfrow = c(2, 3))
par(mar = c(4,5,1,1))
layout(matrix(c(1,2,3, 4,5,6), 2, 3, byrow = TRUE))
fig_ef(1e-6*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,objective="N0", x_cex = 0.7)
mtext("a", adj=-0.3, line=-0.5)
fig_ef_alternative(1e-6*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,
                   objective="Nacc", x_cex = 0.7,parms = pars_fig_alt, level = "N0")
mtext("b", adj=-0.3, line=-0.5)
fig_ef_alternative(1e-6*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,
                   objective="Ncrit", x_cex = 0.7,parms = pars_fig_alt, level = "N0")
mtext("c", adj=-0.3, line=-0.5)
plot.new()
fig_ef(1e-6*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE, objective = "Nacc")
mtext("d", adj=-0.3, line=-0.5)
fig_ef_alternative(1e-6*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,
                   objective="Ncrit", x_cex = 0.7,parms = pars_fig_alt, level = "Nacc")
mtext("e", adj=-0.3, line=-0.5)
dev.off()

# new Figure 3:
pdf("New_Fig3_top_row.pdf", width=7, height=4)
par(mfrow = c(1, 2))
constant_treat(c1=0.67,c2=1,c3=1.5,func=MonroGaffney,logged=TRUE, x_cex = 1)
constant_treat(c1=0.67,c2=1,c3=1.5,func=MonroGaffney,logged=TRUE, x_cex = 1, delayed = TRUE)
dev.off()

# new Figure S3:
pdf("New_FigS3.pdf", width=5.5, height=4)
constant_treat(c1=0.74,c2=1.09,c3=2,func=MonroGaffney,logged=TRUE, x_cex = 1)
dev.off()

# time to progression for ideal intermittent therapy (added during revisions, Oct 2020):
pdf("intermittent_prog.pdf", width=4, height=3.5)
time_prog_var_nmin()
dev.off()

# with mutations (added during revisions, Oct 2020):
cairo_pdf("suppl_dynamics_with_mutations.pdf", width=8, height=6)
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
dev.off()

df_list <- time_intermittent(0.8*pars_fig["N0"], num_steps = 1e4, return_df = TRUE)
df_M <- df_list[["data"]]
x=seq(from=0,to=2200,length=1e5)
pars1 <- pars_fig
yini<-init_yini(pars1)
yini_untreat<-init_yini(pars1)
pars1["Cmax"] <- 2
yini["C"]<-0
untreat<-ode(y=yini_untreat,parms = pars1, times = x, func = MonroGaffney)
yini_untreat["C"]<-0
parms_alt <- pars1
parms_alt["Cmax"] <- 2
N0<-continuous(pars1["N0"],x,parms_alt,MonroGaffneyCont,untreat)
pdf("int_versus_contain.pdf", width=7, height=4)
figure_int_versus_contain(df_list, N0)
dev.off()

