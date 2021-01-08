library(plyr)
library(lattice) # for heatmaps
library(deSolve) # for solving ODEs
library(ggplot2) # for plotting
library(reshape2) # for reshaping dataframes
library(gridExtra) # for plotting
library(viridisLite)
library(dplyr)
library(rasterVis)

source("best_outcome.R")
source("constant_RN.R")
source("continuous_containment_RN.R")
source("delayed_RN.R")
source("figure_plotting_functions.R")
source("GrowthCurvesDifferentModels.R")
source("heatmap_new_RN.R")
source("inter_containment_RN.R")
source("mode_func_RN.R")

###############

pars_fig <- init_pars_fig()
yini_MG <- init_yini(pars_fig)
x<-seq(from=0,to=2200,length=1000)
nref<-seq(from=pars_fig["r0"]*0.8,to=3e12,length=10000)

Gompertz_pars <- pars_fig

logistic_pars <- pars_fig
logistic_pars["K"] <- 6.4e11
logistic_pars["rho"] <- 2.4e-2

vonB_pars <- pars_fig
vonB_pars["K"] <- 5e13
vonB_pars["rho"] <- 90

power_pars <- pars_fig
power_pars["rho"] <- 4.5e-6

exp_pars <- pars_fig
exp_pars["rho"] <- 0.0175

###############

pdf("Fig1ab.pdf", width=6.75, height=2)
continuous_treat_alt(c1=0.9,c2=1.1,c3=2,logged=TRUE, x_cex = 0.7)
dev.off()

pdf("Fig1ce.pdf", width=7.25, height=2)
par(mfrow = c(1, 3))
suppl_constant_treat(pars_fig, TRUE, 2, 2.3e5, 
                     whichlines = c("idMTD", "MTD_R", "contN0", "contN0_R"), title = "", add_legend=TRUE, x_cex = 0.7)
suppl_constant_treat(pars_fig, TRUE, 2, 2.3e5, title = "", 
                     whichlines = c("idMTD", "idcontNacc", "contNacc"), add_legend=TRUE, x_cex = 0.7, 
                     Cmax_in_legend = FALSE, add_dose_plot = TRUE)
dev.off()

pdf("Fig1fgh.pdf", width=7, height=2)
Figure1fgh(x,500,2,logged=TRUE)
dev.off()

png("Fig2a.png", width=850, height=730, res = 200)
Figure2heatmaps(pars_fig, heat_prog_id_cont, heat_prog_id_agg, gap = 0.2, min_z = 1, max_z = 2.8)
dev.off()

png("Fig2b.png", width=850, height=730, res = 200)
Figure2heatmaps(pars_fig, heat_fail_id_cont, heat_fail_id_agg, gap = 0.2, min_z = 1, max_z = 2.8)
dev.off()

pdf("Fig2c.pdf", width=3.5, height=4)
ThreeModelsGains(Gompertz_pars, logistic_pars, vonB_pars, NA, NA)
dev.off()

# takes a long time!
pars_fig_alt <- pars_fig
pars_fig_alt["Kr"] <- pars_fig_alt["K"]
pars_fig_alt["Ks"] <- pars_fig_alt["K"]
pars_fig_alt["alpha"] <- 1
pars_fig_alt["beta"] <- 1
pdf("Fig2defhi.pdf", width=8, height=5)
par(mfrow = c(2, 3))
par(mar = c(4,5,1,1))
layout(matrix(c(1,2,3, 4,5,6), 2, 3, byrow = TRUE))
fig_ef(1e-6*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,objective="N0", x_cex = 0.7)
mtext("d", adj=-0.3, line=-0.5)
fig_ef_alternative(1e-6*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,
                   objective="Nacc", x_cex = 0.7,parms = pars_fig_alt, level = "N0")
mtext("e", adj=-0.3, line=-0.5)
fig_ef_alternative(1e-6*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,
                   objective="Ncrit", x_cex = 0.7,parms = pars_fig_alt, level = "N0")
mtext("f", adj=-0.3, line=-0.5)
plot.new()
fig_ef(1e-6*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE, objective = "Nacc")
mtext("h", adj=-0.3, line=-0.5)
fig_ef_alternative(1e-6*pars_fig["N0"],pars_fig["N0"],length=25,2,x,logged=TRUE,
                   objective="Ncrit", x_cex = 0.7,parms = pars_fig_alt, level = "Nacc")
mtext("i", adj=-0.3, line=-0.5)
dev.off()

png("Fig2g.png", width=750, height=800, res = 200)
FigureFromSupplTable5(150)
dev.off()

pdf("Fig3ab.pdf", width=7, height=4)
par(mfrow = c(1, 2))
constant_treat(c1=0.67,c2=1,c3=1.5,func=MonroGaffney,logged=TRUE, x_cex = 1)
constant_treat(c1=0.67,c2=1,c3=1.5,func=MonroGaffney,logged=TRUE, x_cex = 1, delayed = TRUE)
dev.off()

# takes a very long time!
time_const <- time_dose(0,4,0.01,func=MonroGaffney)
time_const_delayed <- time_dose(0,4,0.01,func=MonroGaffney, delayed = TRUE)
png("Fig3cde.png", width=1500, height=800 * 0.6, res = 200)
par(mfrow = c(1, 3))
plot_time_dose(pars_fig, time_const, "prog")
plot_time_dose(pars_fig, time_const, "fail")
plot_time_dose(pars_fig, time_const_delayed, "fail", delayed = TRUE)
dev.off()

data_approx_long <- heatmap_data(pars_fig, approx = TRUE, length = 200, logZ = 2)
data_approx_long$Z <- pmin(data_approx_long$Z, 8)
data_approx_long$Z <- pmax(data_approx_long$Z, 0)
png("Fig4a.png", width=3*75*7 * 0.6, height=3*75*2*3 * 0.6, res = 200)
Figure4heatmap(data_approx_long)
dev.off()

png("Fig4b.png", width=1100, height=600, res = 200)
plot_best_outcome_orientation2(pars_fig)
dev.off()

pdf("FigS1.pdf", width=4, height=3.5)
time_prog_var_nmin()
dev.off()

pdf("FigS2.pdf", width=8, height=6)
par(mfrow = c(1, 2))
plot_three_model_curves(pars_fig["N0"], Gompertz_pars, logistic_pars, vonB_pars, power_pars, exp_pars)
mtext("a", adj=-0.3, line=-0.5, cex = 1.5)
ThreeModelsGains(Gompertz_pars, logistic_pars, vonB_pars, power_pars, exp_pars)
mtext("b", adj=-0.3, line=-0.5, cex = 1.5)
dev.off()

pdf("FigS3.pdf", width=7, height=3.5)
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
pdf("FigS4.pdf", width=7, height=4)
figure_int_versus_contain(df_list, N0)
dev.off()

cairo_pdf("FigS5.pdf", width=8, height=6)
suppl_dynamics_with_mutations()
dev.off()

pdf("FigS6.pdf", width=7, height=4.5)
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

pdf("FigS7.pdf", width=5.5, height=4)
constant_treat(c1=0.74,c2=1.09,c3=2,func=MonroGaffney,logged=TRUE, x_cex = 1)
dev.off()

# takes a very long time!
data_res3 <- heatmap_data(pars_fig, approx = FALSE, logZ = 2)
pdf("FigS8.pdf", width=5.5, height=4)
Figure4heatmap(data_res3)
dev.off()

