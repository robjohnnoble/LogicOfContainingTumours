library(deSolve)

vonB <- function(rho, K, N0, times) {
  vonB_ode <- function(t, state, parameters) {
    with(
      as.list(c(state, parameters)),
      {
        dN <- rho * (N^(-1/3) - K^(-1/3)) * N
        list(c(dN))
      })
  }
  res <- ode(y = c(N = N0), times = times, func = vonB_ode, parms = c(rho = rho, K = K))
  return(res)
}

Gompertz <- function(K, rho, N0, times) {
  K * exp(log(N0 / K) * exp(-rho * times))
}

logistic <- function(K, rho, N0, times) {
  K * N0*exp(rho * times) / (K + N0*(exp(rho * times) - 1))
}

powerlaw1 <- function(rho, gamma, N0, times) {
  powerlaw_ode <- function(t, state, parameters) {
    with(
      as.list(c(state, parameters)),
      {
        dN <- rho * N^(1-gamma)
        list(c(dN))
      })
  }
  res <- ode(y = c(N = N0), times = times, func = powerlaw_ode, parms = c(rho = rho, gamma = gamma))
  return(res)
}
powerlaw2 <- function(rho, gamma, N0, times) {
  (N0^gamma + gamma*rho*times)^(1/gamma)
}
exponential <- function(rho, N0, times) {
  N0 * exp(rho * times)
}

#########

plot_three_model_curves <- function(N0, Gompertz_pars, logistic_pars, vonB_pars, power_pars, exp_pars) {
  
  tv <- seq(0, 300, length = 100)
  
  # Gompertz:
  plot(Gompertz(Gompertz_pars["K"], Gompertz_pars["rhoG"], N0, tv) ~ tv, type = "l", 
       ylim = c(N0, 1e12), log = "y", 
       lwd=2, 
       xaxt = "n", xlab = "", 
       yaxt = "n", ylab = "number of tumor cells")
  
  # von B:
  res <- vonB(vonB_pars["rho"], vonB_pars["K"], as.numeric(N0), tv)
  lines(res[,"N"] ~ res[,"time"], col = "blue", lwd = 2)
  
  # logistic:
  lines(logistic(logistic_pars["K"], logistic_pars["rho"], N0, tv) ~ tv, col = "red", lwd = 2)
  
  # exponential:
  lines(exponential(exp_pars["rho"], N0, tv) ~ tv, col = "gold", lwd = 2)
  
  # superexponential:
  # res <- powerlaw1(power_pars["rho"], -1/3, as.numeric(N0), tv)
  # lines(res[,"N"] ~ res[,"time"], col = "grey", lwd = 2)
  lines(powerlaw2(power_pars["rho"], -1/3, N0, tv) ~ tv, col = "grey", lwd = 2)
  
  legend("bottomright", legend = c("von Bertalanffy", "Gompertzian", "logisitic", "exponential", "superexponential"), 
         col = c("blue", "black", "red", "gold", "grey"),
         lwd = 2, bty = "n")
  
  #abline(h = N0, lty = 3, col = "red", lwd = 2)
  abline(h = pars_fig["Nacc"], lty = 3, col = "red", lwd = 2)
  abline(h = pars_fig["Ncrit"], lty = 3, col = "red", lwd = 2)
  
  axis(2, pars_fig["Ncrit"], labels = expression(paste(italic("N" [crit]))), las = 2, col = "red", col.axis = "red")
  axis(2, pars_fig["Nacc"], labels = expression(paste(italic("N" [tol]))), las = 2, col = "red", col.axis = "red")
  #axis(2, N0, labels = expression(paste(italic("N")[0], " = 10"^10)), las = 2, col = "red", col.axis = "red")
  axis(2, c(10^10, 10^11, 10^12), labels = parse(text=c("10^10", "10^11", "10^12")), las = 2)
  axis(1, 100*(0:5))
  
  mtext("time (days)", 1, 2, cex = 1)
}

#### gains for three different models
ThreeModelsGains <- function(Gompertz_pars, logistic_pars, vonB_pars, power_pars, exp_pars) {
  pars_fig["Nref"]<-pars_fig["Nacc"]
  ratio <- 10^seq(-6, 0, length = 150)
  
  t1 <- sapply(ratio, heat_fail_id_cont, N0 = pars_fig["N0"], parms = pars_fig)
  t2 <- sapply(ratio, heat_fail_id_agg, N0 = pars_fig["N0"], parms = pars_fig)
  g1 <- t1/t2
  t1 <- sapply(ratio, heat_fail_logistic_id_cont, N0 = pars_fig["N0"], parms = logistic_pars)
  t2 <- sapply(ratio, heat_fail_logistic_id_agg, N0 = pars_fig["N0"], parms = logistic_pars)
  l1 <- t1/t2
  t1 <- sapply(ratio, heat_fail_vonB_id_cont, N0 = pars_fig["N0"], parms = vonB_pars)
  t2 <- sapply(ratio, heat_fail_vonB_id_agg, N0 = pars_fig["N0"], parms = vonB_pars)
  v1 <- t1/t2
  t1 <- sapply(ratio, heat_fail_power_id_cont, N0 = pars_fig["N0"], parms = power_pars, gamma = -1/3)
  t2 <- sapply(ratio, heat_fail_power_id_agg, N0 = pars_fig["N0"], parms = power_pars, gamma = -1/3)
  p1 <- t1/t2
  t1 <- sapply(ratio, heat_fail_exp_id_cont, N0 = pars_fig["N0"], parms = exp_pars)
  t2 <- sapply(ratio, heat_fail_exp_id_agg, N0 = pars_fig["N0"], parms = exp_pars)
  e1 <- t1/t2
  
  plot(g1 ~ ratio, type = "l", ylim = c(0, 5.5), log = "x", lwd = 2, 
       xlab = "",
       ylab = "relative benefit",
       xaxt = "n", yaxt = "n")
  lines(v1 ~ ratio, col = "blue", lwd = 2)
  lines(l1 ~ ratio, col = "red", lwd = 2)
  lines(e1 ~ ratio, col = "gold", lwd = 2)
  lines(p1 ~ ratio, col = "grey", lwd = 2)
  axis(1, 10^(-6:0), labels = parse(text=c("10^-6", "10^-5", "10^-4", "10^-3", "10^-2", "10^-1", "1")))
  axis(2, 0:9, labels = 0:9, las = 2)
  mtext(expression(paste("initial frequency of resistance (", italic("R") [0], " / ", italic("N") [0], ")")), 1, 2, cex = 1)
  
  legend("topright", legend = c("von Bertalanffy", "Gompertzian", "logisitic", "exponential", "superexponential"), 
         col = c("blue", "black", "red", "gold", "grey"),
         lwd = 2, bty = "n")
}

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

pdf("ThreeModelsCurves.pdf", width=3.5, height=4)
par(mfrow = c(1, 1))
plot_three_model_curves(pars_fig["N0"], Gompertz_pars, logistic_pars, vonB_pars, power_pars, exp_pars)
dev.off()

pdf("FiveModels.pdf", width=8, height=6)
par(mfrow = c(1, 2))
plot_three_model_curves(pars_fig["N0"], Gompertz_pars, logistic_pars, vonB_pars, power_pars, exp_pars)
mtext("a", adj=-0.3, line=-0.5, cex = 1.5)
ThreeModelsGains(Gompertz_pars, logistic_pars, vonB_pars, power_pars, exp_pars)
mtext("b", adj=-0.3, line=-0.5, cex = 1.5)
dev.off()

