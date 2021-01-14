# solve a model for fixed parameter values and format output ready for plotting:
solve_odes <- function(state, times, func, parms) {
  parms["rhoG"] <- with(as.list(parms), rhoG(fr0, N0, K))
  parms["rhoLV"] <- with(as.list(parms), rhoLV(fr0, N0, K))
  parms["rhoPL"] <- with(as.list(parms), rhoPL(fr0, N0, gamma))
  # solve the equations:
  out <- ode(y = state, times = times, func = func, parms = parms)
  # convert output to tidy data frame:
  df_model <- melt(as.data.frame(out), id.vars = "time")
  return(df_model)
}

# Solve a model with discontinuous changes in parameter values, 
# assuming that the sum of the variables will alternately reach 
# upper and lower thresholds, starting with the upper threshold. 
# If the thresholds are close together then you might need to 
# increase the value of num_steps.
seq_solve <- function(state, times, func, parms, update_rule, num_steps = 400) {
  if("fr0" %in% names(parms)) {
    parms["rhoG"] <- with(as.list(parms), rhoG(fr0, N0, K))
    parms["rhoLV"] <- with(as.list(parms), rhoLV(fr0, N0, K))
    parms["rhoPL"] <- with(as.list(parms), rhoPL(fr0, N0, gamma))
  }
  # calculate deviation between model and target at time t:
  test_func <- function(t, t_start, parms, target) {
    out <- ode(y = state, times = c(t_start, t), func = func, parms = parms)
    return(sum(as.data.frame(out)[2, -1]) - target) # -1 is to exclude the "time" column
  }
  df <- data.frame() # dataframe that will store all output
  tmin <- min(times) # minimum time value
  tmax <- max(times) # maximum time value
  names_vec <- names(state)
  
  target <- as.numeric(parms["target_max"]) # assume that target_max will be reached before target_min
  
  max_loops <- 1e6 # just to ensure we can't get stuck in an infinite loop
  for(i in 1:max_loops) {
    if(!is.na(parms["target_max"])) {
      testmin <- test_func(tmin, tmin - 1e-6, parms, target) # deviation between model and target at time tmin
      
      # find the earliest time (between tmin and tmax) at which model reaches target
      # or, if there is no such time, then return tmax:
      for(t in seq(tmin, tmax, length = num_steps)) {
        testmax <- test_func(t, tmin - 1e-6, parms, target) # deviation between model and target at time t
        if(sign(testmin) != sign(testmax)) { # if model reaches target between times tmin and t
          soln <- uniroot(test_func, c(tmin - 1e-6, t), t_start = tmin, parms = parms, target = target)
          t <- min(soln$root, t) # precise time at which model reaches target
          break
        }
      }
    } else t <- max(times)
    
    # solve the model for the current time frame:
    out <- ode(y = state, seq(tmin, t, length = 1000), func = func, parms = parms)
    out <- as.data.frame(out)
    df_model <- out
    if(dim(df_model)[2] > 2) df_model$n <- rowSums(df_model[ , -1])
    else df_model$n <- df_model[ , 2]
    df_model <- melt(df_model, id.vars = "time")
    df <- rbind(df, df_model) # apend new results to cumulative results
    
    # break if model has been solved for the entire time period:
    if(t >= max(times)) break
    
    # start at previous end point and end time:
    state <- as.numeric(out[dim(out)[1], -1]) # -1 is to exclude the "time" column
    names(state) <- names_vec
    tmin <- t
    
    # change parameter values and target:
    state <- update_rule(parms, state, target)
    
    if(target == parms["target_max"] & !is.na(parms["target_min"])) target <- as.numeric(parms["target_min"])
    else target <- as.numeric(parms["target_max"])
  }
  return(df)
}

# draw the plot using ggplot2:
plot_mod <- function(df, xvar, labels, my_cols, logged = FALSE, parms = NULL) {
  gg <- ggplot(df, aes_string(xvar, "value", group = "variable", color = "variable", linetype = "variable")) + 
    geom_line() + 
    scale_colour_manual(values = my_cols, labels = labels) + 
    #scale_linetype_manual(labels = labels) + 
    theme_bw()
  
  if(logged) gg <- gg + scale_y_log10()
  
  if(!is.null(parms)) if("target_min" %in% names(parms)  & !is.na(parms["target_min"])) {
    gg <- gg +
      geom_hline(yintercept = pars_MG["target_min"], linetype = "dashed") +
      geom_hline(yintercept = pars_MG["target_max"], linetype = "dashed")
  }
  return(gg)
}

plot_mod_base <- function(df, labels, logged = FALSE, parms = NULL) {
  plot(value ~ time, data = filter(df, variable == "r0"), 
       log = "y", 
       ylim=c(parms["N0"]/2, 2e12),
       lty = 2, 
       lwd = 2, 
       xaxt = "n", 
       yaxt = "n",
       ylab = "number of tumour cells",
       col = "white")
  mtext("time (days)", 1, 2, cex = 0.7)
  abline(h = pars_fig["Ncrit"], col = "red", lty = 3)
  axis(2, pars_fig["Ncrit"], labels = expression(paste(italic("N" [crit]))), las = 2, col = "red", col.axis = "red")
  axis(2, c(10^10, 10^11, 10^12), labels = parse(text=c("10^10", "10^11", "10^12")), las = 2)
  axis(1, 500*(0:5))
  lines(value ~ time, data = filter(df, variable == "r0"), lwd = 2, lty = 2, col = "black")
  lines(value ~ time, data = filter(df, variable == "N0"), lwd = 2, col = "cyan3")
  lines(value ~ time, data = filter(df, variable == "Nacc"), lwd = 2, col = "magenta")
  lines(value ~ time, data = filter(df, variable == "Ncrit"), lwd = 2, col = "orange")
  lines(value ~ time, data = filter(df, variable == "More"), lwd = 2)
  legend("bottomright", col=c("black","black","cyan3","magenta","orange"), 
         lty = c(1, 2, rep(1, 3)),
         legend = c("untreated",
                    "ideal MTD",
                    expression(paste("contain at ", italic("N") [0])),
                    expression(paste("contain at ", italic("N" [tol]))),
                    expression(paste("contain at ", italic("N" [crit])))
         ),
         lwd=2, 
         bty = "n")
}

comgr <- function(fr0) log(1 / fr0) / 100
rhoG <- function(fr0, N0, K) comgr(fr0) / (log(K / N0))
Gompgr <- function(x, fr0, N0, K) rhoG(fr0, N0, K) * (log(K) - x * log(10))
dose_update_MG <- function(parms, state, target) {
  state["C"] <- parms["Cmax"] - state["C"]
  return(state)
}
#NortonSimon kill rate
MonroGaffney <- function(t, state, parameters, tau1 = 0, tau2 = 0) {
  # if either subpopulation becomes very small then it is assumed to stop decreasing
  # (prevents errors in the numerical solver, which can't cope with very small numbers):
  extinct <- function(x) ifelse(x < 1e-10, 0, 1)
  with(
    as.list(c(state, parameters)),
    {
      dy1 <- rhoG * log(K / (y1 + y2)) * (y1 * extinct(y1) * (1 - C - tau1) + tau2 * y2) # susceptible cells (W)
      dy2 <- rhoG * log(K / (y1 + y2)) * (y2 * extinct(y2) * (1 - D - tau2) + tau1 * y1) # resistant cells (R)
      dC <- 0
      dD <- 0
      list(c(dy1, dy2, dC, dD))
    })
}
MonroGaffneyCont<- function(t, state, parameters, tau1 = 0, tau2 = 0) {
  # if either subpopulation becomes very small then it is assumed to stop decreasing
  # (prevents errors in the numerical solver, which can't cope with very small numbers):
  extinct <- function(x) ifelse(x < 1e-10, 0, 1)
  with(
    as.list(c(state, parameters)),
    {
      dy1 <- rhoG * log(K / (y1 + y2)) * (y1 * extinct(y1) * (1 - C - tau1) + tau2 * y2) # susceptible cells (W)
      dy2 <- rhoG * log(K / (y1 + y2)) * (y2 * extinct(y2) * (1 - tau2) + tau1 * y1) # resistant cells (R)
      #   print(sprintf("C=%s,Cmax=%s",C,Cmax))
      dC <- if (C<Cmax) (dy2*y1-dy1*y2)/(y1^2) else 0#Variation of C
      list(c(dy1, dy2, dC, 0))
    })
}
# with additional parameters to account for resistance costs:
MonroGaffneyExtended <- function(t, state, parameters) {
  # if either subpopulation becomes very small then it is assumed to stop decreasing
  # (prevents errors in the numerical solver, which can't cope with very small numbers):
  extinct <- function(x) ifelse(x < 1e-10, 0, 1)
  with(
    as.list(c(state, parameters)),
    {
      dy1 <- y1 * rhoG * log(Ks / (y1 + alpha * y2)) * (1 - C) * extinct(y1) # susceptible cells (W)
      dy2 <- y2 * rhoG * log(Kr / (beta * y1 + y2)) * extinct(y2) # resistant cells (R)
      dC <- 0
      dD <- 0
      list(c(dy1, dy2, dC, dD))
    })
}
dose <- function(y1, y2, Kr, Ks, alpha, beta, level, Cmax) {
  if(y1 + y2 < level) 0
  else min(Cmax, 1 + y2 * log(Kr / (y2 + beta * y1)) / (y1 * log(Ks / (y1 + alpha * y2))))
}
# with additional parameters to account for resistance costs:
MonroGaffneyExtendedCont <- function(t, state, parameters, level = NA) {
  if(is.na(level)) level <- parameters["Nacc"]
  # if either subpopulation becomes very small then it is assumed to stop decreasing
  # (prevents errors in the numerical solver, which can't cope with very small numbers):
  extinct <- function(x) ifelse(x < 1e-10, 0, 1)
  with(
    as.list(c(state, parameters)),
    {
      dy1 <- y1 * rhoG * log(Ks / (y1 + alpha * y2)) * (1 - dose(y1, y2, Kr, Ks, alpha, beta, level, Cmax)) * extinct(y1) # susceptible cells (W)
      dy2 <- y2 * rhoG * log(Kr / (beta * y1 + y2)) * extinct(y2) # resistant cells (R)
      #dC <- if (C<Cmax) (dy2*y1-dy1*y2)/(y1^2) else 0#Variation of C
      dC <- 0
      dD <- 0
      list(c(dy1, dy2, dC, dD))
    })
}
#log kill rate  
MonroGaffneyMu <- function(t, state, parameters){
  extinct <- function(x) ifelse(x < 1e-10, 0, 1)
  with(
    as.list(c(state, parameters)),
    {
      dy1 <- y1*(rhoG * log(K / (y1 + y2)) - rhoG*log(K/1e10)*C) # susceptible cells (W)
      dy2 <- y2 * rhoG * log(K / (y1 + y2)) * extinct(y2) # resistant cells (R)
      dC <- 0
      dD <- 0
      list(c(dy1, dy2, dC, dD))
    })
  
}
times_MG <- c(0, 1000)
pars_MG <- c(K = 2e12, # carrying capacity; denoted Ninfinity in Monro and Gaffney
             rhoG = 0.0204, # basic growth rate
             Cmax = 1, # max treatment dose
             Dmax = 1, # max treatment dose
             target_min = 5e11, 
             target_max = 6e11)
#to use after modification of r0 or N0 in the parameters
init_yini<-function(parms){
  y<-NULL
  y["y1"] <- parms["N0"]-parms["r0"]
  y["y2"] <- parms["r0"]
  y["C"]<-0
  y["D"]<-0
  return(y)
}
init_pars_fig<-function(){
  pars<-c(K = 2e12, ##denoted Ninfinity in Monro and Gaffney
          rhoG = 2.47e-4*24, # basic growth rate
          #rhoG = 0.0204, # basic growth rate
          N0 = 1e10, ##true size at initiation of therapy, varied from 10^9 to 3*10^(11), 
          Nacc = 7e10,##idk if it's the right size
          Ncrit = 5e11)
  pars["r0"]<-pars["N0"]*2.3e-5
  return(pars)
}
#return the time (or row) at which point we cross Nref and then we can begin stabilising it
#if returns i then res[i-1]<=nref, res[i]>nref
find_time<- function(res, nref){
  i<-1
  for(x in (res[,2]+res[,3])){
    if(x>nref){
      return(i)
    }
    i<-i+1
  }
  return(i)
}
###time for untreated growth
time_from_to<-function(a,b,parms= pars_fig){
  a[which(a<parms["r0"])]<-parms["r0"]
  b[which(b>parms["Ncrit"])]<-parms["Ncrit"]
  res <- log(log(parms["K"]/a)/log(parms["K"]/b))
  return(res/parms["rhoG"])
}
#same as find_time
dic_find<-function(value, vector,begin, end){
  n<-(end+begin)
  if((end-begin) == 1){
    return(end)
  }
  div2<-ceiling(n/2)
  a <- vector[div2]
  if(is.na(a)) stop(paste("a is NA"))
  if(is.na(value)) stop("value is NA")
  if(a==value){
    return(div2)
  }else if(a<value){
    return(dic_find(value,vector,div2,end))
  }else{
    return(dic_find(value,vector,begin,div2))
  }
}
