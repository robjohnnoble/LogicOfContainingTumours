# Impact of resistance costs on the best possible outcome

plot_best_outcome_orientation1 <- function(parms) {
  K_s <- parms["K"]
  ylim <- 35
  xlim <- 15
  best_df <- data.frame(beta_vec = seq(0, xlim, length = 1000))
  best_df$line1 <- ifelse(best_df$beta_vec <= 1, 
                          K_s / parms["Nacc"], 
                          K_s / (parms["Nacc"] * best_df$beta_vec))
  best_df$line2 <- ifelse(best_df$beta_vec <= 1, 
                          K_s / parms["Ncrit"], 
                          K_s / (parms["Ncrit"] * best_df$beta_vec))
  best_df$line1bounded <- ifelse(best_df$beta_vec <= 1, 
                                 NA, pmax(best_df$line1, K_s / parms["Ncrit"]))
  best_df$line1bounded2 <- ifelse(best_df$beta_vec <= parms["Ncrit"] / parms["Nacc"], 
                                  NA, best_df$line1)
  best_df$line1bounded3 <- ifelse(best_df$beta_vec >= parms["Ncrit"] / parms["Nacc"], 
                                  NA, pmin(best_df$line1, K_s / parms["Nacc"]))
  best_df$line1bounded4 <- ifelse(best_df$beta_vec <= 1, 
                                  NA, pmin(best_df$line1, K_s / parms["Ncrit"]))
  xlab <- expression(paste(italic(beta)))
  ylab <- expression(paste(italic("K" [s]), " / ", italic("K" [r])))
  ybreak1 <- expression(paste(italic("K" [s]), " / ", italic("N" [crit])))
  ybreak2 <- expression(paste(italic("K" [s]), " / ", italic("N" [tol])))
  xbreak1 <- expression(paste(italic("N" [crit]), " / ", italic("N" [tol])))
  ggplot() + 
    geom_ribbon(aes(x = beta_vec, ymin = K_s / parms["Nacc"], ymax = ylim, fill = "#bbff00"), data = best_df) + 
    geom_ribbon(aes(x = beta_vec, ymin = line1bounded, ymax = K_s / parms["Nacc"], fill = "#ddff00"), data = best_df) + 
    geom_ribbon(aes(x = beta_vec, ymin = line1bounded2, ymax = K_s / parms["Ncrit"], fill = "#ffff00"), data = best_df) + 
    geom_ribbon(aes(x = beta_vec, ymin = K_s / parms["Ncrit"], ymax = line1bounded3, fill = "gold"), data = best_df) + 
    geom_ribbon(aes(x = beta_vec, ymin = line2, ymax = line1bounded4, fill = "orange"), data = best_df) + 
    geom_ribbon(aes(x = beta_vec, ymin = 0, ymax = line2, fill = "red"), data = best_df) + 
    geom_line(aes(x = beta_vec, y = line1), data = best_df) + 
    geom_line(aes(x = beta_vec, y = line2), data = best_df) + 
    geom_hline(yintercept = K_s / parms["Nacc"], linetype = "dashed") + 
    geom_segment(aes(x = 1, y = K_s / parms["Ncrit"], 
                     xend = xlim, yend = K_s / parms["Ncrit"]), linetype = "dashed") + 
    scale_y_continuous(limits = c(0, ylim), name = ylab, expand = c(0, 0), 
                       breaks = c(0, K_s / parms["Ncrit"], K_s / parms["Nacc"], ylim), 
                       labels = c(0, ybreak1, ybreak2, ylim)) +
    scale_x_continuous(name = xlab, expand = c(0, 0), 
                       breaks = c(0, 1, parms["Ncrit"] / parms["Nacc"], xlim), 
                       labels = c(0, 1, xbreak1, xlim)) +
    theme_classic() + 
    theme(legend.text = element_text(margin = margin(t = 5, b = 5, unit = "pt"))) + 
    scale_fill_manual(values=c("#bbff00", "#ddff00", "#ffff00", "gold", "orange", "red"),
                      name="eventual outcome",
                      labels=c("tolerable\n ", 
                               "tolerable idCont\nintolerable idMTD", 
                               "tolerable idCont\nlethal idMTD", 
                               "intolerable\n ", 
                               "intolerable idCont\nlethal idMTD", 
                               "lethal\n "))
}

plot_best_outcome_orientation2 <- function(parms) {
  K_s <- parms["K"]
  ylim <- 15
  xlim <- 32
  best_df <- data.frame(gamma_vec = seq(0, xlim, length = 10000))
  best_df$line1 <- ifelse(best_df$gamma_vec > K_s / parms["Nacc"], 
                          NA, 
                          K_s / (parms["Nacc"] * best_df$gamma_vec))
  best_df$line2 <- ifelse(best_df$gamma_vec > K_s / parms["Ncrit"], 
                          NA, 
                          K_s / (parms["Ncrit"] * best_df$gamma_vec))
  best_df$line1bounded <- ifelse(best_df$gamma_vec < K_s / parms["Ncrit"], 
                                 NA, best_df$line1)
  best_df$line1bounded2 <- ifelse(best_df$gamma_vec > K_s / parms["Ncrit"], 
                                  NA, best_df$line1)
  best_df$line1bounded3 <- ifelse(best_df$gamma_vec > K_s / parms["Ncrit"], 
                                  NA, pmin(best_df$line1, ylim))
  best_df$line2bounded <- ifelse(best_df$gamma_vec >= K_s / parms["Ncrit"], 
                                 NA, best_df$line2)
  best_df$line2bounded2 <- pmin(best_df$line2, ylim)
  ylab <- expression(paste(italic(beta)))
  xlab <- expression(paste(italic("K" [s]), " / ", italic("K" [r])))
  xbreak1 <- expression(paste(italic("K" [s]), " / ", italic("N" [crit])))
  xbreak2 <- expression(paste(italic("K" [s]), " / ", italic("N" [tol])))
  ybreak1 <- expression(paste(italic("N" [crit]), " / ", italic("N" [tol])))
  ggplot() + 
    geom_ribbon(aes(x = c(K_s / parms["Nacc"], xlim), ymin = 0, ymax = ylim, fill = "#bbff00")) + 
    geom_ribbon(aes(x = gamma_vec, ymin = line1bounded, ymax = ylim, fill = "#ddff00"), data = best_df) + 
    geom_ribbon(aes(x = gamma_vec, ymin = line1bounded2, ymax = ylim, fill = "#ffff00"), data = best_df) + 
    geom_ribbon(aes(x = gamma_vec, ymin = 0, ymax = line1bounded, fill = "gold"), data = best_df) + 
    geom_ribbon(aes(x = gamma_vec, ymin = line2bounded, ymax = line1bounded3, fill = "orange"), data = best_df) + 
    geom_ribbon(aes(x = gamma_vec, ymin = 0, ymax = line2bounded2, fill = "red"), data = best_df) + 
    geom_line(aes(x = gamma_vec, y = line1), data = best_df) + 
    geom_line(aes(x = gamma_vec, y = line2), data = best_df) + 
    geom_vline(xintercept = K_s / parms["Nacc"], linetype = "dashed") + 
    geom_vline(xintercept = K_s / parms["Ncrit"], linetype = "dashed") + 
    # geom_vline(xintercept = 10, linetype = "dashed") + 
    # geom_hline(yintercept = 4, linetype = "dashed") + 
    geom_segment(aes(x = K_s / parms["Ncrit"], y = 0, 
                     xend = K_s / parms["Ncrit"], yend = 1)) + 
    geom_segment(aes(x = K_s / parms["Nacc"], y = 0, 
                     xend = K_s / parms["Nacc"], yend = 1)) + 
    scale_x_continuous(name = xlab, expand = c(0, 0), 
                       breaks = c(0, K_s / parms["Ncrit"], 10, 20, K_s / parms["Nacc"]), 
                       labels = c(0, xbreak1, 10, 20, xbreak2)) +
    scale_y_continuous(limits = c(0, ylim), name = ylab, expand = c(0, 0), 
                       breaks = c(0, 1, 5, 10, parms["Ncrit"] / parms["Nacc"], ylim), 
                       labels = c(0, 1, 5, 10, ybreak1, ylim)) +
    theme_classic() + 
    theme(legend.text = element_text(margin = margin(t = 5, b = 5, unit = "pt"))) + 
    scale_fill_manual(values=c("#bbff00", "#ddff00", "#ffff00", "gold", "orange", "red"),
                      name="eventual outcome",
                      labels=c("tolerable\n ", 
                               "tolerable idCont\nintolerable idMTD", 
                               "tolerable idCont\nlethal idMTD", 
                               "intolerable\n ", 
                               "intolerable idCont\nlethal idMTD", 
                               "lethal\n "))
}



