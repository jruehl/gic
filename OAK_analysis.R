# This script includes the functions for the analysis of the OAK study.

# analyze data #################################################################

# compute CIs for Nelson-Aalen & Kaplan-Meier estimators in the standard (two-state) model
#   using Efron's non-parametric bootstrap
EBS_two_state <- function(data_original, survfit_original, times, BS_iter = 1000){
  
  # input:
  ## data_original: original data to base bootstrap on
  ## survfit_original: Nelson-Aalen & Kaplan-Meier information based on original data
  ## times: time points at which CIs should be evaluated
  ## BS_iter: number of bootstrap samples
  
  # output: matrix containing bootstrapped CI limits 
  #   (1st column: lower limit Nelson-Aalen estimator, 
  #    2nd column: upper limit Nelson-Aalen estimator,
  #    3rd column: lower limit Kaplan-Meier estimator,
  #    4th column: upper limit Kaplan-Meier estimator)
  
  # compute EBS NA & KM estimators
  NA_KM_EBS <- sapply(1:BS_iter, FUN = function(j){
    # sample observations
    sample <- sample(1:dim(data_original)[1], 
                     size = dim(data_original)[1], replace = TRUE)
    # retrieve NA & KM estimators
    survfit_sample <- survfit(Surv(time, status) ~ 1, data_original[sample,])
    survfit_sample <- data.frame(time = c(0, survfit_sample$time),
                                 NAE = c(0, survfit_sample$cumhaz),
                                 KME = c(1, survfit_sample$surv))
    unlist(c(survfit_sample[findInterval(times, survfit_sample$time),-1]))
  })
  
  # retrieve original NA & KM estimators
  NA_KM_original <- c(c(0, survfit_original$cumhaz)[findInterval(times, vec = c(0, survfit_original$time))],
                      c(1, survfit_original$surv)[findInterval(times, vec = c(0, survfit_original$time))])
  # compute EBS variances
  Var_EBS <- apply(NA_KM_EBS, MARGIN = 1, FUN = var)
  # compute EBS quantiles
  q_EBS <- apply(1/sqrt(Var_EBS) * (NA_KM_EBS - NA_KM_original), MARGIN = 1, 
                 FUN = quantile, probs = 0.975, na.rm = TRUE)
  # determine EBS CIs for NA & KM estimators
  CI_EBS <- cbind(
    matrix(
      pmax(0, 
           rep(NA_KM_original[1:length(times)], 2) * 
             exp(c(rep(-1, length(times)), rep(1, length(times))) * 
                   rep(q_EBS[1:length(times)], 2) * 
                   sqrt(rep(Var_EBS[1:length(times)], 2)) / 
                   rep(NA_KM_original[1:length(times)], 2))),
      ncol = 2
    ), 
    matrix(
      pmax(0,
           rep(NA_KM_original[(length(times)+1):(2*length(times))], 2)^(
             exp(c(rep(-1, length(times)), rep(1, length(times))) * 
                   rep(q_EBS[(length(times)+1):(2*length(times))], 2) * 
                   sqrt(rep(Var_EBS[(length(times)+1):(2*length(times))], 2)) / 
                   rep(log(NA_KM_original[(length(times)+1):(2*length(times))]), 2)
             ))),
      ncol = 2
    )
  )
  # adjust invalid CI values
  CI_EBS[which(is.na(CI_EBS[,1])),1:2] <- 0
  for(i in which(is.na(CI_EBS[,3]))){
    if(i < min(findInterval(survfit_original$time[survfit_original$n.event == 1], times))) 
      CI_EBS[i,3:4] <- 1
    else 
      CI_EBS[i,3:4] <- 0
  }
  CI_EBS[,3] <- pmin(1, CI_EBS[,3])
  CI_EBS[,4] <- pmin(1, CI_EBS[,4])
  
  return(CI_EBS)
  
}

# compute CIs for Nelson-Aalen & Kaplan-Meier estimators in the standard (two-state) model
#   using wild bootstrap
WBS_two_state <- function(survfit_original, times, BS_iter = 1000){
  
  # input:
  ## survfit_original: Nelson-Aalen & Kaplan-Meier information based on original data
  ## times: time points at which CIs should be evaluated
  ## BS_iter: number of bootstrap samples
  
  # output: matrix containing bootstrapped CI limits 
  #   (1st column: lower limit Nelson-Aalen estimator, 
  #    2nd column: upper limit Nelson-Aalen estimator,
  #    3rd column: lower limit Kaplan-Meier estimator,
  #    4th column: upper limit Kaplan-Meier estimator)
  
  # compute WBS Z values
  NA_KM_WBS <- sapply(1:BS_iter, FUN = function(j){
    G = ifelse(survfit_original$n.event == 0, 0, 
               sapply(survfit_original$n.event, FUN = function(l){sum(rnorm(l, mean = 0, sd = 1))}))
    Z <- data.frame(time = c(0, survfit_original$time),
                    Z_NA = c(0, sqrt(survfit_original$n) * 
                               cumsum(ifelse(survfit_original$n.risk == 0, 0, 1 / survfit_original$n.risk * G))),
                    Z_KM = c(0, -survfit_original$surv * sqrt(survfit_original$n) * 
                               cumsum(ifelse(survfit_original$n.risk == 0, 0, 1 / survfit_original$n.risk * G))))
    unlist(c(Z[findInterval(times, Z$time),-1]))
  })
  
  # retrieve original NA & KM estimators
  NA_KM_original <- c(c(0, survfit_original$cumhaz)[findInterval(times, vec = c(0, survfit_original$time))],
                      c(1, survfit_original$surv)[findInterval(times, vec = c(0, survfit_original$time))])
  # compute WBS variances
  Var_WBS <- apply(NA_KM_WBS, MARGIN = 1, FUN = var)
  # compute WBS quantiles
  q_WBS <- apply(1/sqrt(Var_WBS) * NA_KM_WBS, MARGIN = 1, 
                 FUN = quantile, probs = 0.975, na.rm = TRUE)
  # determine WBS CI for NA & KM estimators
  CI_WBS <- cbind(
    matrix(
      pmax(0,
           rep(NA_KM_original[1:length(times)], 2) * 
             exp(c(rep(-1, length(times)), rep(1, length(times))) * 
                   rep(q_WBS[1:length(times)], 2) * 
                   sqrt(rep(Var_WBS[1:length(times)], 2) / survfit_original$n) / 
                   rep(NA_KM_original[1:length(times)], 2))),
      ncol = 2
    ), 
    matrix(
      pmax(0, 
           rep(NA_KM_original[(length(times)+1):(2*length(times))], 2)^(
             exp(c(rep(-1, length(times)), rep(1, length(times))) * 
                   rep(q_WBS[(length(times)+1):(2*length(times))], 2) * 
                   sqrt(rep(Var_WBS[(length(times)+1):(2*length(times))], 2) / 
                          survfit_original$n) / 
                   rep(log(NA_KM_original[(length(times)+1):(2*length(times))]), 2)
             ))),
      ncol = 2
    )
  )
  # adjust invalid CI values
  CI_WBS[which(is.na(CI_WBS[,1])),1:2] <- 0
  for(i in which(is.na(CI_WBS[,3]))){
    if(i < min(findInterval(survfit_original$time[survfit_original$n.event == 1], times))) 
      CI_WBS[i,3:4] <- 1
    else 
      CI_WBS[i,3:4] <- 0
  }
  CI_WBS[,3] <- pmin(1, CI_WBS[,3])
  CI_WBS[,4] <- pmin(1, CI_WBS[,4])
  
  return(CI_WBS)
  
}

# functions to derive Nelson-Aalen & Aalen-Johansen estimators and compute 
#   non-parametric/wild bootstrap confidence intervals in the multi-state model
source("Bootstrap_functions.R")


# execution ####################################################################

library(readxl)
library(survival)
library(etm)

# load data
OAK <- read_excel("[path]/OAK_study.xlsx", 
                   sheet = "OAK_Clinical_Data")
OAK <- as.data.frame(OAK[,c("PtID", "TRT01P", "PFS", "PFS.CNSR", "OS", "OS.CNSR")])
colnames(OAK) <- c("id", "treatment", "PFS", "PFS_cnsr", "OS", "OS_cnsr")


## standard (two-state) model ####

OAK_atezolizumab <- data.frame(
  time = OAK[OAK$treatment == "MPDL3280A","OS"],
  status = 1 - OAK[OAK$treatment == "MPDL3280A","OS_cnsr"]
)
OAK_docetaxel <- data.frame(
  time = OAK[OAK$treatment == "Docetaxel","OS"],
  status = 1 - OAK[OAK$treatment == "Docetaxel","OS_cnsr"]
)

# derive Nelson-Aalen & Kaplan-Meier estimators
survfit_atezolizumab <- survfit(Surv(time, status) ~ 1, OAK_atezolizumab)
survfit_docetaxel <- survfit(Surv(time, status) ~ 1, OAK_docetaxel)

set.seed(6467249)

# compute EBS CI for Nelson-Aalen & Kaplan-Meier estimators
CI_EBS_atezolizumab <- EBS_two_state(OAK_atezolizumab, survfit_atezolizumab, sort(unique(OAK$OS)))
CI_EBS_docetaxel <- EBS_two_state(OAK_docetaxel, survfit_docetaxel, sort(unique(OAK$OS)))

# compute WBS CI for Nelson-Aalen & Kaplan-Meier estimators
CI_WBS_atezolizumab <- WBS_two_state(survfit_atezolizumab, sort(unique(OAK$OS)))
CI_WBS_docetaxel <- WBS_two_state(survfit_docetaxel, sort(unique(OAK$OS)))


## multi-state model ####

# prepare data for analysis with 'etm_short'
data <- rbind(
  # first transition
  data.frame(
    OAK[,c("id", "treatment")],
    entry = 0, exit = OAK$PFS,
    from = 0, to = ifelse(OAK$PFS_cnsr == 0 & (OAK$PFS < OAK$OS | OAK$OS_cnsr == 1), 1, 
                          ifelse(OAK$OS_cnsr == 0, 2, "cens"))
  ),
  # second transition (if applicable)
  data.frame(
    OAK[OAK$PFS < OAK$OS & OAK$PFS_cnsr == 0, c("id", "treatment")],
    entry = OAK[OAK$PFS < OAK$OS & OAK$PFS_cnsr == 0,]$PFS, exit = OAK[OAK$PFS < OAK$OS & OAK$PFS_cnsr == 0,]$OS,
    from = 1, to = ifelse(OAK$OS_cnsr[OAK$PFS < OAK$OS & OAK$PFS_cnsr == 0] == 0, 2, "cens")
  )
)
data <- data[order(data$id),]

# derive Nelson-Aalen & Aalen-Johansen estimators
etm_atezolizumab <- etm_short(data[data$treatment == "MPDL3280A",], s = 0, t = max(data$exit), TRUE, TRUE)
etm_docetaxel <- etm_short(data[data$treatment == "Docetaxel",], s = 0, t = max(data$exit), TRUE, TRUE)

set.seed(5215477)

# compute EBS CI for Nelson-Aalen & Aalen-Johansen estimators (progression -> death)
CI_EBS_atezolizumab_12 <- EBS(data[data$treatment == "MPDL3280A",], 
                              etm_atezolizumab, 
                              unique(sort(data$exit)), 
                              "1", "2", TRUE, TRUE, 1000)
CI_EBS_docetaxel_12 <- EBS(data[data$treatment == "Docetaxel",], 
                           etm_docetaxel, 
                           unique(sort(data$exit)), 
                           "1", "2", TRUE, TRUE, 1000)

# compute EBS CI for Nelson-Aalen & Aalen-Johansen estimators (progression -> death)
CI_WBS_atezolizumab_12 <- WBS(data[data$treatment == "MPDL3280A",], 
                           etm_atezolizumab, 
                           unique(sort(data$exit)), 
                           "1", "2", TRUE, TRUE, 1000)
CI_WBS_docetaxel_12 <- WBS(data[data$treatment == "Docetaxel",], 
                           etm_docetaxel,
                           unique(sort(data$exit)), 
                           "1", "2", TRUE, TRUE, 1000)


# create plots #################################################################

# create plot
plots <- function(time, type, etm_atezolizumab, etm_docetaxel,
                  CI_EBS_atezolizumab, CI_EBS_docetaxel, CI_WBS_atezolizumab, CI_WBS_docetaxel){
  
  # input:
  ## time: time points at which the confidence intervals are evaluated
  ## type: "NAE" (Nelson-Aalen estimator), 
  ##       "KME" (Kaplan-Meier estimator), 
  ##       "NAE_mult" (Nelson-Aalen estimator in the multi-state model), 
  ##       "AJE" (Aalen-Johansen estimator)
  ## etm_treatment (treatment = atezolizumab/docetaxel): 
  ##   Nelson-Aalen/Kaplan-Meier/Aalen-Johansen information
  ## CI_EBS_[treatment] (treatment = atezolizumab/docetaxel): 
  ##   matrix with lower/upper non-parametric bootstrap confidence limits as columns
  ## CI_WBS_treatment (treatment = atezolizumab/docetaxel): 
  ##   matrix with lower/upper wild bootstrap confidence limits as columns
  
  # derive y axis limit
  ymax <- ifelse(type %in% c("KME", "AJE"), 1, 
                 ceiling(max(CI_EBS_atezolizumab[,2], CI_EBS_docetaxel[,2],
                             CI_WBS_atezolizumab[,2], CI_WBS_docetaxel[,2]) * 10) / 10)
  
  plot <- function(time, CI_EBS, CI_WBS, etm, type, treatment){
    ggplot() +
      geom_stepribbon(aes(x = time, ymin = CI_EBS[,1], ymax = CI_EBS[,2],
                          linetype = "Non-parametric bootstrap"), 
                      colour = "darkturquoise", fill = "darkturquoise", alpha = 0.25, size = 1.3) +
      geom_stepribbon(aes(x = time, ymin = CI_WBS[,1], ymax = CI_WBS[,2],
                          linetype = "Wild bootstrap"), 
                      colour = "chocolate1", fill = "chocolate1", alpha = 0.25, size = 1.3) +
      geom_step(data = data.frame(time = time, 
                                  estimator = switch(type, 
                                                     NAE = etm$cumhaz[findInterval(time, etm$time)],
                                                     KME = etm$surv[findInterval(time, etm$time)],
                                                     NAE_mult = etm$NA_estimator[findInterval(time, etm$time), "1 2"],
                                                     AJE = etm$AJ_estimator["1", "2", findInterval(time, etm$time)])), 
                aes(x = time, y = estimator), size = 1.1) +
      scale_x_continuous(expand = expansion(mult = c(0.03,0)), limits = c(0, 28), breaks = seq(0,28,5), labels = c("0","","10","","20","")) +
      scale_y_continuous(expand = expansion(mult = c(0.03,0)), limits = c(0, ymax), 
                         breaks = seq(0, ymax, ifelse(type %in% c("KME", "AJE"), 0.25, 1)), 
                         labels = if(type %in% c("KME", "AJE")) c("0.0","","0.5","","1.0")
                                     else c(rbind(seq(0,round(ymax), 2), ""))[1:length(seq(0, ymax, 1))]) +
      ggtitle(treatment) + labs(x = "Time [months]", y = ifelse(treatment == "Atezolizumab", 
                                                                switch(type,
                                                                       NAE = "Cumulative hazard",
                                                                       KME = "Survival probability",
                                                                       NAE_mult = "Cumulative hazard",
                                                                       AJE = "Cumulative incidence"), 
                                                                "")) +
      #scale_fill_manual(name = "", values = c("darkturquoise", "chocolate1"), 
      #                  labels = c("Non-parametric bootstrap", "Wild bootstrap")) +
      scale_linetype_manual("", values = c(1, 2),
                            guide = guide_legend(nrow = 1,
                                                 override.aes = 
                                                   list(color = c("darkturquoise", "chocolate1"),
                                                        fill = c("darkturquoise", "chocolate1"),
                                                        size = 1.3))) +
      theme(panel.background = element_blank(), panel.grid.major = element_line(colour = "gray95"),
            axis.text = element_text(colour = "black"), text = element_text(size = 23),
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(margin = margin(t=8, r=0, b=0, l=0)),
            axis.title.y = element_text(margin = margin(t=0, r=8, b=0, l=0)),
            plot.title = element_text(hjust = 0.5, size = 26),
            plot.margin =  unit(c(1,1,1,1), "lines"))
  }
  
  p1 <- plot(time, CI_EBS_atezolizumab, CI_WBS_atezolizumab, etm_atezolizumab, type, "Atezolizumab")
  p2 <- plot(time, CI_EBS_docetaxel, CI_WBS_docetaxel, etm_docetaxel, type, "Docetaxel")
  
  grid.arrange(grobs = list(p1 + theme(legend.position = "none"), 
                            p2 + theme(legend.position = "none"),
                            get_legend(p1)), 
               layout_matrix = matrix(c(1,2, 3,3), byrow = TRUE, nrow = 2),
               heights = c(0.8, 0.2))
  
}


# plots ########################################################################

library(ggplot2)
library(gridExtra)
library(cowplot)
library(pammtools)

# cumulative hazard
plots(unique(sort(OAK$OS)), "NAE", survfit_atezolizumab, survfit_docetaxel,
      CI_EBS_atezolizumab[,1:2], CI_EBS_docetaxel[,1:2], CI_WBS_atezolizumab[,1:2], CI_WBS_docetaxel[,1:2])

# survival probability
plots(unique(sort(OAK$OS)), "KME", survfit_atezolizumab, survfit_docetaxel,
      CI_EBS_atezolizumab[,3:4], CI_EBS_docetaxel[,3:4], CI_WBS_atezolizumab[,3:4], CI_WBS_docetaxel[,3:4])

# cumulative hazard (multi-state model)
plots(unique(sort(data$exit)), "NAE_mult", etm_atezolizumab, etm_docetaxel,
      CI_EBS_atezolizumab_12[,1:2], CI_EBS_docetaxel_12[,1:2], CI_WBS_atezolizumab_12[,1:2], CI_WBS_docetaxel_12[,1:2])

# cumulative incidence (multi-state model)
plots(unique(sort(data$exit)), "AJE", etm_atezolizumab, etm_docetaxel,
      CI_EBS_atezolizumab_12[,3:4], CI_EBS_docetaxel_12[,3:4], CI_WBS_atezolizumab_12[,3:4], CI_WBS_docetaxel_12[,3:4])
