# This script includes the functions that produce results of the first simulation study.

# generate data ################################################################

# generate data with exponentially distributed survival times
data_exponential <- function(n, m, HR, lambda){
  
  # input:
  ## n: number of subjects
  ## m: number of events to observe before censoring is imposed
  ## HR: hazard ratio (treatment vs. control group)
  ## lambda: distribution parameter (in the control group)
  
  # output: data set with the subsequent columns:
  ## treatment: treatment group ('Treatment' vs. 'Control')
  ## entry: study entry [calendar time]
  ## event: event [calendar time]
  ## censoring: last observed time point [calendar time]
  ## study_time: event time [study time]
  ## status: status indicator (0 = censored, 1 = observed event)
  ## recruited_at_entry: number of previously recruited subjects at entry
  
  data <- tibble(
    treatment = c(rep("Control", n/2), rep("Treatment", n/2)),
    # limit recruitment period to the m/n quantile of the event distribution
    #   to reduce probability of subjects entering after the observation period
    entry = runif(n, 0, qexp(m/n, lambda)),
    event = entry + c(rexp(n/2, lambda), rexp(n/2, lambda*HR)),
    censoring = ifelse(event <= sort(event)[m], event, sort(event)[m]),
    study_time = censoring - entry,
    status = ifelse(event <= sort(event)[m], 1, 0),
    recruited_at_entry = rank(entry) - 1
  )
  
}

# generate data with Weibull distributed survival times
data_weibull <- function(n, m, HR, shape, scale){
  
  # input:
  ## n: number of subjects
  ## m: number of events to observe before censoring is imposed
  ## HR: hazard ratio (treatment vs. control group)
  ## shape: distribution shape parameter
  ## scale: distribution scale parameter (in the control group)
  
  # output: data set with the subsequent columns:
  ## treatment: treatment group ('Treatment' vs. 'Control')
  ## entry: study entry [calendar time]
  ## event: event [calendar time]
  ## censoring: last observed time point [calendar time]
  ## study_time: event time [study time]
  ## status: status indicator (0 = censored, 1 = observed event)
  ## recruited_at_entry: number of previously recruited subjects at entry
  
  data <- tibble(
    treatment = c(rep("Control", n/2), rep("Treatment", n/2)),
    # limit recruitment period to the m/n quantile of the event distribution
    #   to reduce probability of subjects entering after the observation period
    entry = runif(n, 0, qweibull(m/n, shape, scale)),
    event = entry + c(rweibull(n/2, shape, scale), 
                      rweibull(n/2, shape, scale/HR^(1/shape))),
    censoring = ifelse(event <= sort(event)[m], event, sort(event)[m]),
    study_time = censoring - entry,
    status = ifelse(event <= sort(event)[m], 1, 0),
    recruited_at_entry = rank(entry) - 1
  )
  
}

# generate data and derive relevant parameters
prepare_data <- function(i, scenario, n, m, HR, lambda, shape, scale){
  
  # input:
  ## i: parameter for parallel processing
  ## scenario: scenario number (1 to 6)
  ## n: number of subjects
  ## m: number of events to observe before censoring is imposed
  ## HR: hazard ratio (treatment vs. control group)
  ## lambda: distribution parameter (in the control group) (exponential scenario)
  ## shape: distribution parameter (Weibull scenario)
  ## scale: distribution parameter (in the control group) (Weibull scenario)
  
  # output:
  ## indicator: indicator implying whether there are at least two observed events 
  ##   in each treatment group
  ## HR_[parameter]_[model] 
  ##  (parameter = treatment/entry/recruited, model = standard/model1/model2): 
  ##  estimated hazard ratio for the specified parameter in the specified group
  ## CI_[parameter]_[model]
  ##   (parameter = treatment/entry/recruited, model = standard/model1/model2): 
  ##   estimated confidence interval for the hazard ratio
  ## breslow_[model] (model = standard/model1/model2): 
  ##   Breslow estimator of the cumulative baseline hazard in the specified model
  
  if(scenario %in% 1:3){
    data <- data_exponential(n, m, HR, lambda)
  }else if(scenario %in% 4:6){
    data <- data_weibull(n, m, HR, shape, scale)
  }
  
  indicator <- 2 <= dim(subset(data, status == 1 & treatment == "Treatment"))[1] & dim(subset(data, status == 1 & treatment == "Treatment"))[1] <= m - 2
  HR_treatment_standard <- summary(coxph(Surv(study_time, status) ~ treatment, data))$coefficients[2]
  HR_treatment_model1 <- summary(coxph(Surv(study_time, status) ~ treatment + entry, data))$coefficients[1,2]
  HR_treatment_model2 <- summary(coxph(Surv(study_time, status) ~ treatment + entry + recruited_at_entry, data))$coefficients[1,2]
  HR_entry_model1 <- summary(coxph(Surv(study_time, status) ~ treatment + entry, data))$coefficients[2,2]
  HR_entry_model2 <- summary(coxph(Surv(study_time, status) ~ treatment + entry + recruited_at_entry, data))$coefficients[2,2]
  HR_recruited_model2 <- summary(coxph(Surv(study_time, status) ~ treatment + entry + recruited_at_entry, data))$coefficients[3,2]
  CI_treatment_standard <- summary(coxph(Surv(study_time, status) ~ treatment, data))$conf.int[1,c(3,4)]
  CI_treatment_model1 <- summary(coxph(Surv(study_time, status) ~ treatment + entry, data))$conf.int[1,c(3,4)]
  CI_treatment_model2 <- summary(coxph(Surv(study_time, status) ~ treatment + entry + recruited_at_entry, data))$conf.int[1,c(3,4)]
  CI_entry_model1 <- summary(coxph(Surv(study_time, status) ~ treatment + entry, data))$conf.int[2,c(3,4)]
  CI_entry_model2 <- summary(coxph(Surv(study_time, status) ~ treatment + entry + recruited_at_entry, data))$conf.int[2,c(3,4)]
  CI_recruited_model2 <- summary(coxph(Surv(study_time, status) ~ treatment + entry + recruited_at_entry, data))$conf.int[3,c(3,4)]
  breslow_standard <- basehaz(coxph(Surv(study_time, status) ~ treatment, data), centered = FALSE)
  breslow_model1 <- basehaz(coxph(Surv(study_time, status) ~ treatment + entry, data), centered = FALSE)
  breslow_model2 <- basehaz(coxph(Surv(study_time, status) ~ treatment + entry + recruited_at_entry, data), centered = FALSE)
  
  return(list(indicator, 
              HR_treatment_standard, HR_treatment_model1, HR_treatment_model2,
              HR_entry_model1, HR_entry_model2, HR_recruited_model2,
              CI_treatment_standard, CI_treatment_model1, CI_treatment_model2,
              CI_entry_model1, CI_entry_model2, CI_recruited_model2,
              breslow_standard, breslow_model1, breslow_model2))
  
}

# run simulations ##############################################################

# run simulations and save summarized results within global variables:
## indicator_[scenario] (scenario = 1 - 6): 
##   indicator implying in which iterations there are at least two observed events in each treatment group 
## HR_[parameter]_[model]_mean_(indicator_)[scenario]
##   (parameter = treatment/entry/recruited, model = standard/model1/model2, scenario = 1 - 6):
##   mean of the estimated hazard ratios for the specified parameter in the specified model and scenario 
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## HR_[parameter]_[model]_median_(indicator_)[scenario]
##   (parameter = treatment/entry/recruited, model = standard/model1/model2, scenario = 1 - 6):
##   median of the estimated hazard ratios for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## CI_[parameter]_[model]_median_width_(indicator_)[scenario]
##   (parameter = treatment/entry/recruited, model = standard/model1/model2, scenario = 1 - 6)
##   median with of the confidence intervals of the estimated hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## empirical_CI_[parameter]_[model]_(indicator_)[scenario]
##   (parameter = treatment/entry/recruited, model = standard/model1/model2, scenario = 1 - 6)
##   empirical confidence interval of the estimated hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## coverage_[parameter]_[model]_(indicator_)[scenario]
##   (parameter = treatment/entry/recruited, model = standard/model1/model2, scenario = 1 - 6)
##   coverage of the confidence intervals of the estimated hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## bias_[parameter]_[model]_mean_(indicator_)[scenario]
##   (parameter = treatment/entry/recruited, model = standard/model1/model2, scenario = 1 - 6)
##   mean bias of the estimated hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## bias_[parameter]_[model]_median_(indicator_)[scenario]
##   (parameter = treatment/entry/recruited, model = standard/model1/model2, scenario = 1 - 6)
##   median bias of the estimated hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## RMSE_[parameter]_[model]_(indicator_)[scenario]
##   (parameter = treatment/entry/recruited, model = standard/model1/model2, scenario = 1 - 6)
##   root mean squared error of the estimated hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## breslow_[model]_[scenario] (model = standard/model1/model2, scenario = 1 - 6):
##   Breslow estimator of the cumulative baseline hazard in the specified model and scenario
simulate <- function(seed, scenario, HR, lambda = 1, shape = 0.5, scale = 1, iter = 100000){
    
  # input: 
  ## seed: seed to use for simulations
  ## scenario: scenario number (1 to 6)
  ## HR: hazard ratio (treatment vs. control group)
  ## lambda: distribution parameter (in the control group) (exponential scenario)
  ## shape: distribution parameter (Weibull scenario)
  ## scale: distribution parameter (in the control group) (Weibull scenario)
  ## iter: number of iterations
  
  n <- switch(scenario,
              50,50,26,50,50,26)
  m <- switch(scenario,
              25,10,13,25,10,13)
  
  # TODO 
  # adapt the number of cores according to your needs
  # (use 8 cores to reproduce the results presented in the manuscript)
  cl <- makeCluster(detectCores() - 1, type="SOCK")
  clusterEvalQ(cl, library(tibble))
  clusterEvalQ(cl, library(survival))
  clusterExport(cl, "data_exponential")
  clusterExport(cl, "data_weibull")
  clusterSetRNGStream(cl, iseed = seed)
  sim <- pblapply(X = 1:iter, FUN = prepare_data, scenario, n, m, HR, lambda, shape, scale, cl = cl)
  stopCluster(cl)
  
  
  indicator <- as.numeric(unlist(lapply(sim, `[[`, 1)))
  HR_treatment_standard <- unlist(lapply(sim, `[[`, 2))
  HR_treatment_model1 <- unlist(lapply(sim, `[[`, 3))
  HR_treatment_model2 <- unlist(lapply(sim, `[[`, 4))
  HR_entry_model1 <- unlist(lapply(sim, `[[`, 5))
  HR_entry_model2 <- unlist(lapply(sim, `[[`, 6))
  HR_recruited_model2 <- unlist(lapply(sim, `[[`, 7))
  
  CI_treatment_standard_lower <- unlist(lapply(lapply(sim, `[[`, 8),`[[`, 1))
  CI_treatment_standard_upper <- unlist(lapply(lapply(sim, `[[`, 8),`[[`, 2))
  CI_treatment_model1_lower <- unlist(lapply(lapply(sim, `[[`, 9),`[[`, 1))
  CI_treatment_model1_upper <- unlist(lapply(lapply(sim, `[[`, 9),`[[`, 2))
  CI_treatment_model2_lower <- unlist(lapply(lapply(sim, `[[`, 10),`[[`, 1))
  CI_treatment_model2_upper <- unlist(lapply(lapply(sim, `[[`, 10),`[[`, 2))
  CI_entry_model1_lower <- unlist(lapply(lapply(sim, `[[`, 11),`[[`, 1))
  CI_entry_model1_upper <- unlist(lapply(lapply(sim, `[[`, 11),`[[`, 2))
  CI_entry_model2_lower <- unlist(lapply(lapply(sim, `[[`, 12),`[[`, 1))
  CI_entry_model2_upper <- unlist(lapply(lapply(sim, `[[`, 12),`[[`, 2))
  CI_recruited_model2_lower <- unlist(lapply(lapply(sim, `[[`, 13),`[[`, 1))
  CI_recruited_model2_upper <- unlist(lapply(lapply(sim, `[[`, 13),`[[`, 2))
  
  
  assign(paste0("breslow_standard_", scenario), lapply(sim, `[[`, 14), envir = .GlobalEnv)
  assign(paste0("breslow_model1_", scenario), lapply(sim, `[[`, 15), envir = .GlobalEnv)
  assign(paste0("breslow_model2_", scenario), lapply(sim, `[[`, 16), envir = .GlobalEnv)
  
  
  assign(paste0("HR_treatment_standard_mean_", scenario), mean(HR_treatment_standard), envir = .GlobalEnv)
  assign(paste0("HR_treatment_model1_mean_", scenario), mean(HR_treatment_model1), envir = .GlobalEnv)
  assign(paste0("HR_treatment_model2_mean_", scenario), mean(HR_treatment_model2), envir = .GlobalEnv)
  assign(paste0("HR_entry_model1_mean_", scenario), mean(HR_entry_model1), envir = .GlobalEnv)
  assign(paste0("HR_entry_model2_mean_", scenario), mean(HR_entry_model2), envir = .GlobalEnv)
  assign(paste0("HR_recruited_model2_mean_", scenario), mean(HR_recruited_model2), envir = .GlobalEnv)
  
  assign(paste0("HR_treatment_standard_median_", scenario), median(HR_treatment_standard), envir = .GlobalEnv)
  assign(paste0("HR_treatment_model1_median_", scenario), median(HR_treatment_model1), envir = .GlobalEnv)
  assign(paste0("HR_treatment_model2_median_", scenario), median(HR_treatment_model2), envir = .GlobalEnv)
  assign(paste0("HR_entry_model1_median_", scenario), median(HR_entry_model1), envir = .GlobalEnv)
  assign(paste0("HR_entry_model2_median_", scenario), median(HR_entry_model2), envir = .GlobalEnv)
  assign(paste0("HR_recruited_model2_median_", scenario), median(HR_recruited_model2), envir = .GlobalEnv)
  
  assign(paste0("CI_treatment_standard_median_width_", scenario), median(CI_treatment_standard_upper - CI_treatment_standard_lower), envir = .GlobalEnv)
  assign(paste0("CI_treatment_model1_median_width_", scenario), median(CI_treatment_model1_upper - CI_treatment_model1_lower), envir = .GlobalEnv)
  assign(paste0("CI_treatment_model2_median_width_", scenario), median(CI_treatment_model2_upper - CI_treatment_model2_lower), envir = .GlobalEnv)
  assign(paste0("CI_entry_model1_median_width_", scenario), median(CI_entry_model1_upper - CI_entry_model1_lower), envir = .GlobalEnv)
  assign(paste0("CI_entry_model2_median_width_", scenario), median(CI_entry_model2_upper - CI_entry_model2_lower), envir = .GlobalEnv)
  assign(paste0("CI_recruited_model2_median_width_", scenario), median(CI_recruited_model2_upper - CI_recruited_model2_lower), envir = .GlobalEnv)
  
  assign(paste0("empirical_CI_treatment_standard_", scenario), quantile(HR_treatment_standard, c(0.025, 0.975)), envir = .GlobalEnv)
  assign(paste0("empirical_CI_treatment_model1_", scenario), quantile(HR_treatment_model1, c(0.025, 0.975)), envir = .GlobalEnv)
  assign(paste0("empirical_CI_treatment_model2_", scenario), quantile(HR_treatment_model2, c(0.025, 0.975)), envir = .GlobalEnv)
  assign(paste0("empirical_CI_entry_model1_", scenario), quantile(HR_entry_model1, c(0.025, 0.975)), envir = .GlobalEnv)
  assign(paste0("empirical_CI_entry_model2_", scenario), quantile(HR_entry_model2, c(0.025, 0.975)), envir = .GlobalEnv)
  assign(paste0("empirical_CI_recruited_model2_", scenario), quantile(HR_recruited_model2, c(0.025, 0.975)), envir = .GlobalEnv)
  
  assign(paste0("coverage_treatment_standard_", scenario), mean(CI_treatment_standard_lower <= HR & HR <= CI_treatment_standard_upper), envir = .GlobalEnv)
  assign(paste0("coverage_treatment_model1_", scenario), mean(CI_treatment_model1_lower <= HR & HR <= CI_treatment_model1_upper), envir = .GlobalEnv)
  assign(paste0("coverage_treatment_model2_", scenario), mean(CI_treatment_model2_lower <= HR & HR <= CI_treatment_model2_upper), envir = .GlobalEnv)
  assign(paste0("coverage_entry_model1_", scenario), mean(CI_entry_model1_lower <= 1 & 1 <= CI_entry_model1_upper), envir = .GlobalEnv)
  assign(paste0("coverage_entry_model2_", scenario), mean(CI_entry_model2_lower <= 1 & 1 <= CI_entry_model2_upper), envir = .GlobalEnv)
  assign(paste0("coverage_recruited_model2_", scenario), mean(CI_recruited_model2_lower <= 1 & 1 <= CI_recruited_model2_upper), envir = .GlobalEnv)
  
  assign(paste0("bias_treatment_standard_mean_", scenario), mean(HR_treatment_standard - HR), envir = .GlobalEnv)
  assign(paste0("bias_treatment_model1_mean_", scenario), mean(HR_treatment_model1 - HR), envir = .GlobalEnv)
  assign(paste0("bias_treatment_model2_mean_", scenario), mean(HR_treatment_model2 - HR), envir = .GlobalEnv)
  assign(paste0("bias_entry_model1_mean_", scenario), mean(HR_entry_model1 - 1), envir = .GlobalEnv)
  assign(paste0("bias_entry_model2_mean_", scenario), mean(HR_entry_model2 - 1), envir = .GlobalEnv)
  assign(paste0("bias_recruited_model2_mean_", scenario), mean(HR_recruited_model2 - 1), envir = .GlobalEnv)
  
  assign(paste0("bias_treatment_standard_median_", scenario), median(HR_treatment_standard - HR), envir = .GlobalEnv)
  assign(paste0("bias_treatment_model1_median_", scenario), median(HR_treatment_model1 - HR), envir = .GlobalEnv)
  assign(paste0("bias_treatment_model2_median_", scenario), median(HR_treatment_model2 - HR), envir = .GlobalEnv)
  assign(paste0("bias_entry_model1_median_", scenario), median(HR_entry_model1 - 1), envir = .GlobalEnv)
  assign(paste0("bias_entry_model2_median_", scenario), median(HR_entry_model2 - 1), envir = .GlobalEnv)
  assign(paste0("bias_recruited_model2_median_", scenario), median(HR_recruited_model2 - 1), envir = .GlobalEnv)
  
  assign(paste0("RMSE_treatment_standard_", scenario), sqrt(mean((HR_treatment_standard - HR)^2)), envir = .GlobalEnv)
  assign(paste0("RMSE_treatment_model1_", scenario), sqrt(mean((HR_treatment_model1 - HR)^2)), envir = .GlobalEnv)
  assign(paste0("RMSE_treatment_model2_", scenario), sqrt(mean((HR_treatment_model2 - HR)^2)), envir = .GlobalEnv)
  assign(paste0("RMSE_entry_model1_", scenario), sqrt(mean((HR_entry_model1 - 1)^2)), envir = .GlobalEnv)
  assign(paste0("RMSE_entry_model2_", scenario), sqrt(mean((HR_entry_model2 - 1)^2)), envir = .GlobalEnv)
  assign(paste0("RMSE_recruited_model2_", scenario), sqrt(mean((HR_recruited_model2 - 1)^2)), envir = .GlobalEnv)
  
  
  if(scenario %in% c(2,3,5,6)){
    
    assign(paste0("indicator_", scenario), indicator, envir = .GlobalEnv)
    
    assign(paste0("HR_treatment_standard_mean_indicator_", scenario), mean(HR_treatment_standard[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("HR_treatment_model1_mean_indicator_", scenario), mean(HR_treatment_model1[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("HR_treatment_model2_mean_indicator_", scenario), mean(HR_treatment_model2[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("HR_entry_model1_mean_indicator_", scenario), mean(HR_entry_model1[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("HR_entry_model2_mean_indicator_", scenario), mean(HR_entry_model2[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("HR_recruited_model2_mean_indicator_", scenario), mean(HR_recruited_model2[indicator == 1]), envir = .GlobalEnv)
    
    assign(paste0("HR_treatment_standard_median_indicator_", scenario), median(HR_treatment_standard[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("HR_treatment_model1_median_indicator_", scenario), median(HR_treatment_model1[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("HR_treatment_model2_median_indicator_", scenario), median(HR_treatment_model2[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("HR_entry_model1_median_indicator_", scenario), median(HR_entry_model1[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("HR_entry_model2_median_indicator_", scenario), median(HR_entry_model2[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("HR_recruited_model2_median_indicator_", scenario), median(HR_recruited_model2[indicator == 1]), envir = .GlobalEnv)
    
    assign(paste0("CI_treatment_standard_median_width_indicator_", scenario), median(CI_treatment_standard_upper[indicator == 1] - CI_treatment_standard_lower[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("CI_treatment_model1_median_width_indicator_", scenario), median(CI_treatment_model1_upper[indicator == 1] - CI_treatment_model1_lower[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("CI_treatment_model2_median_width_indicator_", scenario), median(CI_treatment_model2_upper[indicator == 1] - CI_treatment_model2_lower[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("CI_entry_model1_median_width_indicator_", scenario), median(CI_entry_model1_upper[indicator == 1] - CI_entry_model1_lower[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("CI_entry_model2_median_width_indicator_", scenario), median(CI_entry_model2_upper[indicator == 1] - CI_entry_model2_lower[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("CI_recruited_model2_median_width_indicator_", scenario), median(CI_recruited_model2_upper[indicator == 1] - CI_recruited_model2_lower[indicator == 1]), envir = .GlobalEnv)
    
    assign(paste0("empirical_CI_treatment_standard_indicator_", scenario), quantile(HR_treatment_standard[indicator == 1], c(0.025, 0.975)), envir = .GlobalEnv)
    assign(paste0("empirical_CI_treatment_model1_indicator_", scenario), quantile(HR_treatment_model1[indicator == 1], c(0.025, 0.975)), envir = .GlobalEnv)
    assign(paste0("empirical_CI_treatment_model2_indicator_", scenario), quantile(HR_treatment_model2[indicator == 1], c(0.025, 0.975)), envir = .GlobalEnv)
    assign(paste0("empirical_CI_entry_model1_indicator_", scenario), quantile(HR_entry_model1[indicator == 1], c(0.025, 0.975)), envir = .GlobalEnv)
    assign(paste0("empirical_CI_entry_model2_indicator_", scenario), quantile(HR_entry_model2[indicator == 1], c(0.025, 0.975)), envir = .GlobalEnv)
    assign(paste0("empirical_CI_recruited_model2_indicator_", scenario), quantile(HR_recruited_model2[indicator == 1], c(0.025, 0.975)), envir = .GlobalEnv)
    
    assign(paste0("coverage_treatment_standard_indicator_", scenario), mean(CI_treatment_standard_lower[indicator == 1] <= HR & HR <= CI_treatment_standard_upper[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("coverage_treatment_model1_indicator_", scenario), mean(CI_treatment_model1_lower[indicator == 1] <= HR & HR <= CI_treatment_model1_upper[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("coverage_treatment_model2_indicator_", scenario), mean(CI_treatment_model2_lower[indicator == 1] <= HR & HR <= CI_treatment_model2_upper[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("coverage_entry_model1_indicator_", scenario), mean(CI_entry_model1_lower[indicator == 1] <= 1 & 1 <= CI_entry_model1_upper[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("coverage_entry_model2_indicator_", scenario), mean(CI_entry_model2_lower[indicator == 1] <= 1 & 1 <= CI_entry_model2_upper[indicator == 1]), envir = .GlobalEnv)
    assign(paste0("coverage_recruited_model2_indicator_", scenario), mean(CI_recruited_model2_lower[indicator == 1] <= 1 & 1 <= CI_recruited_model2_upper[indicator == 1]), envir = .GlobalEnv)
    
    assign(paste0("bias_treatment_standard_mean_indicator_", scenario), mean(HR_treatment_standard[indicator == 1] - HR), envir = .GlobalEnv)
    assign(paste0("bias_treatment_model1_mean_indicator_", scenario), mean(HR_treatment_model1[indicator == 1] - HR), envir = .GlobalEnv)
    assign(paste0("bias_treatment_model2_mean_indicator_", scenario), mean(HR_treatment_model2[indicator == 1] - HR), envir = .GlobalEnv)
    assign(paste0("bias_entry_model1_mean_indicator_", scenario), mean(HR_entry_model1[indicator == 1] - 1), envir = .GlobalEnv)
    assign(paste0("bias_entry_model2_mean_indicator_", scenario), mean(HR_entry_model2[indicator == 1] - 1), envir = .GlobalEnv)
    assign(paste0("bias_recruited_model2_mean_indicator_", scenario), mean(HR_recruited_model2[indicator == 1] - 1), envir = .GlobalEnv)
    
    assign(paste0("bias_treatment_standard_median_indicator_", scenario), median(HR_treatment_standard[indicator == 1] - HR), envir = .GlobalEnv)
    assign(paste0("bias_treatment_model1_median_indicator_", scenario), median(HR_treatment_model1[indicator == 1] - HR), envir = .GlobalEnv)
    assign(paste0("bias_treatment_model2_median_indicator_", scenario), median(HR_treatment_model2[indicator == 1] - HR), envir = .GlobalEnv)
    assign(paste0("bias_entry_model1_median_indicator_", scenario), median(HR_entry_model1[indicator == 1] - 1), envir = .GlobalEnv)
    assign(paste0("bias_entry_model2_median_indicator_", scenario), median(HR_entry_model2[indicator == 1] - 1), envir = .GlobalEnv)
    assign(paste0("bias_recruited_model2_median_indicator_", scenario), median(HR_recruited_model2[indicator == 1] - 1), envir = .GlobalEnv)
    
    assign(paste0("RMSE_treatment_standard_indicator_", scenario), sqrt(mean((HR_treatment_standard[indicator == 1] - HR)^2)), envir = .GlobalEnv)
    assign(paste0("RMSE_treatment_model1_indicator_", scenario), sqrt(mean((HR_treatment_model1[indicator == 1] - HR)^2)), envir = .GlobalEnv)
    assign(paste0("RMSE_treatment_model2_indicator_", scenario), sqrt(mean((HR_treatment_model2[indicator == 1] - HR)^2)), envir = .GlobalEnv)
    assign(paste0("RMSE_entry_model1_indicator_", scenario), sqrt(mean((HR_entry_model1[indicator == 1] - 1)^2)), envir = .GlobalEnv)
    assign(paste0("RMSE_entry_model2_indicator_", scenario), sqrt(mean((HR_entry_model2[indicator == 1] - 1)^2)), envir = .GlobalEnv)
    assign(paste0("RMSE_recruited_model2_indicator_", scenario), sqrt(mean((HR_recruited_model2[indicator == 1] - 1)^2)), envir = .GlobalEnv)
  }
  
}

# output results ###############################################################

# output results summarizing the bias of the estimated treatment hazard ratio
bias_treatment <- function(){
  
  data.frame(
    Distribution = c(rep("Exponential",15), rep("Weibull",15)),
    n = rep(c(rep(50,9), rep(26,6)),2),
    m = rep(c(rep(25,3), rep(10,6), rep(13,6)),2),
    Model = rep(c(c("Standard Model", "Model 1", "Model 2"), rep(rep(c("Standard Model", "Model 1", "Model 2"), each = 2), 2)), 2),
    mean_bias = c(bias_treatment_standard_mean_1, bias_treatment_model1_mean_1, bias_treatment_model2_mean_1,
                  bias_treatment_standard_mean_2, bias_treatment_standard_mean_indicator_2, 
                  bias_treatment_model1_mean_2, bias_treatment_model1_mean_indicator_2, 
                  bias_treatment_model2_mean_2, bias_treatment_model2_mean_indicator_2,
                  bias_treatment_standard_mean_3, bias_treatment_standard_mean_indicator_3, 
                  bias_treatment_model1_mean_3, bias_treatment_model1_mean_indicator_3, 
                  bias_treatment_model2_mean_3, bias_treatment_model2_mean_indicator_3,
                  bias_treatment_standard_mean_4, bias_treatment_model1_mean_4, bias_treatment_model2_mean_4,
                  bias_treatment_standard_mean_5, bias_treatment_standard_mean_indicator_5, 
                  bias_treatment_model1_mean_5, bias_treatment_model1_mean_indicator_5, 
                  bias_treatment_model2_mean_5, bias_treatment_model2_mean_indicator_5,
                  bias_treatment_standard_mean_6, bias_treatment_standard_mean_indicator_6, 
                  bias_treatment_model1_mean_6, bias_treatment_model1_mean_indicator_6, 
                  bias_treatment_model2_mean_6, bias_treatment_model2_mean_indicator_6),
    median_bias = c(bias_treatment_standard_median_1, bias_treatment_model1_median_1, bias_treatment_model2_median_1,
                    bias_treatment_standard_median_2, bias_treatment_standard_median_indicator_2, 
                    bias_treatment_model1_median_2, bias_treatment_model1_median_indicator_2, 
                    bias_treatment_model2_median_2, bias_treatment_model2_median_indicator_2,
                    bias_treatment_standard_median_3, bias_treatment_standard_median_indicator_3, 
                    bias_treatment_model1_median_3, bias_treatment_model1_median_indicator_3, 
                    bias_treatment_model2_median_3, bias_treatment_model2_median_indicator_3,
                    bias_treatment_standard_median_4, bias_treatment_model1_median_4, bias_treatment_model2_median_4,
                    bias_treatment_standard_median_5, bias_treatment_standard_median_indicator_5, 
                    bias_treatment_model1_median_5, bias_treatment_model1_median_indicator_5, 
                    bias_treatment_model2_median_5, bias_treatment_model2_median_indicator_5,
                    bias_treatment_standard_median_6, bias_treatment_standard_median_indicator_6, 
                    bias_treatment_model1_median_6, bias_treatment_model1_median_indicator_6, 
                    bias_treatment_model2_median_6, bias_treatment_model2_median_indicator_6),
    rmse = c(RMSE_treatment_standard_1, RMSE_treatment_model1_1, RMSE_treatment_model2_1,
             RMSE_treatment_standard_2, RMSE_treatment_standard_indicator_2,
             RMSE_treatment_model1_2, RMSE_treatment_model1_indicator_2,
             RMSE_treatment_model2_2, RMSE_treatment_model2_indicator_2,
             RMSE_treatment_standard_3, RMSE_treatment_standard_indicator_3,
             RMSE_treatment_model1_3, RMSE_treatment_model1_indicator_3,
             RMSE_treatment_model2_3, RMSE_treatment_model2_indicator_3,
             RMSE_treatment_standard_4, RMSE_treatment_model1_4, RMSE_treatment_model2_4,
             RMSE_treatment_standard_5, RMSE_treatment_standard_indicator_5,
             RMSE_treatment_model1_5, RMSE_treatment_model1_indicator_5,
             RMSE_treatment_model2_5, RMSE_treatment_model2_indicator_5,
             RMSE_treatment_standard_6, RMSE_treatment_standard_indicator_6,
             RMSE_treatment_model1_6, RMSE_treatment_model1_indicator_6,
             RMSE_treatment_model2_6, RMSE_treatment_model2_indicator_6),
    coverage = c(coverage_treatment_standard_1, coverage_treatment_model1_1, coverage_treatment_model2_1,
                 coverage_treatment_standard_2, coverage_treatment_standard_indicator_2,
                 coverage_treatment_model1_2, coverage_treatment_model1_indicator_2,
                 coverage_treatment_model2_2, coverage_treatment_model2_indicator_2,
                 coverage_treatment_standard_3, coverage_treatment_standard_indicator_3,
                 coverage_treatment_model1_3, coverage_treatment_model1_indicator_3,
                 coverage_treatment_model2_3, coverage_treatment_model2_indicator_3,
                 coverage_treatment_standard_4, coverage_treatment_model1_4, coverage_treatment_model2_4,
                 coverage_treatment_standard_5, coverage_treatment_standard_indicator_5,
                 coverage_treatment_model1_5, coverage_treatment_model1_indicator_5,
                 coverage_treatment_model2_5, coverage_treatment_model2_indicator_5,
                 coverage_treatment_standard_6, coverage_treatment_standard_indicator_6,
                 coverage_treatment_model1_6, coverage_treatment_model1_indicator_6,
                 coverage_treatment_model2_6, coverage_treatment_model2_indicator_6),
    Excluded = c(rep(0, 3), 
                 rep(c(0, 100000 - sum(indicator_2)), 3),
                 rep(c(0, 100000 - sum(indicator_3)), 3),
                 rep(0, 3),
                 rep(c(0, 100000 - sum(indicator_5)), 3),
                 rep(c(0, 100000 - sum(indicator_6)), 3))
  )
  
}

# output results summarizing the bias of the estimated hazard ratios for the additional covariates
bias_additional <- function(){
  
  data.frame(
    Distribution = c(rep("Exponential",15), rep("Weibull",15)),
    n = rep(c(rep(50,9), rep(26,6)),2),
    m = rep(c(rep(25,3), rep(10,6), rep(13,6)),2),
    Model = rep(c(c("Model 1", "Model 2", "Model 2"), rep(c(rep("Model 1", 2), rep("Model 2", 4)), 2)), 2),
    covariate = rep(c(c("Entry time", "Entry time", "Recruited subjects"), rep(c(rep("Entry time", 4), rep("Recruited subjects", 2)), 2)), 2),
    mean_bias = c(bias_entry_model1_mean_1, bias_entry_model2_mean_1, bias_recruited_model2_mean_1,
                  bias_entry_model1_mean_2, bias_entry_model1_mean_indicator_2, 
                  bias_entry_model2_mean_2, bias_entry_model2_mean_indicator_2, 
                  bias_recruited_model2_mean_2, bias_recruited_model2_mean_indicator_2,
                  bias_entry_model1_mean_3, bias_entry_model1_mean_indicator_3, 
                  bias_entry_model2_mean_3, bias_entry_model2_mean_indicator_3, 
                  bias_recruited_model2_mean_3, bias_recruited_model2_mean_indicator_3,
                  bias_entry_model1_mean_4, bias_entry_model2_mean_4, bias_recruited_model2_mean_4,
                  bias_entry_model1_mean_5, bias_entry_model1_mean_indicator_5, 
                  bias_entry_model2_mean_5, bias_entry_model2_mean_indicator_5, 
                  bias_recruited_model2_mean_5, bias_recruited_model2_mean_indicator_5,
                  bias_entry_model1_mean_6, bias_entry_model1_mean_indicator_6, 
                  bias_entry_model2_mean_6, bias_entry_model2_mean_indicator_6, 
                  bias_recruited_model2_mean_6, bias_recruited_model2_mean_indicator_6),
    median_bias = c(bias_entry_model1_median_1, bias_entry_model2_median_1, bias_recruited_model2_median_1,
                    bias_entry_model1_median_2, bias_entry_model1_median_indicator_2, 
                    bias_entry_model2_median_2, bias_entry_model2_median_indicator_2, 
                    bias_recruited_model2_median_2, bias_recruited_model2_median_indicator_2,
                    bias_entry_model1_median_3, bias_entry_model1_median_indicator_3, 
                    bias_entry_model2_median_3, bias_entry_model2_median_indicator_3, 
                    bias_recruited_model2_median_3, bias_recruited_model2_median_indicator_3,
                    bias_entry_model1_median_4, bias_entry_model2_median_4, bias_recruited_model2_median_4,
                    bias_entry_model1_median_5, bias_entry_model1_median_indicator_5, 
                    bias_entry_model2_median_5, bias_entry_model2_median_indicator_5, 
                    bias_recruited_model2_median_5, bias_recruited_model2_median_indicator_5,
                    bias_entry_model1_median_6, bias_entry_model1_median_indicator_6, 
                    bias_entry_model2_median_6, bias_entry_model2_median_indicator_6, 
                    bias_recruited_model2_median_6, bias_recruited_model2_median_indicator_6),
    rmse = c(RMSE_entry_model1_1, RMSE_entry_model2_1, RMSE_recruited_model2_1,
             RMSE_entry_model1_2, RMSE_entry_model1_indicator_2, 
             RMSE_entry_model2_2, RMSE_entry_model2_indicator_2, 
             RMSE_recruited_model2_2, RMSE_recruited_model2_indicator_2,
             RMSE_entry_model1_3, RMSE_entry_model1_indicator_3, 
             RMSE_entry_model2_3, RMSE_entry_model2_indicator_3, 
             RMSE_recruited_model2_3, RMSE_recruited_model2_indicator_3,
             RMSE_entry_model1_4, RMSE_entry_model2_4, RMSE_recruited_model2_4,
             RMSE_entry_model1_5, RMSE_entry_model1_indicator_5, 
             RMSE_entry_model2_5, RMSE_entry_model2_indicator_5, 
             RMSE_recruited_model2_5, RMSE_recruited_model2_indicator_5,
             RMSE_entry_model1_6, RMSE_entry_model1_indicator_6, 
             RMSE_entry_model2_6, RMSE_entry_model2_indicator_6, 
             RMSE_recruited_model2_6, RMSE_recruited_model2_indicator_6),
    coverage = c(coverage_entry_model1_1, coverage_entry_model2_1, coverage_recruited_model2_1,
                 coverage_entry_model1_2, coverage_entry_model1_indicator_2, 
                 coverage_entry_model2_2, coverage_entry_model2_indicator_2, 
                 coverage_recruited_model2_2, coverage_recruited_model2_indicator_2,
                 coverage_entry_model1_3, coverage_entry_model1_indicator_3, 
                 coverage_entry_model2_3, coverage_entry_model2_indicator_3, 
                 coverage_recruited_model2_3, coverage_recruited_model2_indicator_3,
                 coverage_entry_model1_4, coverage_entry_model2_4, coverage_recruited_model2_4,
                 coverage_entry_model1_5, coverage_entry_model1_indicator_5, 
                 coverage_entry_model2_5, coverage_entry_model2_indicator_5, 
                 coverage_recruited_model2_5, coverage_recruited_model2_indicator_5,
                 coverage_entry_model1_6, coverage_entry_model1_indicator_6, 
                 coverage_entry_model2_6, coverage_entry_model2_indicator_6, 
                 coverage_recruited_model2_6, coverage_recruited_model2_indicator_6),
    Excluded = c(rep(0, 3), 
                 rep(c(0, 100000 - sum(indicator_2)), 3),
                 rep(c(0, 100000 - sum(indicator_3)), 3),
                 rep(0, 3),
                 rep(c(0, 100000 - sum(indicator_5)), 3),
                 rep(c(0, 100000 - sum(indicator_6)), 3))
  )
  
}


# perform simulations ##########################################################

library(parallel)
library(pbapply)

## HR = 1 ####
simulate(seed = 8151188, scenario = 1, HR = 1)
simulate(seed = 9259582, scenario = 2, HR = 1)
simulate(seed = 5940754, scenario = 3, HR = 1)
simulate(seed = 3101981, scenario = 4, HR = 1)
simulate(seed = 9332161, scenario = 5, HR = 1)
simulate(seed = 7647627, scenario = 6, HR = 1)

bias_treatment()
bias_additional()

# note: this overwrites the variables in the workspace;
# save results using 'save.image()'
## HR = 0.8 ####
simulate(seed = 8151188, scenario = 1, HR = 0.8)
simulate(seed = 9259582, scenario = 2, HR = 0.8)
simulate(seed = 5940754, scenario = 3, HR = 0.8)
simulate(seed = 3101981, scenario = 4, HR = 0.8)
simulate(seed = 9332161, scenario = 5, HR = 0.8)
simulate(seed = 7647627, scenario = 6, HR = 0.8)

bias_treatment()
bias_additional()

# note: this overwrites the variables in the workspace;
# save results using 'save.image()'
## HR = 1.25 ####
simulate(seed = 8151188, scenario = 1, HR = 1.25)
simulate(seed = 9259582, scenario = 2, HR = 1.25)
simulate(seed = 5940754, scenario = 3, HR = 1.25)
simulate(seed = 3101981, scenario = 4, HR = 1.25)
simulate(seed = 9332161, scenario = 5, HR = 1.25)
simulate(seed = 7647627, scenario = 6, HR = 1.25)

bias_treatment()
bias_additional()
