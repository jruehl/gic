# This script includes the functions that produce the results of the first simulation study.

# generate data ################################################################

# generate event-driven data with exponentially distributed survival times
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

# generate event-driven data with Weibull distributed survival times
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
    event = entry + c(rweibull(n/2, shape, scale), rweibull(n/2, shape, scale/HR^(1/shape))),
    censoring = ifelse(event <= sort(event)[m], event, sort(event)[m]),
    study_time = censoring - entry,
    status = ifelse(event <= sort(event)[m], 1, 0),
    recruited_at_entry = rank(entry) - 1
  )
  
}

# generate randomly censored data with exponentially distributed survival times
#   (for comparison)
data_exponential_randomCensoring <- function(n, m, HR, lambda){
  
  # input:
  ## n: (for comparibility to event-driven scenarios)
  ## m: number of events to observe
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
  
  data <- data.frame()
  # ensure that m events are observed 
  #   (or else, scenarios with few observed events will lead to very extreme results that affect the mean outcomes)
  while(sum(data$status) < m){
    data <- rbind(
      data, 
      tibble(
        treatment = sample(c("Control", "Treatment"), 1),
        entry = runif(1, 0, qexp(m/n, lambda)),
        event = entry + ifelse(treatment == "Control", 
                               rexp(1, lambda), 
                               rexp(1, lambda*HR)),
        censoring = entry + rexp(1, ((1+HR)*(n-2*m)*lambda + 
                                       sqrt(((1+HR)*(n-2*m)*lambda)^2 + 
                                              16*HR*m*(n-m)*lambda^2))/(4*m)),
        study_time = min(event, censoring) - entry,
        status = ifelse(event <= censoring, 1, 0)
      )
    )
  }
  data$recruited_at_entry = rank(data$entry) - 1
  
  return(data)
  
}

# generate data and derive relevant parameters
prepare_data <- function(i, scenario, n, m, HR, lambda, shape, scale){
  
  # input:
  ## i: parameter for parallel processing
  ## scenario: scenario number (1 to 8)
  ## n: number of subjects
  ## m: number of events to observe before censoring is imposed
  ## HR: hazard ratio (treatment vs. control group)
  ## lambda: distribution parameter (in the control group) (exponential scenario)
  ## shape: distribution parameter (Weibull scenario)
  ## scale: distribution parameter (in the control group) (Weibull scenario)
  
  # output:
  ## indicator: indicator implying whether there are at least two observed events 
  ##   in each treatment group
  ## logHR_[model]_[parameter] 
  ##  (model = standard/model1/model2, parameter = treatment/entry/recruited):
  ##  estimated log-hazard ratio for the specified parameter in the specified group
  ## CI_[model]_[parameter]
  ##   (model = standard/model1/model2, parameter = treatment/entry/recruited):
  ##   estimated confidence interval for the log-hazard ratio
  ## breslow_[model] (model = standard/model1/model2): 
  ##   Breslow estimator of the cumulative baseline hazard in the specified model
  
  if(scenario %in% 1:5)
    data <- data_exponential(n, m, HR, lambda)
  else if(scenario %in% 6:10)
    data <- data_weibull(n, m, HR, shape, scale)
  else if(scenario %in% 11:13)
    data <- data_exponential_randomCensoring(n, m, HR, lambda)
  
  indicator <- sum(data$status == 1 & data$treatment == "Treatment") %in% 2:(m-2)
  
  models <- c("standard", "model1", "model2")
  for(model in models){
    
    vars <- switch(model, 
                   standard = c("treatment"),
                   model1 = c("treatment", "entry"),
                   model2 = c("treatment", "entry", "recruited"))
    modelvars <- paste(ifelse(vars == "recruited", "recruited_at_entry", vars), 
                       collapse = " + ")
    for(var in vars){
      
      # retrieve log-HR estimates
      assign(paste0("logHR_", model, "_", var), 
             eval(parse(text = paste0("summary(coxph(Surv(study_time, status) ~ ", 
                                      modelvars, ", data))$coefficients[", 
                                      switch(var, 
                                             treatment = ifelse(model == "standard", "", "1,"), 
                                             entry = "2,", 
                                             recruited = "3,"), 
                                      "1]"))))
      # retrieve confidence intervals
      assign(paste0("CI_", model, "_", var), 
             eval(parse(text = paste0("confint(coxph(Surv(study_time, status) ~ ", 
                                      modelvars, ", data))[", 
                                      switch(var, 
                                             treatment = "1", 
                                             entry = "2", 
                                             recruited = "3"), 
                                      ",]"))))
      
    }
    
    # retrieve breslow estimates
    assign(paste0("breslow_", model), 
           eval(parse(text = paste0("basehaz(coxph(Surv(study_time, status) ~ ", 
                                    modelvars, ", data), centered = FALSE)"))))
    
  }
  
  eval(parse(text = paste0("return(list(indicator, ",
                           paste0("logHR_", models[c(1,2,3,2,3,3)], "_", 
                                  rep(c("treatment", "entry", "recruited"), c(3,2,1)),
                                  sep = ", ", collapse = ""),
                           paste0("CI_", models[c(1,2,3,2,3,3)], "_", 
                                  rep(c("treatment", "entry", "recruited"), c(3,2,1)),
                                  sep = ", ", collapse = ""),
                           paste0("breslow_", models, collapse = ", "), 
                           "))")))
  
}

# run simulations ##############################################################

# run simulations and save summarized results within global variables:
## indicator_HR[HR]_scenario[scenario] (scenario = 1 - 10): 
##   indicator implying in which iterations there are at least two observed events in each treatment group 
## meanlogHR_HR[HR]_scenario[scenario]_[model]_[parameter]_(indicator)
##   (scenario = 1 - 10, model = standard/model1/model2, parameter = treatment/entry/recruited):
##   mean of the estimated log-hazard ratios for the specified parameter in the specified model and scenario 
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## medianlogHR_HR[HR]_scenario[scenario]_[model]_[parameter]_(indicator)
##   (scenario = 1 - 10, model = standard/model1/model2, parameter = treatment/entry/recruited):
##   median of the estimated log-hazard ratios for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## medianCIwidth_HR[HR]_scenario[scenario]_[model]_[parameter]_(indicator)
##   (scenario = 1 - 10, model = standard/model1/model2, parameter = treatment/entry/recruited):
##   median with of the confidence intervals of the estimated log-hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## empiricalCI_HR[HR]_scenario[scenario]_[model]_[parameter]_(indicator)
##   (scenario = 1 - 10, model = standard/model1/model2, parameter = treatment/entry/recruited):
##   empirical confidence interval of the estimated log-hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## coverage_HR[HR]_scenario[scenario]_[model]_[parameter]_(indicator)
##   (scenario = 1 - 10, model = standard/model1/model2, parameter = treatment/entry/recruited):
##   coverage of the confidence intervals of the estimated log-hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## meanbias_HR[HR]_scenario[scenario]_[model]_[parameter]_(indicator)
##   (scenario = 1 - 10, model = standard/model1/model2, parameter = treatment/entry/recruited):
##   mean bias of the estimated log-hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## medianbias_HR[HR]_scenario[scenario]_[model]_[parameter]_(indicator)
##   (scenario = 1 - 10, model = standard/model1/model2, parameter = treatment/entry/recruited):
##   median bias of the estimated log-hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## RMSE_HR[HR]_scenario[scenario]_[model]_[parameter]_(indicator)
##   (scenario = 1 - 10, model = standard/model1/model2, parameter = treatment/entry/recruited):
##   root mean squared error of the estimated log-hazard ratio for the specified parameter in the specified model and scenario
##   (indicator: restricted to iterations with at least two observed events in each treatment group)
## breslow_HR[HR]_scenario[scenario]_[model] (scenario = 1 - 10, model = standard/model1/model2):
##   Breslow estimator of the cumulative baseline hazard in the specified model and scenario
simulate <- function(seed, scenario, HR, lambda = 1, shape = 0.5, scale = 1, 
                     iter = 100000, cores = detectCores() - 1){
    
  # input: 
  ## seed: seed to use for simulations
  ## scenario: scenario number (1 to 10)
  ## HR: hazard ratio (treatment vs. control group)
  ## lambda: distribution parameter (in the control group) (exponential scenario)
  ## shape: distribution parameter (Weibull scenario)
  ## scale: distribution parameter (in the control group) (Weibull scenario)
  ## iter: number of iterations
  ## cores: number of cores for parallel computations 
  ##   (use 8 cores to reproduce the results from the manuscript)
  
  n <- switch(scenario,
              600,300,50,50,26,600,300,50,50,26,50,50,26)
  m <- switch(scenario,
              300,150,25,10,13,300,150,25,10,13,25,10,13)
  
  cl <- makeCluster(cores, type="SOCK")
  clusterEvalQ(cl, library(tibble))
  clusterEvalQ(cl, library(survival))
  clusterExport(cl, "data_exponential")
  clusterExport(cl, "data_weibull")
  clusterExport(cl, "data_exponential_randomCensoring")
  clusterSetRNGStream(cl, iseed = seed)
  sim <- pblapply(X = 1:iter, FUN = prepare_data, scenario, n, m, HR, lambda, shape, scale, cl = cl)
  stopCluster(cl)
  
  
  indicator <- as.numeric(unlist(lapply(sim, `[[`, 1)))
  assign(paste0("indicator_HR", HR, "_scenario", scenario), indicator, envir = .GlobalEnv)
  
  for(model in c("standard", "model1", "model2")){
    
    vars <- switch(model, 
                   standard = c("treatment"),
                   model1 = c("treatment", "entry"),
                   model2 = c("treatment", "entry", "recruited"))
    for(var in vars){
      
      # retrieve log-HR estimates
      assign(paste0("logHR_", model, "_", var), 
             eval(parse(text = paste0("unlist(lapply(sim, `[[`, ",
                                      switch(model,
                                             standard = 0,
                                             model1 = 1,
                                             model2 = 2) + 
                                        switch(var,
                                               treatment = 2,
                                               entry = 4,
                                               recruited = 5), 
                                      "))"))))
      
      # retrieve confidence intervals
      assign(paste0("CI_", model, "_", var, "_lower"), 
             eval(parse(text = paste0("unlist(lapply(lapply(sim, `[[`, ",
                                      switch(model,
                                             standard = 6,
                                             model1 = 7,
                                             model2 = 8) + 
                                        switch(var,
                                               treatment = 2,
                                               entry = 4,
                                               recruited = 5), 
                                      "), `[[`, 1))"))))
      assign(paste0("CI_", model, "_", var, "_upper"), 
             eval(parse(text = paste0("unlist(lapply(lapply(sim, `[[`, ",
                                      switch(model,
                                             standard = 6,
                                             model1 = 7,
                                             model2 = 8) + 
                                        switch(var,
                                               treatment = 2,
                                               entry = 4,
                                               recruited = 5), 
                                      "), `[[`, 2))"))))
      
    }
    
    assign(paste0("breslow_HR", HR, "_scenario", scenario, "_", model), 
           lapply(sim, `[[`, switch(model,
                                    standard = 14,
                                    model1 = 15,
                                    model2 = 16)), envir = .GlobalEnv)
    
  }
  
  
  for(model in c("standard", "model1", "model2")){
    
    vars <- switch(model, 
                   standard = c("treatment"),
                   model1 = c("treatment", "entry"),
                   model2 = c("treatment", "entry", "recruited"))
    for(var in vars){
      
      assign(paste0("meanlogHR_HR", HR, "_scenario", scenario, "_", model, "_", var),
             eval(parse(text = paste0("mean(logHR_", model, "_", var, ")"))), 
             envir = .GlobalEnv)
      assign(paste0("medianlogHR_HR", HR, "_scenario", scenario, "_", model, "_", var),
             eval(parse(text = paste0("median(logHR_", model, "_", var, ")"))), 
             envir = .GlobalEnv)
      assign(paste0("medianCIwidth_HR", HR, "_scenario", scenario, "_", model, "_", var),
             eval(parse(text = paste0("median(CI_", model, "_", var, "_upper - CI_", model, "_", var, "_lower)"))), 
             envir = .GlobalEnv)
      assign(paste0("empiricalCI_HR", HR, "_scenario", scenario, "_", model, "_", var),
             eval(parse(text = paste0("quantile(logHR_", model, "_", var, ", c(0.025, 0.975))"))), 
             envir = .GlobalEnv)
      assign(paste0("coverage_HR", HR, "_scenario", scenario, "_", model, "_", var),
             eval(parse(text = paste0("mean(CI_", model, "_", var, "_lower <= ", ifelse(var == "treatment", log(HR), 0), 
                                      " & ", ifelse(var == "treatment", log(HR), 0), " <= CI_", model, "_", var, "_upper)"))), 
             envir = .GlobalEnv)
      assign(paste0("meanbias_HR", HR, "_scenario", scenario, "_", model, "_", var),
             eval(parse(text = paste0("mean(logHR_", model, "_", var, " - ", ifelse(var == "treatment", log(HR), 0), ")"))), 
             envir = .GlobalEnv)
      assign(paste0("medianbias_HR", HR, "_scenario", scenario, "_", model, "_", var),
             eval(parse(text = paste0("median(logHR_", model, "_", var, " - ", ifelse(var == "treatment", log(HR), 0), ")"))), 
             envir = .GlobalEnv)
      assign(paste0("RMSE_HR", HR, "_scenario", scenario, "_", model, "_", var),
             eval(parse(text = paste0("sqrt(mean((logHR_", model, "_", var, " - ", ifelse(var == "treatment", log(HR), 0), ")^2))"))), 
             envir = .GlobalEnv)
      
      assign(paste0("MC_log_HR", HR, "_scenario", scenario, "_", model, "_", var),
             eval(parse(text = paste0("round(sqrt(var(logHR_", model, "_", var, ") / iter), 5)"))),
             envir = .GlobalEnv)
      
      
      if(!all(indicator == 1)){
        
        assign(paste0("meanlogHR_HR", HR, "_scenario", scenario, "_", model, "_", var, "_indicator"),
               eval(parse(text = paste0("mean(logHR_", model, "_", var, "[indicator == 1])"))), 
               envir = .GlobalEnv)
        assign(paste0("medianlogHR_HR", HR, "_scenario", scenario, "_", model, "_", var, "_indicator"),
               eval(parse(text = paste0("median(logHR_", model, "_", var, "[indicator == 1])"))), 
               envir = .GlobalEnv)
        assign(paste0("medianCIwidth_HR", HR, "_scenario", scenario, "_", model, "_", var, "_indicator"),
               eval(parse(text = paste0("median(CI_", model, "_", var, "_upper[indicator == 1] - CI_", model, "_", var, "_lower[indicator == 1])"))), 
               envir = .GlobalEnv)
        assign(paste0("empiricalCI_HR", HR, "_scenario", scenario, "_", model, "_", var, "_indicator"),
               eval(parse(text = paste0("quantile(logHR_", model, "_", var, "[indicator == 1], c(0.025, 0.975))"))), 
               envir = .GlobalEnv)
        assign(paste0("coverage_HR", HR, "_scenario", scenario, "_", model, "_", var, "_indicator"),
               eval(parse(text = paste0("mean(CI_", model, "_", var, "_lower[indicator == 1] <=", ifelse(var == "treatment", log(HR), 0), 
                                        " & ", ifelse(var == "treatment", log(HR), 0), "<= CI_", model, "_", var, "_upper[indicator == 1])"))), 
               envir = .GlobalEnv)
        assign(paste0("meanbias_HR", HR, "_scenario", scenario, "_", model, "_", var, "_indicator"),
               eval(parse(text = paste0("mean(logHR_", model, "_", var, "[indicator == 1] - ", ifelse(var == "treatment", log(HR), 0), ")"))), 
               envir = .GlobalEnv)
        assign(paste0("medianbias_HR", HR, "_scenario", scenario, "_", model, "_", var, "_indicator"),
               eval(parse(text = paste0("median(logHR_", model, "_", var, "[indicator == 1] - ", ifelse(var == "treatment", log(HR), 0), ")"))), 
               envir = .GlobalEnv)
        assign(paste0("RMSE_HR", HR, "_scenario", scenario, "_", model, "_", var, "_indicator"),
               eval(parse(text = paste0("sqrt(mean((logHR_", model, "_", var, "[indicator == 1] - ", ifelse(var == "treatment", log(HR), 0), ")^2))"))), 
               envir = .GlobalEnv)
        
      }
      
    }
    
  }
  
}

# output results ###############################################################

# output results summarizing the bias of the estimated treatment log-hazard ratio
bias_treatment <- function(HR, iter = 100000){
  
  ## HR: hazard ratio (treatment vs. control group)
  ## iter: number of iterations
  
  models <- c("standard", "model1", "model2")
  data.frame(
    Distribution = rep(c("Exponential", "Weibull"), each=21),
    n = rep(c(rep(c(600,300), each=3), rep(50,9), rep(26,6)), 2),
    m = rep(c(rep(c(300,150,25), each=3), rep(c(10,13), each=6)), 2),
    Model = rep(c("Standard Model", "Model 1", "Model 2")[c(rep(c(1,2,3), 3), rep(rep(c(1,2,3), each=2), 2))], 2),
    mean_bias = eval(parse(text = 
                             paste0("c(",
                                    paste0("meanbias_HR", HR,
                                           "_scenario", c(rep(c(1,2,3), each=3), rep(c(4,5), each=6), rep(c(6,7,8), each=3), rep(c(9,10), each=6)),
                                           "_", rep(c(rep(models, 3), rep(rep(models, each=2), 2)), 2),
                                           "_treatment",
                                           rep(c(rep("", 9), rep(c("", "_indicator"), 6)), 2),
                                           collapse = ", "),
                                    ")"))),
    median_bias = eval(parse(text = 
                               paste0("c(", 
                                      paste0("medianbias_HR", HR,
                                             "_scenario", c(rep(c(1,2,3), each=3), rep(c(4,5), each=6), rep(c(6,7,8), each=3), rep(c(9,10), each=6)),
                                             "_", rep(c(rep(models, 3), rep(rep(models, each=2), 2)), 2),
                                             "_treatment",
                                             rep(c(rep("", 9), rep(c("", "_indicator"), 6)), 2),
                                             collapse = ", "),
                                      ")"))),
    rmse = eval(parse(text = 
                        paste0("c(", 
                               paste0("RMSE_HR", HR,
                                      "_scenario", c(rep(c(1,2,3), each=3), rep(c(4,5), each=6), rep(c(6,7,8), each=3), rep(c(9,10), each=6)),
                                      "_", rep(c(rep(models, 3), rep(rep(models, each=2), 2)), 2),
                                      "_treatment",
                                      rep(c(rep("", 9), rep(c("", "_indicator"), 6)), 2),
                                      collapse = ", "),
                               ")"))),
    coverage = eval(parse(text = 
                            paste0("c(", 
                                   paste0("coverage_HR", HR,
                                          "_scenario", c(rep(c(1,2,3), each=3), rep(c(4,5), each=6), rep(c(6,7,8), each=3), rep(c(9,10), each=6)),
                                          "_", rep(c(rep(models, 3), rep(rep(models, each=2), 2)), 2),
                                          "_treatment",
                                          rep(c(rep("", 9), rep(c("", "_indicator"), 6)), 2),
                                          collapse = ", "),
                                   ")"))),
    Excluded = c(rep(0, 9), 
                 rep(c(0, iter - sum(eval(parse(text = paste0("indicator_HR", HR, "_scenario2"))))), 3),
                 rep(c(0, iter - sum(eval(parse(text = paste0("indicator_HR", HR, "_scenario3"))))), 3),
                 rep(0, 9),
                 rep(c(0, iter - sum(eval(parse(text = paste0("indicator_HR", HR, "_scenario5"))))), 3),
                 rep(c(0, iter - sum(eval(parse(text = paste0("indicator_HR", HR, "_scenario6"))))), 3))
  )
  
}

# output results summarizing the bias of the estimated log-hazard ratios for the additional covariates
bias_additional <- function(HR, iter = 100000){
  
  ## HR: hazard ratio (treatment vs. control group)
  ## iter: number of iterations
  
  models <- c("model1", "model2", "model2")
  vars <- c("entry", "entry", "recruited")
  data.frame(
    Distribution = rep(c("Exponential", "Weibull"), each=21),
    n = rep(c(rep(c(600,300), each=3), rep(50,9), rep(26,6)), 2),
    m = rep(c(rep(c(300,150,25), each=3), rep(c(10,13), each=6)), 2),
    Model = rep(c("Model 1", "Model 2", "Model 2")[c(rep(c(1,2,3), 3), rep(rep(c(1,2,3), each=2), 2))], 2),
    covariate = rep(c("Entry time", "Entry time", "Recruited subjects")[c(rep(c(1,2,3), 3), rep(rep(c(1,2,3), each=2), 2))], 2),
    mean_bias = eval(parse(text = 
                             paste0("c(",
                                    paste0("meanbias_HR", HR,
                                           "_scenario", c(rep(c(1,2,3), each=3), rep(c(4,5), each=6), rep(c(6,7,8), each=3), rep(c(9,10), each=6)),
                                           "_", rep(c(rep(models, 3), rep(rep(models, each=2), 2)), 2),
                                           "_", rep(c(rep(vars, 3), rep(rep(vars, each=2), 2)), 2), 
                                           rep(c(rep("", 9), rep(c("", "_indicator"), 6)), 2),
                                           collapse = ", "),
                                    ")"))),
    median_bias = eval(parse(text = 
                               paste0("c(",
                                      paste0("medianbias_HR", HR,
                                             "_scenario", c(rep(c(1,2,3), each=3), rep(c(4,5), each=6), rep(c(6,7,8), each=3), rep(c(9,10), each=6)),
                                             "_", rep(c(rep(models, 3), rep(rep(models, each=2), 2)), 2),
                                             "_", rep(c(rep(vars, 3), rep(rep(vars, each=2), 2)), 2), 
                                             rep(c(rep("", 9), rep(c("", "_indicator"), 6)), 2),
                                             collapse = ", "),
                                      ")"))),
    rmse = eval(parse(text = 
                        paste0("c(",
                               paste0("RMSE_HR", HR,
                                      "_scenario", c(rep(c(1,2,3), each=3), rep(c(4,5), each=6), rep(c(6,7,8), each=3), rep(c(9,10), each=6)),
                                      "_", rep(c(rep(models, 3), rep(rep(models, each=2), 2)), 2),
                                      "_", rep(c(rep(vars, 3), rep(rep(vars, each=2), 2)), 2), 
                                      rep(c(rep("", 9), rep(c("", "_indicator"), 6)), 2),
                                      collapse = ", "),
                               ")"))),
    coverage = eval(parse(text = 
                            paste0("c(",
                                   paste0("coverage_HR", HR,
                                          "_scenario", c(rep(c(1,2,3), each=3), rep(c(4,5), each=6), rep(c(6,7,8), each=3), rep(c(9,10), each=6)),
                                          "_", rep(c(rep(models, 3), rep(rep(models, each=2), 2)), 2),
                                          "_", rep(c(rep(vars, 3), rep(rep(vars, each=2), 2)), 2), 
                                          rep(c(rep("", 9), rep(c("", "_indicator"), 6)), 2),
                                          collapse = ", "),
                                   ")"))),
    Excluded = c(rep(0, 9), 
                 rep(c(0, iter - sum(eval(parse(text = paste0("indicator_HR", HR, "_scenario2"))))), 3),
                 rep(c(0, iter - sum(eval(parse(text = paste0("indicator_HR", HR, "_scenario3"))))), 3),
                 rep(0, 9),
                 rep(c(0, iter - sum(eval(parse(text = paste0("indicator_5"))))), 3),
                 rep(c(0, iter - sum(eval(parse(text = paste0("indicator_6"))))), 3))
  )
  
}


# execution ####################################################################

library(parallel)
library(pbapply)

## event-driven censoring: HR = 1 ####
simulate(seed = 1469597, scenario = 1, HR = 1, cores = 8)
simulate(seed = 9016046, scenario = 2, HR = 1, cores = 8)
simulate(seed = 8151188, scenario = 3, HR = 1, cores = 8)
simulate(seed = 9259582, scenario = 4, HR = 1, cores = 8)
simulate(seed = 5940754, scenario = 5, HR = 1, cores = 8)
simulate(seed = 5438504, scenario = 6, HR = 1, cores = 8)
simulate(seed = 4106576, scenario = 7, HR = 1, cores = 8)
simulate(seed = 3101981, scenario = 8, HR = 1, cores = 8)
simulate(seed = 9332161, scenario = 9, HR = 1, cores = 8)
simulate(seed = 7647627, scenario = 10, HR = 1, cores = 8)

bias_treatment(HR = 1)
bias_additional(HR = 1)

## event-driven censoring: HR = 0.8 ####
simulate(seed = 1469597, scenario = 1, HR = 0.8, cores = 8)
simulate(seed = 9016046, scenario = 2, HR = 0.8, cores = 8)
simulate(seed = 8151188, scenario = 3, HR = 0.8, cores = 8)
simulate(seed = 9259582, scenario = 4, HR = 0.8, cores = 8)
simulate(seed = 5940754, scenario = 5, HR = 0.8, cores = 8)
simulate(seed = 5438504, scenario = 6, HR = 0.8, cores = 8)
simulate(seed = 4106576, scenario = 7, HR = 0.8, cores = 8)
simulate(seed = 3101981, scenario = 8, HR = 0.8, cores = 8)
simulate(seed = 9332161, scenario = 9, HR = 0.8, cores = 8)
simulate(seed = 7647627, scenario = 10, HR = 0.8, cores = 8)

bias_treatment(HR = 0.8)
bias_additional(HR = 0.8)

## event-driven censoring: HR = 1.25 ####
simulate(seed = 1469597, scenario = 1, HR = 1.25, cores = 8)
simulate(seed = 9016046, scenario = 2, HR = 1.25, cores = 8)
simulate(seed = 8151188, scenario = 3, HR = 1.25, cores = 8)
simulate(seed = 9259582, scenario = 4, HR = 1.25, cores = 8)
simulate(seed = 5940754, scenario = 5, HR = 1.25, cores = 8)
simulate(seed = 5438504, scenario = 6, HR = 1.25, cores = 8)
simulate(seed = 4106576, scenario = 7, HR = 1.25, cores = 8)
simulate(seed = 3101981, scenario = 8, HR = 1.25, cores = 8)
simulate(seed = 9332161, scenario = 9, HR = 1.25, cores = 8)
simulate(seed = 7647627, scenario = 10, HR = 1.25, cores = 8)

bias_treatment(HR = 1.25)
bias_additional(HR = 1.25)

## random censoring (for comparison) ####
simulate(seed = 8151188, scenario = 11, HR = 1, cores = 8)
simulate(seed = 9259582, scenario = 12, HR = 1, cores = 8)
simulate(seed = 5940754, scenario = 13, HR = 1, cores = 8)
