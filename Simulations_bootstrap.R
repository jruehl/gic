# This script includes the functions that produce the results of the second simulation study.

# generate data ################################################################

# generate data from an illness-death model without recovery
generate_data <- function(n, lambda01, lambda02, lambda12, staggered, censoring, m){
  
  # input:
  ## n: number of subjects
  ## lambda01: distribution parameter for "0 -> 1" transitions
  ## lambda02: distribution parameter for "0 -> 2" transitions
  ## lambda12: distribution parameter for "1 -> 2" transitions
  ## staggered: logical value indicating whether entry should be staggered
  ## censoring: logical value indicating whether type II censoring should be applied
  ## m: number of events to observe before censoring is imposed
  
  # output: data set with (one row per transition and) the subsequent columns:
  ## id: subject id
  ## entry: (study) time of entry into the respective state
  ## exit: (study) time of exit from the respective state
  ## from: state from where the transition occurs
  ## to: state to which the transition occurs
  
  # generate "0 -> 1" & "0 -> 2" transitions
  data <- as.data.frame(tibble(
    id = 1:n,
    entry = sapply(1:n, FUN = function(i){ifelse(staggered == T, runif(1, 0, 60), 0)}),
    exit = entry + rexp(n, lambda01 + lambda02),
    from = 0, 
    to = rbinom(n, 1, lambda02/(lambda01+lambda02)) + 1
  ))
  # add "1 -> 2" transitions
  try({
    data <- rbind(
      data,
      data.frame(
        id = data[data$to == 1, "id"],
        entry = data[data$to == 1, "exit"],
        exit = data[data$to == 1, "exit"] + rexp(dim(data[data$to == 1,])[1], lambda12),
        from = 1,
        to = 2
      )
    )
  }, silent = TRUE)
  if(censoring){
    # implement type II censoring
    max_time <- sort(data[data$to == 2, "exit"])[m]
    data <- data[data$entry <= max_time,]
    data$to <- ifelse(data$exit > max_time, "cens", data$to)
    data$exit[data$to == "cens"] <- max_time
  }
  if(staggered){
    # transform entry and exit times to study time scale
    data <- merge(data, 
                  data[data$from == 0, c("id", "entry")],
                  by = "id")
    data$entry <- data$entry.x - data$entry.y
    data$exit <- data$exit - data$entry.y
    data <- data[,c("id", "entry", "exit", "from", "to")]
  }
  
  return(data)
  
}

# analyze data #################################################################

# functions to derive Nelson-Aalen/Aalen-Johansen estimator and compute 
#   non-parametric/wild bootstrap confidence intervals
source("Bootstrap_functions.R")

# generate data, apply both bootstrap approaches, and summarize results
bootstrap <- function(i, n, m, lambda01, lambda02, lambda12, staggered, times, 
                      NAE, AJE, BS_iter, NA_AJ_true){
  
  # input:
  ## i: parameter for parallel processing
  ## n: number of subjects
  ## m: number of events to observe before censoring is imposed
  ## lambda01: distribution parameter for "0 -> 1" transitions
  ## lambda02: distribution parameter for "0 -> 2" transitions
  ## lambda12: distribution parameter for "1 -> 2" transitions
  ## staggered: logical value indicating whether entry should be staggered
  ## times: time points at which the CIs should be evaluated
  ## NAE: logical value indicating whether CIs for Nelson-Aalen estimator should be computed
  ## AJE: logical value indicating whether CIs for Aalen-Johansen estimator should be computed
  ## BS_iter: number of bootstrap samples
  ## NA_AJ_true: true Nelson-Aalen/Aalen-Johansen estimators
  
  # output: vector with CI widths & coverages at given time points for both EBS & WBS
  
  ## generate data
  data <- generate_data(n, lambda01, lambda02, lambda12, staggered, censoring = TRUE, m)
  
  # retrieve NA/AJ info
  etm <- etm_short(data, 0, max(times), NAE, AJE)
  
  # compute CIs for NA/AJ estimators for "1 -> 2" transitions 
  #   by means of EBS
  CI_EBS <- EBS(data, etm, times, "1", "2", NAE, AJE, BS_iter)
  
  # compute CIs for NA/AJ estimators for "1 -> 2 transitions 
  #   by means of WBS
  CI_WBS <- WBS(data, etm, times, "1", "2", NAE, AJE, BS_iter)
  
  return(c(
    if(NAE){
      c(
        width_EBS_NA = CI_EBS[,2] - CI_EBS[,1], 
        coverage_EBS_NA = CI_EBS[,1] <= NA_AJ_true[1:length(times)] & 
          NA_AJ_true[1:length(times)] <= CI_EBS[,2], 
        width_WBS_NA = CI_WBS[,2] - CI_WBS[,1], 
        coverage_WBS_NA = CI_WBS[,1] <= NA_AJ_true[1:length(times)] & 
          NA_AJ_true[1:length(times)] <= CI_WBS[,2])
    }, 
    if(AJE){
      c(
        width_EBS_AJ = CI_EBS[,2*NAE+2] - CI_EBS[,2*NAE+1], 
        coverage_EBS_AJ = CI_EBS[,2*NAE+1] <= NA_AJ_true[(NAE*length(times)+1):length(NA_AJ_true)] & 
          NA_AJ_true[(NAE*length(times)+1):length(NA_AJ_true)] <= CI_EBS[,2*NAE+2], 
        width_WBS_AJ = CI_WBS[,2*NAE+2] - CI_WBS[,2*NAE+1], 
        coverage_WBS_AJ = CI_WBS[,2*NAE+1] <= NA_AJ_true[(NAE*length(times)+1):length(NA_AJ_true)] & 
          NA_AJ_true[(NAE*length(times)+1):length(NA_AJ_true)] <= CI_WBS[,2*NAE+2]
      )
    }
    ))
  
}

# run simulations ##############################################################

# run simulations and save summarized results
simulate <- function(seed, n, m, lambda01 = 0.01, lambda02 = 0.03, lambda12 = 0.1, 
                     staggered, times = c(16,18,20), NAE = TRUE, AJE = FALSE, 
                     BS_iter = 1000, cores = detectCores() - 1){
  
  # input:
  ## seed: seed to use for simulations
  ## n: number of subjects
  ## m: number of events to observe before censoring is imposed
  ## lambda01: distribution parameter for "0 -> 1" transitions
  ## lambda02: distribution parameter for "0 -> 2" transitions
  ## lambda12: distribution parameter for "1 -> 2" transitions
  ## staggered: logical value indicating whether entry should be staggered
  ## times: time points at which the CIs should be evaluated
  ## NAE: logical value indicating whether CIs for Nelson-Aalen estimator should be computed
  ## AJE: logical value indicating whether CIs for Aalen-Johansen estimator should be computed
  ## BS_iter: number of bootstrap samples
  ## cores: number of cores for parallel computations 
  ##   (use 8 cores to reproduce the results from the manuscript)
  
  # output: vector with mean CI widths & coverages at given time points for both EBS & WBS, 
  #   for Nelson-Aalen & Aalen-Johansen estimators
  
  set.seed(seed)
  
  # compute true NA/AJ estimators (based on uncensored data) for comparison
  NA_AJ_true <- rowMeans(sapply(1:10000, FUN = function(i){
    etm <- etm_short(generate_data(n, lambda01, lambda02, lambda12, 
                                 staggered, censoring = FALSE),
                     0, max(times), NAE, AJE)
    c(etm$NA_estimator[findInterval(times, vec = etm$time),"1 2"],
      etm$AJ_estimator["1","2", findInterval(times, vec = etm$time)])
  }))
  
  cl <- makeCluster(cores, type="SOCK")
  clusterEvalQ(cl, library(tibble))
  clusterEvalQ(cl, library(etm))
  clusterExport(cl, "generate_data")
  clusterExport(cl, "etm_short")
  clusterExport(cl, "EBS")
  clusterExport(cl, "WBS")
  clusterSetRNGStream(cl, iseed = seed)
  res <- pblapply(X = 1:1000, FUN = bootstrap, n, m, lambda01, lambda02, lambda12, 
                  staggered, times, NAE, AJE, BS_iter, NA_AJ_true, cl = cl)
  stopCluster(cl)
  
  result <- colMeans(matrix(unlist(res), ncol = sum(NAE,AJE)*4*length(times), byrow = TRUE), 
                     na.rm = TRUE)
  names(result) <- c(
    if(NAE){
      c(
        paste0("width EBS NA (t = ", times, ")"),
        paste0("coverage EBS NA (t = ", times, ")"),
        paste0("width WBS NA (t = ", times, ")"),
        paste0("coverage WBS NA (t = ", times, ")")
      )
    },
    if(AJE){
      c(
        paste0("width EBS AJ (t = ", times, ")"),
        paste0("coverage EBS AJ (t = ", times, ")"),
        paste0("width WBS AJ (t = ", times, ")"),
        paste0("coverage WBS AJ (t = ", times, ")")
      )
    }
  )
  
  return(result)
  
}


# execution ####################################################################

library(tibble)
library(etm)
library(parallel)
library(pbapply)

result1 <- simulate(seed = 8570144, n = 600, m = 300, staggered = FALSE, cores = 8)
result1_staggered <- simulate(seed = 8570144, n = 600, m = 300, staggered = TRUE, cores = 8)
result2 <- simulate(seed = 5346732, n = 400, m = 200, staggered = FALSE, cores = 8)
result2_staggered <- simulate(seed = 5346732, n = 400, m = 200, staggered = TRUE, cores = 8)
result3 <- simulate(seed = 7755153, n = 200, m = 100, staggered = FALSE, cores = 8)
result3_staggered <- simulate(seed = 7755153, n = 200, m = 100, staggered = TRUE, cores = 8)
result4 <- simulate(seed = 8466261, n = 100, m = 50, staggered = FALSE, cores = 8)
result4_staggered <- simulate(seed = 8466261, n = 100, m = 50, staggered = TRUE, cores = 8)
result5 <- simulate(seed = 4055037, n = 80, m = 40, staggered = FALSE, cores = 8)
result5_staggered <- simulate(seed = 4055037, n = 80, m = 40, staggered = TRUE, cores = 8)
result6 <- simulate(seed = 7416887, n = 50, m = 25, staggered = FALSE, cores = 8)
result6_staggered <- simulate(seed = 7416887, n = 50, m = 25, staggered = TRUE, cores = 8)
