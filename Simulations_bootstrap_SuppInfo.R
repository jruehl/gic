# This script includes the functions that produce the results of the second simulation study.

# generate data ################################################################

# generate data from an illness-death model without recovery
generate_data <- function(n, lambda01, lambda02, lambda12, staggered = FALSE, censoring, m){
  
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
  
  # generate "0 -> 1" and "0 -> 2" transitions
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

# shorter version of etm() including Nelson-Aalen estimator
etm_short <- function(data, s, t){
  
  # input:
  ## data: data to derive Nelson-Aalen/Aalen-Johansen estimator from
  ## s: starting time for Aalen-Johansen estimator
  ## t: end time for Aalen-Johansen estimator
  
  # output: list with the subsequent entries:
  ## time: event times at which the transition probabilities are computed
  ## NA_est: Nelson-Aalen estimator
  ## est: Aalen-Johansen estimator
  ## n.risk: number of individuals at risk
  ## n.event: number of transitions
  
  # recode states (including censoring)
  data$from <- as.integer(as.character(factor(data$from, 
                                              levels = c("0", "1"), 
                                              labels = c("1", "2"), 
                                              ordered = TRUE)))
  data$to <- as.integer(as.character(factor(data$to, 
                                            levels = c("cens", "0", "1", "2"), 
                                            labels = c("0", "1", "2", "3"), 
                                            ordered = TRUE)))
  # restrict to transition times of interest
  transition_times <- sort(unique(data$exit[data$to != 0]))
  transition_times <- transition_times[s < transition_times & transition_times <= t]
  time <- NULL
  NA_estimator <- NULL
  AJ_estimator <- array(diag(1, nrow = 3, ncol = 3), dim = c(3, 3, 1))
  n.risk <- NULL
  n.event <- NULL
  if(length(transition_times) > 0){
    # estimate transition probabilities
    trans <- .Call("gen_msm", 
                   transition_times, data$entry, data$exit, 
                   data$from, data$to, 3, 
                   matrix(0, nrow = length(transition_times), ncol = 3))
    # clean data 
    if(s == 0){
      # add entries for t = 0 to enable search with 'findInterval'
      time <- c(0, trans$time)
      n.risk <- rbind(c(length(unique(data$id)), 0), trans$n.risk[, c(1,2)])
      n.event <- array(c(rep(0, 9), trans$n.event), dim = c(3, 3, length(time)))
      AJ_estimator <- array(c(diag(1, nrow = 3, ncol = 3), trans$est), 
                            dim = c(3, 3, length(time)))
    }else{
      time <- trans$time
      n.risk <- trans$n.risk[, c(1,2), drop = FALSE]
      n.event <- trans$n.event
      AJ_estimator <- trans$est
    }
    colnames(n.risk) <- c("0", "1")
    dimnames(n.event) <- list(c("0", "1", "2"), c("0", "1", "2"), time)
    dimnames(AJ_estimator) <- list(c("0", "1", "2"), c("0", "1", "2"), time)
    # compute NA estimator
    NA_estimator <- cbind(cumsum(replace(n.event["0","1",] / n.risk[,"0"], 
                                         list = n.risk[,"0"] == 0, values = 0)),
                          cumsum(replace(n.event["0","2",] / n.risk[,"0"], 
                                         list = n.risk[,"0"] == 0, values = 0)),
                          cumsum(replace(n.event["1","2",] / n.risk[,"1"], 
                                         list = n.risk[,"1"] == 0, values = 0)))
    colnames(NA_estimator) <- c("0 1", "0 2", "1 2")
  }
  
  list(time = time,
       NA_estimator = NA_estimator, 
       AJ_estimator= AJ_estimator,
       n.risk = n.risk,
       n.event = n.event)
  
}

# apply bootstrap ##############################################################

# compute CIs for the Nelson-Aalen & Aalen-Johansen estimators for "1 -> 2" transitions 
#   using Efron's non-parametric bootstrap
EBS <- function(data_original, NA_AJ_original, times, BS_iter){
  
  # input:
  ## data_original: original data to base bootstrap on
  ## NA_AJ_original: Nelson-Aalen & Aalen-Johansen estimators for "1 -> 2" transitions 
  ##   based on original data
  ## times: time points at which the CIs should be evaluated
  ## BS_iter: number of bootstrap samples
  
  # output: matrix containing bootstrapped CI limits 
  #   (1st column: lower limit Nelson-Aalen estimator, 
  #    2nd column: upper limit Nelson-Aalen estimator,
  #    3rd column: lower limit Aalen-Johansen estimator,
  #    4th column: upper limit Aalen-Johansen estimator)
  
  # compute EBS NA & AJ estimators
  NA_AJ_EBS <- sapply(1:BS_iter, FUN = function(j){
    # sample observations
    sample <- sample(unique(data_original$id), 
                     size = length(unique(data_original$id)), replace = TRUE)
    sample_data <- merge(data.frame(id = sample), data_original, 
                         by = "id", all.x = TRUE, sort = TRUE)
    # determine new subject IDs to distinguish between observations
    sample_data[sample_data$from == 0, "id2"] <- 
      seq(1:dim(sample_data[sample_data$from == 0,])[1])
    sample_data$id2[is.na(sample_data$id2)] <- unlist(tapply(
      sample_data$id2, INDEX = sample_data$id, 
      FUN = function(x){if(any(is.na(x))){x[!is.na(x)]}}
    ))
    # clean data
    sample_data <- sample_data[,c("id2", "entry", "exit", "from", "to")]
    names(sample_data)[1] <- "id"
    # retrieve NA/AJ estimators
    etm_sample <- etm_short(sample_data, s = 0, t = max(times))

    c(etm_sample$NA_estimator[findInterval(times, vec = etm_sample$time),"1 2"],
      etm_sample$AJ_estimator["1","2", findInterval(times, vec = etm_sample$time)])
    
  })
  # compute EBS variances
  Var_EBS <- apply(NA_AJ_EBS, MARGIN = 1, FUN = var)
  # compute EBS quantiles
  q_EBS <- apply(1/sqrt(Var_EBS) * (NA_AJ_EBS - NA_AJ_original), MARGIN = 1, 
                 FUN = quantile, probs = 0.975, na.rm = TRUE)
  # determine EBS CI for Nelson-Aalen & Aalen-Johansen estimators
  CI_EBS <- matrix(
    pmax(0,
         rep(NA_AJ_original, 2) * 
           exp(c(rep(-1, 2*length(times)), rep(1, 2*length(times))) * 
                 rep(q_EBS, 2) * sqrt(rep(Var_EBS, 2)) / rep(NA_AJ_original, 2))),
    ncol = 4)[,c(1,3,2,4)]
  CI_EBS[,4] <- pmin(1, CI_EBS[,4])
  
  return(CI_EBS)
  
}

# compute CIs for the Nelson-Aalen & Aalen-Johansen estimators for "1 -> 2" transitions 
#  using wild bootstrap
WBS <- function(data_original, etm_original, times, BS_iter){
  
  # input:
  ## data_original: original data
  ## etm_original: Nelson-Aalen & Aalen-Johansen information based on original data
  ## times: time points at which the CIs should be evaluated
  ## BS_iter: number of bootstrap samples
  
  # output: matrix containing bootstrapped CI limits 
  #  (1st column: lower limit Nelson-Aalen estimator, 
  #   2nd column: upper limit Nelson-Aalen estimator,
  #   3rd column: lower limit Aalen-Johansen estimator,
  #   4th column: upper limit Aalen-Johansen estimator)
  
  # determine sample size (after omitting unobserved events)
  n <- length(unique(data_original$id))
  
  # compute WBS Z values (for NA residuals)
  delta_N <- etm_original$n.event["1","2",][etm_original$n.risk[,"1"] > 0]
  Y <- etm_original$n.risk[,"1"][etm_original$n.risk[,"1"] > 0]
  index <- findInterval(times, vec = etm_original$time[etm_original$n.risk[,"1"] > 0])
  Z <- sapply(1:BS_iter, FUN = function(i){
    G = ifelse(delta_N == 0, 0,
               sapply(delta_N, FUN = function(l){sum(rnorm(l, mean = 0, sd = 1))}))
    Z <- rep(0, length(times))
    Z[which(index != 0)] <- sqrt(n) * cumsum(G / Y)[index]
    return(Z)
  })
  
  # compute WBS zeta values (for AJ residuals)
  P.u.v.all <- 
    lapply(etm_original$time[findInterval(times, vec = etm_original$time)], 
           FUN = function(v){
             P.u.v.it <- lapply(etm_original$time[etm_original$time <= v], 
                                FUN = function(u){
                                  if (u != v){
                                    P.u.v.etm <- etm_short(data_original,
                                                           s = u, t = v)$AJ_estimator
                                    P.u.v <- P.u.v.etm[, , dim(P.u.v.etm)[3]]
                                  }else{
                                    P.u.v <- diag(1, 3)
                                    rownames(P.u.v) <- c("0", "1", "2")
                                    colnames(P.u.v) <- c("0", "1", "2")
                                  }
                                  return(P.u.v)
                                })
             names(P.u.v.it) <- as.character(etm_original$time[etm_original$time <= v])
             return(P.u.v.it)
           })
  names(P.u.v.all) <- as.character(
    etm_original$time[findInterval(times, vec = etm_original$time)]
  )
  zeta <- sapply(1:BS_iter, FUN = function(i){
    dxi.it <- lapply(1:3, FUN = function(j){
      sqrt(n) * 
        ifelse(etm_original$n.risk[, c("0", "0", "1")[j]] == 0, 0,
               sapply(etm_original$n.event[c("0", "0", "1")[j], c("1", "2", "2")[j], ], 
                      FUN = function(l){sum(rnorm(l, mean = 0, sd = 1))}) / 
                 etm_original$n.risk[, c("0", "0", "1")[j]])
    })
    names(dxi.it) <- c("0 1", "0 2", "1 2")
    dxi.u.all <- lapply(etm_original$time, FUN = function(u){
      dxi.u <- matrix(0, nrow = 3, ncol = 3)
      colnames(dxi.u) <- c("0", "1", "2")
      rownames(dxi.u) <- c("0", "1", "2")
      for(j in 1:3){
        dxi.u[c("0", "0", "1")[j], c("1", "2", "2")[j]] <- 
          dxi.it[[c("0 1", "0 2", "1 2")[j]]][[which(etm_original$time == u)]]
      }
      diag(dxi.u) <- -rowSums(dxi.u)
      return(dxi.u)
    })
    names(dxi.u.all) <- as.character(etm_original$time)
    sapply(etm_original$time[findInterval(times, vec = etm_original$time)], 
           FUN = function(v){
             sum(
               sapply(etm_original$time[etm_original$time <= v], FUN = function(u){
                 etm_original$AJ_estimator["1", , as.character(u)] %*% 
                   dxi.u.all[[as.character(u)]] %*% 
                   P.u.v.all[[as.character(v)]][[as.character(u)]][,"2"]
               })
             )
           })
  })
  
  # compute WBS variances
  Var_WBS <- apply(rbind(Z, zeta), MARGIN = 1, FUN = var)
  # compute WBS quantiles
  q_WBS <- apply(1/sqrt(Var_WBS) * rbind(Z, zeta), MARGIN = 1, 
                 FUN = quantile, probs = 0.975, na.rm = TRUE)
  # determine WBS CI for Nelson-Aalen & Aalen-Johansen estimators
  CI_WBS <- matrix(
    pmax(0, 
         rep(c(
           etm_original$NA_estimator[findInterval(times, vec = etm_original$time), "1 2"],
           etm_original$AJ_estimator["1","2", findInterval(times, vec = etm_original$time)]
         ), 2) * 
           exp(c(rep(-1, 2*length(times)), rep(1, 2*length(times))) * 
                 rep(q_WBS, 2) * sqrt(rep(Var_WBS, 2) / n) / 
                 rep(c(
                   etm_original$NA_estimator[
                     findInterval(times, vec = etm_original$time), "1 2"
                     ],
                   etm_original$AJ_estimator[
                     "1","2", findInterval(times, vec = etm_original$time)
                     ]
                 ), 2)
           )
    ),
    ncol = 4)[,c(1,3,2,4)]
  CI_WBS[,4] <- pmin(1, CI_WBS[,4])
  
  return(CI_WBS)
  
}

# combine ######################################################################

# generate data, apply both bootstrap approaches, and summarize results
bootstrap <- function(i, n, m, lambda01, lambda02, lambda12, staggered, 
                      times, BS_iter, NA_AJ_true){
  
  # input:
  ## i: parameter for parallel processing
  ## n: number of subjects
  ## m: number of events to observe before censoring is imposed
  ## lambda01: distribution parameter for "0 -> 1" transitions
  ## lambda02: distribution parameter for "0 -> 2" transitions
  ## lambda12: distribution parameter for "1 -> 2" transitions
  ## staggered: logical value indicating whether entry should be staggered
  ## times: time points at which the CIs should be evaluated
  ## BS_iter: number of bootstrap samples
  ## NA_AJ_true: true Nelson-Aalen & Aalen-Johansen estimators
  
  # output: vector with CI widths & coverages at given time points for both EBS & WBS
  
  ## generate data
  data <- generate_data(n, lambda01, lambda02, lambda12, staggered, censoring = TRUE, m)
  
  # retrieve NA/AJ info
  etm <- etm_short(data, s = 0, t = max(times))
  
  # compute CIs for the Nelson-Aalen & Aalen-Johansen estimators for "1 -> 2" transitions 
  #   by means of EBS
  CI_EBS <- EBS(data, c(etm$NA_estimator[findInterval(times, vec = etm$time),"1 2"],
                        etm$AJ_estimator["1","2", findInterval(times, vec = etm$time)]), 
                times, BS_iter)
  
  # compute CIs for the Nelson-Aalen & Aalen-Johansen estimators for "1 -> 2 transitions 
  #   by means of WBS
  CI_WBS <- WBS(data, etm, times, BS_iter)
  
  return(c(
    # NA estimator
    width_EBS_NA = CI_EBS[,2] - CI_EBS[,1], 
    coverage_EBS_NA = CI_EBS[,1] <= NA_AJ_true[1:length(times)] & 
      NA_AJ_true[1:length(times)] <= CI_EBS[,2], 
    width_WBS_NA = CI_WBS[,2] - CI_WBS[,1], 
    coverage_WBS_NA = CI_WBS[,1] <= NA_AJ_true[1:length(times)] & 
      NA_AJ_true[1:length(times)] <= CI_WBS[,2],
    # AJ estimator
    width_EBS_AJ = CI_EBS[,4] - CI_EBS[,3], 
    coverage_EBS_AJ = CI_EBS[,3] <= NA_AJ_true[(length(times)+1):length(NA_AJ_true)] & 
      NA_AJ_true[(length(times)+1):length(NA_AJ_true)] <= CI_EBS[,4], 
    width_WBS_AJ = CI_WBS[,4] - CI_WBS[,3], 
    coverage_WBS_AJ = CI_WBS[,3] <= NA_AJ_true[(length(times)+1):length(NA_AJ_true)] & 
      NA_AJ_true[(length(times)+1):length(NA_AJ_true)] <= CI_WBS[,4]
    ))
  
}

# run simulations ##############################################################

# run simulations and save summarized results
simulate <- function(seed, n, m, lambda01 = 0.01, lambda02 = 0.03, lambda12 = 0.1, 
                     staggered = FALSE, times = c(8,12,16), BS_iter = 1000){
  
  # input:
  ## seed: seed to use for simulations
  ## n: number of subjects
  ## m: number of events to observe before censoring is imposed
  ## lambda01: distribution parameter for "0 -> 1" transitions
  ## lambda02: distribution parameter for "0 -> 2" transitions
  ## lambda12: distribution parameter for "1 -> 2" transitions
  ## staggered: logical value indicating whether entry should be staggered
  ## times: time points at which the CIs should be evaluated
  ## BS_iter: number of bootstrap samples
  
  # output: vector with mean CI widths & coverages at given time points for both EBS & WBS, 
  #   for Nelson-Aalen & Aalen-Johansen estimators
  
  
  set.seed(seed)
  # compute true NA & AJ estimators (based on uncensored data) for comparison
  NA_AJ_true <- rowMeans(sapply(1:10000, FUN = function(i){
    etm <- etm_short(generate_data(n, lambda01, lambda02, lambda12, 
                                 staggered, censoring = FALSE),
                     s = 0, t = max(times))
    c(etm$NA_estimator[findInterval(times, vec = etm$time),"1 2"],
      etm$AJ_estimator["1","2", findInterval(times, vec = etm$time)])
  }))
  
  
  # TODO 
  # adapt the number of cores according to your needs
  # (use 7 cores to reproduce the results presented in the manuscript)
  cl <- makeCluster(detectCores() - 1, type="SOCK")
  clusterEvalQ(cl, library(tibble))
  clusterEvalQ(cl, library(etm))
  clusterExport(cl, "generate_data")
  clusterExport(cl, "etm_short")
  clusterExport(cl, "EBS")
  clusterExport(cl, "WBS")
  clusterSetRNGStream(cl, iseed = seed)
  res <- pblapply(X = 1:1000, FUN = bootstrap, n, m, lambda01, lambda02, lambda12, 
                  staggered, times, BS_iter, NA_AJ_true, cl = cl)
  stopCluster(cl)
  
  result <- colMeans(matrix(unlist(res), ncol = 8*length(times), byrow = TRUE), na.rm = TRUE)
  names(result) <- c(paste0("width EBS NA (t = ", times, ")"),
                     paste0("coverage EBS NA (t = ", times, ")"),
                     paste0("width WBS NA (t = ", times, ")"),
                     paste0("coverage WBS NA (t = ", times, ")"),
                     paste0("width EBS AJ (t = ", times, ")"),
                     paste0("coverage EBS AJ (t = ", times, ")"),
                     paste0("width WBS AJ (t = ", times, ")"),
                     paste0("coverage WBS AJ (t = ", times, ")"))
  return(result)
  
}


# perform simulations ##########################################################

library(tibble)
library(etm)
library(parallel)
library(pbapply)

result1 <- simulate(seed = 12345, n = 200, m = 100)
result2 <- simulate(seed = 12345, n = 200, m = 100, staggered = TRUE)
result3 <- simulate(seed = 12346, n = 100, m = 50)
result4 <- simulate(seed = 12346, n = 100, m = 50, staggered = TRUE)
result5 <- simulate(seed = 12347, n = 80, m = 40)
result6 <- simulate(seed = 12347, n = 80, m = 40, staggered = TRUE)
result7 <- simulate(seed = 12348, n = 50, m = 25)
result8 <- simulate(seed = 12348, n = 50, m = 25, staggered = TRUE)
