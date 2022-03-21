# This script includes the functions to derive confidence intervals for the
#   Nelson-Aalen/Aalen-Johansen estimator for illness-death data using both the
#   non-parametric and the wild bootstrap

# analyze data #################################################################

# derive Nelson-Aalen/Aalen-Johansen estimator
#   (shorter version of etm() including Nelson-Aalen estimator)
etm_short <- function(data, s, t, NAE, AJE){
  
  # input:
  ## data: data to derive Nelson-Aalen/Aalen-Johansen estimator from
  ## s: starting time for Aalen-Johansen estimator
  ## t: end time for Aalen-Johansen estimator
  ## NAE: logical value indicating whether Nelson-Aalen estimator should be computed
  ## AJE: logical value indicating whether Aalen-Johansen estimator should be computed
  
  # output: list with the subsequent entries:
  ## time: event times at which the transition probabilities are computed
  ## NA_estimator: Nelson-Aalen estimator (if applicable)
  ## AJ_estimator: Aalen-Johansen estimator (if applicable)
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
  if(AJE) 
    AJ_estimator <- array(diag(1, nrow = 3, ncol = 3), dim = c(3, 3, 1)) 
  else 
    AJ_estimator <- NULL
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
      if(AJE) AJ_estimator <- array(c(diag(1, nrow = 3, ncol = 3), trans$est),
                                    dim = c(3, 3, length(time)))
    }else{
      time <- trans$time
      n.risk <- trans$n.risk[, c(1,2), drop = FALSE]
      n.event <- trans$n.event
      if(AJE) AJ_estimator <- trans$est
    }
    colnames(n.risk) <- c("0", "1")
    dimnames(n.event) <- list(c("0", "1", "2"), c("0", "1", "2"), time)
    if(AJE) dimnames(AJ_estimator) <- list(c("0", "1", "2"), c("0", "1", "2"), time)
    if(NAE){
      NA_estimator <- cbind(cumsum(replace(n.event["0","1",] / n.risk[,"0"], 
                                           list = n.risk[,"0"] == 0, values = 0)),
                            cumsum(replace(n.event["0","2",] / n.risk[,"0"], 
                                           list = n.risk[,"0"] == 0, values = 0)),
                            cumsum(replace(n.event["1","2",] / n.risk[,"1"], 
                                           list = n.risk[,"1"] == 0, values = 0)))
      colnames(NA_estimator) <- c("0 1", "0 2", "1 2")
    }
  }
  
  list(time = time,
       NA_estimator = NA_estimator, 
       AJ_estimator= AJ_estimator,
       n.risk = n.risk,
       n.event = n.event)
  
}

# apply bootstrap ##############################################################

# compute CIs for the Nelson-Aalen/Aalen-Johansen estimators for the transition of interest
#   using Efron's non-parametric bootstrap
EBS <- function(data_original, etm_original, times, from = "1", to = "2", NAE, AJE, BS_iter){
  
  # input:
  ## data_original: original data to base bootstrap on
  ## etm_original: Nelson-Aalen/Aalen-Johansen information based on original data
  ## times: time points at which CIs should be evaluated
  ## from: initial state of the transition of interest
  ## to: target state of the transition of interest
  ## NAE: logical value indicating whether CIs for Nelson-Aalen estimator should be computed
  ## AJE: logical value indicating whether CIs for Aalen-Johansen estimator should be computed
  ## BS_iter: number of bootstrap samples
  
  # output: matrix containing bootstrapped CI limits 
  #   (if NAE = TRUE & AJE = TRUE: 
  #    1st column: lower limit Nelson-Aalen estimator, 
  #    2nd column: upper limit Nelson-Aalen estimator,
  #    3rd column: lower limit Aalen-Johansen estimator,
  #    4th column: upper limit Aalen-Johansen estimator)
  
  # compute EBS NA/AJ estimators
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
    etm_sample <- etm_short(sample_data, s = 0, t = max(times), NAE = NAE, AJE = AJE)
    
    c(etm_sample$NA_estimator[findInterval(times, vec = etm_sample$time), 
                              paste0(from, " ", to)],
      etm_sample$AJ_estimator[as.character(from), as.character(to), 
                              findInterval(times, vec = etm_sample$time)])
    
  })
  
  # retrieve original NA/AJ estimators
  NA_AJ_original <- c(etm_original$NA_estimator[findInterval(times, vec = etm_original$time), 
                                                paste0(from, " ", to)],
                      etm_original$AJ_estimator[as.character(from), as.character(to), 
                                                findInterval(times, vec = etm_original$time)])
  # compute EBS variances
  Var_EBS <- apply(NA_AJ_EBS, MARGIN = 1, FUN = var)
  # compute EBS quantiles
  q_EBS <- apply(1/sqrt(Var_EBS) * (NA_AJ_EBS - NA_AJ_original), MARGIN = 1, 
                 FUN = quantile, probs = 0.975, na.rm = TRUE)
  # determine EBS CIs for NA/AJ estimators
  CI_EBS <- matrix(
    pmax(0,
         rep(NA_AJ_original, 2) * 
           exp(c(rep(-1, sum(NAE,AJE)*length(times)), rep(1, sum(NAE,AJE)*length(times))) * 
                 rep(q_EBS, 2) * sqrt(rep(Var_EBS, 2)) / 
                 rep(NA_AJ_original, 2))),
    ncol = 2*(sum(NAE,AJE)))
  if(sum(NAE,AJE) == 2) CI_EBS <- CI_EBS[,c(1,3,2,4)]
  # adjust invalid CI values
  if(NAE) CI_EBS[which(is.na(CI_EBS[,1])), 1:2] <- 0
  if(AJE){
    CI_EBS[which(is.na(CI_EBS[,2*NAE+1])), c(2*NAE+1,2*NAE+2)] <- 0
    CI_EBS[,2*NAE+2] <- pmin(1, CI_EBS[,2*NAE+2])
  }
  
  return(CI_EBS)
  
}

# compute CIs for the Nelson-Aalen/Aalen-Johansen estimators for the transition of interest
#  using wild bootstrap
WBS <- function(data_original, etm_original, times, from = "1", to = "2", NAE, AJE, BS_iter){
  
  # input:
  ## data_original: original data
  ## etm_original: Nelson-Aalen & Aalen-Johansen information based on original data
  ## times: time points at which the CIs should be evaluated
  ## from: initial state of the transition of interest
  ## to: target state of the transition of interest
  ## NAE: logical value indicating whether CIs for Nelson-Aalen estimator should be computed
  ## AJE: logical value indicating whether CIs for Aalen-Johansen estimator should be computed
  ## BS_iter: number of bootstrap samples
  
  # output: matrix containing bootstrapped CI limits 
  #   (if NAE = TRUE & AJE = TRUE: 
  #    1st column: lower limit Nelson-Aalen estimator, 
  #    2nd column: upper limit Nelson-Aalen estimator,
  #    3rd column: lower limit Aalen-Johansen estimator,
  #    4th column: upper limit Aalen-Johansen estimator)
  
  # determine sample size
  n <- length(unique(data_original$id))
  
  if(NAE){
    # compute WBS Z values (for NA residuals)
    delta_N <- etm_original$n.event[as.character(from), as.character(to),
                                    etm_original$n.risk[,as.character(from)] > 0]
    Y <- etm_original$n.risk[etm_original$n.risk[,as.character(from)] > 0, as.character(from)]
    index <- findInterval(times, vec = etm_original$time[etm_original$n.risk[,as.character(from)] > 0])
    Z <- sapply(1:BS_iter, FUN = function(i){
      G = ifelse(delta_N == 0, 0,
                 sapply(delta_N, FUN = function(l){sum(rnorm(l, mean = 0, sd = 1))}))
      replace(rep(0, length(times)), which(index != 0), sqrt(n) * cumsum(G / Y)[index])
    })
  }
  
  if(AJE){
    # compute WBS zeta values (for AJ residuals)
    P.u.v.all <- 
      lapply(etm_original$time[findInterval(times, vec = etm_original$time)], 
             FUN = function(v){
               P.u.v.it <- lapply(etm_original$time[etm_original$time <= v], 
                                  FUN = function(u){
                                    if (u != v){
                                      P.u.v.etm <- etm_short(data_original,
                                                             s = u, t = v, 
                                                             NAE = FALSE, AJE = TRUE)$AJ_estimator
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
                     P.u.v.all[[as.character(v)]][[as.character(u)]][,as.character(to)]
                 })
               )
             })
    })
  }
  
  # retrieve original NA/AJ estimators
  NA_AJ_original <- c(etm_original$NA_estimator[findInterval(times, vec = etm_original$time), 
                                                paste0(from, " ", to)],
                      etm_original$AJ_estimator[as.character(from), as.character(to), 
                                                findInterval(times, vec = etm_original$time)])
  # compute WBS variances
  Var_WBS <- apply(rbind(if(NAE) Z, 
                         if(AJE) zeta), 
                   MARGIN = 1, FUN = var)
  # compute WBS quantiles
  q_WBS <- apply(1/sqrt(Var_WBS) * rbind(if(NAE) Z, 
                                         if(AJE) zeta), 
                 MARGIN = 1, 
                 FUN = quantile, probs = 0.975, na.rm = TRUE)
  # determine WBS CIs for NA/AJ estimators
  CI_WBS <- matrix(
    pmax(0, 
         rep(NA_AJ_original, 2) * 
           exp(c(rep(-1, sum(NAE,AJE)*length(times)), rep(1, sum(NAE,AJE)*length(times))) * 
                 rep(q_WBS, 2) * sqrt(rep(Var_WBS, 2) / n) / 
                 rep(NA_AJ_original, 2)
           )
    ),
    ncol = 2*sum(NAE,AJE))
  if(sum(NAE,AJE) == 2) CI_WBS <- CI_WBS[,c(1,3,2,4)]
  # adjust invalid CI values
  if(NAE) CI_WBS[which(is.na(CI_WBS[,1])), 1:2] <- 0
  if(AJE){
    CI_WBS[which(is.na(CI_WBS[,2*NAE+1])), c(2*NAE+1,2*NAE+2)] <- 0
    CI_WBS[,2*NAE+2] <- pmin(1, CI_WBS[,2*NAE+2])
  }
  
  return(CI_WBS)
  
}
