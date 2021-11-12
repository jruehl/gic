The R scripts in this repository can be used to reproduce the numerical results in the manuscript
'General Independent Censoring in Event-Driven Trials with Staggered Entry' (Rühl, Beyersmann, Friedrich):

- The script 'Simulations_calendar_times.R' includes the R code to reproduce the 
  results of the first simulation study.
- The sript 'Breslow_estimators.R' includes the R code to reproduce the shadow plots 
  of the Breslow estimators based on the output of 'Simulations_calendar_times.R'.
- The script 'Simulations_bootstrap.R' includes the R code to reproduce the results of
  the second simulation study.

---

R version 3.5.1 (2018-07-02)

Platform: x86_64-w64-mingw32/x64 (64-bit)

Running under: Windows >= 8 x64 (build 9200)

Packages:
- parallel (v3.5.1)
- pbapply (v1.4.2)
- tibble (v2.1.3)
- survival (v3.1.12)
- dplyr (v0.8.3)
- ggplot2 (v3.3.2.9000)
- gridExtra (v2.3)
- grid (v3.5.1)
- cowplot (v1.0.0)
- etm (v1.1)
