The R scripts in this repository can be used to reproduce the numerical results in the manuscript
'General Independent Censoring in Event-Driven Trials with Staggered Entry' (Rühl, Beyersmann, Friedrich):

- The script 'Simulations_calendar_times.R' includes the R code to reproduce the 
  results of the first simulation study.
- The script 'Breslow_estimators.R' includes the R code to reproduce the shadow plots 
  of the Breslow estimators based on the output of 'Simulations_calendar_times.R'.
- The script 'Bootstrap_functions.R' includes the R code to derive bootstrap confidence intervals
  for illness-death data without recovery (as used for the second simulation study and the real data analysis).
- The script 'Simulations_bootstrap.R' includes the R code to reproduce the results of
  the second simulation study.
- The script 'OAK_analysis.R' includes the R code for the analysis of the OAK study.

---

R version 4.1.2 (2021-11-01)

Platform: x86_64-w64-mingw32/x64 (64-bit)

Running under: Windows 10 x64 (build 19042)

Packages:
- parallel (v4.1.2)
- pbapply (v1.5.0)
- tibble (v3.1.6)
- survival (v3.2.13)
- dplyr (v1.0.7)
- ggplot2 (v3.3.5)
- gridExtra (v2.3)
- grid (v4.1.2)
- cowplot (v1.1.1)
- etm (v1.1)
- readxl (v1.3.1)
- pammtools (v0.5.8)
