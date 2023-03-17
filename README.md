# Longitudinal Incremental Propensity Score Interventions for Limited Resource Settings
# Introduction
Here we provide the code to reproduce the analysis described in: Sarvet AL, Wanis KN, Young J, Hernandez-Alejandro R, Stensrud MJ. Longitudinal Incremental Propensity Score Interventions for Limited Resource Settings, published in Biometrics, 2023.
> 

# Abstract
 Many real-life treatments are of limited supply and cannot be provided to all individuals in the population. For example, patients on the liver transplant waiting list usually cannot be assigned a liver transplant immediately at the time they reach highest priority because a suitable organ is not immediately available. In settings with limited supply, investigators are often interested in the effects of treatment strategies in which a limited proportion of patients receive an organ at a given time, that is, treatment regimes satisfying resource constraints. Here, we describe an estimand that allows us to define causal effects of treatment strategies that satisfy resource constraints:  Incremental Propensity Score Interventions for limited resources (IPSIs). IPSIs flexibly constrain time-varying resource utilization through proportional scaling of patients' natural propensities for treatment, thereby preserving existing propensity rank ordering compared to the status quo. We derive a simple class of inverse probability weighted estimators, and we apply one such estimator to evaluate the effect of restricting or expanding utilization of ‘increased risk’ liver organs to treat patients with end-stage liver disease.

Data are available from the Scientific Registry of Transplant Recipients by application

# Organization
- `Main.R` — R file which contains further details describing the analysis, and which allows users to reproduce the main analyses by running the code contained in `scripts`.
- `scripts`  — Scripts to reproduce the main analyses in the manuscript.

# Correspondence
If you have any questions, comments, or discover an error, please contact Aaron Sarvet at aaron.sarvet@epfl.ch
