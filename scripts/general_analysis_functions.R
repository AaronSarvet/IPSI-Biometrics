


estimate <- function(data, baseline_covariates,
										 time_varying_covariates_transplant, time_varying_covariates_censoring, 
										 increased_risk_usage_factor, standard_risk_usage_factor, 
										 use_existing_models=F, existing_model=NULL, K=136,
										 boot=F, seed=NULL, regime=NULL){
	
	if(boot){
		set.seed(seed)
		bootsamp<-sample(unique(data$patient_id), replace=T)
		bootsamp<-data.frame(patient_id=bootsamp, bootID=1:(length(bootsamp)))
		data<-bootsamp %>% left_join(data) %>% select(-patient_id) %>% rename(patient_id=bootID)
	}
	NPresults<-vector(2, mode="list")
	names(NPresults)<-c("Input", "Results")
	NPresults[["Input"]]<-vector(3, mode="list")
	names(NPresults[["Input"]])<-c("Model covariates", "Usage factors", "Quantiles")
	NPresults[["Input"]][["Model covariates"]]<-vector(3, mode="list")
	names(NPresults[["Input"]][["Model covariates"]])<-c("Baseline", "Time-varying (transplant)", "Time-varying (censoring)")
	if(use_existing_models==F){
		NPresults[["Input"]][["Model covariates"]][["Baseline"]]<-baseline_covariates
		NPresults[["Input"]][["Model covariates"]][["Time-varying (transplant)"]]<-time_varying_covariates_transplant
		NPresults[["Input"]][["Model covariates"]][["Time-varying (censoring)"]]<-time_varying_covariates_censoring
	}
	NPresults[["Input"]][["Usage factors"]]<-c(increased_risk_usage_factor, standard_risk_usage_factor)
	names(NPresults[["Input"]][["Usage factors"]])<-c("Increased risk_usage factor", "Standard_risk usage factor")
	NPresults[["Input"]][["Quantiles"]]<-vector(10, mode="list")
	if(use_existing_models==F){
		NPresults[["Input"]][["Quantiles"]][[1]]<- quantile(data$baseline_MELD,probs=c(0.35,0.65))
		NPresults[["Input"]][["Quantiles"]][[2]]<- quantile(data$baseline_MELD, probs=c(0.05,0.95))
		NPresults[["Input"]][["Quantiles"]][[3]]<- quantile(data$age,probs=c(0.35,0.65))
		NPresults[["Input"]][["Quantiles"]][[4]]<- quantile(data$age, probs=c(0.05,0.95))
		NPresults[["Input"]][["Quantiles"]][[5]]<- quantile(data$height,probs=c(0.35,0.65))
		NPresults[["Input"]][["Quantiles"]][[6]]<- quantile(data$height, probs=c(0.05,0.95))
		NPresults[["Input"]][["Quantiles"]][[7]]<- quantile(data$weight,probs=c(0.35,0.65))
		NPresults[["Input"]][["Quantiles"]][[8]]<- quantile(data$weight, probs=c(0.05,0.95))
		NPresults[["Input"]][["Quantiles"]][[9]]<- quantile(data$MELD,probs=c(0.35,0.65))
		NPresults[["Input"]][["Quantiles"]][[10]]<-quantile(data$MELD, probs=c(0.05,0.95))
	}
	NPresults[["Results"]]<-vector(3, mode="list")
	names(NPresults[["Results"]])<-c("Models", "Weights", "Non-parametric summaries (by WL)")
	NPresults[["Results"]][["Models"]]<-vector(5, mode="list")
	names(NPresults[["Results"]][["Models"]])<-c("Censoring (Admin)", "Censoring (LTFU)", "Censoring (transplant)", "Standard Risk", "Increased Risk")
	NPresults[["Results"]][["Weights"]]<-data.frame( 
		#Administrative variables 
		patient_id = data$patient_id,
		interval_WL = data$interval_WL,
		#Events 
		post_transplant = data$post_transplant,
		transplant_event = data$transplant_event,
		censoring_event_admin = data$censoring_event_admin,
		censoring_event = data$censoring_event,
		transplant_and_standard_risk = data$transplant_and_standard_risk,
		transplant_and_increased_risk = data$transplant_and_increased_risk,
		transplant_and_censored_at_tx = data$transplant_and_censored_at_tx,
		death_event = data$death_event,
		#Predicted treatment probabilities
		prob_censored_admin = NA,
		prob_censored = NA,
		prob_censored_at_transplant = NA,
		prob_increased_risk_transplant = NA,
		prob_standard_risk_transplant = NA,
		#Weights for censoring
		Wt_censoring_admin= NA,
		Wt_censoring= NA,
		Wt_censored_at_transplant=NA,
		Wt_no_int=NA,
		Wt_naive=NA,
		#Weights for resource constraints (wait-list entry time-scale)
		Wt_d_WL=1,
		Wt_dd_WL=1,
		Wt_A_1_nat_WL=1,
		Wt_A_2_nat_WL=1,
		Wt_g_WL=1 ) 
	NPresults[["Results"]][["Non-parametric summaries (by WL)"]]<-matrix(NA, K+1, 30)
	colnames(NPresults[["Results"]][["Non-parametric summaries (by WL)"]])<-c(
		"Lambda_A_1",  "Lambda_A_1_g_prime", "Lambda_A_1_g", "Lambda_A_1_g_naive",  "chi_1", "rho_1", "delta",   "E_A_1_KM", "E_A_1_g_KM", "E_A_1_gplus_KM", "E_A_1_gplus_naive_KM",
		"Lambda_A_2",  "Lambda_A_2_g_prime", "Lambda_A_2_g","chi_2", "rho_2", "gamma",  "E_A_2_KM", "E_A_2_g_KM",  "E_A_2_gplus_KM", "E_A_2_gplus_naive",
		"Lambda_Y_sub", "Lambda_Y_sub_g", "Lambda_Y_sub_naive", "Lambda_Y", "Lambda_Y_g", "Lambda_Y_naive",  "E_Y_KM",  "E_Y_g_KM", "E_Y_naive_KM"
	)
	rownames(NPresults[["Results"]][["Non-parametric summaries (by WL)"]])<-0:K

	
	
	
	
	if(use_existing_models==F){      
		#Model formulae      
		transplant_var_list <- list("ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))",
																paste(baseline_covariates, collapse=" + "),
																paste(time_varying_covariates_transplant, collapse=" + "))
		censoring_var_list <- list("ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))",
															 paste(baseline_covariates, collapse=" + "),
															 paste(time_varying_covariates_censoring, collapse=" + "))
		
		transplant_formula_RHS <- paste(
			transplant_var_list[transplant_var_list != ""], 
			collapse=" + "
		)
		censoring_formula_RHS <- paste(
			censoring_var_list[censoring_var_list != ""], 
			collapse=" + "
		)
		
		#Saving model results
		
		Begin<-Sys.time()
		
		NPresults[["Results"]][["Models"]][["Censoring (Admin)"]]<- (data %>% {table(.$interval_WL, .$censoring_event_admin)} %>% prop.table(margin=1) )[,2]
		
		gc()
		NPresults[["Results"]][["Models"]][["Censoring (LTFU)"]]<-  bigglm(formula=as.formula(paste("censoring_event", censoring_formula_RHS, sep=" ~ ")), 
																																			 data=data %>% filter(censoring_event_admin==0), 
																																			 family=binomial(), maxit=20)
		gc()
		NPresults[["Results"]][["Models"]][["Censoring (transplant)"]]<-  bigglm(formula=as.formula(paste("transplant_and_censored_at_tx", transplant_formula_RHS, sep=" ~ ")), 
																																						 data=data %>% filter(censoring_event_admin==0 & censoring_event==0 & post_transplant==0), 
																																						 family=binomial(), chunksize=(data %>% filter(censoring_event_admin==0 & censoring_event==0 & post_transplant==0) %>% nrow)/2, maxit=20)
		gc()
		NPresults[["Results"]][["Models"]][["Standard Risk"]]<-  bigglm(formula=as.formula(paste("transplant_and_standard_risk", transplant_formula_RHS, sep=" ~ ")), 
																																		data=data %>% filter(censoring_event_admin==0 & censoring_event==0 & post_transplant==0 & transplant_and_censored_at_tx == 0), 
																																		family=binomial(), chunksize=(data %>% filter(censoring_event_admin==0 & censoring_event==0 & post_transplant==0 & transplant_and_censored_at_tx == 0) %>% nrow)/3, maxit=20)
		gc()
		NPresults[["Results"]][["Models"]][["Increased Risk"]]<-  bigglm(formula=as.formula(paste("transplant_and_increased_risk", transplant_formula_RHS, sep=" ~ ")), 
																																		 data=data %>% filter(censoring_event_admin==0 & censoring_event==0 & post_transplant==0 & transplant_and_censored_at_tx == 0 & transplant_and_standard_risk==0), 
																																		 family=binomial(), chunksize=(data %>% filter(censoring_event_admin==0 & censoring_event==0 & post_transplant==0 & transplant_and_censored_at_tx == 0 & transplant_and_standard_risk==0) %>% nrow)/5, maxit=20)
		gc()
		#Saving probabilities of intervention events
		NPresults[["Results"]][["Weights"]]<-NPresults[["Results"]][["Weights"]] %>% mutate(
			prob_censored_admin = NPresults[["Results"]][["Models"]][["Censoring (Admin)"]][interval_WL+1],
			prob_censored = predict(NPresults[["Results"]][["Models"]][["Censoring (LTFU)"]], 
															newdata=data, type="response"),
			prob_censored_at_transplant = 
				predict(NPresults[["Results"]][["Models"]][["Censoring (transplant)"]], 
								newdata=data, type="response"),
			prob_standard_risk_transplant =
				predict(NPresults[["Results"]][["Models"]][["Standard Risk"]], 
								newdata=data, type="response"),
			prob_increased_risk_transplant = 
				predict(NPresults[["Results"]][["Models"]][["Increased Risk"]], 
								newdata=data, type="response"),
			Wt_censoring_admin = case_when(censoring_event_admin == 1 ~ 0,
																		 TRUE ~ 1/(1-prob_censored_admin)),
			Wt_censoring = case_when(censoring_event == 1 ~ 0,
															 TRUE ~ 1/(1-prob_censored)),
			Wt_censored_at_transplant = case_when(post_transplant == 1 ~ 1,
																						transplant_and_censored_at_tx == 1 ~ 0,
																						TRUE ~ 1/(1-prob_censored_at_transplant)),
			Wt_censored_high_risk = case_when(post_transplant == 1 ~ 1,
																				transplant_and_increased_risk == 1 ~ 0,
																				TRUE ~ 1/(1-prob_increased_risk_transplant)),
			Wt_d_WL = 1,
			Wt_dd_WL = 1
		) %>% group_by(patient_id) %>% mutate(
			Wt_no_int = cumprod(Wt_censoring_admin*Wt_censoring*Wt_censored_at_transplant),
			Wt_naive  = cumprod(Wt_censoring_admin*Wt_censoring*Wt_censored_at_transplant*Wt_censored_high_risk)
		) %>%  ungroup()
		
		gc()  
	}
	if(use_existing_models==T){       
		NPresults[["Results"]][["Weights"]]<-NPresults[["Results"]][["Weights"]] %>% mutate(
			prob_censored_admin = existing_model[["Censoring (Admin)"]][interval_WL+1],
			prob_censored = predict(existing_model[["Censoring (LTFU)"]], 
															newdata=data, type="response"),
			prob_censored_at_transplant = 
				predict(existing_model[["Censoring (transplant)"]], 
								newdata=data, type="response"),
			prob_standard_risk_transplant =
				predict(existing_model[["Standard Risk"]], 
								newdata=data, type="response"),
			prob_increased_risk_transplant = 
				predict(existing_model[["Increased Risk"]], 
								newdata=data, type="response"),
			Wt_censoring_admin = case_when(censoring_event_admin == 1 ~ 0,
																		 TRUE ~ 1/(1-prob_censored_admin)),
			Wt_censoring = case_when(censoring_event == 1 ~ 0,
															 TRUE ~ 1/(1-prob_censored)),
			Wt_censored_at_transplant = case_when(post_transplant == 1 ~ 1,
																						transplant_and_censored_at_tx == 1 ~ 0,
																						TRUE ~ 1/(1-prob_censored_at_transplant)),
			Wt_censored_high_risk = case_when(post_transplant == 1 ~ 1,
																				transplant_and_increased_risk == 1 ~ 0,
																				TRUE ~ 1/(1-prob_increased_risk_transplant)),
			Wt_d_WL = 1,
			Wt_dd_WL = 1
			) %>% group_by(patient_id) %>% mutate(
			Wt_no_int = cumprod(Wt_censoring_admin*Wt_censoring*Wt_censored_at_transplant),
			Wt_naive = cumprod(Wt_censoring_admin*Wt_censoring*Wt_censored_at_transplant*Wt_censored_high_risk)
		) %>%  ungroup()
		
		gc()          
	}     
	
	NPresults[["Results"]][2:3]<- adjust_wts(
		data_with_weights = NPresults[["Results"]][["Weights"]], 
		increased_risk_usage_factor = increased_risk_usage_factor, 
		standard_risk_usage_factor = standard_risk_usage_factor,
		NPsumm_WL  = NPresults[["Results"]][["Non-parametric summaries (by WL)"]],
		maxtime = K
	) 
	
	NPresults[["Results"]][["Weights"]]<-NPresults[["Results"]][["Weights"]] %>% select(Wt_g_WL)
	if(boot==F){
		return(NPresults)
	}
	if(boot==T){
		bootstrap_results[[1]][seed,1:(K+1) , ,regime]<<-NPresults[["Results"]][[3]]
		if(regime=="g_1" & use_existing_models==F){
			bootstrap_results[[2]][[seed]]<<-NPresults[["Results"]][[1]]
		}
		return(NULL)
		gc()
	}
	
}

adjust_wts <- function(data_with_weights, increased_risk_usage_factor, standard_risk_usage_factor, NPsumm_WL,  maxtime) {
		#Initialize cumulative treatment weights
		data_with_weights$Wt_g_WL<-data_with_weights$Wt_no_int
		n<-length(unique(data_with_weights$patient_id))
		#Iterate over time-scale
		for (j in 1:(maxtime+1)) {
			#Update weights in interval j-1, prior to intervention on transplants of either type 
			data_with_weights <- data_with_weights %>%  
				mutate(
					Wt_A_1_nat_WL  =  case_when(
						interval_WL==0 ~   Wt_no_int,
						interval_WL==j-1 ~ lag(Wt_g_WL)*Wt_censoring_admin*Wt_censoring*Wt_censored_at_transplant,
						TRUE             ~ Wt_A_1_nat_WL
					)
				) 
			#Compute law-dependent parameters for IPSI intervention on A1 in interval j-1
			data_j<-data_with_weights %>% filter(interval_WL == j-1)
			NPsumm_WL[j, "Lambda_A_1"]   <- data_j %>% {sum(.$transplant_and_standard_risk*.$Wt_no_int)     / sum((1-.$post_transplant)*.$Wt_no_int)}
			NPsumm_WL[j, "Lambda_A_1_g_prime"] <- data_j %>% {sum(.$transplant_and_standard_risk*.$Wt_A_1_nat_WL) / sum((1-.$post_transplant)*.$Wt_A_1_nat_WL)}
			NPsumm_WL[j, "Lambda_A_1_g_naive"]   <- data_j %>% {sum(.$transplant_and_standard_risk*.$Wt_naive)     / sum((1-.$post_transplant)*.$Wt_naive)}
			NPsumm_WL[j, "E_A_1_KM"]      <- ifelse(j==1, NPsumm_WL[j, "Lambda_A_1"],
																							NPsumm_WL[j, "Lambda_A_1"]*
																								prod(1-NPsumm_WL[1:(j-1), "Lambda_A_1"])*
																								prod(1-NPsumm_WL[1:(j-1), "Lambda_A_2"])*
																								prod(1-NPsumm_WL[1:(j-1), "Lambda_Y_sub"])          
			)
			NPsumm_WL[j, "E_A_1_g_KM"]<-ifelse(j==1, NPsumm_WL[j, "Lambda_A_1_g_prime"],
																				 NPsumm_WL[j, "Lambda_A_1_g_prime"]*
																				 	prod(1-NPsumm_WL[1:(j-1), "Lambda_A_1_g"])*
																				 	prod(1-NPsumm_WL[1:(j-1), "Lambda_A_2_g"])*
																				 	prod(1-NPsumm_WL[1:(j-1), "Lambda_Y_sub_g"])          
			)
			NPsumm_WL[j, "chi_1"]      <-NPsumm_WL[j, "Lambda_A_1"] <= 1
			NPsumm_WL[j, "rho_1"]      <-ifelse(NPsumm_WL[j, "Lambda_A_1"]==0, NPsumm_WL[j-1, "rho_1"],
																					(standard_risk_usage_factor*NPsumm_WL[j, "E_A_1_KM"] <= NPsumm_WL[j, "E_A_1_g_KM"])
			)
			
			NPsumm_WL[j, "delta"]    <- ifelse(!NPsumm_WL[j, "Lambda_A_1"], 1,
																					 ( ( (   standard_risk_usage_factor*NPsumm_WL[j, "E_A_1_KM"]) / (   NPsumm_WL[j, "E_A_1_g_KM"]) )^ (  NPsumm_WL[j, "rho_1"]) *
																					 		( (1- standard_risk_usage_factor*NPsumm_WL[j, "E_A_1_KM"]*NPsumm_WL[j, "Lambda_A_1_g_prime"]/NPsumm_WL[j, "E_A_1_g_KM"]) / (1- NPsumm_WL[j, "Lambda_A_1_g_prime"]) )^ (1-NPsumm_WL[j, "rho_1"])
																					 )*NPsumm_WL[j, "chi_1"]
			)
			
			
			#Update weights for IPSI intervention on A1 in interval j-1  
			data_with_weights <- data_with_weights %>% mutate(
				Wt_d_WL = case_when(interval_WL                  != j-1 ~ Wt_d_WL,
														post_transplant              == 1   ~ 1,
														transplant_and_standard_risk == 1   ~ 
															(   NPsumm_WL[j, "delta"]                                                                     )^(  NPsumm_WL[j, "rho_1"])*
															((1-NPsumm_WL[j, "delta"]*(1-prob_standard_risk_transplant))/   prob_standard_risk_transplant )^(1-NPsumm_WL[j, "rho_1"]),
															((1-NPsumm_WL[j, "delta"]*   prob_standard_risk_transplant )/(1-prob_standard_risk_transplant))^(  NPsumm_WL[j, "rho_1"])*
															(   NPsumm_WL[j, "delta"]                                                                     )^(1-NPsumm_WL[j, "rho_1"])
				),
				Wt_A_2_nat_WL =  case_when(
					interval_WL                  == j-1 ~ Wt_A_1_nat_WL*Wt_d_WL,
					TRUE                                ~ Wt_A_2_nat_WL
				)
				
			)
			#Compute law-dependent parameters for IPSI intervention on A2 in interval j-1
			data_j<-data_with_weights %>% filter(interval_WL == j-1)
			
			
			NPsumm_WL[j, "Lambda_A_1_g"]<-data_j %>% {sum(.$transplant_and_standard_risk*.$Wt_A_2_nat_WL) / sum((1-.$post_transplant)*.$Wt_A_2_nat_WL)}
			NPsumm_WL[j, "E_A_1_gplus_KM"]<-ifelse(j==1, NPsumm_WL[j, "Lambda_A_1_g"],
																						 NPsumm_WL[j, "Lambda_A_1_g"]*
																						 	prod(1-NPsumm_WL[1:(j-1), "Lambda_A_1_g"])*
																						 	prod(1-NPsumm_WL[1:(j-1), "Lambda_A_2_g"])*
																						 	prod(1-NPsumm_WL[1:(j-1), "Lambda_Y_sub_g"])          
			)
			NPsumm_WL[j, "E_A_1_gplus_naive_KM"]<-ifelse(j==1, NPsumm_WL[j, "Lambda_A_1_g_naive"],
																									 NPsumm_WL[j, "Lambda_A_1_g_naive"]*
																									 	prod(1-NPsumm_WL[1:(j-1), "Lambda_A_1_g_naive"])*
																									 	prod(1-NPsumm_WL[1:(j-1), "Lambda_Y_sub_naive"])          
			)
			NPsumm_WL[j, "Lambda_A_2"]   <- data_j %>% {sum(.$transplant_and_increased_risk*.$Wt_no_int)     / sum((1-.$post_transplant)*(1-.$transplant_and_standard_risk)*.$Wt_no_int)}
			NPsumm_WL[j, "Lambda_A_2_g_prime"] <- data_j %>% {sum(.$transplant_and_increased_risk*.$Wt_A_2_nat_WL) / sum((1-.$post_transplant)*(1-.$transplant_and_standard_risk)*.$Wt_A_2_nat_WL)}
			NPsumm_WL[j, "E_A_2_KM"]      <- ifelse(j==1, NPsumm_WL[j, "Lambda_A_2"],
																							NPsumm_WL[j, "Lambda_A_2"]*
																								prod(1-NPsumm_WL[1:(j)  , "Lambda_A_1"])*
																								prod(1-NPsumm_WL[1:(j-1), "Lambda_A_2"])*
																								prod(1-NPsumm_WL[1:(j-1), "Lambda_Y_sub"])         
			)
			NPsumm_WL[j, "E_A_2_g_KM"]<-ifelse(j==1, NPsumm_WL[j, "Lambda_A_2_g_prime"],
																				 NPsumm_WL[j, "Lambda_A_2_g_prime"]*
																				 	prod(1-NPsumm_WL[1:(j), "Lambda_A_1_g"])*
																				 	prod(1-NPsumm_WL[1:(j-1), "Lambda_A_2_g"])*
																				 	prod(1-NPsumm_WL[1:(j-1), "Lambda_Y_sub_g"])          
			)
			NPsumm_WL[j, "chi_2"]      <-NPsumm_WL[j, "Lambda_A_2"] <= 1
			NPsumm_WL[j, "rho_2"]      <-ifelse(NPsumm_WL[j, "Lambda_A_2"]==0, NPsumm_WL[j-1, "rho_2"],
																					(increased_risk_usage_factor*NPsumm_WL[j, "E_A_2_KM"] <= NPsumm_WL[j, "E_A_2_g_KM"])
			)
			
			
			NPsumm_WL[j, "gamma"]    <- ifelse(!NPsumm_WL[j, "Lambda_A_2"], 1,
																					 ( ( (   increased_risk_usage_factor*NPsumm_WL[j, "E_A_2_KM"]) / (   NPsumm_WL[j, "E_A_2_g_KM"]) )^ (  NPsumm_WL[j, "rho_2"]) *
																					 		( (1- increased_risk_usage_factor*NPsumm_WL[j, "E_A_2_KM"]*NPsumm_WL[j, "Lambda_A_2_g_prime"]/NPsumm_WL[j, "E_A_2_g_KM"]) / (1- NPsumm_WL[j, "Lambda_A_2_g_prime"]) )^ (1-NPsumm_WL[j, "rho_2"])
																					 )*NPsumm_WL[j, "chi_2"]
			)      
			#Update weights for IPSI intervention on A2 in interval j-1  
			data_with_weights <- data_with_weights %>% mutate(
				Wt_dd_WL = case_when(interval_WL                  != j-1 ~ Wt_dd_WL,
														 post_transplant              == 1 |
														 	transplant_and_standard_risk == 1    ~ 1,
														 transplant_and_increased_risk == 1   ~ 
														 	(   NPsumm_WL[j, "gamma"]                                                                       )^(  NPsumm_WL[j, "rho_2"])*
														 	((1-NPsumm_WL[j, "gamma"]*(1-prob_increased_risk_transplant))/   prob_increased_risk_transplant )^(1-NPsumm_WL[j, "rho_2"]),
														 TRUE                                ~
														 	((1-NPsumm_WL[j, "gamma"]*   prob_increased_risk_transplant )/(1-prob_increased_risk_transplant))^(  NPsumm_WL[j, "rho_2"])*
														 	(   NPsumm_WL[j, "gamma"]                                                                       )^(1-NPsumm_WL[j, "rho_2"])
				),
				Wt_g_WL =  case_when(
					interval_WL                  == j-1 ~ Wt_A_2_nat_WL*Wt_dd_WL,
					TRUE                                ~ Wt_g_WL
				)
			)
			
			#Compute final law-dependent parameters for interval j-1
			data_j<-data_with_weights %>% filter(interval_WL == j-1)
			NPsumm_WL[j, "Lambda_A_2_g"]<-data_j %>% {sum(.$transplant_and_increased_risk*.$Wt_g_WL) / sum((1-.$post_transplant)*(1-.$transplant_and_standard_risk)*.$Wt_g_WL)}
			NPsumm_WL[j, "E_A_2_gplus_KM"]<-ifelse(j==1, NPsumm_WL[j, "Lambda_A_2_g"],
																						 NPsumm_WL[j, "Lambda_A_2_g"]*
																						 	prod(1-NPsumm_WL[1:(j), "Lambda_A_1_g"])*
																						 	prod(1-NPsumm_WL[1:(j-1), "Lambda_A_2_g"])*
																						 	prod(1-NPsumm_WL[1:(j-1), "Lambda_Y_sub_g"])          
			)
			NPsumm_WL[j, "E_A_2_gplus_naive"]<-0
			NPsumm_WL[j, "Lambda_Y_sub_g"]  <-data_j %>% {sum(.$death_event*(1-.$post_transplant)*(1-.$transplant_and_standard_risk)*(1-.$transplant_and_increased_risk)*.$Wt_g_WL)   / sum((1-.$post_transplant)*(1-.$transplant_and_standard_risk)*(1-.$transplant_and_increased_risk)*.$Wt_g_WL)}        
			NPsumm_WL[j, "Lambda_Y_sub"]    <-data_j %>% {sum(.$death_event*(1-.$post_transplant)*(1-.$transplant_and_standard_risk)*(1-.$transplant_and_increased_risk)*.$Wt_no_int) / sum((1-.$post_transplant)*(1-.$transplant_and_standard_risk)*(1-.$transplant_and_increased_risk)*.$Wt_no_int)}
			NPsumm_WL[j, "Lambda_Y_sub_naive"]    <-data_j %>% {sum(.$death_event*(1-.$post_transplant)*(1-.$transplant_and_standard_risk)*(1-.$transplant_and_increased_risk)*.$Wt_naive) / sum((1-.$post_transplant)*(1-.$transplant_and_standard_risk)*(1-.$transplant_and_increased_risk)*.$Wt_naive)}
			NPsumm_WL[j, "Lambda_Y_g"]    <-data_j %>% {sum(.$death_event*.$Wt_g_WL)   / sum(.$Wt_g_WL)}        
			NPsumm_WL[j, "Lambda_Y"]      <-data_j %>% {sum(.$death_event*.$Wt_no_int)   / sum(.$Wt_no_int)}
			NPsumm_WL[j, "Lambda_Y_naive"]      <-data_j %>% {sum(.$death_event*.$Wt_naive)   / sum(.$Wt_naive)}
			NPsumm_WL[j, "E_Y_g_KM"]      <-1- prod(1-NPsumm_WL[1:j, "Lambda_Y_g"])    
			NPsumm_WL[j, "E_Y_KM"]        <-1- prod(1-NPsumm_WL[1:j, "Lambda_Y"])       
			NPsumm_WL[j, "E_Y_naive_KM"]  <-1- prod(1-NPsumm_WL[1:j, "Lambda_Y_naive"])          
			
			print(j)
			gc()
		}    
	
	return(list(data_with_weights, NPsumm_WL))
	
}




