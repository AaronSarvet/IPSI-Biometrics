# First load and install the following packages which we will need for data cleaning and analysis.

requiredPackages = c(
	"plyr",
	"tidyverse",
	"data.table",
	"pbapply",
	"ggthemes",
	"splines",
	"parallel",
	"haven",
	"Hmisc",
	"biglmm",
	"showtext"
)

for (p in requiredPackages) {
	if (!require(p, character.only = TRUE))
		install.packages(p)
	library(p, character.only = TRUE)
}

# We will need several functions which can be found in \scripts\General analysis functions.R

source("./scripts/general_analysis_functions.R")

# The data for this analysis can be obtained by application to the Scientific Registry of Transplant Recipients
# We need to clean the data using the code in \scripts\data cleaning.R
# The data cleaning step will save a file called candidates_analysis.RData in the current directory in \data folder 
# which can be loaded later instead of repeating the data cleaning 

source("./scripts/data_cleaning.R")





################
################
#Point estimates
################
################


results<-vector(3, mode="list")
names(results)<-c("1, 0", "1, 1.25", "1, 1.5")

results[[1]]<-
	estimate(
		data=candidates_analysis,
		baseline_covariates = list("ns(baseline_MELD,knots=NPresults[[1]][[3]][[1]], Boundary.knots=NPresults[[1]][[3]][[2]])", 
															 "baseline_MELD_exception", "status1", "gender", "race", "year_of_listing",
															 "ns(age,knots=NPresults[[1]][[3]][[3]], Boundary.knots=NPresults[[1]][[3]][[4]])", 
															 "ns(height,knots=NPresults[[1]][[3]][[5]], Boundary.knots=NPresults[[1]][[3]][[6]])", 
															 "ns(weight,knots=NPresults[[1]][[3]][[7]], Boundary.knots=NPresults[[1]][[3]][[8]])",
															 "accept_incompatible_blood_type", "accept_extra_corporeal_liver", "accept_liver_segment",
															 "accept_HBV_positive_donor", "accept_HCV_positive_donor",
															 "patient_on_life_support", "functional_status", "primary_diagnosis", "spontaneous_bacterial_peritonitis",
															 "history_of_PV_thrombosis", "history_of_TIPSS"),
		time_varying_covariates_transplant = list("ns(MELD,knots=NPresults[[1]][[3]][[9]], Boundary.knots=NPresults[[1]][[3]][[10]])*
																							MELD_exception*
																							ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))"),
		time_varying_covariates_censoring = list("ns(MELD,knots=NPresults[[1]][[3]][[9]], Boundary.knots=NPresults[[1]][[3]][[10]])*
																						 MELD_exception*
																						 ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))",
																						 "post_transplant*
																						 ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))"),
		increased_risk_usage_factor = 0, #theta_1
		standard_risk_usage_factor = 1  #theta_2
		)



results[[2]]<-
	estimate(
		data=candidates_analysis,
		baseline_covariates = list("ns(baseline_MELD,knots=NPresults[[1]][[3]][[1]], Boundary.knots=NPresults[[1]][[3]][[2]])", 
															 "baseline_MELD_exception", "status1", "gender", "race", "year_of_listing",
															 "ns(age,knots=NPresults[[1]][[3]][[3]], Boundary.knots=NPresults[[1]][[3]][[4]])", 
															 "ns(height,knots=NPresults[[1]][[3]][[5]], Boundary.knots=NPresults[[1]][[3]][[6]])", 
															 "ns(weight,knots=NPresults[[1]][[3]][[7]], Boundary.knots=NPresults[[1]][[3]][[8]])",
															 "accept_incompatible_blood_type", "accept_extra_corporeal_liver", "accept_liver_segment",
															 "accept_HBV_positive_donor", "accept_HCV_positive_donor",
															 "patient_on_life_support", "functional_status", "primary_diagnosis", "spontaneous_bacterial_peritonitis",
															 "history_of_PV_thrombosis", "history_of_TIPSS"),
		time_varying_covariates_transplant = list("ns(MELD,knots=NPresults[[1]][[3]][[9]], Boundary.knots=NPresults[[1]][[3]][[10]])*
																							MELD_exception*
																							ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))"),
		time_varying_covariates_censoring = list("ns(MELD,knots=NPresults[[1]][[3]][[9]], Boundary.knots=NPresults[[1]][[3]][[10]])*
																						 MELD_exception*
																						 ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))",
																						 "post_transplant*
																						 ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))"),
		increased_risk_usage_factor = 1.25, #theta_1
		standard_risk_usage_factor = 1,  #theta_2
		use_existing_models=T, existing_model=results[[1]][["Results"]][["Models"]] #Weights for all regimes are computed using the same model
		)

results[[3]]<-
	estimate(
		data=candidates_analysis,
		baseline_covariates = list("ns(baseline_MELD,knots=NPresults[[1]][[3]][[1]], Boundary.knots=NPresults[[1]][[3]][[2]])", 
															 "baseline_MELD_exception", "status1", "gender", "race", "year_of_listing",
															 "ns(age,knots=NPresults[[1]][[3]][[3]], Boundary.knots=NPresults[[1]][[3]][[4]])", 
															 "ns(height,knots=NPresults[[1]][[3]][[5]], Boundary.knots=NPresults[[1]][[3]][[6]])", 
															 "ns(weight,knots=NPresults[[1]][[3]][[7]], Boundary.knots=NPresults[[1]][[3]][[8]])",
															 "accept_incompatible_blood_type", "accept_extra_corporeal_liver", "accept_liver_segment",
															 "accept_HBV_positive_donor", "accept_HCV_positive_donor",
															 "patient_on_life_support", "functional_status", "primary_diagnosis", "spontaneous_bacterial_peritonitis",
															 "history_of_PV_thrombosis", "history_of_TIPSS"),
		time_varying_covariates_transplant = list("ns(MELD,knots=NPresults[[1]][[3]][[9]], Boundary.knots=NPresults[[1]][[3]][[10]])*
																							MELD_exception*
																							ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))"),
		time_varying_covariates_censoring = list("ns(MELD,knots=NPresults[[1]][[3]][[9]], Boundary.knots=NPresults[[1]][[3]][[10]])*
																						 MELD_exception*
																						 ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))",
																						 "post_transplant*
																						 ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))"),
		increased_risk_usage_factor = 1.5, #theta_1
		standard_risk_usage_factor = 1,  #theta_2
		use_existing_models=T, existing_model=results[[1]][["Results"]][["Models"]] #Weights for all regimes are computed using the same model
		)



################
################
#Bootstraps
################
################

bootstrap_results<-vector(2, mode="list")
names(bootstrap_results)<-c("Estimate Array", "Models")
bootstrap_results[[1]]<-array(data=NA, dim=c(500, 137, 30, 3), dimnames=list(1:500, 0:136, c(
	"Lambda_A_1",  "Lambda_A_1_g", "Lambda_A_1_gplus", "Lambda_A_1_gplus_naive",  "chi_1", "rho_1", "delta_1",   "E_A_1_KM", "E_A_1_g_KM",  "E_A_1_gplus_KM",  "E_A_1_gplus_naive_KM",
	"Lambda_A_2",  "Lambda_A_2_g", "Lambda_A_2_gplus","chi_2", "rho_2", "delta_2", "E_A_2_KM", "E_A_2_g_KM",  "E_A_2_gplus_KM", "E_A_2_gplus_naive",
	"Lambda_Y_r", "Lambda_Y_r_g", "Lambda_Y_r_naive", "Lambda_Y", "Lambda_Y_g", "Lambda_Y_naive",  "E_Y_KM",  "E_Y_g_KM",  "E_Y_naive_KM"
), c("g_1", "g_2", "g_3")))
bootstrap_results[[2]]<-vector(500, mode="list")


for(i in 1:500){
	bootstrap_results[[2]][[i]]<-vector(5, mode="list")
	names(bootstrap_results[[2]][[i]])<-c("Censoring (Admin)", "Censoring (LTFU)", "Censoring (transplant)", "Standard Risk", "Increased Risk")
}


#We run the bootstraps in batches to save progress at intermediate points.
bootin<-c(seq(1,326,25), 351:375, seq(376,476,25))
bootout<-c(seq(25,350,25), 351:375, seq(400,500,25))
regimes<-c("g_2", "g_3")
prevs<-c(1.25, 1.5)

Begin<-Sys.time()
for(i in bootstart:bootend){
	load(paste0("./Bootstraps/bootresults_",bootin[i],"_",bootout[i],".RData"))
	for(j in bootin[i]:bootout[i]){
		for(k in 1:2){
			estimate( 
				seed=j,
				data=candidates_analysis,
				baseline_covariates = list("ns(baseline_MELD,knots=NPresults[[1]][[3]][[1]], Boundary.knots=NPresults[[1]][[3]][[2]])", 
																	 "baseline_MELD_exception", "status1", "gender", "race", "year_of_listing",
																	 "ns(age,knots=NPresults[[1]][[3]][[3]], Boundary.knots=NPresults[[1]][[3]][[4]])", 
																	 "ns(height,knots=NPresults[[1]][[3]][[5]], Boundary.knots=NPresults[[1]][[3]][[6]])", 
																	 "ns(weight,knots=NPresults[[1]][[3]][[7]], Boundary.knots=NPresults[[1]][[3]][[8]])",
																	 "accept_incompatible_blood_type", "accept_extra_corporeal_liver", "accept_liver_segment",
																	 "accept_HBV_positive_donor", "accept_HCV_positive_donor",
																	 "patient_on_life_support", "functional_status", "primary_diagnosis", "spontaneous_bacterial_peritonitis",
																	 "history_of_PV_thrombosis", "history_of_TIPSS"),
				time_varying_covariates_transplant = list("ns(MELD,knots=NPresults[[1]][[3]][[9]], Boundary.knots=NPresults[[1]][[3]][[10]])*
																									MELD_exception*
																									ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))"),
				time_varying_covariates_censoring = list("ns(MELD,knots=NPresults[[1]][[3]][[9]], Boundary.knots=NPresults[[1]][[3]][[10]])*
																								 MELD_exception*
																								 ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))",
																								 "post_transplant*
																								 ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))"),
				increased_risk_usage_factor = 0,
				standard_risk_usage_factor = 1,  
				K=89, boot=T,  regime="g_1", use_existing_models = F
			)
			gc()
			estimate( 
				seed=j,
				data=candidates_analysis,
				baseline_covariates = list("ns(baseline_MELD,knots=NPresults[[1]][[3]][[1]], Boundary.knots=NPresults[[1]][[3]][[2]])", 
																	 "baseline_MELD_exception", "status1", "gender", "race", "year_of_listing",
																	 "ns(age,knots=NPresults[[1]][[3]][[3]], Boundary.knots=NPresults[[1]][[3]][[4]])", 
																	 "ns(height,knots=NPresults[[1]][[3]][[5]], Boundary.knots=NPresults[[1]][[3]][[6]])", 
																	 "ns(weight,knots=NPresults[[1]][[3]][[7]], Boundary.knots=NPresults[[1]][[3]][[8]])",
																	 "accept_incompatible_blood_type", "accept_extra_corporeal_liver", "accept_liver_segment",
																	 "accept_HBV_positive_donor", "accept_HCV_positive_donor",
																	 "patient_on_life_support", "functional_status", "primary_diagnosis", "spontaneous_bacterial_peritonitis",
																	 "history_of_PV_thrombosis", "history_of_TIPSS"),
				time_varying_covariates_transplant = list("ns(MELD,knots=NPresults[[1]][[3]][[9]], Boundary.knots=NPresults[[1]][[3]][[10]])*
																									MELD_exception*
																									ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))"),
				time_varying_covariates_censoring = list("ns(MELD,knots=NPresults[[1]][[3]][[9]], Boundary.knots=NPresults[[1]][[3]][[10]])*
																								 MELD_exception*
																								 ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))",
																								 "post_transplant*
																								 ns(days_since_start, knots = c(60,120,360,720,1620), Boundary.knots = c(30,3600))"),
				increased_risk_usage_factor = prevs[k], 
				standard_risk_usage_factor = 1,  
				K=89, boot=T,  regime=regimes[k], use_existing_models = T, existing_model=bootstrap_results[[2]][[j]]
				)
			gc()
		}
	}
	bootstrap_results[[2]]<-NULL
	save(bootstrap_results, file=paste0("./Bootstraps/bootresults_",bootin[i],"_",bootout[i],".RData"))
}



bootstrap_results_all<-bootstrap_results
for(i in 1:44){
	load(paste0("./Bootstraps/bootresults_",bootin[i],"_",bootout[i],".RData"))
	bootstrap_results_all[[1]][bootin[i]:bootout[i], , ,]<-bootstrap_results[[1]][bootin[i]:bootout[i], , ,]
	gc()
}



quantile(bootstrap_results_all[[1]][1:500, 90, "E_Y_KM", "g_1"], c(0.025, 0.975))
quantile(bootstrap_results_all[[1]][1:500, 90, "E_Y_g_KM", "g_1"], c(0.025, 0.975))
quantile(bootstrap_results_all[[1]][1:500, 90, "E_Y_g_KM", "g_1"] - bootstrap_results_all[[1]][1:500, 90, "E_Y_KM", "g_1"], c(0.025, 0.975))
quantile(bootstrap_results_all[[1]][1:500, 90, "E_Y_naive_KM", "g_1"], c(0.025, 0.975))

quantile(bootstrap_results_all[[1]][1:500, 90, "E_Y_g_KM", "g_2"], c(0.025, 0.975))
quantile(bootstrap_results_all[[1]][1:500, 90, "E_Y_g_KM", "g_2"] - bootstrap_results_all[[1]][1:500, 90, "E_Y_KM", "g_2"], c(0.025, 0.975))
quantile(bootstrap_results_all[[1]][1:500, 90, "E_Y_g_KM", "g_3"], c(0.025, 0.975))
quantile(bootstrap_results_all[[1]][1:500, 90, "E_Y_g_KM", "g_3"] - bootstrap_results_all[[1]][1:500, 90, "E_Y_KM", "g_3"], c(0.025, 0.975))



################
################
#Figures
################
################


################
#Colors
################
lightblue<-rgb( 0+ (255-0)*.9 ,0 + (255-0)*.9, 255 + (255-255)*.9 , max = 255)
lightpink<-rgb(255 + (255-255)*.9,0 + (255-0)*.9, 0+ (255-0)*.9  , max = 255)
myblue<-rgb(0 + (255-0)*.6,0 + (255-0)*.6, 255+ (255-255)*.6, max=255)
mygreenlight<-"#8ed4a1"
mygreen<-"#618f6e"
font_add(family = "ArialMSuni", regular = "/System/Library/Fonts/Supplemental/Arial Unicode.ttf")
showtext_auto()

################
#General plot function
################

plot_surv<-function(regime1, regime2, regime3, regime4=NULL, regime5=NULL, linetypes=c(1,1,1), linewidths=c(2,2,2), 
										background_color, line_colors, regimes, ylabel, years, AddLegend=T, defaultX=T, defaultY=T){
	par(family="ArialMSuni")
	regimes<-list(regime1[1:(years*12)], regime2[1:(years*12)],regime3[1:(years*12)],regime4[1:(years*12)],regime5[1:(years*12)])
	plot(x=c(0,  rep(seq(30, years*12*30, 30), each=2), years*12*30), type="n",
			 xlim=c(0,years*12*30),
			 ylim=c(0, max(unlist(regimes))),
			 xlab="Days since waitlisted", ylab=ylabel, xaxt = "n"
	)
	rect(par("usr")[1], par("usr")[3],
			 par("usr")[2], par("usr")[4],
			 col = background_color)
	for(i in 1:5){
		lines(x=c(0,  rep(seq(30, years*12*30, 30), each=2), years*12*30),
					y=c(0, 0, rep(regimes[[i]], each=2)),
					lwd=linewidths[i], lty=linetypes[i], col=line_colors[i]
		)
	}
	if(defaultX){
		axis(1, at = seq(0, years*12*30, 500))
	}
	if(defaultY){
		axis(2, at = seq(0, round(max(unlist(regimes)), digits = 1), round(max(unlist(regimes)), digits = 1)/5))
	}
	if(AddLegend){
		legend(legend=regimes,
					 lty=linetypes, lwd=linewidths-0.5,
					 col=line_colors,
					 x=years*12*30,
					 xjust=1,
					 y=max(unlist(regimes))*.1,
					 yjust=0
		)
	}
}



################
#Treatment plot
################
pdf(file="trx.pdf",width=7.5,height=10)
par(mfcol=c(2,1), oma = c(5,2.5,1,1))
{
	par(mar=c(.5,1,0,0))
	plot_surv(regime1 = cumsum(results[[1]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_A_1_gplus_KM"])      ,
						regime2 = cumsum(results[[1]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_A_1_gplus_naive_KM"]),
						regime3 = cumsum(results[[2]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_A_1_gplus_KM"])      ,
						regime4 = cumsum(results[[3]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_A_1_gplus_KM"])      ,
						regime5 = cumsum(results[[1]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_A_1_KM"])            ,
						background_color = lightblue,
						line_colors = c(myblue,  "red", mygreenlight, mygreen, "black"),
						linetypes = c(1,2,2,3, 3),
						linewidths = c(4,2,2,2,2),
						regimes = c(expression("Eliminating Increased Risk Organs, IPSI (g=g"[1]*")"),
												expression("Eliminating Increased Risk Organs, naive (g=g"[0]*")"),
												expression(paste("Natural Course (g=" , "\u2205", ")", sep = ""))),
						ylabel = expression("Cumulative incidence of transplant (standard risk), cumsum(E["*{A["1,k"]^"g+"}*"])"),
						years=7.5, AddLegend=F,  defaultX=F, defaultY=F
	)
	text(expression("A. Standard Risk, "*{A["1,k"]^"g+"}), x=par("usr")[1]*.5, y=par("usr")[4]*.98, adj=c(0,1))
	par(mar=c(1,1,0,0))
	plot_surv(regime1 = cumsum(results[[1]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_A_2_gplus_KM"])      ,
						regime2 = cumsum(results[[1]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_A_2_gplus_naive"]),
						regime3 = cumsum(results[[2]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_A_2_gplus_KM"])      ,
						regime4 = cumsum(results[[3]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_A_2_gplus_KM"])      ,
						regime5 = cumsum(results[[1]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_A_2_KM"])            ,
						background_color = lightblue,
						line_colors = c(myblue,  "red", mygreenlight, mygreen, "black"),
						linetypes = c(1,2,2,3, 3),
						linewidths = c(4,2,2,2,2),
						regimes = c(expression("Eliminating Increased Risk Organs, IPSI (g=g"[1]*")"),
												expression("Eliminating Increased Risk Organs, naive (g=g"[0]*")"),
												expression(paste("Natural Course (g=" , "\u2205", ")", sep = ""))),
						ylabel = expression("Cumulative incidence of transplant (increased risk), cumsum(E["*{A["2,k"]^"g+"}*"])"),
						years=7.5, AddLegend=F,  defaultX=T, defaultY=F
	)
	text(expression("B. Increased Risk, "*{A["2,k"]^"g+"}), x=par("usr")[1]*.5, y=par("usr")[4]*.98, adj=c(0,1))
	mtext("Days since waitlisted", side=1, line=2)
	mtext(expression("Cumulative incidence of transplant"), side=2, line=1,  outer = TRUE)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
	legend('bottom',
				 legend = c(
				 	expression("Eliminating Increased Risk Organs, IPSI (g=g"[1]*")"),
				 	expression("Increasing Increased Risk Organs +25%, IPSI (g=g"[2]*")"),
				 	expression("Increasing Increased Risk Organs +50%, IPSI (g=g"[3]*")"),
				 	expression(paste("Natural Course (g=" , "\u2205", ")", sep = "")),
				 	expression("Eliminating Increased Risk Organs, naive (g=g"[0]*")")
				 ),
				 col = c(myblue,   mygreenlight, mygreen,"black", "red" ), lwd = 3, lty=c(1,2, 3, 3,2), xpd = TRUE, ncol = 2, seg.len=3, bty = 'n', cex=0.7, text.width=0.8, x.intersp=0.25
	)
}
dev.off()
################
#Survival plot
################
pdf(file="surv.pdf",width=8,height=6.25)
par(mfcol=c(1,1), oma = c(5,2.5,1,1))
{
	par(mar=c(1,1,0,0))
	plot_surv(regime1 = results[[1]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_Y_g_KM"]          ,
						regime2 = results[[1]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_Y_naive_KM"],
						regime3 = results[[2]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_Y_g_KM"]          ,
						regime4 = results[[3]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_Y_g_KM"]          ,
						regime5 = results[[1]][["Results"]][["Non-parametric summaries (by WL)"]][, "E_Y_KM"]                  ,
						background_color = lightpink,
						line_colors = c(myblue,  "red", mygreenlight, mygreen, "black"),
						linetypes = c(1,2,1,1,1),
						linewidths = c(2,1,2,2,2),
						regimes = c(expression("Eliminating Increased Risk Organs, IPSI (g=g"[1]*")"),
												expression("Eliminating Increased Risk Organs, naive (g=g"[0]*")"),
												expression(paste("Natural Course (g=" , "\u2205", ")", sep = ""))),
						ylabel = expression("Cumulative incidence of transplant (increased risk), cumsum(E["*{A["2,k"]^"g+"}*"])"),
						years=7.5, AddLegend=F,  defaultX=T, defaultY=F
	)
	mtext("Days since waitlisted", side=1, line=2)
	mtext(expression("Cumulative incidence of death"), side=2, line=1,  outer = TRUE)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
	legend('bottom',
				 legend = c(
				 	expression("Eliminating Increased Risk Organs, IPSI (g=g"[1]*")"),
				 	expression("Increasing Increased Risk Organs +25%, IPSI (g=g"[2]*")"),
				 	expression("Increasing Increased Risk Organs +50%, IPSI (g=g"[3]*")"),
				 	expression(paste("Natural Course (g=" , "\u2205", ")", sep = "")),
				 	expression("Eliminating Increased Risk Organs, naive (g=g"[0]*")")
				 ),
				 col = c(myblue,   mygreenlight, mygreen,"black", "red" ), lwd = 3, lty=c(1,2,3,3,2), xpd = TRUE, ncol = 2, seg.len=3, bty = 'n', cex=0.7, text.width=0.8, x.intersp=0.25
	)
}
dev.off()
