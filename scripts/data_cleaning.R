#read in the raw data
cand_liin <- as.data.frame(read_sas("./data/cand_liin.sas7bdat"))
stathist_liin <- as.data.frame(read_sas("./data/stathist_liin.sas7bdat"))
tx_li <- as.data.frame(read_sas("./data/tx_li.sas7bdat"))
donor_deceased <- as.data.frame(read_sas("./data/donor_deceased.sas7bdat"))


#clean the 3 databases that we need to link
cand_liin <- cand_liin %>%
	mutate(
		status_id = PX_ID,
		patient_id = PERS_ID,
		death_date = case_when(
			is.na(PERS_SSA_DEATH_DT) == F ~ as.Date(PERS_SSA_DEATH_DT, format =
																								"%m/%d/%y"),
			is.na(PERS_OPTN_DEATH_DT) == F ~ as.Date(PERS_OPTN_DEATH_DT, format =
																							 	"%m/%d/%y")
		),
		age = as.numeric(CAN_AGE_IN_MONTHS_AT_LISTING) / 12,
		race = as.factor(
			case_when(
				CAN_RACE == 128 ~ "Pacific Islander",
				CAN_RACE == 2000 ~ "Hispanic/Latino",
				CAN_RACE == 64 ~ "Asian",
				CAN_RACE == 16 ~ "Black/African-American",
				CAN_RACE == 32 ~ "Native American",
				CAN_RACE == 8 ~ "White",
				CAN_RACE_SRTR == "MULTI" ~ "Multi-Racial",
				TRUE ~ NA_character_
			)
		),
		gender = as.factor(
			case_when(
				CAN_GENDER == "F" ~ "F",
				CAN_GENDER == "M" ~ "M",
				TRUE ~ NA_character_
			)
		),
		height = as.numeric(CAN_HGT_CM),
		weight = as.numeric(CAN_WGT_KG),
		accept_incompatible_blood_type = as.factor(
			case_when(
				CAN_ACPT_ABO_INCOMP == "Y" ~ "Y",
				CAN_ACPT_ABO_INCOMP == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		accept_a2_donor = as.factor(
			case_when(
				CAN_ACPT_A2_DON == "Y" ~ "Y",
				CAN_ACPT_A2_DON == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		accept_extra_corporeal_liver = as.factor(
			case_when(
				CAN_ACPT_EXTRACORP_LI == "Y" ~ "Y",
				CAN_ACPT_EXTRACORP_LI == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		accept_liver_segment = as.factor(
			case_when(
				CAN_ACPT_LI_SEG == "Y" ~ "Y",
				CAN_ACPT_LI_SEG == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		accept_HBV_positive_donor = as.factor(
			case_when(
				CAN_ACPT_HBC_POS == "Y" ~ "Y",
				CAN_ACPT_HBC_POS == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		accept_HCV_positive_donor = as.factor(
			case_when(
				CAN_ACPT_HCV_POS == "Y" ~ "Y",
				CAN_ACPT_HCV_POS == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		patient_state = as.factor(
			case_when(
				CAN_PERM_STATE == "ZZ: UNKNOWN" ~ NA_character_,
				CAN_PERM_STATE %in% c(
					"AS",
					"GU",
					"MP",
					"PR",
					"VI"
				) == T ~ "U.S territory",
				CAN_PERM_STATE == "NA" ~ "Non-U.S",
				TRUE ~ as.character(CAN_PERM_STATE)
			)
		),
		patient_educational_status = as.factor(
			case_when(
				CAN_EDUCATION == 1 ~ "none",
				CAN_EDUCATION == 2 ~ "0-8",
				CAN_EDUCATION == 3 ~ "9-12",
				CAN_EDUCATION == 4 ~ "some college / technical school",
				CAN_EDUCATION == 5 ~ "associate/bachelor degree",
				CAN_EDUCATION == 6 ~ "post-college / graduate degree",
				TRUE ~ NA_character_
			)
		),
		medical_condition = as.factor(
			case_when(
				CAN_MED_COND == 1 ~ "in ICU",
				CAN_MED_COND == 2 ~ "hospitalized not in ICU",
				CAN_MED_COND == 3 ~ "not hospitalized",
				TRUE ~ NA_character_
			)
		),
		patient_on_life_support = as.factor(
			case_when(
				CAN_LIFE_SUPPORT == "N" ~ "N",
				CAN_LIFE_SUPPORT == "Y" ~ "Y",
				TRUE ~ NA_character_
			)
		),
		functional_status = as.factor(
			case_when(
				CAN_FUNCTN_STAT    ==  1    |
					CAN_FUNCTN_STAT ==  2100 |
					CAN_FUNCTN_STAT ==  2090 |
					CAN_FUNCTN_STAT ==  2080 |
					CAN_FUNCTN_STAT ==  2070 ~ "no assistance",
				CAN_FUNCTN_STAT ==     2    |
					CAN_FUNCTN_STAT ==  2060 |
					CAN_FUNCTN_STAT ==  2050 |
					CAN_FUNCTN_STAT ==  2040 |
					CAN_FUNCTN_STAT ==  2030 ~ "some assistance",
				CAN_FUNCTN_STAT ==     3   |
					CAN_FUNCTN_STAT ==  2020|
					CAN_FUNCTN_STAT ==  2010~ "total assistance",
				TRUE ~ NA_character_
			)
		),
		primary_diagnosis = as.factor(
			case_when(
				CAN_DGN %in% c(
					999 ,
					4593,
					4592,
					4520,
					4510,
					4500,
					4455,
					4451,
					4450,
					4285,
					4290,
					4597
				) == T ~ "Other",
				CAN_DGN %in% c(
					4430,
					4420,
					4410,
					4405,
					4404,
					4403,
					4402,
					4401,
					4400
				) == T ~ "Malignant neoplasm",
				CAN_DGN %in% c(
					4315,
					4308,
					4307,
					4306,
					4305,
					4304,
					4303,
					4302,
					4301,
					4300
				) == T ~ "Metabolic",
				CAN_DGN %in% c(
					4275,
					4272,
					4271,
					4270,  
					4265,
					4264,
					4260,
					4255,
					4250,
					4245,
					4242,
					4241,
					4240,
					4235,
					4231,
					4230,
					4220
				) == T ~ "Cholestatic"                                                        ,
				CAN_DGN %in% c(
					4217,
					4216,
					4215,
					4214,
					4213,
					4212,
					4210,
					4209,
					4208,
					4207,
					4206,
					4205,
					4204,
					4203,
					4202,
					4201,
					4200,
					4280
				) == T ~ "Non-cholestatic",
				CAN_DGN %in% c(
					4110,
					4108,
					4107,
					4106,
					4105,
					4104,
					4102,
					4101,
					4100
				) == T ~ "Fulminant hepatic failure",
				TRUE ~ NA_character_
			)
		),
		diabetes = as.factor(
			case_when(
				CAN_DIAB %in% c(2:4, 998) == T ~ "Y",
				CAN_DIAB == 1 ~ "N",
				TRUE ~ NA_character_
			)
		),
		history_of_malignancy = as.factor(
			case_when(
				CAN_MALIG == "Y" ~ "Y",
				CAN_MALIG == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		coronary_artery_disease = as.factor(
			case_when(
				CAN_ANGINA %in% c(
					30,
					6 ,
					7 
				) == T ~ "Y",
				CAN_ANGINA == 1 ~ "N",
				TRUE ~ NA_character_
			)
		),
		hypertension = as.factor(
			case_when(
				CAN_DRUG_TREAT_HYPERTEN == "Y" ~ "Y",
				CAN_DRUG_TREAT_HYPERTEN == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		COPD = as.factor(
			case_when(
				CAN_DRUG_TREAT_COPD == "Y" ~ "Y",
				CAN_DRUG_TREAT_COPD == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		symptomatic_cerebrovascular_disease = as.factor(
			case_when(
				CAN_CEREB_VASC == "Y" ~ "Y",
				CAN_CEREB_VASC == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		symptomatic_peripheral_vascular_disease = as.factor(
			case_when(
				CAN_PERIPH_VASC == "Y" ~ "Y",
				CAN_PERIPH_VASC == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		pulmonary_embolism = as.factor(
			case_when(
				CAN_PULM_EMBOL == "Y" ~ "Y",
				CAN_PULM_EMBOL == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		spontaneous_bacterial_peritonitis = as.factor(
			case_when(
				CAN_BACTERIA_PERIT == "Y" ~ "Y",
				CAN_BACTERIA_PERIT == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		history_of_PV_thrombosis = as.factor(
			case_when(
				CAN_PORTAL_VEIN == "Y" ~ "Y",
				CAN_PORTAL_VEIN == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		history_of_TIPSS = as.factor(
			case_when(
				CAN_TIPSS == "Y" ~ "Y",
				CAN_TIPSS == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		variceal_bleeding = as.factor(
			case_when(
				CAN_VARICEAL_BLEEDING == "Y" ~ "Y",
				CAN_VARICEAL_BLEEDING == "N" ~ "N",
				TRUE ~ NA_character_
			)
		)
	) %>% mutate(death = as.numeric(case_when(is.na(death_date) == F  ~ 1,
																						TRUE ~ 0))) %>%
	filter(age >= 18) %>%
	filter(WL_ORG == "LI") %>%
	filter(CAN_PREV_LI == 0)

stathist_liin <- stathist_liin %>%
	mutate(
		status_id = PX_ID,
		MELDcat = as.factor(
			case_when(
				CANHX_SRTR_LAB_MELD >= 6235                              ~ ">= 35",
				CANHX_SRTR_LAB_MELD >= 6230 & CANHX_SRTR_LAB_MELD < 6235 ~ "30-35",
				CANHX_SRTR_LAB_MELD >= 6225 & CANHX_SRTR_LAB_MELD < 6230 ~ "25-29",
				CANHX_SRTR_LAB_MELD >= 6220 & CANHX_SRTR_LAB_MELD < 6225 ~ "20-24",
				CANHX_SRTR_LAB_MELD >= 6215 & CANHX_SRTR_LAB_MELD < 6220 ~ "15-19",
				CANHX_SRTR_LAB_MELD <  6215                              ~ "< 15" ,
				TRUE ~ NA_character_
			)
		),
		MELD = CANHX_SRTR_LAB_MELD - 6200,
		MELD_exception = as.factor(
			case_when(
				CANHX_EXC_FLG == 1 ~ "Y",
				CANHX_EXC_FLG == 0 ~ "N",
				TRUE ~ NA_character_
			)
		), 
		status1 = as.factor(
			case_when(
				CANHX_STAT_CD %in% c(6010:6012, 3010) == T ~ "status1",
				TRUE ~ "other status"
			)
		),
		INR     = CANHX_INR,
		albumin = CANHX_ALBUMIN,
		ascites = as.factor(
			case_when(
				CANHX_ASCITES == 1 ~ "absent",
				CANHX_ASCITES == 2 ~ "slight",
				CANHX_ASCITES == 3 ~ "moderate",
				TRUE ~ NA_character_
			)
		),
		bilirubin  = CANHX_BILI,
		creatinine = CANHX_SERUM_CREAT,
		dialysis_past_week = as.factor(
			case_when(
				CAN_LAST_DIAL_PRIOR_WEEK == "Y" ~ "Y",
				CAN_LAST_DIAL_PRIOR_WEEK == "N" ~ "N",
				TRUE ~ NA_character_
			)
		),
		encephalopathy = as.factor(
			case_when(
				CAN_LAST_ENCEPH == 1 ~ "none",
				CAN_LAST_ENCEPH == 2 ~ "1-2",
				CAN_LAST_ENCEPH == 3 ~ "3-4",
				TRUE ~ NA_character_
			)
		),
		sodium = CANHX_SERUM_SODIUM,
		height = CANHX_HGT_CM,
		weight = CANHX_WGT_KG,
		record_start = as.Date(CANHX_BEGIN_DT_TM, format = "%d%b%y"),
		record_end   = as.Date(CANHX_END_DT_TM  , format = "%d%b%y"),
		listing_date = as.Date(CAN_LISTING_DT   , format = "%m/%d/%y")
	) %>%
	filter(WL_ORG == "LI") %>%
	filter(status_id %in% cand_liin$status_id)

tx_li <-
	join(
		tx_li %>% filter(CAN_PREV_LI == 0) %>% dplyr::select(-DON_ANTI_HCV),
		donor_deceased %>% dplyr::select(
			c(
				DONOR_ID,
				DON_ANTI_HIV,
				DON_HBV_SURF_ANTIGEN,
				DON_ANTI_HCV
			)
		),
		by = "DONOR_ID"
	) %>%
	mutate(
		patient_id = PERS_ID,
		transplant_date = as.Date(REC_TX_DT, format = "%m/%d/%y"),
		last_graft_follow_up_date = as.Date(TFL_LAFUDATE, format =
																					"%m/%d/%y"),
		non_heart_beating_donor = as.factor(
			case_when(
				DON_NON_HR_BEAT == "N" ~ "N",
				DON_NON_HR_BEAT == "Y" ~ "Y",
				TRUE ~ NA_character_
			)
		),
		deceased_donor = as.factor(
			case_when(
				DON_TY == "C" ~ "deceased",
				DON_TY == "L" ~ "living",
				TRUE ~ NA_character_
			)
		),
		CDC_increased_risk =
			case_when(
				DON_MEET_CDC_HIGH_RISK == "Y" ~ 1,
				DON_MEET_CDC_HIGH_RISK == "N" ~ 0,
				TRUE ~ NA_real_
			), 
		HIV_status_of_donor = as.factor(
			case_when(
				DON_ANTI_HIV == "N" ~ "negative",
				DON_ANTI_HIV %in% c(
					"C",
					"I",
					"ND",
					"P",
					"U"
				) ~ "positive or indeterminate",
				TRUE ~ NA_character_
			)
		),
		HBV_status_of_donor = as.factor(
			case_when(
				DON_HBV_SURF_ANTIGEN == "N" ~ "negative",
				DON_HBV_SURF_ANTIGEN %in% c(
					"C",
					"I",
					"ND",
					"P",
					"U"
				) ~ "positive or indeterminate",
				TRUE ~ NA_character_
			)
		),
		HCV_status_of_donor = as.factor(
			case_when(
				DON_ANTI_HCV == "N" ~ "negative",
				DON_ANTI_HCV %in% c(
					"C",
					"I",
					"ND",
					"PD",
					"P",
					"U"
				) ~ "positive or indeterminate",
				TRUE ~ NA_character_
			)
		)
	) %>%
	filter(patient_id %in% cand_liin$patient_id)

#join the 3 databases into a combined database
candidates <- join(stathist_liin %>% dplyr::select(
	c(
		MELD,
		MELD_exception,
		status1,
		record_start,
		record_end,
		listing_date,
		status_id
	)
),
cand_liin %>% dplyr::select(
	c(
		age,
		race,
		gender,
		height,
		weight,
		accept_incompatible_blood_type,
		accept_extra_corporeal_liver,
		accept_liver_segment,
		accept_HBV_positive_donor,
		accept_HCV_positive_donor,
		patient_on_life_support,
		functional_status,
		primary_diagnosis,
		spontaneous_bacterial_peritonitis,
		history_of_PV_thrombosis,
		history_of_TIPSS,
		CAN_REM_CD,
		death,
		death_date,
		status_id,
		patient_id
	)
),
by = "status_id") %>% join(., tx_li %>% dplyr::select(
	c(
		non_heart_beating_donor,
		deceased_donor,
		CDC_increased_risk,
		HIV_status_of_donor,
		HBV_status_of_donor,
		HCV_status_of_donor,
		patient_id,
		transplant_date,
		last_graft_follow_up_date
	)
), by = "patient_id") %>% group_by(patient_id) %>%
	mutate(listing_date = min(listing_date)) %>% 
	filter(  year(listing_date) >= 2005 &
					 	year(listing_date) <= 2015) %>%
	filter(record_start <= as.Date("05/31/16", format = "%m/%d/%y")) %>%
	mutate(record_end = case_when(
		is.na(record_end) == T ~ as.Date("05/31/16", format = "%m/%d/%y"),
		TRUE ~ pmin(record_end, as.Date("05/31/16", format = "%m/%d/%y"))
	)) %>% 
	mutate(
		age = min(age),
		max_record_end = max(record_end),
		reason_record_end = CAN_REM_CD[which(record_end == max(record_end))[1]],
		transplanted = as.numeric(case_when(is.na(transplant_date) == F ~ 1, TRUE ~ 0)),
		censored_at_transplant = case_when((HIV_status_of_donor == "negative" &
																					HCV_status_of_donor == "negative" &
																					HBV_status_of_donor == "negative" &
																					non_heart_beating_donor == "N" &
																					deceased_donor == "deceased" &
																					is.na(CDC_increased_risk) == F) ~ 0,
																			 is.na(transplant_date) == F ~ 1),
		censored_at_transplant_forDCD = case_when((HIV_status_of_donor == "negative" &
																							 	HCV_status_of_donor == "negative" &
																							 	HBV_status_of_donor == "negative" &
																							 	deceased_donor == "deceased" &
																							 	is.na(non_heart_beating_donor) == F) ~ 0,
																							is.na(transplant_date) == F ~ 1)
	) %>%
	mutate(
		last_follow_up = case_when(
			death == 1 ~ min(death_date, as.Date("05/31/16", format = "%m/%d/%y")),
			transplanted == 1 &
				death == 0 ~ min(last_graft_follow_up_date, as.Date("05/31/16", format = "%m/%d/%y")),
			transplanted == 0 &
				death == 0 ~ max_record_end
		)
	) %>%
	mutate(
		censored = case_when(death != 1 & last_follow_up < as.Date("05/31/16", format = "%m/%d/%y") ~ 1,
												 TRUE ~ 0)
	) %>% 
	{
		.[order(.$patient_id, .$record_start), ]
	} %>% ungroup() %>%
	mutate(year_of_listing    = as.factor(year(listing_date   )),
				 year_of_transplant = as.factor(year(transplant_date)))

#remove duplicated records
candidates <- bind_rows(
	candidates[candidates$patient_id %in% unique((candidates %>% group_by(patient_id) %>% filter(duplicated(record_start) ==
																																															 	T))$patient_id),] %>%
		group_by(patient_id, record_start) %>% filter((
			is.na(MELD) + is.na(status1) + is.na(MELD_exception)
		)  ==
			min((
				is.na(MELD) + is.na(status1)  + is.na(MELD_exception)
			))) %>%
		filter(duplicated(record_start) == F),
	candidates[!candidates$patient_id %in% unique((candidates %>% group_by(patient_id) %>% filter(duplicated(record_start) ==
																																																	T))$patient_id), ]
) %>% filter(record_start <= transplant_date |
						 	is.na(transplant_date) == T) %>%
						 	{
						 		.[order(.$patient_id, .$record_start), ]
						 	} %>% ungroup()

#carry forward covariate information
candidates <-
	candidates %>% group_by(patient_id) %>% fill() %>% ungroup()

#save candidates file
save(candidates, file="/Users/sarvet/Dropbox/EPFL Postdoc/Limited Resources/Paper C/Code/SRTR data/pubsaf1909/candidates.RData")
	#load("/Users/sarvet/Dropbox/EPFL Postdoc/Limited Resources/Paper C/Code/SRTR data/pubsaf1909/candidates.RData")



#create expanded file
expand<-function(x) {
	temp <- candidates %>% filter(patient_id == x)
	temp %>% mutate(next_date = lead(record_start, 1)) %>%
		mutate(next_date = case_when(is.na(next_date) == T ~ last_follow_up,
																 TRUE ~ next_date)) %>%
		mutate(days = map2(record_start, next_date, `:`)) %>%
		unnest() %>%  group_by(days) %>% filter(row_number() == n()) %>% ungroup() %>%
		mutate(days_since_start = as.numeric(seq_along(1:(n()))) - 1) %>% filter(
			days_since_start == 0 |
				days_since_start == max(days_since_start) |
				days_since_start %% 30 == 0
		) %>% mutate(
			date = min(record_start) + days_since_start,
			#Approximate everyone's record_start by setting to the max date less than their true record_start in a sequence of dates 30 days apart starting from "01/01/05"
			record_start_approx = as.Date(
				as.numeric(record_start)  - (
					(as.numeric(record_start) - as.numeric(as.Date("01/01/05", format = "%m/%d/%y"))) %% 30
				)
				, origin = "1970-01-01"
			) ,
			days_since_start_cal = (as.numeric(min(record_start_approx)) - as.numeric(as.Date("01/01/05", format = "%m/%d/%y"))) + days_since_start,
			date_approx = as.Date("01/01/05", format = "%m/%d/%y") + days_since_start_cal,
			interval_WL  = ceiling(days_since_start/30),
			interval_cal = ceiling(days_since_start_cal/30)
		)
}



#create expanded file
candidates_expanded <-
	bind_rows(pblapply(unique(candidates$patient_id), expand))

#save expanded candidates file
save(candidates_expanded, file="./data/candidates_expanded.RData")
	#load(file="./data/candidates_expanded.RData")

#place censoring and death variables
candidates_analysis <-
	candidates_expanded %>% group_by(patient_id) %>%  mutate(
		death_event = case_when(
			death_date <= last_follow_up &
				(
					(death_date == lead(date) & lead(date) == max(date) & (lead(date)-date != 30)) |
						(death_date == date & date == max(date) & (date-lag(date) == 30)) |
						(death_date >= date & death_date < lead(date)) |
						(death_date < min(date) & date == min(date)) |
						(n() == 1 & death == 1)
				) ~ 1,
			TRUE ~ 0
		),
		censoring_event = case_when(censored == 1 &
																	((last_follow_up == lead(date) & lead(date) == max(date) & (lead(date)-date != 30)) |
																	 	(last_follow_up == date & date == max(date) & (date-lag(date) == 30))) ~ 1,
																TRUE ~ 0),
		transplant_event = case_when(
			transplant_date <= last_follow_up &
				(
					(transplant_date == lead(date) & lead(date) == max(date) & (lead(date)-date != 30)) |
						(transplant_date == date & date == max(date) & (date-lag(date) == 30)) |
						(transplant_date >= date &
						 	transplant_date < lead(date)) |
						transplant_date < min(date) &
						date == min(date) | (n() == 1 & transplanted == 1)
				) ~ 1,
			TRUE ~ 0
		),
		post_transplant = case_when(transplant_date < date ~ 1,
																TRUE ~ 0),
		time_since_transplant = as.numeric(case_when(
			is.na(transplant_date) == F ~ pmax(date - transplant_date, 0),
			TRUE ~ 0
		))
	) %>% ungroup() %>% mutate(
		transplant_and_increased_risk = case_when(transplant_event == 1 & censored_at_transplant == 0 & CDC_increased_risk == 1 ~ 1,
																							TRUE ~ 0),
		transplant_and_standard_risk = case_when(transplant_event == 1 & censored_at_transplant == 0 & CDC_increased_risk == 0 ~ 1,
																						 TRUE ~ 0),
		transplant_and_DCD = case_when(transplant_event == 1 & censored_at_transplant == 0 & non_heart_beating_donor == 1 ~ 1,
																	 TRUE ~ 0),
		transplant_and_nonDCD = case_when(transplant_event == 1 & censored_at_transplant == 0 & non_heart_beating_donor == 0 ~ 1,
																			TRUE ~ 0),
		transplant_and_censored_at_tx = case_when(transplant_event == 1 & censored_at_transplant == 1 ~ 1,
																							TRUE ~ 0),
		transplant_and_censored_at_tx_forDCD = case_when(transplant_event == 1 & censored_at_transplant_forDCD == 1 ~ 1,
																										 TRUE ~ 0))

#remove rows after death/censoring]
candidates_analysis <-
	candidates_analysis %>% group_by(patient_id) %>%
	mutate(
		filtervardeath = case_when(death_event == 1 ~ 0, TRUE ~ cumsum(death_event)),
		filtervarcensor = case_when(censoring_event == 1 ~ 0, TRUE ~ cumsum(censoring_event)),
		filtervartransplantcensor = case_when(transplant_and_censored_at_tx == 1 ~ 0, TRUE ~ cumsum(transplant_and_censored_at_tx))
	) %>%
	filter(filtervardeath == 0 &
				 	filtervarcensor == 0 &
				 	filtervartransplantcensor == 0) %>% ungroup() %>% select(-c(filtervardeath, filtervarcensor, filtervartransplantcensor))

#filter person-visits which are administratively censored and persons added and removed from waitlist on the same day (n=76)
candidates_analysis <-
	candidates_analysis %>% 
	filter((date + 30) <= as.Date("05/31/16", format = "%m/%d/%y") & interval_cal<=136 & !(listing_date==last_follow_up)) 

#select variables for analysis file
candidates_analysis <-
	candidates_analysis %>% dplyr::select(
		c(
			patient_id,
			status_id,
			interval_WL, 
			interval_cal,
			year_of_listing,
			date,
			days_since_start,
			last_follow_up,
			reason_record_end,
			censored,
			censoring_event,
			death,
			death_event,
			death_date,
			transplanted,
			censored_at_transplant,
			transplant_event,
			transplant_and_increased_risk,
			transplant_and_standard_risk,
			transplant_and_censored_at_tx,
			transplant_date,
			post_transplant,
			time_since_transplant,
			CDC_increased_risk,
			year_of_transplant,
			MELD,
			MELD_exception,
			status1,
			gender,
			age,
			race,
			height,
			weight,
			accept_incompatible_blood_type,
			accept_extra_corporeal_liver,
			accept_liver_segment,
			accept_HBV_positive_donor,
			accept_HCV_positive_donor,
			patient_on_life_support,
			functional_status,
			primary_diagnosis,
			spontaneous_bacterial_peritonitis,
			history_of_PV_thrombosis,
			history_of_TIPSS
		)
	)


ids_with_missing_baseline <- 
	candidates_analysis %>% 
	filter(days_since_start==0) %>% 
	{.[which(
		is.na(.$MELD) | 
			is.na(.$MELD_exception) |
			is.na(.$status1) |
			is.na(.$gender) |
			is.na(.$age) |
			is.na(.$race) |
			is.na(.$height) |
			is.na(.$weight) |
			is.na(.$accept_incompatible_blood_type) |
			is.na(.$accept_extra_corporeal_liver) |
			is.na(.$accept_liver_segment) |
			is.na(.$accept_HBV_positive_donor) |
			is.na(.$accept_HCV_positive_donor) |
			is.na(.$patient_on_life_support) |
			is.na(.$functional_status) |
			is.na(.$primary_diagnosis) |
			is.na(.$spontaneous_bacterial_peritonitis) |
			is.na(.$history_of_PV_thrombosis) |
			is.na(.$history_of_TIPSS)
	),]$patient_id}

candidates_analysis <- 
	candidates_analysis %>% filter(!patient_id %in% ids_with_missing_baseline)

#fix non-time-varying covariates at baseline level
candidates_analysis <-
	candidates_analysis %>% group_by(patient_id) %>%
	mutate(
		baseline_MELD = first(MELD),
		baseline_MELD_exception = first(MELD_exception),
		status1 = first(status1),
		gender = first(gender),
		age = first(age),
		race = first(race),
		height = first(height),
		weight = first(weight),
		accept_incompatible_blood_type = first(accept_incompatible_blood_type),
		accept_extra_corporeal_liver = first(accept_extra_corporeal_liver),
		accept_liver_segment = first(accept_liver_segment),
		accept_HBV_positive_donor = first(accept_HBV_positive_donor),
		accept_HCV_positive_donor = first(accept_HCV_positive_donor),
		patient_on_life_support = first(patient_on_life_support),
		functional_status = first(functional_status),
		primary_diagnosis = first(primary_diagnosis),
		spontaneous_bacterial_peritonitis = first(spontaneous_bacterial_peritonitis),
		history_of_PV_thrombosis = first(history_of_PV_thrombosis),
		history_of_TIPSS = first(history_of_TIPSS)
	) %>% ungroup()


#Patient-level variables for event failure during follow-up
candidates_analysis <-
	candidates_analysis %>% mutate(
		censoring_event_admin = case_when(interval_cal==max(candidates_analysis$interval_cal) ~ 1,
																			TRUE ~0),
		censoring_event = case_when(censoring_event_admin == 1  ~ 0,
																TRUE ~ censoring_event),
		transplant_and_censored_at_tx = case_when(censoring_event_admin == 1 | censoring_event == 1  ~ 0,
																							TRUE ~ transplant_and_censored_at_tx),
		transplant_and_standard_risk= case_when(censoring_event_admin == 1 | censoring_event == 1  ~ 0,
																						TRUE ~ transplant_and_standard_risk),
		transplant_and_increased_risk= case_when(censoring_event_admin == 1 | censoring_event == 1  ~ 0,
																						 TRUE ~ transplant_and_increased_risk),
		death_event = case_when(censoring_event_admin == 1 | censoring_event == 1 | transplant_and_censored_at_tx==1 ~ 0,
														TRUE ~ death_event)
	) %>% group_by(patient_id) %>%
	mutate(
		censored_admin = sum(censoring_event_admin),
		censored=sum(censoring_event),
		censored_at_transplant=sum(transplant_and_censored_at_tx),
		transplanted_and_increased_risk=sum(transplant_and_increased_risk),
		transplanted_and_standard_risk=sum(transplant_and_standard_risk), 
		transplanted = transplanted_and_increased_risk + transplanted_and_standard_risk,
		death=sum(death_event)
	) %>% ungroup


#free up memory after cleaning
rm(stathist_liin, cand_liin, tx_li, donor_deceased, candidates, candidates_expanded, ids_with_missing_baseline)
gc()

#save analysis file
save(candidates_analysis, file="./data/candidates_analysis.RData")
	#load("/Users/sarvet/Dropbox/EPFL Postdoc/Limited Resources/Paper C/Code/SRTR data/pubsaf1909/candidates_analysis.RData")