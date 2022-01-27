********************************************************************************
*
*	Do-file:		an_cox_imputed.do
*
*	Project:		SGTF Omicron
*
*	Programmed by:	Daniel Grint
*
*	Data used:		output/main_imputed.dta
*
*	Data created:	output/an_imputed_eth2
*
*	Other output:	an_cox_imputed.log
*
*
********************************************************************************
*
*	Purpose:		This do-file imputes missing ethnicity data
*  
********************************************************************************

* Open a log file
cap log close
log using ./logs/an_cox_imputed, replace t

clear


use ./output/main_imputed.dta

* DROP MISSING UTLA
noi di "DROPPING MISSING UTLA DATA"
drop if utla_group==""

* DROP IF NO DATA ON SGTF
noi di "DROPPING NO SGTF DATA" 
drop if has_sgtf==0

noi di "SUBSETTING ON COX CENSORED POPULATION"
keep if cox_pop==1

noi di "SUBSETTING ON START/STOP DATES"
keep if inrange(start_week,49,52)

tab sgtf cox_ae, row
tab sgtf cox_ae if ae_surv_d > study_start, row

* Declare survival data
mi stset ae_surv_d, origin(study_start) fail(cox_ae) scale(1) id(patient_id)


* Stratified by region
mi estimate, eform: stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 i.home_bin ///
			 , strata(utla_group)
			 
estimates save ./output/an_imputed_eth2, replace


* N (events)
tab sgtf cox_ae if e(sample)



log close
