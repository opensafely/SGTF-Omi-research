********************************************************************************
*
*	Do-file:		an_cox.do
*
*	Project:		SGTF Omicron
*
*	Programmed by:	Daniel Grint
*
*	Data used:		output/main.dta
*
*	Data created:	
*
*	Other output:	an_cox.log containing:
*					1-Unadjusted Cox models
*					2-Adjusted Cox models
*
********************************************************************************
*
*	Purpose:		This do-file runs Cox PH models, calculating HR for Omicron 
*					vs. Delta
*  
********************************************************************************

* Open a log file
cap log close
log using ./logs/an_cox, replace t

clear

use ./output/main.dta

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

summ cox_ae_time, d
summ cox_ae_time if cox_ae_time > 0, d


*Set up output file
cap file close tablecontent

file open tablecontent using ./output/table2.txt, write text replace

file write tablecontent ("Table 2: Hazard ratios for S-Fail vs. S-Pos") _n _n

file write tablecontent ("Estimate")	_tab ///
						("HR (95% CI)")	_tab ///
						("P-value")		_n


*******************
/* Unadjusted HR */
*******************

stset ae_surv_d, origin(study_start) fail(cox_ae) scale(1) id(patient_id)

stcox i.sgtf

* Stratified by UTLA
stcox i.sgtf, strata(utla_group)

* N (events)
tab sgtf cox_ae if e(sample)

* Output unadjusted
lincom 1.sgtf, eform
file write tablecontent _n ("Region stratified") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n


* Interaction with time
stcox i.sgtf, tvc(i.sgtf) strata(utla_group)


*******************************
/* Causal minimum adjustment */
*******************************

* Stratified by UTLA

stcox i.sgtf ib2.vax ///
			 , strata(utla_group)
			 
* N (events)
tab sgtf cox_ae if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("Vax adj.") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n


*********************************************************************
/* Demographically adjusted HR - age as spline, cat hh size 	   */
/* Not adjusting for comorbidities, obesity, or smoking	status	   */
*********************************************************************

* Stratified by UTLA
* Excluding missing ethnicity

stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.hh_total_cat i.home_bin ///
			 ib1.rural_urban5 ib49.start_week age1 age2 age3 ///
			 if eth2 != 6 ///
			 , strata(utla_group)
			 
* N (events)
tab sgtf cox_ae if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("Demographically adj.") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n


			 
*********************************************************************
/* Demographically adjusted HR - age grouped, cat hh size		   */
/* Not adjusting for comorbidities, obesity, or smoking	status	   */
*********************************************************************

* Stratified by UTLA
* Excluding missing ethnicity

stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.hh_total_cat i.home_bin ///
			 ib1.rural_urban5 ib49.start_week ib1.agegroup6 ///
			 if eth2 != 6 ///
			 , strata(utla_group)


			 
****************************************************
/* Fully adjusted HR - age as spline, cat hh size */
****************************************************

* Stratified by UTLA
* Excluding missing ethnicity

stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)
			 
est store e_no_int

* N (events)
tab sgtf cox_ae if e(sample)
bysort vax: tab sgtf cox_ae if e(sample)
bysort start_week: tab sgtf cox_ae if e(sample)
bysort comorb_cat: tab sgtf cox_ae if e(sample)
bysort eth2: tab sgtf cox_ae if e(sample)
bysort imd: tab sgtf cox_ae if e(sample)
bysort agegroup6: tab sgtf cox_ae if e(sample)

estat phtest, d


lincom 1.sgtf, eform
file write tablecontent ("Fully adj.") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n


* Plot scaled schoenfeld residuals
estat phtest, plot(1.sgtf)
graph export ./output/cox_shoen.svg, as(svg) replace


* KM plot
sts graph,	surv by(sgtf) ci risktable(, order(1 "S-Pos" 2 "S-Fail") size(small)) ///
			ylabel(0.994(0.001)1, format(%5.3f)) ///
			legend(order(2 4) label(2 "S-Pos") label(4 "S-Fail") rows(1))
graph export ./output/cox_km.svg, as(svg) replace


* Cumulative hazard plot
sts graph,	cumhaz by(sgtf) ci ///
			ylabel(minmax, format(%5.3f)) ///
			legend(order(2 4) label(2 "S-Pos") label(4 "S-Fail") rows(1))
graph export ./output/cox_cumhaz.svg, as(svg) replace

		
* Smoothed hazard plot
sts graph,	haz by(sgtf) ///
			legend(label(1 "S-Pos") label(2 "S-Fail"))
graph export ./output/cox_haz.svg, as(svg) replace


***********************
/* Subgroup analyses */
***********************

* Epi week
stcox i.sgtf##ib49.start_week ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)

est store e_weekX

* Test for interaction
lrtest e_no_int e_weekX

file write tablecontent _n ("Subgroup analyses") _n 

file write tablecontent _n ("Epi. week") _tab _tab %6.4f (r(p)) _n

* Epi week Omicron vs. Delta HR
lincom 1.sgtf, eform	// week 49
file write tablecontent ("05Dec-11Dec") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#50.start_week, eform	// week 50
file write tablecontent ("12Dec-18Dec") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#51.start_week, eform	// week 51
file write tablecontent ("19Dec-25Dec") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#52.start_week, eform	// week 52
file write tablecontent ("26Dec-01Jan") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n



* Vaccination status
stcox i.sgtf##ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)

est store e_vaxX

* Test for interaction
lrtest e_no_int e_vaxX

file write tablecontent _n ("Vaccination") _tab _tab %6.4f (r(p)) _n

* Vax status Omicron vs. Delta HR
lincom 1.sgtf + 1.sgtf#0.vax, eform			// unvax
file write tablecontent ("Unvax") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#1.vax, eform			// 1 dose
file write tablecontent ("1 dose") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf, eform						// 2 doses
file write tablecontent ("2 doses") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#3.vax, eform			// booster
file write tablecontent ("Booster") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n



* Comorbidities
stcox i.sgtf##ib0.comorb_cat ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)

est store e_comorbX

* Test for interaction
lrtest e_no_int e_comorbX

file write tablecontent _n ("Comorbidities") _tab _tab %6.4f (r(p)) _n

* Comorbidities VOC vs. non-VOC HR
lincom 1.sgtf, eform						// no comorbs
file write tablecontent ("None") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#1.comorb_cat, eform	// 1 comorb
file write tablecontent ("1") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#2.comorb_cat, eform	// 2+ comorbs
file write tablecontent ("2+") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n


/* Specific comorbidities */

* Renal
stcox i.sgtf ib2.vax age1 age2 age3 ///
			if renal_flag == 1 ///
			, strata(utla_group)
			
tab sgtf cox_ae if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("Renal v.adj") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n
			
stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 & renal_flag == 1 ///
			 , strata(utla_group)
			 
tab sgtf cox_ae if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("Renal f.adj") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

* Diabetes
stcox i.sgtf ib2.vax age1 age2 age3 ///
			if dm == 1 ///
			, strata(utla_group)
			
tab sgtf cox_ae if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("DM v.adj") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n
			
stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 & dm == 1 ///
			 , strata(utla_group)
			 
tab sgtf cox_ae if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("DM f.adj") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n
			
* CVD
stcox i.sgtf ib2.vax age1 age2 age3 ///
			if chronic_cardiac_disease == 1 ///
			, strata(utla_group)
			
tab sgtf cox_ae if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("CVD v.adj") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n
			
stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 & chronic_cardiac_disease == 1 ///
			 , strata(utla_group)
			 
tab sgtf cox_ae if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("CVD f.adj") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

/*
* Test for trend
stcox i.sgtf i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 c.comorb_cat ib1.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)

est store e_linco

stcox i.sgtf##c.comorb_cat i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib1.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)

est store e_lincoX

lrtest e_linco e_lincoX
local lin_lr_p = r(p)

lincom 1.sgtf#c.comorb_cat, eform
file write tablecontent ("Per unit increase") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (`lin_lr_p') _n

*/


* Ethnicity
stcox i.sgtf##ib1.eth2 ib2.vax i.male ib1.imd ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)

est store e_eth2X

* Test for interaction
lrtest e_no_int e_eth2X

file write tablecontent _n ("Ethnicity") _tab _tab %6.4f (r(p)) _n

* Ethnicity VOC vs. non-VOC HR
lincom 1.sgtf, eform					// White
file write tablecontent ("White") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#5.eth2, eform	// Other
file write tablecontent ("Not white") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n


/*

* IMD
stcox i.sgtf##ib1.imd i.male ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib1.start_week ib0.comorb_cat age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)

est store e_imdX

* Test for interaction
lrtest e_no_int e_imdX

file write tablecontent _n ("IMD") _tab _tab %6.4f (r(p)) _n

* IMD VOC vs. non-VOC HR
lincom 1.sgtf, eform						// 1
file write tablecontent ("1 Least deprived") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#2.imd, eform	// 2
file write tablecontent ("2") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#3.imd, eform	// 3
file write tablecontent ("3") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#4.imd, eform	// 4
file write tablecontent ("4") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#5.imd, eform	// 5
file write tablecontent ("5 Most deprived") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

*/

* Age group
stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week ib0.agegroup6 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)
			 
est store e_age

* N (events)
tab sgtf agegroup6 if e(sample)

estat phtest, d


stcox i.sgtf##ib0.agegroup6 ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)
			 
est store e_ageX

* Test for interaction
lrtest e_age e_ageX

file write tablecontent _n ("Age group") _tab _tab %6.4f (r(p)) _n

* Age group VOC vs. non-VOC HR
lincom 1.sgtf, eform						// 0-39
file write tablecontent ("0-39") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#1.agegroup6, eform	// 40-54
file write tablecontent ("40-54") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#2.agegroup6, eform	// 55-64
file write tablecontent ("55-64") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#3.agegroup6, eform	// 65-74
file write tablecontent ("65-74") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#4.agegroup6, eform	// 75-84
file write tablecontent ("75-84") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

lincom 1.sgtf + 1.sgtf#5.agegroup6, eform	// 85+
file write tablecontent ("85+") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

/*

* Test for trend
stcox i.sgtf c.agegroupA i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib1.start_week i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)
			 
est store e_cage

stcox i.sgtf##c.agegroupA i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib1.start_week i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)
			 
est store e_cageX

lrtest e_cage e_cageX
local lin_age_p = r(p)

lincom 1.sgtf#c.agegroupA, eform
file write tablecontent ("Per group increase") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (`lin_age_p') _n

*/


**************************
/* Sensitivity analyses */
**************************

file write tablecontent _n ("Sensitivity analyses") _n

* Excluding care home status

stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 ///
			 if eth2 != 6 ///
			 , strata(utla_group)
			 
* N (events)
tab sgtf cox_ae if e(sample)
			 
lincom 1.sgtf, eform
file write tablecontent ("Excluding care home") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

* Excluding ethnicity

stcox i.sgtf ib2.vax i.male ib1.imd ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 ///
			 , strata(utla_group)
			 
* N (events)
tab sgtf cox_ae if e(sample)
			 
lincom 1.sgtf, eform
file write tablecontent ("Excluding ethnicity") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

* Min 14-days pre-EC data censor

stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 i.home_bin ///
			 if ae_14_pop == 1 & eth2 != 6 ///
			 , strata(utla_group)
			 
* N (events)
tab sgtf cox_ae if e(sample)
			 
lincom 1.sgtf, eform
file write tablecontent ("Min 14-days FU") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

* Add 1 to FU

stset ae_surv_d1, origin(study_start) fail(cox_ae) scale(1) id(patient_id)

stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)
			 
* N (events)
tab sgtf cox_ae if e(sample)
			 
lincom 1.sgtf, eform
file write tablecontent ("FU plus 1") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n
			 

**************************
/* Admission as outcome */
**************************

stset ae_surv_d, origin(study_start) fail(cox_admit) scale(1) id(patient_id)

* Stratified by UTLA
* Excluding missing ethnicity

stcox i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			 , strata(utla_group)

* N (events)
tab sgtf cox_admit if e(sample)
bysort vax: tab sgtf cox_admit if e(sample)
bysort start_week: tab sgtf cox_admit if e(sample)
bysort comorb_cat: tab sgtf cox_admit if e(sample)
bysort eth2: tab sgtf cox_admit if e(sample)
bysort imd: tab sgtf cox_admit if e(sample)
bysort agegroup6: tab sgtf cox_admit if e(sample)

estat phtest, d


lincom 1.sgtf, eform
file write tablecontent ("Fully adj. Admission") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n


/*

* 28-days follow-up censor
stset stime_death28, origin(study_start) fail(cox_death28) scale(1) id(patient_id)


* Stratified by region
stcox i.sgtf i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib1.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 & risk_pop==1 ///
			 , strata(utla_group)
			 
* N (events)
tab sgtf cox_death if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("28-day follow-up censor") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n

tab risk_pop risk_pop_40


stset stime_death, origin(study_start) fail(cox_death) scale(1) id(patient_id)

* Include with 40-days follow-up
* Stratified by region
stcox i.sgtf i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib1.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 & risk_pop_40==1 ///
			 , strata(utla_group)
			 
* N (events)
tab sgtf cox_death if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("Min. 40-days follow-up") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n


* No adjustment for care home
* Stratified by region
stcox i.sgtf i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib2.start_week age1 age2 age3 ///
			 if eth2 != 6 ///
			 , strata(utla_group)
			 
* N (events)
tab sgtf cox_death if e(sample)

lincom 1.sgtf, eform
file write tablecontent ("No care home adj.") _tab 
file write tablecontent %4.2f (r(estimate)) (" (") %4.2f (r(lb)) ("-") %4.2f (r(ub)) (")") _tab %6.4f (r(p)) _n



*/


* Close output table
file write tablecontent _n _n
file close tablecontent
			 

			 
log close


clear

insheet using ./output/table2.txt, clear
export excel using ./output/table2.xlsx, replace

