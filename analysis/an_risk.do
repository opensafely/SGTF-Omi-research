********************************************************************************
*
*	Do-file:		an_risk.do
*
*	Project:		SGTF Omicron
*
*	Programmed by:	Daniel Grint
*
*	Data used:		output/cr_main.dta
*
*	Data created:	
*
*	Other output:	an_risk.log containing:
*					1-Subgroup analyses of OR and absolute risk
*
********************************************************************************
*
*	Purpose:		This do-file runs logit models, calculating relative (glm) and
*					absolute odds (margins)
*  
********************************************************************************

* Open a log file
cap log close
log using ./logs/an_risk, replace t

clear


use ./output/main.dta

* DROP MISSING UTLA
noi di "DROPPING MISSING UTLA DATA"
drop if utla_group==""

* DROP IF NO DATA ON SGTF
noi di "DROPPING NO SGTF DATA" 
drop if has_sgtf==0

noi di "SUBSETTING ON COX 14-DAY CENSORED POPULATION"
keep if ae_14_pop==1

noi di "SUBSETTING ON START/STOP DATES"
keep if inrange(start_week,49,52)

tab sgtf cox_ae14, row



*******************
/* Unadjusted OR */
*******************

glm cox_ae14 i.sgtf, family(bin) link(logit) eform

* Absolute odds
margins sgtf



***************************************
/* Fully adjusted OR - age as spline */
***************************************

glm cox_ae14 i.sgtf ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			, family(bin) link(logit) eform
			
			
* Vaccine/SGTF interaction

glm cox_ae14 i.sgtf##ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week age1 age2 age3 i.home_bin ///
			 if eth2 != 6 ///
			, family(bin) link(logit) eform
		
est store vax_int

margins comorb_cat#vax if sgtf==0, asobserved

est restore vax_int

margins comorb_cat#vax if sgtf==1, asobserved

est restore vax_int

margins comorb_cat#vax, dydx(sgtf) asobserved


*************************************
/* Fully adjusted OR - age grouped */
*************************************

glm cox_ae14 i.sgtf##ib2.vax i.male ib1.imd ib1.eth2 ib1.smoke_nomiss2 ib1.obese4cat ib1.hh_total_cat ///
			 ib1.rural_urban5 ib0.comorb_cat ib49.start_week ib1.agegroup6 i.home_bin ///
			 if eth2 != 6 ///
			, family(bin) link(logit) eform

est store fully


* Adjusted absolute risks
margins comorb_cat#vax#agegroup6 if sgtf==0, post asobserved

* Save risk estimates
matrix est = e(b)
matrix inv_est = est'
svmat inv_est

* Save SE estimates
matrix var = e(V)
matrix diag_var = vecdiag(var)
matrix inv_var = diag_var'
svmat inv_var
gen sq_var = sqrt(inv_var1)

noi disp "CHECK MARGINS ARE CORRECTLY CALCULATED TO MATCH ABOVE"
list inv_est1 sq_var in 1/72

* Re-Calculate CI
gen risk0 = inv_est1*100
gen lb0 = (inv_est1 - invnormal(0.975)*sq_var)*100
gen ub0 = (inv_est1 + invnormal(0.975)*sq_var)*100

order lb ub, after(risk0)


est restore fully

* Adjusted absolute risks
margins comorb_cat#vax#agegroup6 if sgtf==1, post asobserved

* Save risk estimates
matrix est1 = e(b)
matrix inv_estx = est1'
svmat inv_estx

* Save SE estimates
matrix var1 = e(V)
matrix diag_var1 = vecdiag(var1)
matrix inv_varx = diag_var1'
svmat inv_varx
gen sq_varx = sqrt(inv_varx1)

noi disp "CHECK MARGINS ARE CORRECTLY CALCULATED TO MATCH ABOVE"
list inv_estx1 sq_varx in 1/72

* Re-Calculate CI
gen risk1 = inv_estx1*100
gen lb1 = (inv_estx1 - invnormal(0.975)*sq_varx)*100
gen ub1 = (inv_estx1 + invnormal(0.975)*sq_varx)*100

order lb1 ub1, after(risk1)


* Adjusted risk difference
est restore fully

margins comorb_cat#vax#agegroup6, dydx(sgtf) post asobserved

* Save risk estimates
matrix diff = e(b)
matrix diff = diff[1, 73..144]
matrix inv_diff = diff'
svmat inv_diff

* Save SE estimates
matrix dvar = e(V)
matrix diag_dvar = vecdiag(dvar)
matrix inv_dvar = diag_dvar'
matrix inv_dvar = inv_dvar[73..144,1]
svmat inv_dvar
gen sq_dvar = sqrt(inv_dvar1)

noi disp "CHECK MARGINS ARE CORRECTLY CALCULATED TO MATCH ABOVE"
list inv_diff1 sq_dvar in 1/72

* Re-Calculate CI
gen diff1 = inv_diff1*100
gen dlb = (inv_diff1 - invnormal(0.975)*sq_dvar)*100
gen dub = (inv_diff1 + invnormal(0.975)*sq_dvar)*100

order dlb dub, after(diff1)


gen risk_labels = "Unvax: 0-<40" in 1
replace risk_labels = "40-<55" in 2
replace risk_labels = "55-<65" in 3
replace risk_labels = "65-<75" in 4
replace risk_labels = "75-<85" in 5
replace risk_labels = "85+" in 6

replace risk_labels = "First dose: 0-<40" in 7
replace risk_labels = "40-<55" in 8
replace risk_labels = "55-<65" in 9
replace risk_labels = "65-<75" in 10
replace risk_labels = "75-<85" in 11
replace risk_labels = "85+" in 12

replace risk_labels = "Second dose: 0-<40" in 13
replace risk_labels = "40-<55" in 14
replace risk_labels = "55-<65" in 15
replace risk_labels = "65-<75" in 16
replace risk_labels = "75-<85" in 17
replace risk_labels = "85+" in 18

replace risk_labels = "Booster: 0-<40" in 19
replace risk_labels = "40-<55" in 20
replace risk_labels = "55-<65" in 21
replace risk_labels = "65-<75" in 22
replace risk_labels = "75-<85" in 23
replace risk_labels = "85+" in 24

replace risk_labels = "Unvax: 0-<40" in 25
replace risk_labels = "40-<55" in 26
replace risk_labels = "55-<65" in 27
replace risk_labels = "65-<75" in 28
replace risk_labels = "75-<85" in 29
replace risk_labels = "85+" in 30

replace risk_labels = "First dose: 0-<40" in 31
replace risk_labels = "40-<55" in 32
replace risk_labels = "55-<65" in 33
replace risk_labels = "65-<75" in 34
replace risk_labels = "75-<85" in 35
replace risk_labels = "85+" in 36

replace risk_labels = "Second dose: 0-<40" in 37
replace risk_labels = "40-<55" in 38
replace risk_labels = "55-<65" in 39
replace risk_labels = "65-<75" in 40
replace risk_labels = "75-<85" in 41
replace risk_labels = "85+" in 42

replace risk_labels = "Booster: 0-<40" in 43
replace risk_labels = "40-<55" in 44
replace risk_labels = "55-<65" in 45
replace risk_labels = "65-<75" in 46
replace risk_labels = "75-<85" in 47
replace risk_labels = "85+" in 48

replace risk_labels = "Unvax: 0-<40" in 49
replace risk_labels = "40-<55" in 50
replace risk_labels = "55-<65" in 51
replace risk_labels = "65-<75" in 52
replace risk_labels = "75-<85" in 53
replace risk_labels = "85+" in 54

replace risk_labels = "First dose: 0-<40" in 55
replace risk_labels = "40-<55" in 56
replace risk_labels = "55-<65" in 57
replace risk_labels = "65-<75" in 58
replace risk_labels = "75-<85" in 59
replace risk_labels = "85+" in 60

replace risk_labels = "Second dose: 0-<40" in 61
replace risk_labels = "40-<55" in 62
replace risk_labels = "55-<65" in 63
replace risk_labels = "65-<75" in 64
replace risk_labels = "75-<85" in 65
replace risk_labels = "85+" in 66

replace risk_labels = "Booster: 0-<40" in 67
replace risk_labels = "40-<55" in 68
replace risk_labels = "55-<65" in 69
replace risk_labels = "65-<75" in 70
replace risk_labels = "75-<85" in 71
replace risk_labels = "85+" in 72




***********************************
/* Output table of absolute risk */
***********************************

cap file close tablecontent

file open tablecontent using ./output/table3.txt, write text replace

file write tablecontent ("Table 3: Absolute risk of AE") _n _n

file write tablecontent ("Comorbidities/Vax/Age group")		_tab ///
						("S-Pos (95% CI)")					_tab ///
						("S-Fail (95% CI)")					_tab ///
						("Diff (95% CI)")					_n

forvalues i=1/72 {
	
	preserve
		keep if _n == `i'
		if inlist(`i',7,13,19,31,37,43,55,61,67) {
			file write tablecontent _n
		}
		if `i'==1 {
			file write tablecontent _n ("No Comorbidities") _n
		}
		if `i'==25 {
			file write tablecontent _n ("1 Comorbidity") _n
		}
		if `i'==49 {
			file write tablecontent _n ("2+ Comorbidities") _n
		}
		file write tablecontent %9s (risk_labels) _tab %4.2f (risk0) (" (") %4.2f (lb0) ("-") %4.2f (ub0) (")") _tab %4.2f (risk1) (" (") %4.2f (lb1) ("-") %4.2f (ub1) (")") _tab %4.2f (diff1) (" (") %4.2f (dlb) ("-") %4.2f (dub) (")") _n
	restore

}

file close tablecontent



log close



insheet using ./output/table3.txt, clear
export excel using ./output/table3.xlsx, replace
