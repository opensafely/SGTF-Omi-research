********************************************************************************
*
*	Do-file:		cr_imputed.do
*
*	Project:		SGTF Omicron
*
*	Programmed by:	Daniel Grint
*
*	Data used:		output/main.dta
*
*	Data created:	output/main_imputed.dta
*
*	Other output:	cr_imputed.log
*
*
********************************************************************************
*
*	Purpose:		This do-file imputes missing ethnicity data
*  
********************************************************************************

* Open a log file
cap log close
log using ./logs/cr_imputed, replace t

clear

use ./output/main.dta

tab eth2
recode eth2 6=. 5=0
tab eth2, m


egen inc = rowmiss(age1 age2 age3 male obese4cat smoke_nomiss imd comorb_cat region ///
					vax rural_urban hh_total_cat home_bin sgtf start_week cox_ae)
					
keep if inc==0


mi set wide
mi register imputed eth2

mi impute logit eth2				///
			age1 age2 age3 			///
			i.male 					///
			i.obese4cat				///
			i.smoke_nomiss			///
			i.imd 					///
			i.comorb_cat			///
			i.region				///
			i.vax					///
			i.rural_urban			///
			i.hh_total_cat			///
			i.home_bin				///
			i.sgtf					///
			i.start_week			///
			cox_ae, add(10) rseed(06012022) noisily iter(20)
			
			


label data "SGTF OMICRON IMPUTED DATASET: $S_DATE"

save ./output/main_imputed.dta, replace

log close
