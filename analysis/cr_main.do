********************************************************************************
*
*	Do-file:		cr_main.do
*
*	Project:		SGTF-Omi
*
*	Programmed by:	Daniel Grint
*					Adapted from SGTF-CFR
*
*	Data used:		output/input.csv
*					lookups/MSOA_lookup.dta
*
*	Data created:	output/main.dta  (main analysis dataset)
*
*	Other output:	None
*
********************************************************************************
*
*	Purpose:		This do-file creates the variables required for the 
*					main analysis and saves into Stata datasets.
*  
********************************************************************************

* Open a log file
cap log close
log using ./logs/cr_main, replace t

clear


import delimited ./output/input.csv
merge m:1 msoa using ./lookups/MSOA_lookup
drop if _merge==2
drop _merge


rename hiv hiv_code

****************************
*  Create required cohort  *
****************************

gen study_start = date(sgss_pos_inrange, "YMD")
summ study_start
noi disp "MINIMUM START DATE: " %td r(min)
noi disp "MAXIMUM START DATE: " %td r(max)

gen study_end = date("01jan2022", "DMY")
format %td study_start study_end

* DROP IF NO POSITIVE PCR TEST IN SGSS DURING STUDY PERIOD
noi di "NO POSITIVE TEST IN STUDY PERIOD"
drop if sgss_pos_inrange == ""

* DROP IF DIED ON/BEFORE STUDY START DATE
noi di "DIED ON/BEFORE STUDY START DATE:" 
drop if date(died_date_ons, "YMD") <= study_start

* Age: Exclude those with implausible ages
assert age<.
noi di "DROPPING AGE>105:" 
drop if age>105

* Sex: Exclude categories other than M and F
assert inlist(sex, "M", "F", "I", "U")
noi di "DROPPING GENDER NOT M/F:" 
drop if inlist(sex, "I", "U")


**************************
*  CREATE SGTF VARIABLE  *
**************************

desc sgtf

replace sgtf=99 if sgtf==.

gen has_sgtf=0
replace has_sgtf=1 if inrange(sgtf,0,1)

label define sgtfLab 0 "S-Pos" 1 "S-Fail" 9 "Unclassified" 99 "Blank"
label values sgtf sgtfLab


******************************
*  Convert strings to dates  *
******************************

* Covariates
foreach var of varlist 	bp_sys_date 					///
						bp_dias_date 					///
						hba1c_percentage_date			///
						hba1c_mmol_per_mol_date			///
						hypertension					///
						bmi_date_measured				///
						chronic_respiratory_disease 	///
						chronic_cardiac_disease 		///
						diabetes 						///
						lung_cancer 					///
						haem_cancer						///
						other_cancer 					///
						chronic_liver_disease 			///
						stroke							///
						dementia		 				///
						other_neuro 					///
						organ_transplant 				///	
						dysplenia						///
						sickle_cell 					///
						aplastic_anaemia 				///
						hiv_date						///
						permanent_immunodeficiency 		///
						temporary_immunodeficiency		///
						ra_sle_psoriasis  dialysis 	{
	confirm string variable `var'
	replace `var' = `var' + "-15"
	rename `var' `var'_dstr
	replace `var'_dstr = " " if `var'_dstr == "-15"
	gen `var'_date = date(`var'_dstr, "YMD") 
	order `var'_date, after(`var'_dstr)
	drop `var'_dstr
	format `var'_date %td
}

rename bmi_date_measured_date      bmi_date_measured
rename bp_dias_date_measured_date  bp_dias_date
rename bp_sys_date_measured_date   bp_sys_date
rename hba1c_percentage_date_date  hba1c_percentage_date
rename hba1c_mmol_per_mol_date_date  hba1c_mmol_per_mol_date
rename hiv_date_date hiv_date

* Dates of: covid tests, ONS death, vaccinations, etc.
foreach var of varlist 	dereg_date died_date_ons ae_covid_date ae_any_date vaxdate1 ///
						vaxdate2 vaxdate3 last_covid_tpp_probable last_pos_test_sgss sgss_pos_inrange {
		confirm string variable `var'
		rename `var' _tmp
		gen `var' = date(_tmp, "YMD")
		drop _tmp
		format %d `var'
}


***************************************************
*  Define vaccination and prior infection status  *
***************************************************

gen vax = 0
replace vax = 1 if (vaxdate1+14) < study_start	// First dose >14-days prior to infection
replace vax = 2 if (vaxdate2+14) < study_start	// Second dose >14-days prior to infection
replace vax = 3 if (vaxdate3+14) < study_start	// Booster >14-days prior to infection

label define vaxLab 0 "Unvax" 1 "First dose" 2 "Second dose" 3 "Booster"
label values vax vaxLab

tab vax, m
tab vax sgtf if has_sgtf == 1, m col

gen prev_inf = 0
replace prev_inf = 1 if last_covid_tpp_probable < study_start
replace prev_inf = 1 if last_pos_test_sgss < study_start

label define prev_infLab 0 "None" 1 "Prior inf."
label values prev_inf prev_infLab

tab prev_inf, m
tab prev_inf sgtf if has_sgtf == 1, m col
 
*******************************
*  Recode implausible values  *
*******************************

* BMI 

* Only keep if within certain time period? using bmi_date_measured ?
* NB: Some BMI dates in future or after cohort entry

* Set implausible BMIs to missing:
replace bmi = . if !inrange(bmi, 15, 50)

* Sex
assert inlist(sex, "M", "F")
gen male = (sex=="M")
drop sex

label define maleLab 0 "F" 1 "M"
label values male maleLab


* Smoking
label define smoke_nomissLab 1 "Never" 2 "Former" 3 "Current" 

gen     smoke = 1  if smoking_status=="N"
replace smoke = 2  if smoking_status=="E"
replace smoke = 3  if smoking_status=="S"
replace smoke = . if smoking_status=="M"
label values smoke smoke_nomissLab
drop smoking_status

* Ethnicity (5 category)
tab ethnicity, m

replace ethnicity = . if ethnicity==.
label define ethnicityLab 	1 "White"  					///
							2 "Mixed" 					///
							3 "Asian or Asian British"	///
							4 "Black"  					///
							5 "Other"					
						
label values ethnicity ethnicityLab

* Re-order ethnicity
gen eth5=1 if ethnicity==1
replace eth5=2 if ethnicity==3
replace eth5=3 if ethnicity==4
replace eth5=4 if ethnicity==2
replace eth5=5 if ethnicity==5
replace eth5=6 if ethnicity==.

label define eth5Lab	1 "White"  					///
						2 "South Asian"				///
						3 "Black"  					///
						4 "Mixed" 					///
						5 "Other"					///
						6 "Missing"
 
label values eth5 eth5Lab

recode eth5 2/4=5, gen(eth2)
order eth2, after(eth5)

label values eth2 eth5Lab

tab eth2, m


* Ethnicity (16 category)
replace ethnicity_16 = . if ethnicity==.
label define ethnicity_16 									///
						1 "British or Mixed British" 		///
						2 "Irish" 							///
						3 "Other White" 					///
						4 "White + Black Caribbean" 		///
						5 "White + Black African"			///
						6 "White + Asian" 					///
 						7 "Other mixed" 					///
						8 "Indian or British Indian" 		///
						9 "Pakistani or British Pakistani" 	///
						10 "Bangladeshi or British Bangladeshi" ///
						11 "Other Asian" 					///
						12 "Caribbean" 						///
						13 "African" 						///
						14 "Other Black" 					///
						15 "Chinese" 						///
						16 "Other" 							
						
label values ethnicity_16 ethnicity_16

* Ethnicity (16 category grouped further)
* Generate a version of the full breakdown with mixed in one group
gen ethnicity_16_combinemixed = ethnicity_16
recode ethnicity_16_combinemixed 4/7 = 4
label define ethnicity_16_combinemixed 	///
						1 "British or Mixed British" ///
						2 "Irish" ///
						3 "Other White" ///
						4 "All mixed" ///
						8 "Indian or British Indian" ///
						9 "Pakistani or British Pakistani" ///
						10 "Bangladeshi or British Bangladeshi" ///
						11 "Other Asian" ///
						12 "Caribbean" ///
						13 "African" ///
						14 "Other Black" ///
						15 "Chinese" ///
						16 "Other" 
						
label values ethnicity_16_combinemixed ethnicity_16_combinemixed

* STP 
tab stp

rename stp stp_old
bysort stp_old: gen stp = 1 if _n==1
replace stp = sum(stp)
drop stp_old

* MSOA/UTLA

egen n_msoa = tag(msoa)
count if n_msoa

bysort msoa: gen count1 = _N
summ count1, d

egen n_utla = tag(utla)
count if n_utla

bysort utla: gen count2 = _N
summ count2, d

* Regroup UTLAs with small case numbers

gen utla_group = utla_name
tab utla_group

replace utla_group = "Redbridge, Barking and Dagenham" if utla_name == "Barking and Dagenham"
replace utla_group = "Redbridge, Barking and Dagenham" if utla_name == "Redbridge"

replace utla_group = "Bucks/Ox/West. Berks/Swindon" if utla_name == "Buckinghamshire"
replace utla_group = "Bucks/Ox/West. Berks/Swindon" if utla_name == "Oxfordshire"
replace utla_group = "Bucks/Ox/West. Berks/Swindon" if utla_name == "Swindon"
replace utla_group = "Bucks/Ox/West. Berks/Swindon" if utla_name == "West Berkshire"

replace utla_group = "Camden and Westminster" if utla_name == "Camden"
replace utla_group = "Camden and Westminster" if utla_name == "Westminster"

replace utla_group = "" if utla_name == "Isles of Scilly"

replace utla_group = "Richmond and Hounslow" if utla_name == "Richmond upon Thames"
replace utla_group = "Richmond and Hounslow" if utla_name == "Hounslow"

replace utla_group = "Rutland and Lincoln" if utla_name == "Rutland"
replace utla_group = "Rutland and Lincoln" if utla_name == "Lincolnshire"

replace utla_group = "Bolton and Tameside" if utla_name == "Bolton"
replace utla_group = "Bolton and Tameside" if utla_name == "Tameside"

tab utla_group, m



* NHS England regions
tab region, m

gen region2=.
replace region2=0 if region=="East"
replace region2=1 if region=="East Midlands"
replace region2=2 if region=="London"
replace region2=3 if region=="North East"
replace region2=4 if region=="North West"
replace region2=5 if region=="South East"
replace region2=6 if region=="South West"
replace region2=7 if region=="West Midlands"
replace region2=8 if region=="Yorkshire and The Humber"

drop region
rename region2 region

label define regionLab	0 "East" ///
						1 "East Midlands"  ///
						2 "London" ///
						3 "North East" ///
						4 "North West" ///
						5 "South East" ///
						6 "South West" ///
						7 "West Midlands" ///
						8 "Yorkshire and the Humber"
						
label values region regionLab


**************************
*  Categorise variables  *
**************************

/*  Age variables  */

noi di "DROPPING IF NO AGE" 
drop if age>=.

* Create categorised age
recode age 	0/17.9999=0 ///
			18/29.9999 = 1 /// 
		    30/39.9999 = 2 /// 
			40/49.9999 = 3 ///
			50/59.9999 = 4 ///
			60/69.9999 = 5 ///
			70/79.9999 = 6 ///
			80/max = 7, gen(agegroup) 

label define agegroupLab 	0 "0-<18" ///
							1 "18-<30" ///
							2 "30-<40" ///
							3 "40-<50" ///
							4 "50-<60" ///
							5 "60-<70" ///
							6 "70-<80" ///
							7 "80+"
						
label values agegroup agegroupLab

* For subgroup analysis
recode age 	0/64.9999=1 ///
			65/74.9999=2 ///
			75/84.9999=3 ///
			85/max=4, gen(agegroupA) 

label define agegroupALab 	1 "0-<65" ///
							2 "65-<75" ///
							3 "75-<85" ///
							4 "85+"
							
label values agegroupA agegroupALab


recode agegroupA 4=3, gen(agegroupB)

label define agegroupBLab 	1 "0-<65" ///
							2 "65-<75" ///
							3 "75+"
							
label values agegroupB agegroupBLab


* More age categories
recode age 	0/39.9999=0 ///
			40/54.9999 = 1 /// 
		    55/64.9999 = 2 /// 
			65/74.9999 = 3 ///
			75/84.9999 = 4 ///
			85/max = 5, gen(agegroup6) 

label define agegroup6Lab 	0 "0-39" ///
							1 "40-54" ///
							2 "55-64" ///
							3 "65-74" ///
							4 "75-84" ///
							5 "85+"
							
label values agegroup6 agegroup6Lab

recode age 	0/39.9999=0 ///
			40/69.9999 = 1 /// 
		    70/max = 2, gen(agegroup3)
			
label define agegroup3Lab 	0 "0-39" ///
							1 "40-69" ///
							2 "70+"
							
label values agegroup3 agegroup3Lab


* Create binary age
recode age min/69.999=0 70/max=1, gen(age70)

* Check there are no missing ages
assert age<.
assert agegroup<.
assert age70<.

* Create restricted cubic splines for age centred on 65
gen age65 = age - 65
mkspline age = age65, cubic nknots(4)


/*  Body Mass Index  */

* BMI (NB: watch for missingness)
gen 	bmicat = .
recode  bmicat . = 1 if bmi<18.5
recode  bmicat . = 2 if bmi<25
recode  bmicat . = 3 if bmi<30
recode  bmicat . = 4 if bmi<35
recode  bmicat . = 5 if bmi<40
recode  bmicat . = 6 if bmi<.
replace bmicat = . if bmi>=.

label define bmicat 1 "Underweight (<18.5)" 	///
					2 "Normal (18.5-24.9)"		///
					3 "Overweight (25-29.9)"	///
					4 "Obese I (30-34.9)"		///
					5 "Obese II (35-39.9)"		///
					6 "Obese III (40+)"			
					
label values bmicat bmicat

* Create more granular categorisation
recode bmicat 1/3 . = 1 4=2 5=3 6=4, gen(obese4cat)

label define obese4catLab 	1 "No record of obesity" 	///
							2 "Obese I (30-34.9)"		///
							3 "Obese II (35-39.9)"		///
							4 "Obese III (40+)"		
label values obese4cat obese4catLab
order obese4cat, after(bmicat)

tab obese4cat, m


/*  Smoking  */

* Create non-missing 3-category variable for current smoking
recode smoke .=1, gen(smoke_nomiss)
order smoke_nomiss, after(smoke)
label values smoke_nomiss smoke_nomissLab

recode smoke_nomiss 3=2, gen(smoke_nomiss2)
order smoke_nomiss2, after(smoke_nomiss)

label values smoke_nomiss2 smoke_nomissLab

tab smoke_nomiss2, m

/*  Asthma  */

* Asthma  (coded: 0 No, 1 Yes no OCS, 2 Yes with OCS)
rename asthma asthmacat
recode asthmacat 0=1 1=2 2=3 .=1
label define asthmacat 1 "No" 2 "Yes, no OCS" 3 "Yes with OCS"
label values asthmacat asthmacat

gen asthma = (asthmacat==2|asthmacat==3)


/*  Blood pressure   */

* Categorise
gen     bpcat = 1 if bp_sys < 120 &  bp_dias < 80
replace bpcat = 2 if inrange(bp_sys, 120, 130) & bp_dias<80
replace bpcat = 3 if inrange(bp_sys, 130, 140) | inrange(bp_dias, 80, 90)
replace bpcat = 4 if (bp_sys>=140 & bp_sys<.) | (bp_dias>=90 & bp_dias<.) 
replace bpcat = . if bp_sys>=. | bp_dias>=. | bp_sys==0 | bp_dias==0

label define bpcat 1 "Normal" 2 "Elevated" 3 "High, stage I"	///
					4 "High, stage II" 
label values bpcat bpcat

recode bpcat .=1, gen(bpcat_nomiss)
label values bpcat_nomiss bpcat

* Create non-missing indicator of known high blood pressure
gen bphigh = (bpcat==4)
order bpcat bphigh, after(bp_dias_date)


/*  IMD  */

* Group into 5 groups
rename imd imd_o
egen imd = cut(imd_o), group(5) icodes
replace imd = imd + 1
replace imd = . if imd_o==-1
drop imd_o
tab imd

* Reverse the order (so high is more deprived)
recode imd 5=1 4=2 3=3 2=4 1=5 .=.
tab imd

label define imdLab 1 "1 least deprived" 2 "2" 3 "3" 4 "4" 5 "5 most deprived" 
label values imd imdLab 

noi di "DROPPING IF NO IMD" 
drop if imd>=.


/*  HOUSEHOLD SIZE  */

gen hh_total_cat=.
replace hh_total_cat=1 if household_size >=1 & household_size<=2
replace hh_total_cat=2 if household_size >=3 & household_size<=5
replace hh_total_cat=3 if household_size >=6 & household_size<=10
replace hh_total_cat=4 if household_size >=11 & household_size !=.

label define hh_total_catLab	1 "1-2" ///
								2 "3-5" ///
								3 "6-10" ///
								4 "11+"
											
label values hh_total_cat hh_total_catLab

tab hh_total_cat, m


/*  RURAL OR URBAN  */

* Create a 5 category rural urban variable based upon meeting with Roz 21st October
gen rural_urban5=.
replace rural_urban5=1 if rural_urban==1
replace rural_urban5=2 if rural_urban==2
replace rural_urban5=3 if rural_urban==3|rural_urban==4
replace rural_urban5=4 if rural_urban==5|rural_urban==6
replace rural_urban5=5 if rural_urban==7|rural_urban==8

label define rural_urban5Lab	1 "Urban major conurbation" ///
								2 "Urban minor conurbation" ///
								3 "Urban city and town" ///
								4 "Rural town and fringe" ///
								5 "Rural village and dispersed"
								
label values rural_urban5 rural_urban5Lab
tab rural_urban5, m


/*  CARE HOME TYPE  */

tab care_home_type, m

gen home_bin=0
replace home_bin=1 if care_home_type=="PC" | care_home_type=="PN" | care_home_type=="PS"

tab care_home_type home_bin, m

label define home_binLab	0 "Private home" ///
							1 "Care home"

label values home_bin home_binLab


/*  Centred age, sex, IMD, ethnicity (for adjusted KM plots)  */ 

* Centre age (linear)
summ age
gen c_age = age-r(mean)



**************************************************
*  Create binary comorbidity indices from dates  *
**************************************************

* Comorbidities ever before study_start
foreach var of varlist	chronic_respiratory_disease_date 	///
						chronic_cardiac_disease_date 		///
						diabetes_date 						///
						chronic_liver_disease_date 			///
						stroke_date							///
						dementia_date						///
						other_neuro_date					///
						organ_transplant_date 				///
						aplastic_anaemia_date				///
						hypertension_date					///
						dysplenia_date 						///
						sickle_cell_date 					///
						hiv_date							///
						permanent_immunodeficiency_date		///
						temporary_immunodeficiency_date		///
						ra_sle_psoriasis_date dialysis_date {
	local newvar =  substr("`var'", 1, length("`var'") - 5)
	gen `newvar' = (`var'< study_start)
	order `newvar', after(`var')
}


**************************
*  Epidemiological week  *
**************************

gen start_week = 52 if study_start <= date("01jan2022", "DMY")
replace start_week = 51 if study_start <= date("25dec2021", "DMY")
replace start_week = 50 if study_start <= date("18dec2021", "DMY")
replace start_week = 49 if study_start <= date("11dec2021", "DMY")
replace start_week = 48 if study_start <= date("04dec2021", "DMY")
replace start_week = 47 if study_start <= date("27nov2021", "DMY")
replace start_week = 46 if study_start <= date("20nov2021", "DMY")
replace start_week = 45 if study_start <= date("13nov2021", "DMY")
replace start_week = 44 if study_start <= date("06nov2021", "DMY")
replace start_week = 43 if study_start <= date("30Oct2021", "DMY")
replace start_week = 42 if study_start <= date("23Oct2021", "DMY")
replace start_week = 41 if study_start <= date("16Oct2021", "DMY")
replace start_week = 40 if study_start <= date("09Oct2021", "DMY")

tab start_week, m

label define start_weekLab	40 "03Oct-09Oct" ///
							41 "10Oct-16Oct" ///
							42 "17Oct-23Oct" ///
							43 "24Oct-30Oct" ///
							44 "31Oct-06Nov" ///
							45 "07Nov-13Nov" ///
							46 "14Nov-20Nov" ///
							47 "21Nov-27Nov" ///
							48 "28Nov-04Dec" ///
							49 "05Dec-11Dec" ///
							50 "12Dec-18Dec" ///
							51 "19Dec-25Dec" ///
							52 "26Dec-01Jan"

							
label values start_week start_weekLab

tab start_week, m

* Recode small epi weeks 1 and 2 for epi week interaction
/*
recode start_week 1=2, gen(start_weekA)
label define start_weekLabA	2 "16Nov-29Nov"	///
							3 "30Nov-06Dec" ///
							4 "07Dec-13Dec" ///
							5 "14Dec-20Dec" ///
							6 "21Dec-27Dec" ///
							7 "28Dec-03Jan" ///
							8 "04Jan-11Jan"
							
label values start_weekA start_weekLabA

tab start_week start_weekA
*/


***************************
*  Grouped comorbidities  *
***************************

/*  Neurological  */

* Stroke and dementia
egen stroke_dementia = rowmax(stroke dementia)
order stroke_dementia, after(dementia_date)


/*  Spleen  */

* Spleen problems (dysplenia/splenectomy/etc and sickle cell disease)   
egen spleen = rowmax(dysplenia sickle_cell) 
order spleen, after(sickle_cell)


/*  Cancer  */

label define cancer 1 "Never" 2 "Last year" 3 "2-5 years ago" 4 "5+ years"

local fiveybefore = study_start-5*365.25
local oneybefore = study_start-365.25

* Haematological malignancies
gen     cancer_haem_cat = 4 if inrange(haem_cancer_date, d(1/1/1900), `fiveybefore')
replace cancer_haem_cat = 3 if inrange(haem_cancer_date, `fiveybefore', `oneybefore')
replace cancer_haem_cat = 2 if inrange(haem_cancer_date, `oneybefore', study_start)
recode  cancer_haem_cat . = 1
label values cancer_haem_cat cancer

* All other cancers
gen     cancer_exhaem_cat = 4 if inrange(lung_cancer_date,  d(1/1/1900), `fiveybefore') | ///
								 inrange(other_cancer_date, d(1/1/1900), `fiveybefore') 
replace cancer_exhaem_cat = 3 if inrange(lung_cancer_date,  `fiveybefore', `oneybefore') | ///
								 inrange(other_cancer_date, `fiveybefore', `oneybefore') 
replace cancer_exhaem_cat = 2 if inrange(lung_cancer_date,  `oneybefore', study_start) | ///
								 inrange(other_cancer_date, `oneybefore', study_start)
recode  cancer_exhaem_cat . = 1
label values cancer_exhaem_cat cancer

* Put variables together
order cancer_exhaem_cat cancer_haem_cat, after(other_cancer_date)


/*  Immunosuppression  */

* Immunosuppressed:
* HIV, permanent immunodeficiency ever, OR 
* temporary immunodeficiency or aplastic anaemia last year
gen temp1  = max(hiv, permanent_immunodeficiency)
gen temp2  = inrange(temporary_immunodeficiency_date, `oneybefore', study_start)
gen temp3  = inrange(aplastic_anaemia_date, `oneybefore', study_start)

egen other_immunosuppression = rowmax(temp1 temp2 temp3)
drop temp1 temp2 temp3
order other_immunosuppression, after(temporary_immunodeficiency)


/*  Hypertension  */

gen htdiag_or_highbp = bphigh
recode htdiag_or_highbp 0 = 1 if hypertension==1 


************
*   eGFR   *
************

* Set implausible creatinine values to missing (Note: zero changed to missing)
replace creatinine = . if !inrange(creatinine, 20, 3000) 
	
* Divide by 88.4 (to convert umol/l to mg/dl)
gen SCr_adj = creatinine/88.4

gen min=.
replace min = SCr_adj/0.7 if male==0
replace min = SCr_adj/0.9 if male==1
replace min = min^-0.329  if male==0
replace min = min^-0.411  if male==1
replace min = 1 if min<1

gen max=.
replace max=SCr_adj/0.7 if male==0
replace max=SCr_adj/0.9 if male==1
replace max=max^-1.209
replace max=1 if max>1

* egfr calculated using CKD-EPI formula with no eth
gen egfr=min*max*141
replace egfr=egfr*(0.993^age)
replace egfr=egfr*1.018 if male==0

* Categorise into ckd stages
* CKD stage calc without eth
egen egfr_cat = cut(egfr), at(0, 15, 30, 45, 60, 5000)
recode egfr_cat 0=5 15=4 30=3 45=2 60=0, generate(ckd)
* 0 = "No CKD" 	2 "stage 3a" 3 "stage 3b" 4 "stage 4" 5 "stage 5"
label define ckd 0 "No CKD" 1 "CKD"
label values ckd ckd

* Convert into CKD group
*recode ckd 2/5=1, gen(chronic_kidney_disease)
*replace chronic_kidney_disease = 0 if creatinine==. 

recode ckd 0=1 2/3=2 4/5=3, gen(reduced_kidney_function_cat)
replace reduced_kidney_function_cat = 1 if creatinine==. 
label define reduced_kidney_function_catlab ///
	1 "None" 2 "Stage 3a/3b egfr 30-60	" 3 "Stage 4/5 egfr<30"
label values reduced_kidney_function_cat reduced_kidney_function_catlab 

*More detailed version incorporating stage 5 or dialysis as a separate category	
recode ckd 0=1 2/3=2 4=3 5=4, gen(reduced_kidney_function_cat2)
replace reduced_kidney_function_cat2 = 1 if creatinine==. 
replace reduced_kidney_function_cat2 = 4 if dialysis==1 

label define reduced_kidney_function_cat2lab ///
	1 "None" 2 "Stage 3a/3b egfr 30-60	" 3 "Stage 4 egfr 15-<30" 4 "Stage 5 egfr <15 or dialysis"
label values reduced_kidney_function_cat2 reduced_kidney_function_cat2lab 
 
	
***********
*  Hba1c  *
***********

/*  Diabetes severity  */

* Set zero or negative to missing
replace hba1c_percentage   = . if hba1c_percentage<=0
replace hba1c_mmol_per_mol = . if hba1c_mmol_per_mol<=0

local fifteenmbefore = study_start-15*(365.25/12)

* Only consider measurements in last 15 months
replace hba1c_percentage   = . if hba1c_percentage_date   < `fifteenmbefore'
replace hba1c_mmol_per_mol = . if hba1c_mmol_per_mol_date < `fifteenmbefore'


/* Express  HbA1c as percentage  */ 

* Express all values as perecentage 
noi summ hba1c_percentage hba1c_mmol_per_mol 
gen 	hba1c_pct = hba1c_percentage 
replace hba1c_pct = (hba1c_mmol_per_mol/10.929)+2.15 if hba1c_mmol_per_mol<. 

* Valid % range between 0-20  
replace hba1c_pct = . if !inrange(hba1c_pct, 0, 20) 
replace hba1c_pct = round(hba1c_pct, 0.1)


/* Categorise hba1c and diabetes  */

* Group hba1c
gen 	hba1ccat = 0 if hba1c_pct <  6.5
replace hba1ccat = 1 if hba1c_pct >= 6.5  & hba1c_pct < 7.5
replace hba1ccat = 2 if hba1c_pct >= 7.5  & hba1c_pct < 8
replace hba1ccat = 3 if hba1c_pct >= 8    & hba1c_pct < 9
replace hba1ccat = 4 if hba1c_pct >= 9    & hba1c_pct !=.
label define hba1ccat 0 "<6.5%" 1">=6.5-7.4" 2">=7.5-7.9" 3">=8-8.9" 4">=9"
label values hba1ccat hba1ccat
tab hba1ccat

* Create diabetes, split by control/not
gen     diabcat = 1 if diabetes==0
replace diabcat = 2 if diabetes==1 & inlist(hba1ccat, 0, 1)
replace diabcat = 3 if diabetes==1 & inlist(hba1ccat, 2, 3, 4)
replace diabcat = 4 if diabetes==1 & !inlist(hba1ccat, 0, 1, 2, 3, 4)

label define diabcat 	1 "No diabetes" 			///
						2 "Controlled diabetes"		///
						3 "Uncontrolled diabetes" 	///
						4 "Diabetes, no hba1c measure"
label values diabcat diabcat

* Delete unneeded variables
drop hba1c_pct hba1c_percentage hba1c_mmol_per_mol


******************************
*  Aggregated comorbidities  *
******************************

/*
	Create a comorbidities variable based upon Fizz's JCVI work that has 0, 1, 2 or more of 
	the following comorbdities: 
	(1) respiratory disease, (2) severe asthma, (3) chronic cardiac disease, (4) diabetes, 
	(5) non-haematological cancer (diagnosed in last year), (6) haematological cancer (diagnosed within 5 years), 
	(7) liver disease, (8) stroke, (9) dementia, (10) poor kidney function, (11) organ transplant, 
	(12) asplenia, (13) other immunosuppression.
*/


*(1) respiratory disease
tab chronic_respiratory_disease

*think I need level "3" of this, as this is asthma that requires OCS
*(2) severe asthma
tab asthmacat
generate asthma_severe=0
replace asthma_severe=1 if asthmacat==3
tab asthma_severe

*(3) cardiac disease
tab chronic_cardiac_disease

*(4) diabetes
tab diabcat
tab diabcat, nolabel
generate dm=0
replace dm=1 if diabcat>1
tab dm
tab dm diabcat

*(5) non-haem cancer (in previous year)
tab cancer_exhaem
tab cancer_exhaem, nolab
generate cancer_nonhaemPrevYear=0
replace cancer_nonhaemPrevYear=1 if cancer_exhaem_cat==2
tab cancer_nonhaemPrevYear

*(6) haem cancer (within previous 5 years)
tab cancer_haem
tab cancer_haem, nolab
generate cancer_haemPrev5Years=0
replace cancer_haemPrev5Years=1 if inrange(cancer_haem_cat,2,3)
tab cancer_haemPrev5Years

*(7) liver disease
tab chronic_liver_disease

*(8 and 9) stroke or dementia
tab stroke_dementia

*(10) poor kidney function
tab reduced_kidney_function_cat2
tab reduced_kidney_function_cat2, nolabel
gen egfr60 = 0
replace egfr60 = 1 if reduced_kidney_function_cat2 > 1
tab egfr60 reduced_kidney_function_cat2

*(11) organ transplant
tab organ_transplant

*(12) asplenia
tab spleen

*(13) other immunosuppression
tab other_immuno

*create renal flag
gen renal_flag=0
replace renal_flag=1 if inlist(reduced_kidney_function_cat2,2,3,4,5)
tab renal_flag reduced_kidney_function_cat2
replace renal_flag=1 if organ_transplant==1

*create a total comborb var
order chronic_respiratory_disease asthma_severe chronic_cardiac_disease dm cancer_nonhaemPrevYear cancer_haemPrev5Years chronic_liver_disease stroke_dementia egfr60 organ_transplant spleen other_immuno
egen totComorbsOfInterest=rowtotal(chronic_respiratory_disease - other_immuno)
summ totComorbsOfInterest, d
order totComorbsOfInterest

*create the covariate var I need
generate comorb_cat=.
replace comorb_cat=0 if totComorbsOfInterest==0
replace comorb_cat=1 if totComorbsOfInterest==1
replace comorb_cat=2 if totComorbsOfInterest>1

label define comorb_catLab 	0 "No comorbidity" 	///
							1 "1 comorbidity"		///
							2 "2+ comorbidities"			
label values comorb_cat comorb_catLab
tab comorb_cat, m



********************************
*  Outcomes and survival time  *
********************************

/*  Censoring dates  */

noi di "REMEMBER TO UPDATE DATE OF EC DATA UPLOAD"
gen ec_data_date = date("28jan2022", "DMY")
gen ec_data_cens = ec_data_date-7				// Censor AE data 1 week prior to data upload

gen risk_14_days = study_start+14
gen ae_14_pop = (risk_14_days <= ec_data_cens)	// Indicator for has 14-days follow-up

/*

noi di "REMEMBER TO UPDATE DATE OF ONS DATA UPLOAD"
gen ons_data_date = date("5may2021", "DMY")
gen ons_data_cens = ons_data_date-14			// Censor 14 days prior to ONS death data upload
gen risk_28_days = study_start+28
gen risk_40_days = study_start+40

* 28-day risk population
gen risk_pop = (risk_28_days <= ons_data_cens)	// Indicator for has 28-days follow-up
gen time_check_28 = ons_data_cens-study_start
summ time_check_28 if risk_pop == 1, d

* 40-day risk population
gen risk_pop_40 = (risk_40_days <= ons_data_cens)	// Indicator for has 40-days follow-up
gen time_check_40 = ons_data_cens-study_start
summ time_check_40 if risk_pop_40 == 1, d

*/

/*   Outcomes   */

gen ae_pre_cens = (ae_covid_date < ec_data_cens)
gen ae_time = ae_covid_date-study_start if ae_pre_cens == 1
summ ae_time, d

gen any_ae = (ae_covid_date < .)
gen all_ae = (ae_any_date < .)
gen died = (died_date_ons < .)

gen ae_dest = "Home" if inlist(ae_destination,306689006,306691003,306694006)
replace ae_dest = "Police/legal" if inlist(ae_destination,306705005,50861005)
replace ae_dest = "Admitted" if inlist(ae_destination,306706006,1066371000000106,1066391000000105,1066341000000100,1066361000000104,1066331000000109)
replace ae_dest = "Baby unit" if inlist(ae_destination,1066401000000108,1066381000000108)
replace ae_dest = "Transfer/hospice" if inlist(ae_destination,1066351000000102,183919006,19712007)
replace ae_dest = "Mortuary" if inlist(ae_destination,305398007)

tab ae_destination ae_dest if all_ae == 1, m
tab ae_destination ae_dest if any_ae == 1, m


gen ae_admit = 0 if !missing(ae_dest)
replace ae_admit = 1 if ae_dest == "Admitted"

tab ae_admit

/*

* 28-day risk indicator
gen died_pre_cens = (died_date_ons < ons_data_cens)
gen death_time = died_date_ons-study_start if died_pre_cens == 1
summ death_time, d

gen risk_28 = (death_time <= 28)
tab died_pre_cens risk_28, row

* 40-day risk indicator
gen risk_40 = (death_time <= 40)
tab died_pre_cens risk_40, row

*/



/* Survival time */

* Censoring date for Cox
gen cox_pop = (study_start < ec_data_cens)			// Include if data before EC data censor

* Censor at death or EC data censor
gen ae_surv_d = min(ae_covid_date, died_date_ons, ec_data_cens)
gen ae_surv_d1 = ae_surv_d+1
gen ae_surv_d14 = min(ae_covid_date, died_date_ons, ec_data_cens, study_start+14)

gen cox_ae = (ae_covid_date < .)
replace cox_ae = 0 if (ae_covid_date > ae_surv_d)

tab cox_ae, m

* Censor at 14-days
gen cox_ae14 = (ae_covid_date < .)
replace cox_ae14 = 0 if (ae_covid_date > ae_surv_d14)

gen cox_admit = cox_ae
replace cox_admit = 0 if ae_admit != 1

tab cox_admit, m

gen cox_ae_time = ae_surv_d-study_start

tab ae_destination ae_dest if cox_ae == 1, m

/*

* Censor at death, ons data censor, or 7 days prior to vaccine
gen stime_death = min(died_date_ons, ons_data_cens, vacc_cens)
gen stime_death28 = min(died_date_ons, ons_data_cens, vacc_cens, study_start+28)

gen cox_death = (died_date_ons < .)
replace cox_death = 0 if (died_date_ons > stime_death)

gen cox_death28 = (died_date_ons < .)
replace cox_death28 = 0 if (died_date_ons > stime_death28)

gen cox_time = stime_death-study_start
gen cox_time_d = stime_death-study_start if cox_death==1


*/

format ec_data_date ec_data_cens ae_surv_d ae_surv_d1 ae_surv_d14 %td
		
*********************
*  Label variables  *
*********************

* Remove old labels
foreach var of varlist _all {
	label var `var' ""
}

* Demographics
label var patient_id					"Patient ID"
label var start_week					"Epidemiological week of study"

label var age 							"Age (years)"
label var agegroup						"Grouped age"
label var agegroupA						"Age subgroups"
label var agegroupB						"Age subgroups"
label var agegroup6						"Six age groups"
label var agegroup3						"Three age groups"
label var age70 						"70 years and older"
label var age1 							"Age65 spline 1"
label var age2 							"Age65 spline 2"
label var age3 							"Age65 spline 3"
label var male 							"Male"
label var bmi 							"Body Mass Index (BMI, kg/m2)"
label var bmicat 						"Grouped BMI"
label var bmi_date  					"Body Mass Index (BMI, kg/m2), date measured"
label var obese4cat						"Evidence of obesity (missing set to none)"
label var smoke		 					"Smoking status"
label var smoke_nomiss	 				"Smoking status (missing set to never)"
label var smoke_nomiss2	 				"Binary smoking status (missing set to never)"
label var imd 							"Index of Multiple Deprivation (IMD)"
label var eth5							"Ethnicity in 5 categories"
label var eth2							"Ethnicity in 2 catergories"
label var ethnicity_16					"Ethnicity in 16 categories"
label var ethnicity_16_combinemixed		"Ethnicity detailed with mixed groups combined"
label var stp 							"Sustainability and Transformation Partnership"
label var region 						"NHS England region"
label var utla_group					"UTLA"
label var rural_urban					"Rural urban classification"
label var rural_urban5					"Rural Urban in five categories"

*label var household_size				"Household size"
label var hh_total_cat					"Categorical household size"
*label var care_home_type				"Care home status"
label var home_bin						"Binary care home status"

/*
label var hba1ccat						"Categorised hba1c"
label var egfr_cat						"Calculated eGFR"
label var bp_sys 						"Systolic blood pressure"
label var bp_sys_date 					"Systolic blood pressure, date"
label var bp_dias 						"Diastolic blood pressure"
label var bp_dias_date 					"Diastolic blood pressure, date"
label var bpcat 						"Grouped blood pressure"
label var bphigh						"Binary high (stage 1/2) blood pressure"
label var htdiag_or_highbp				"Diagnosed hypertension or high blood pressure"
*/

/*
label var c_age							"Centred age"
label var c_male 						"Centred sex (code: -1/+1)"
label var c_imd							"Centred Index of Multiple Deprivation (values: -2/+2)"
label var c_ethnicity					"Centred ethnicity (values: -2/+2)"
*/

* Comorbidities
/*
label var chronic_respiratory_disease	"Respiratory disease (excl. asthma)"
label var asthmacat						"Asthma, grouped by severity (OCS use)"
label var asthma						"Asthma"
label var chronic_cardiac_disease		"Heart disease"
label var diabetes						"Diabetes"
label var diabcat						"Diabetes, grouped"
label var cancer_exhaem_cat				"Cancer (exc. haematological), grouped by time since diagnosis"
label var cancer_haem_cat				"Haematological malignancy, grouped by time since diagnosis"
label var chronic_liver_disease			"Chronic liver disease"
label var stroke_dementia				"Stroke or dementia"
label var other_neuro					"Neuro condition other than stroke/dementia"	
label var reduced_kidney_function_cat	"Reduced kidney function" 
label var organ_transplant 				"Organ transplant recipient"
label var dysplenia						"Dysplenia (splenectomy, other, not sickle cell)"
label var sickle_cell 					"Sickle cell"
label var spleen						"Spleen problems (dysplenia, sickle cell)"
label var ra_sle_psoriasis				"RA, SLE, Psoriasis (autoimmune disease)"
label var aplastic_anaemia				"Aplastic anaemia"
label var hiv 							"HIV"
label var permanent_immunodeficiency 	"Permanent immunodeficiency"
label var temporary_immunodeficiency 	"Temporary immunosuppression"
label var other_immunosuppression		"Immunosuppressed (combination algorithm)"
label var chronic_respiratory_disease_date	"Respiratory disease (excl. asthma), date"
label var chronic_cardiac_disease_date	"Heart disease, date"
label var diabetes_date					"Diabetes, date"
label var lung_cancer_date				"Lung cancer, date"
label var haem_cancer_date				"Haem. cancer, date"
label var other_cancer_date				"Any cancer, date"
label var chronic_liver_disease_date	"Liver, date"
label var stroke_date					"Stroke, date"
label var dementia_date					"Dementia, date"
label var other_neuro_date				"Neuro condition other than stroke/dementia, date"	
label var organ_transplant_date			"Organ transplant recipient, date"
label var dysplenia_date				"Splenectomy etc, date"
label var sickle_cell_date 				"Sickle cell, date"
label var ra_sle_psoriasis_date			"RA, SLE, Psoriasis (autoimmune disease), date"
label var aplastic_anaemia_date			"Aplastic anaemia, date"
label var hiv_date 						"HIV, date"
label var permanent_immunodeficiency_date "Permanent immunodeficiency, date"
label var temporary_immunodeficiency_date "Temporary immunosuppression, date"
label var dialysis						"Dialysis"
*/
label var comorb_cat					"Categorical number of comorbidites"
label var renal_flag					"Flag for renal disease"

	
* Outcomes and follow-up



label var study_start					"Date of study entry"
label var study_end						"Date of last entry"
label var cox_pop						"1=Population for Cox analysis"
label var ec_data_cens					"EC data censor"
label var ae_covid_date					"Raw date of AE"
label var died_date_ons					"ONS death date"
label var ae_surv_d						"Cox survival time date"
label var ae_surv_d1					"Cox survival time date plus 1 day"
label var ae_surv_d14					"Cox survival time (censored at 14-days)"
label var cox_ae						"AE outcome for Cox"
label var cox_ae14						"AE outcome for Cox (censored at 14-days)"
label var cox_admit						"AE admission outcome for Cox"
label var cox_ae_time					"Cox follow-up time"
label var any_ae						"AE covid any time"
label var ae_pre_cens					"AE covid pre-censor"
label var ae_time						"AE time from start"
label var all_ae						"All cause AE"
label var ae_destination				"AE Destination (numeric)"
label var ae_dest						"AE Destination"
label var ae_admit						"Binary AE admission destination"
label var died							"Died"

label var ae_14_pop						"1=Population for 14-day follow-up"

/*

label var risk_pop						"1=Population for 28-day risk analysis"
label var risk_pop_40					"1=Population for 40-day risk analysis"
label var risk_28						"28-day outcome"
label var risk_40						"40-day outcome"
label var stime_death					"Date of study exit"
label var stime_death28					"Date of study exit (28-day censor)"
label var stime_comp_death				"Date of study exit, discharge as censor"
label var stime_hosp_test				"Date of exit (hospital admission)"
label var stime_icu_test				"Date of exit (icu admission)"
label var cox_death28					"Outcome for Cox (28-day censor)"
label var cox_time						"Follow-up time"
label var cox_time_d					"Time to death"
label var end_death_hosp				"Outcome death|hosp"
label var comp_death_hosp				"Discharge censor outcome death|hosp"
label var end_hosp_test					"Outcome hosp|test"
label var end_icu_test					"Outcome icu|test"
label var end_death_icu					"Outcome death|icu"
label var time_death_hosp				"Follow-up time (death|hosp)"
label var time_death_hosp1				"Follow-up time (death|hosp)"
label var time_comp_death				"Follow-up time (death|hosp) discharge as censor"
label var time_hosp_test				"Follow-up time (hosp|test)"
label var time_icu_test					"Follow-up time (icu|test)"
label var time_death_icu				"Follow-up time (death|icu)"
label var time_death_icu1				"Follow-up time (death|icu)"
label var hosp_28						"28-day hospitalisation outcome"

*/

label var sgtf							"SGTF (exposure)"
label var has_sgtf						"1=Has SGTF data"

label var vax							"Vaccination status"
label var prev_inf						"Prior infection status"

/*
label var covid_admission_date			"Date of hospital admission" 
label var icu_admission_date			"Date of icu admission" 
label var icu_pop						"Population with days spent in ICU"
label var covid_icu_days				"Days spent in ICU"
label var covid_discharge_date			"Date of hospital discharge" 
label var hosp_discharge				"Discharged from hospital"
label var death_inout					"Death with or without hospital admission"


* Deaths before exclusions
tab cox_death utla_group if cox_pop==1 & has_sgtf==1 & missing(utla_group), m
tab cox_death eth2 if cox_pop==1 & has_sgtf==1, m
tab cox_death hh_total_cat if cox_pop==1 & has_sgtf==1, m
tab cox_death comorb_cat if cox_pop==1 & has_sgtf==1, m
tab cox_death start_week if cox_pop==1 & has_sgtf==1, m

*/

********************
*  Missing values  *
********************

*tab smoke sgtf if cox_pop==1 & has_sgtf==1 & utla_group !="", m
*tab bmicat sgtf if cox_pop==1 & has_sgtf==1 & utla_group !="", m
*tab care_home_type sgtf if cox_pop==1 & has_sgtf==1 & utla_group !="", m


***************
*  Tidy data  *
***************

* REDUCE DATASET SIZE TO VARIABLES NEEDED - KEEP THOSE WITH LABELS
ds , has(varl)
keep `r(varlist)'


***************
*  Save data  *
***************

sort patient_id
order patient_id sgtf

label data "SGTF-OMICRON ANALYSIS DATASET: $S_DATE"

save ./output/main.dta, replace

log close

