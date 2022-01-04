# IMPORT STATEMENTS
# This imports the cohort extractor package. This can be downloaded via pip
from cohortextractor import (
    StudyDefinition,
    patients,
    codelist_from_csv,
    codelist,
    combine_codelists,
    filter_codes_by_category,
)

# dictionary of MSOA codes (for dummy data)
from dictionaries import dict_msoa

# IMPORT CODELIST DEFINITIONS FROM CODELIST.PY (WHICH PULLS THEM FROM
# CODELIST FOLDER
from codelists import *


# STUDY DEFINITION
# Defines both the study population and points to the important covariates and outcomes
study = StudyDefinition(
    default_expectations={
        "date": {"earliest": "1970-01-01", "latest": "today"},
        "rate": "uniform",
        "incidence": 0.2,
    },

   # STUDY POPULATION - registered 1 year prior to positive test in SGSS
    population=patients.satisfying(
        "one_practice AND sgss_pos_inrange",
        one_practice=patients.registered_with_one_practice_between(
            "sgss_pos_inrange - 1 year", "sgss_pos_inrange"
        ),
    ),

    dereg_date=patients.date_deregistered_from_all_supported_practices(
        on_or_after="2021-11-01", date_format="YYYY-MM",
    ),


    ### SGSS positive in study period
    sgss_pos_inrange=patients.with_test_result_in_sgss(
       pathogen="SARS-CoV-2",
       test_result="positive",
       find_first_match_in_period=True,
       between=["2021-10-03", "2022-01-01"],
       returning="date",
       date_format="YYYY-MM-DD",
       return_expectations={"date": {"earliest": "2021-10-03", "latest": "today"},
                            "incidence": 0.9
       },
    ), 

    sgtf=patients.with_test_result_in_sgss(
       pathogen="SARS-CoV-2",
       test_result="positive",
       find_first_match_in_period=True,
       between=["2021-10-03", "2022-01-01"],
       returning="s_gene_target_failure",
       return_expectations={
            "rate": "universal",
            "category": {"ratios": {"0": 0.7, "1": 0.1, "9": 0.1, "": 0.1}},
       },
    ), 
    
    
    ### OUTCOMES - ONS death dates
    died_date_ons=patients.died_from_any_cause(
        on_or_after="2020-02-01",
        returning="date_of_death",
        include_month=True,
        include_day=True,
        return_expectations={"date": {"earliest": "2021-10-01"}, "incidence" : 0.1},
    ),
    
    ### AE ATTENDANCE
    ae_covid_date=patients.attended_emergency_care(
        returning= "date_arrived",
        with_these_diagnoses=snomed_covid,
        on_or_after="sgss_pos_inrange",
        find_first_match_in_period=True,  
        date_format="YYYY-MM-DD",  
        return_expectations={"date": {"earliest": "2021-10-02"}, "incidence" : 0.2},
    ),    

    ae_destination=patients.attended_emergency_care(
        returning= "discharge_destination",
        with_these_diagnoses=snomed_covid,
        on_or_after="sgss_pos_inrange",
        find_first_match_in_period=True,  
        return_expectations={
            "rate": "universal",
            "category": 
                {"ratios": {
                    "1": 0.1,
                    "2": 0.1,
                    "3": 0.1,
                    "4": 0.1,
                    "5": 0.1,
                    "6": 0.1,
                    "7": 0.2,
                    "8": 0.2,
                }
            },
        },
    ),
    
    ae_any_date=patients.attended_emergency_care(
        returning= "date_arrived",
        on_or_after="sgss_pos_inrange",
        find_first_match_in_period=True,  
        date_format="YYYY-MM-DD",  
        return_expectations={"date": {"earliest": "2021-10-02"}, "incidence" : 0.3},
    ),


    ### COVID Vaccines

    ## Any covid vaccination
    vaxdate1=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        on_or_after="2020-12-08",
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": "2020-12-08", "latest": "today"},
            "incidence": 0.8
        },
    ),
    vaxdate2=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        on_or_after="vaxdate1 + 1 day",
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": "2021-01-08", "latest" : "today"},
            "incidence": 0.7
        },
    ),
    vaxdate3=patients.with_tpp_vaccination_record(
        target_disease_matches="SARS-2 CORONAVIRUS",
        on_or_after="vaxdate2 + 1 day",
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD",
         return_expectations={
            "date": {"earliest": "2021-02-08", "latest" : "today"},
            "incidence": 0.5
        },
    ),
    
    ### COVID infections
    
    # Primary care
    last_covid_tpp_probable=patients.with_these_clinical_events(
        combine_codelists(covid_identification_in_primary_care_case_codes_clinical,
                          covid_identification_in_primary_care_case_codes_test,
                          covid_identification_in_primary_care_case_codes_seq),
        return_last_date_in_period=True,
        between=["2020-08-01", "sgss_pos_inrange - 7 days"],
        include_day=True,
        return_expectations={"date": {"earliest": "2020-11-01"}, "incidence" : 0.2},
    ),
    
    # SGSS
    last_pos_test_sgss=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="positive",
        find_last_match_in_period=True,
        returning="date",
        between=["2020-08-01", "sgss_pos_inrange - 1 days"],       
        date_format="YYYY-MM-DD",
        return_expectations={"date": {"earliest": "2020-11-01"},
                            "incidence": 0.1
        },
    ), 

    ### DEMOGRAPHIC COVARIATES
    # AGE
    age=patients.age_as_of(
        "sgss_pos_inrange",
        return_expectations={
            "rate": "universal",
            "int": {"distribution": "population_ages"},
        },
    ),

    # SEX
    sex=patients.sex(
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {"M": 0.49, "F": 0.51}},
        }
    ),

    # DEPRIVIATION
    imd=patients.address_as_of(
        "sgss_pos_inrange",
        returning="index_of_multiple_deprivation",
        round_to_nearest=100,
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {"100": 0.1, "200": 0.2, "300": 0.2, "400": 0.2, "500": 0.2, "600": 0.1}},
        },
    ),

    # GEOGRAPHIC REGION CALLED STP
    stp=patients.registered_practice_as_of(
        "sgss_pos_inrange",
        returning="stp_code",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": {
                    "STP1": 0.1,
                    "STP2": 0.1,
                    "STP3": 0.1,
                    "STP4": 0.1,
                    "STP5": 0.1,
                    "STP6": 0.1,
                    "STP7": 0.1,
                    "STP8": 0.1,
                    "STP9": 0.1,
                    "STP10": 0.1,
                }
            },
        },
    ),

    # GEOGRAPHIC REGION MSOA
    msoa=patients.registered_practice_as_of(        
        "sgss_pos_inrange",
        returning="msoa_code",
        return_expectations={
            "rate": "universal",
            "category": {"ratios": dict_msoa},
        },
    ),

    # REGION - one of NHS England 9 regions
    region=patients.registered_practice_as_of(
        "sgss_pos_inrange",
        returning="nuts1_region_name",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": {
                    "North East": 0.1,
                    "North West": 0.05,
                    "Yorkshire and The Humber": 0.1,
                    "East Midlands": 0.1,
                    "West Midlands": 0.05,
                    "East of England": 0.3,
                    "London": 0.15,
                    "South East": 0.1,
                    "South West": 0.05,
                },
            },
        },
    ),

    # RURAL OR URBAN LOCATION
    rural_urban=patients.address_as_of(
        "sgss_pos_inrange",
        returning="rural_urban_classification",
        return_expectations={
            "rate": "universal",
            "category": 
                {"ratios": {
                    "1": 0.1,
                    "2": 0.1,
                    "3": 0.1,
                    "4": 0.1,
                    "5": 0.1,
                    "6": 0.1,
                    "7": 0.2,
                    "8": 0.2,
                }
            },
        },
    ),

    # HOUSEHOLD INFORMATION
    household_id=patients.household_as_of(
        "2020-02-01",
        returning="pseudo_id",
        return_expectations={
            "int": {"distribution": "normal", "mean": 1000, "stddev": 200},
            "incidence": 1,
        },
    ),

    household_size=patients.household_as_of(
        "2020-02-01",
        returning="household_size",
        return_expectations={
            "int": {"distribution": "normal", "mean": 3, "stddev": 1},
            "incidence": 1,
        },
    ),

    care_home_type=patients.care_home_status_as_of(
        "sgss_pos_inrange",
        categorised_as={
            "PC": """
              IsPotentialCareHome
              AND LocationDoesNotRequireNursing='Y'
              AND LocationRequiresNursing='N'
            """,
            "PN": """
              IsPotentialCareHome
              AND LocationDoesNotRequireNursing='N'
              AND LocationRequiresNursing='Y'
            """,
            "PS": "IsPotentialCareHome",
            "U": "DEFAULT",
        },
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {"PC": 0.05, "PN": 0.05, "PS": 0.05, "U": 0.85,},},
        },
    ), 

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/10
    bmi=patients.most_recent_bmi(
        between=["2010-02-01", "sgss_pos_inrange"],
        minimum_age_at_measurement=16,
        include_measurement_date=True,
        include_month=True,
        return_expectations={
            "date": {},
            "float": {"distribution": "normal", "mean": 35, "stddev": 10},
            "incidence": 0.95,
        },
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/6
    smoking_status=patients.categorised_as(
        {
            "S": "most_recent_smoking_code = 'S'",
            "E": """
                 most_recent_smoking_code = 'E' OR (
                   most_recent_smoking_code = 'N' AND ever_smoked
                 )
            """,
            "N": "most_recent_smoking_code = 'N' AND NOT ever_smoked",
            "M": "DEFAULT",
        },
        return_expectations={
            "category": {"ratios": {"S": 0.6, "E": 0.1, "N": 0.2, "M": 0.1}}
        },
        most_recent_smoking_code=patients.with_these_clinical_events(
            clear_smoking_codes,
            find_last_match_in_period=True,
            on_or_before="sgss_pos_inrange",
            returning="category",
        ),
        ever_smoked=patients.with_these_clinical_events(
            filter_codes_by_category(clear_smoking_codes, include=["S", "E"]),
            on_or_before="sgss_pos_inrange",
        ),
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/27
    ethnicity=patients.with_these_clinical_events(
        ethnicity_codes,
        returning="category",
        find_last_match_in_period=True,
        include_date_of_match=True,
        return_expectations={
            "category": {"ratios": {"1": 0.75, "2": 0.05, "3": 0.05, "4": 0.05, "5": 0.1}},
            "incidence": 0.75,
        },
    ),
    ethnicity_16=patients.with_these_clinical_events(
        ethnicity_codes_16,
        returning="category",
        find_last_match_in_period=True,
        include_date_of_match=True,
        return_expectations={
            "category": {"ratios": {"1": 0.8, "5": 0.1, "3": 0.1}},
            "incidence": 0.75,
        },
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/21
    chronic_respiratory_disease=patients.with_these_clinical_events(
        chronic_respiratory_disease_codes,
        return_first_date_in_period=True,
        include_month=True,
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/55
    asthma=patients.categorised_as(
        {
            "0": "DEFAULT",
            "1": """
                (
                  recent_asthma_code OR (
                    asthma_code_ever AND NOT
                    copd_code_ever
                  )
                ) AND (
                  prednisolone_last_year = 0 OR 
                  prednisolone_last_year > 4
                )
            """,
            "2": """
                (
                  recent_asthma_code OR (
                    asthma_code_ever AND NOT
                    copd_code_ever
                  )
                ) AND
                prednisolone_last_year > 0 AND
                prednisolone_last_year < 5
                
            """,
        },
        return_expectations={"category": {"ratios": {"0": 0.8, "1": 0.1, "2": 0.1}},},
        recent_asthma_code=patients.with_these_clinical_events(
            asthma_codes, between=["2018-02-01", "sgss_pos_inrange"],
        ),
        asthma_code_ever=patients.with_these_clinical_events(asthma_codes),
        copd_code_ever=patients.with_these_clinical_events(
            chronic_respiratory_disease_codes
        ),
        prednisolone_last_year=patients.with_these_medications(
            pred_codes,
            between=["2020-11-01", "sgss_pos_inrange"],
            returning="number_of_matches_in_period",
        ),
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/7
    chronic_cardiac_disease=patients.with_these_clinical_events(
        chronic_cardiac_disease_codes,
        return_first_date_in_period=True,
        include_month=True,
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/30
    diabetes=patients.with_these_clinical_events(
        diabetes_codes, return_first_date_in_period=True, include_month=True,
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/32
    lung_cancer=patients.with_these_clinical_events(
        lung_cancer_codes, return_first_date_in_period=True, include_month=True,
    ),
    haem_cancer=patients.with_these_clinical_events(
        haem_cancer_codes, return_first_date_in_period=True, include_month=True,
    ),
    other_cancer=patients.with_these_clinical_events(
        other_cancer_codes, return_first_date_in_period=True, include_month=True,
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/12
    chronic_liver_disease=patients.with_these_clinical_events(
        chronic_liver_disease_codes,
        return_first_date_in_period=True,
        include_month=True,
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/14
    other_neuro=patients.with_these_clinical_events(
        other_neuro, return_first_date_in_period=True, include_month=True,
    ),
    stroke=patients.with_these_clinical_events(
        stroke, return_first_date_in_period=True, include_month=True,
    ),
    dementia=patients.with_these_clinical_events(
        dementia, return_first_date_in_period=True, include_month=True,
    ),

    # Chronic kidney disease
    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/17
    creatinine=patients.with_these_clinical_events(
        creatinine_codes,
        find_last_match_in_period=True,
        on_or_before="sgss_pos_inrange",
        returning="numeric_value",
        include_date_of_match=True,
        include_month=True,
        return_expectations={
            "float": {"distribution": "normal", "mean": 60.0, "stddev": 15},
            "date": {"earliest": "2019-02-28", "latest": "2021-11-16"},
            "incidence": 0.95,
        },
    ),
    dialysis=patients.with_these_clinical_events(
        dialysis_codes, return_first_date_in_period=True, include_month=True,
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/31
    organ_transplant=patients.with_these_clinical_events(
        organ_transplant_codes, return_first_date_in_period=True, include_month=True,
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/13
    dysplenia=patients.with_these_clinical_events(
        spleen_codes, return_first_date_in_period=True, include_month=True,
    ),
    sickle_cell=patients.with_these_clinical_events(
        sickle_cell_codes, return_first_date_in_period=True, include_month=True,
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/36
    aplastic_anaemia=patients.with_these_clinical_events(
        aplastic_codes, return_last_date_in_period=True, include_month=True,
    ),
    hiv=patients.with_these_clinical_events(
        hiv_codes,
        returning="category", 
        find_first_match_in_period=True, 
        include_date_of_match=True,
        include_month=True,
        return_expectations={
            "category": {"ratios": {"43C3.": 0.8, "XaFuL": 0.2}},
            },
    ),   
    permanent_immunodeficiency=patients.with_these_clinical_events(
        permanent_immune_codes, return_first_date_in_period=True, include_month=True,
    ),
    temporary_immunodeficiency=patients.with_these_clinical_events(
        temp_immune_codes, return_last_date_in_period=True, include_month=True,
    ),

    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/23
    # immunosuppressant_med=
    # hypertension
    hypertension=patients.with_these_clinical_events(
        hypertension_codes, return_first_date_in_period=True, include_month=True,
    ),

    # Blood pressure
    # https://github.com/ebmdatalab/tpp-sql-notebook/issues/35
    bp_sys=patients.mean_recorded_value(
        systolic_blood_pressure_codes,
        on_most_recent_day_of_measurement=True,
        on_or_before="sgss_pos_inrange",
        include_measurement_date=True,
        include_month=True,
        return_expectations={
            "float": {"distribution": "normal", "mean": 80, "stddev": 10},
            "date": {"latest": "2021-11-16"},
            "incidence": 0.95,
        },
    ),
    bp_dias=patients.mean_recorded_value(
        diastolic_blood_pressure_codes,
        on_most_recent_day_of_measurement=True,
        on_or_before="sgss_pos_inrange",
        include_measurement_date=True,
        include_month=True,
        return_expectations={
            "float": {"distribution": "normal", "mean": 120, "stddev": 10},
            "date": {"latest": "2021-11-16"},
            "incidence": 0.95,
        },
    ),
    hba1c_mmol_per_mol=patients.with_these_clinical_events(
        hba1c_new_codes,
        find_last_match_in_period=True,
        on_or_before="sgss_pos_inrange",
        returning="numeric_value",
        include_date_of_match=True,
        include_month=True,
        return_expectations={
            "date": {"latest": "2021-11-16"},
            "float": {"distribution": "normal", "mean": 40.0, "stddev": 20},
            "incidence": 0.95,
        },
    ),
    hba1c_percentage=patients.with_these_clinical_events(
        hba1c_old_codes,
        find_last_match_in_period=True,
        on_or_before="sgss_pos_inrange",
        returning="numeric_value",
        include_date_of_match=True,
        include_month=True,
        return_expectations={
            "date": {"latest": "2021-11-16"},
            "float": {"distribution": "normal", "mean": 5, "stddev": 2},
            "incidence": 0.95,
        },
    ),

    # # https://github.com/ebmdatalab/tpp-sql-notebook/issues/49
    ra_sle_psoriasis=patients.with_these_clinical_events(
        ra_sle_psoriasis_codes, return_first_date_in_period=True, include_month=True,
    ),
)
