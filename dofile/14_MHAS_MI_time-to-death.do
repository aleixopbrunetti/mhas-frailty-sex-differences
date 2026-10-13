********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Multiple imputation - terminal decline among decedents (time-to-death scale)
********************************************************************************

* ==============================================================================
* 3. Load Data
* ==============================================================================

use "$cleanmhasdata/wave1to5_wide_deaths_code.dta", clear

* Set seed
set seed 14122025

* ==============================================================================
* 4. Prepare Time-to-Death Variables (before imputation)
* ==============================================================================

* Baseline year for wave 1
local baseline_year 2001

* Interview year for each wave
gen intyear1 = `baseline_year' + 0   // Wave 1: 2001
gen intyear2 = `baseline_year' + 2   // Wave 2: 2003
gen intyear3 = `baseline_year' + 11  // Wave 3: 2012
gen intyear4 = `baseline_year' + 14  // Wave 4: 2015
gen intyear5 = `baseline_year' + 17  // Wave 5: 2018

* Time to death at each wave (scaled: decades)
gen ttd1 = (dyear - intyear1) / 10
gen ttd2 = (dyear - intyear2) / 10
gen ttd3 = (dyear - intyear3) / 10
gen ttd4 = (dyear - intyear4) / 10
gen ttd5 = (dyear - intyear5) / 10

label var ttd1 "Time to death at Wave 1 (years)"
label var ttd2 "Time to death at Wave 2 (years)"
label var ttd3 "Time to death at Wave 3 (years)"
label var ttd4 "Time to death at Wave 4 (years)"
label var ttd5 "Time to death at Wave 5 (years)"

* ==============================================================================
* 5. Structural missingness: FI undefined after death
* ==============================================================================

* Eligible at wave = died after that wave's interview year
gen elig_w1 = (dyear > intyear1)
gen elig_w2 = (dyear > intyear2)
gen elig_w3 = (dyear > intyear3)
gen elig_w4 = (dyear > intyear4)
gen elig_w5 = (dyear > intyear5)

label var elig_w1 "Eligible for Wave 1 FI (alive at wave)"
label var elig_w2 "Eligible for Wave 2 FI (alive at wave)"
label var elig_w3 "Eligible for Wave 3 FI (alive at wave)"
label var elig_w4 "Eligible for Wave 4 FI (alive at wave)"
label var elig_w5 "Eligible for Wave 5 FI (alive at wave)"

* Set post-death FI to missing (structural missingness)
replace w1fi = . if elig_w1==0
replace w2fi = . if elig_w2==0
replace w3fi = . if elig_w3==0
replace w4fi = . if elig_w4==0
replace w5fi = . if elig_w5==0

* ==============================================================================
* 5.1 Missingness diagnostics among decedents (structural vs non-response)
* ==============================================================================

display as text "=========================================================================="
display as text "MISSINGNESS DIAGNOSTICS FOR DECEDENTS"
display as text "=========================================================================="

* Total decedents
count if death==1
local n_decedents = r(N)
display as text "Total decedents: `n_decedents'"
display ""

* By wave: total missing, structural (died before wave), non-response (eligible but missing), observed
foreach wave in 1 2 3 4 5 {
    display as text "---------- Wave `wave' ----------"
    
    * Total missing
    count if death==1 & missing(w`wave'fi)
    local miss_total = r(N)
    local miss_pct = (`miss_total' / `n_decedents') * 100
    
    * Structural missingness (died before wave)
    count if death==1 & missing(w`wave'fi) & elig_w`wave'==0
    local miss_structural = r(N)
    local struct_pct = (`miss_structural' / `n_decedents') * 100
    
    * Non-response (eligible but missing)
    count if death==1 & missing(w`wave'fi) & elig_w`wave'==1
    local miss_nonresp = r(N)
    local nonresp_pct = (`miss_nonresp' / `n_decedents') * 100
    
    * Observed (eligible and observed)
    count if death==1 & !missing(w`wave'fi) & elig_w`wave'==1
    local obs_data = r(N)
    local obs_pct = (`obs_data' / `n_decedents') * 100
    
    display as text "  Total missing: `miss_total' (`miss_pct'%)"
    display as result "    Structural (died before wave): `miss_structural' (`struct_pct'%)"
    display as result "    Non-response (alive but missing): `miss_nonresp' (`nonresp_pct'%)"
    display as text "  Observed data: `obs_data' (`obs_pct'%)"
    display ""
}

display as text "=========================================================================="
display ""

/* Imputation strategy:
- Impute all waves (w1fi-w5fi) to fill non-response missingness
- After imputation, reset post-death values to missing using eligibility flags
- Each participant therefore only has imputed values for waves before death
- Baseline covariates with missing data (educ, wealth_gbp) also imputed
*/ 

* ==============================================================================
* 6. Multiple Imputation (decedents analysed)
* ==============================================================================

*---------------------------
* 6.1 Set MI structure
*---------------------------

mi set wide
mi svyset [pweight=w1respwt]

mi register regular patid sex w1respwt rural death dyear age ///
    intyear1 intyear2 intyear3 intyear4 intyear5 ///
    ttd1 ttd2 ttd3 ttd4 ttd5 ///
    elig_w1 elig_w2 elig_w3 elig_w4 elig_w5

mi register imputed wealth_gbp educ vigact ///
    smoke_status srh drink_status higov cogstat ///
    w1fi w2fi w3fi w4fi w5fi
    
* Passive variables (defined after imputation)
generate wealth_gbp_c = .
generate age_c = .
generate age_cat = .

mi register passive wealth_gbp_c age_c age_cat

*---------------------------
* 6.2 MICE, stratified by sex
*---------------------------
local Mtarget 40
local Burnin 20

mi impute chained ///
	(pmm, knn(10)) 		w1fi w2fi w3fi w4fi w5fi ///
	(pmm, knn(10)) 		wealth_gbp ///
    (mlogit, augment) 	smoke_status educ ///
    (ologit, augment) 	srh cogstat ///
    (logit, augment) 	vigact drink_status higov ///
    = c.age i.rural [pweight=w1respwt], by(sex) ///
    add(`Mtarget') burnin(`Burnin') dots noisily

*---------------------------
* 6.2.1 Post-imputation: remove imputed values after death
*---------------------------

* Number of imputations
mi describe
local M = r(M)

* Apply to m=0 and all imputations 1..M
mi xeq 0/`M': replace w1fi = . if elig_w1==0
mi xeq 0/`M': replace w2fi = . if elig_w2==0
mi xeq 0/`M': replace w3fi = . if elig_w3==0
mi xeq 0/`M': replace w4fi = . if elig_w4==0
mi xeq 0/`M': replace w5fi = . if elig_w5==0

* ------------------------------------------------------------------------------
* 6.3 Passive variables
* ------------------------------------------------------------------------------
* Wealth centred on observed mean
mi xeq 0: quietly summarize wealth_gbp if !missing(wealth_gbp), meanonly
scalar WEALTH_MEAN = r(mean)
mi passive: replace wealth_gbp_c = wealth_gbp - WEALTH_MEAN

* Age centred on observed mean
mi xeq 0: quietly summarize age if !missing(age), meanonly
scalar AGE_MEAN = r(mean)
mi passive: replace age_c= age - AGE_MEAN

* Age categories
mi passive: replace age_cat = 1 if age >= 50 & age < 60
mi passive: replace age_cat = 2 if age >= 60 & age < 70
mi passive: replace age_cat = 3 if age >= 70 & age < 80
mi passive: replace age_cat = 4 if age >= 80 & !missing(age)

label define age_cat_lbl 1 "50-59" 2 "60-69" 3 "70-79" 4 "80+"
label values age_cat age_cat_lbl
label variable age_cat "Age category"

*-------------------------------------------------------------------------------
* 6.4 Convert data from wide to long format
*-------------------------------------------------------------------------------
* Rename for reshape
rename w1fi wfi1
rename w2fi wfi2
rename w3fi wfi3
rename w4fi wfi4
rename w5fi wfi5

rename ttd1 time_to_death1
rename ttd2 time_to_death2
rename ttd3 time_to_death3
rename ttd4 time_to_death4
rename ttd5 time_to_death5

* Time since baseline for each wave
gen time1 = 0   // Wave 1: 2001 - baseline
gen time2 = 2   // Wave 2: 2003
gen time3 = 11  // Wave 3: 2012
gen time4 = 14  // Wave 4: 2015
gen time5 = 17  // Wave 5: 2018

*-- Reshape to long format
capture drop wave
mi reshape long wfi time time_to_death, i(patid) j(wave)

label variable time "Time since baseline (in years)"
label variable time_to_death "Years until death (0=death year)"

* Drop impossible values (interview after death)
drop if time_to_death < 0

* Sex-specific time variables for level-2 heteroskedasticity
gen female = (sex == 2)
gen male = (sex == 1)
gen time_female = female * time
gen time_male = male * time

* Sex-specific time-to-death variables for random effects
gen female_ttd = female
gen male_ttd = male
gen time_female_ttd = female * time_to_death
gen time_male_ttd = male * time_to_death

label var female_ttd "Female indicator"
label var male_ttd "Male indicator"
label var time_female_ttd "Female x time_to_death"
label var time_male_ttd "Male x time_to_death"

*------------------------------------------------------------------------------
* Keep decedents only
*------------------------------------------------------------------------------
keep if death ==1

* Save MI dataset
save "$cleanmhasdata/wave1to5_long_MI_TTD.dta", replace

mi describe

mi varying

* Check: no observed FI after death
mi xeq 1: tab wave if !missing(wfi) & death==1 & time_to_death < 0, missing

*===============================================================================
* 6. Time-to-death model on MI data (decedents)
*===============================================================================

use "$cleanmhasdata/wave1to5_long_MI_TTD.dta", clear

capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots saving(miests_ttd, replace): ///
	 mixed wfi i.sex ///
	 c.time_to_death c.time_to_death#c.time_to_death c.time_to_death#c.time_to_death#c.time_to_death ///
     i.sex#c.time_to_death c.age_c c.age_c#c.time_to_death ///
     c.wealth_gbp_c i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status i.drink_status#c.time_to_death i.higov ///
     i.cogstat i.cogstat#c.time_to_death ///
     || patid: female_ttd time_female_ttd, noconstant covariance(unstructured) ///
     || patid: male_ttd   time_male_ttd,   noconstant covariance(unstructured) ///
     ml pw(w1respwt) nolog stddev difficult iterate(100) tolerance(1e-4)
		
estimates store mi_ttd_model

* Test sex difference in time-to-death slope
mi test 2.sex#c.time_to_death

* ------------------------------------------------------------------------------
* Margins for time-to-death model (decedents)
* -----------------------------------------------------------------------------
* Requirement: net install mimrgns, from(http://fmwww.bc.edu/RePEc/bocode/m)

* time_to_death in decades: 0, 2, 5, 11, 15, 17 years before death
mimrgns sex, vce(unconditional) at(time_to_death=(0 0.2 0.5 1.1 1.5 1.7))

* Extract and transpose results matrix
matrix margins_sex = r(table)
matrix list margins_sex
matrix margins_sex_t = margins_sex'

* Save margins results
preserve
    clear
    svmat margins_sex_t
    save "$working/margins_MI_sex_TTD.dta", replace
restore

*--- Linear combinations of parameters

* Sex difference in the year of death (time_to_death = 0)
lincom 2.sex

* Female vs male slope difference on the time-to-death scale
lincom 2.sex#c.time_to_death - 1.sex#c.time_to_death

* Female–male difference at selected years before death (time_to_death in decades)
lincom 2.sex + 0.2*(2.sex#c.time_to_death)   // 2 years
lincom 2.sex + 0.5*(2.sex#c.time_to_death)   // 5 years
lincom 2.sex + 1.0*(2.sex#c.time_to_death)   // 10 years
lincom 2.sex + 1.1*(2.sex#c.time_to_death)   // 11 years
lincom 2.sex + 1.5*(2.sex#c.time_to_death)   // 15 years
lincom 2.sex + 1.7*(2.sex#c.time_to_death)   // 17 years

* ---------------------
* Predicted FI by sex with shaded CIs (reversed x-axis: years before death)
* -----------------------

preserve
    use "$working/margins_MI_sex_TTD.dta", clear
    
    * Rename variables
    rename margins_sex_t1 pred_fi
    rename margins_sex_t2 se
    rename margins_sex_t3 t
    rename margins_sex_t4 pvalue
    rename margins_sex_t5 ci_lb
    rename margins_sex_t6 ci_ub
    rename margins_sex_t7 df
    rename margins_sex_t8 crit
    rename margins_sex_t9 eform
    
    * time_to_death and sex identifiers (rows ordered by time point, then sex)
    gen obs = _n
    gen time_to_death = .
    gen sex = .
    
    forvalues i = 1/12 {
        local at_num = ceil(`i'/2)
        local sex_num = mod(`i'-1, 2) + 1
        local timevals "0 0.2 0.5 1.1 1.5 1.7"
        local tval : word `at_num' of `timevals'
        replace time_to_death = `tval' in `i'
        replace sex = `sex_num' in `i'
    }
    
    gen sex_label = "Male" if sex == 1
    replace sex_label = "Female" if sex == 2
    
    twoway ///
        (rarea ci_lb ci_ub time_to_death if sex==1, color(blue%20) lwidth(none)) ///
        (rarea ci_lb ci_ub time_to_death if sex==2, color(red%20) lwidth(none)) ///
        (line pred_fi time_to_death if sex==1, lcolor(blue) lwidth(medium)) ///
        (line pred_fi time_to_death if sex==2, lcolor(red) lwidth(medium)), ///
        title("Predicted FI by Sex Over Time") ///
        ytitle("Predicted Frailty Index") ///
        xtitle("Years before death") ///
        ylabel(, nogrid) ///
        xlabel(1.7 "17" 1.5 "15" 1.1 "11" 0.5 "5" 0.2 "2" 0 "0", nogrid) ///
        xscale(reverse) ///
        legend(order(3 "Male" 4 "Female") position(6) ring(1) cols(2)) ///
        name(fi_sex_ttd_custom_reversed, replace)
        
    graph export "$outputgraphs/mi_margins_sex_TTD_custom_reversed.jpg", replace
    graph export "$outputgraphs/mi_margins_sex_TTD_custom_reversed.tif", replace width(4000)
    
restore

*===============================================================================
* 7. Session end
*===============================================================================
display as text "------------------------------------------------------------"
display as text "MHAS_MI_TTD_completed"
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close