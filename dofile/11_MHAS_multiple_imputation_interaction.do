********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Multiple Imputation for MHAS Frailty Index Study
********************************************************************************

* ==============================================================================
* 3. Load Data
* ==============================================================================

use "$cleanmhasdata/wave1to5_wide.dta", clear

* Set seed
set seed 14122025

* ==============================================================================
* 4 Multiple Imputation for missing data
* ==============================================================================

*---------------------------
* 4.1 Set MI structure
*---------------------------
mi set wide
mi svyset [pweight=w1respwt]

mi register regular patid sex w1respwt rural
mi register imputed wealth_gbp educ vigact ///
    smoke_status srh drink_status higov cogstat ///
    w1fi w2fi w3fi w4fi w5fi
	
* Create passive variables (centred wealth and age, age category)
generate wealth_gbp_c = .
generate age_c        = .
generate age_cat      = .

mi register passive wealth_gbp_c age_c age_cat

*------------------------------------------------------------
* 4.2 Multiple imputation by chained equations (MICE), stratified by sex
*---------------------------------------------------------------
// FI imputed in wide format (one variable per wave) to preserve within-person correlation
// PMM for FI and wealth; mlogit for nominal; ologit for ordered; logit for binary variables

* Number of imputations and burn-in
local Mtarget 40
local Burnin 20

mi impute chained ///
    (pmm, knn(10))          w1fi w2fi w3fi w4fi w5fi  ///
    (pmm, knn(10))          wealth_gbp ///
    (mlogit, augment)       smoke_status educ ///
    (ologit, augment)       srh cogstat ///
    (logit, augment)        vigact drink_status higov ///
    =  c.age i.rural [pweight=w1respwt], by(sex) ///
    add(`Mtarget') burnin(`Burnin') dots noisily

* ------------------------------------------------------------------------------
* 4.3 Passive variables
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
* 4.4 Convert data from wide to long format & create time variable
*-------------------------------------------------------------------------------
* Rename FI variables for reshape
rename w1fi wfi1
rename w2fi wfi2
rename w3fi wfi3
rename w4fi wfi4
rename w5fi wfi5

* Time variables for each wave (years since baseline)
gen time1 = 0   // Wave 1: 2001 - baseline
gen time2 = 2   // Wave 2: 2003
gen time3 = 11  // Wave 3: 2012
gen time4 = 14  // Wave 4: 2015
gen time5 = 17  // Wave 5: 2018

*-- Reshape to long format
drop (wave)
mi reshape long wfi time, i(patid) j(wave) string

label variable time "Time since baseline (in years)"

* Sex-specific time variables for level-2 heteroskedasticity
gen female = (sex == 2)
gen male = (sex == 1)
gen time_female = female * c.time
gen time_male = male * c.time

* Save the complete MI dataset
save "$cleanmhasdata/wave1to5_long_MI.dta", replace

mi describe

*===============================================================================
* 6. Mixed model on MI data - continuous age
*===============================================================================

use "$cleanmhasdata/wave1to5_long_MI.dta", clear

capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots saving(miests, replace): ///
		mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
		i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
		i.educ i.rural i.vigact i.smoke_status i.srh ///
		i.drink_status i.drink_status#c.time i.higov ///
		i.cogstat i.cogstat#c.time ///
		|| patid: female time_female, noconstant covariance(unstructured) ///
		|| patid: male time_male, noconstant covariance(unstructured) ///
		ml pw(w1respwt) nolog stddev iterate(300) vce(robust)
		
estimates store mi_model

* ------------------------------------------------------------------------------
* Margins - continuous age
* -----------------------------------------------------------------------------
* Requirement: net install mimrgns, from(http://fmwww.bc.edu/RePEc/bocode/m)

mimrgns sex, vce(unconditional) at(time=(0 2 11 14 17))

* Extract and transpose results matrix
matrix margins_sex = r(table)
matrix list margins_sex
matrix margins_sex_t = margins_sex'

* Save margins results
preserve
    clear
    svmat margins_sex_t
    save "$working/margins_MI_sex_results.dta", replace
restore

*--- Linear combinations of parameters

* Baseline female–male gap (time = 0)
lincom 2.sex

* Difference in slopes (female − male, per year)
lincom 2.sex#c.time - 1.sex#c.time

* Female–male gap at each wave (time = 2, 11, 14, 17)
lincom 2.sex + 2*(2.sex#c.time)    // Wave 2 (time = 2)
lincom 2.sex + 11*(2.sex#c.time)   // Wave 3 (time = 11)
lincom 2.sex + 14*(2.sex#c.time)   // Wave 4 (time = 14)
lincom 2.sex + 17*(2.sex#c.time)   // Wave 5 (time = 17)

* Annual slope in females
lincom c.time + 2.sex#c.time

* Slope ratio: female / male
nlcom (slope_ratio: ( _b[c.time] + _b[2.sex#c.time] ) / _b[c.time])
nlcom (pct_faster: 100 * ( _b[2.sex#c.time] / _b[c.time] ))

* Relative female/male gap at t = 17 (ratio of predicted means)
nlcom (gap_ratio_t17: ( _b[2.sex] + 17*_b[2.sex#c.time] + _b[_cons] + 17*_b[c.time] ) / ///
                      ( _b[_cons] + 17*_b[c.time] ))

*------------------------------------------------------------------------------
* Plotting predicted FI by sex - continuous age
* ----------------------------------------------------------------------------
* Predicted FI by sex (lines with capped CIs)
preserve
    clear
    svmat margins_sex_t
    
    * Rename variables to match r(table) rows
    rename margins_sex_t1 b
    rename margins_sex_t2 se
    rename margins_sex_t3 t
    rename margins_sex_t4 pvalue
    rename margins_sex_t5 ll
    rename margins_sex_t6 ul
    rename margins_sex_t7 df
    rename margins_sex_t8 crit
    rename margins_sex_t9 eform
    
    * Time and sex identifiers
    gen obs = _n
    gen time = .
    gen sex = .
    
    forvalues i = 1/10 {
        local t = mod(`i'-1, 5)
        local timevals "0 2 11 14 17"
        local tval : word `=`t'+1' of `timevals'
        replace time = `tval' in `i'
        replace sex = ceil(`i'/5) in `i'
    }
    
twoway (rcap ll ul time if sex==1, lcolor(blue)) ///
           (rcap ll ul time if sex==2, lcolor(red)) ///
           (connected b time if sex==1, mcolor(blue) lcolor(blue) msymbol(O)) ///
           (connected b time if sex==2, mcolor(red) lcolor(red) msymbol(O)), ///
        title("Predicted FI trajectories by sex") ///
        xtitle("Time since baseline (years)") ytitle("Linear prediction, fixed portion") ///
        xlabel(0 2 11 14 17, nogrid) ylabel(, nogrid) ///
        legend(order(3 "1.Male" 4 "2.Female") pos(4) row(2)) ///
        name(margins_sex, replace)
graph export "$outputgraphs/mi_margins_sex_time.jpg", replace
graph export "$outputgraphs/mi_margins_sex_time.tif", replace width(4000)

restore

* ---------------------
* Predicted FI by sex with shaded CIs
* -----------------------

preserve
    use "$working/margins_MI_sex_results.dta", clear
    
    * Rename variables
    rename margins_sex_t1 b
    rename margins_sex_t2 se
    rename margins_sex_t3 t
    rename margins_sex_t4 pvalue
    rename margins_sex_t5 ll
    rename margins_sex_t6 ul
    rename margins_sex_t7 df
    rename margins_sex_t8 crit
    rename margins_sex_t9 eform
    
    * Time and sex identifiers
    gen obs = _n
    gen time = .
    gen sex = .
    
    forvalues i = 1/10 {
        local at_num = ceil(`i'/2)
        local sex_num = mod(`i'-1, 2) + 1
        local timevals "0 2 11 14 17"
        local tval : word `at_num' of `timevals'
        replace time = `tval' in `i'
        replace sex = `sex_num' in `i'
    }
    
    twoway (rarea ll ul time if sex==1, fcolor(blue%20) lwidth(none)) ///
           (rarea ll ul time if sex==2, fcolor(red%20) lwidth(none)) ///
           (line b time if sex==1, lcolor(blue) lwidth(medium)) ///
           (line b time if sex==2, lcolor(red) lwidth(medium)), ///
        title("Predicted FI by Sex Over Time") ///
        xtitle("Time since baseline (years)") ytitle("Predicted Frailty Index") ///
        xlabel(0 2 11 14 17, nogrid) ylabel(, nogrid) ///
        legend(order(3 "Male" 4 "Female") pos(6) row(1)) ///
        name(margins_sex, replace)
    graph export "$outputgraphs/mi_margins_sex_time_CIover.jpg", replace
    graph export "$outputgraphs/mi_margins_sex_time_CIover.tif", replace width(4000)
restore
	
/* Alternative graph (not used)
preserve
    clear
    svmat margins_sex_t
    
    * Rename variables to match row names
    rename margins_sex_t1 b
    rename margins_sex_t2 se
    rename margins_sex_t3 t
    rename margins_sex_t4 pvalue
    rename margins_sex_t5 ll
    rename margins_sex_t6 ul
    rename margins_sex_t7 df
    rename margins_sex_t8 crit
    rename margins_sex_t9 eform
    
    * Generate observation numbers and group/time variables
    gen obs = _n
    gen time = .
    gen sex = .
    
    * Fill in time and sex values based on mimrgns output structure
    * (Adjust indices based on your actual output)
    forvalues i = 1/10 {
        local t = mod(`i'-1, 5)
        local timevals "0 2 11 14 17"
        local tval : word `=`t'+1' of `timevals'
        replace time = `tval' in `i'
        replace sex = ceil(`i'/5) in `i'
    }
    
* Plot with confidence intervals
twoway (rarea ll ul time if sex==1, fcolor(blue%20) lwidth(none)) ///
       (rarea ll ul time if sex==2, fcolor(red%20) lwidth(none)) ///
       (line b time if sex==1, lcolor(blue) lwidth(medium)) ///
       (line b time if sex==2, lcolor(red) lwidth(medium)), ///
    title("Predicted FI by Sex Over Time") ///
    xtitle("Time since baseline (years)") ytitle("Predicted Frailty Index") ///
    xlabel(0 2 11 14 17, nogrid) ylabel(, nogrid) ///
    legend(order(3 "Male" 4 "Female") pos(6) row(1)) ///
    name(margins_sex, replace)
graph export "$outputgraphs/mi_margins_sex_time_CIover.jpg", replace
* Export as high-resolution TIFF
graph export "$outputgraphs/mi_margins_sex_time_CIover.tif", replace width(4000)
restore
*/

*----------------------------------------------------------------------
* Individual FI trajectories (observed data, m = 0)
*-----------------------------------------------------------------------

preserve
* Males
mi extract 0, clear
twoway ///
    (line wfi time if sex==1, lcolor(blue%1) lwidth(vthin) connect(ascending)) ///
    (lowess wfi time if sex==1, lcolor(black) lwidth(midthick) lpattern(solid)), ///
    title("Individual FI Trajectories (Males)") ///
    ytitle("Frailty Index") ///
    xtitle("Years since baseline") ///
    ylabel(0(0.2) 0.8, nogrid) ///
    xlabel(0 2 11 14 17, nogrid) ///
    legend(off) ///
    name(fi_trajectories_male, replace)
graph export "$outputgraphs/fi_trajectories_male.jpg", replace

* Females
twoway ///
    (line wfi time if sex==2, lcolor(red%1) lwidth(vthin) connect(ascending)) ///
    (lowess wfi time if sex==2, lcolor(black) lwidth(midthick) lpattern(solid)), ///
    title("Individual FI Trajectories (Females)") ///
    ytitle("Frailty Index") ///
    xtitle("Years since baseline") ///
    ylabel(0(0.2) 0.8, nogrid) ///
    xlabel(0 2 11 14 17, nogrid) ///
    legend(off) ///
    name(fi_trajectories_female, replace)
graph export "$outputgraphs/fi_trajectories_female.jpg", replace

* Restore MI data
use "$cleanmhasdata/wave1to5_long_MI.dta", clear
restore

* =======================================================================================
* 7 Mixed model on MI data - age categories
* =======================================================================================

capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots saving(miests, replace): ///
		mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
		i.sex#c.time i.age_cat i.age_cat#c.time c.wealth_gbp_c ///
		i.educ i.rural i.vigact i.smoke_status i.srh ///
		i.drink_status i.drink_status#c.time i.higov ///
		i.cogstat i.cogstat#c.time ///
		|| patid: female time_female, noconstant covariance(unstructured) ///
		|| patid: male time_male, noconstant covariance(unstructured) ///
		ml pw(w1respwt) nolog stddev iterate(300) vce(robust)

* -----------------------------------------------------------------
* Margins - age categories
* ---------------------------------------------------------------
* Predicted FI by age band and sex (waves 1 to 5)
mimrgns sex#age_cat, vce(unconditional) at(time=(0 2 11 14 17))

* Extract, transpose and save results matrix
matrix margins_sex_age = r(table)
matrix list margins_sex_age
matrix margins_sex_age_t = margins_sex_age'

preserve
    clear
    svmat margins_sex_age_t
    save "$working/margins_MI_sex_age_cat_results.dta", replace
restore

*--- Linear combinations of parameters

* Sex difference in baseline FI
lincom _b[2.sex]

* Sex difference in linear time slope
lincom _b[2.sex#c.time]

* Age category differences at baseline (80+ vs 50-59)
lincom _b[4.age_cat] - _b[1.age_cat]

* Age category differences in time slope (80+ vs 50-59)
lincom _b[4.age_cat#c.time] - _b[1.age_cat#c.time]

* Other age groups vs 50-59: baseline
lincom _b[2.age_cat] - _b[1.age_cat]  // 60-69 vs 50-59
lincom _b[3.age_cat] - _b[1.age_cat]  // 70-79 vs 50-59

* Other age groups vs 50-59: time slope
lincom _b[2.age_cat#c.time] - _b[1.age_cat#c.time]  // 60-69 vs 50-59
lincom _b[3.age_cat#c.time] - _b[1.age_cat#c.time]  // 70-79 vs 50-59

* --- Adjacent age group differences
* 70-79 vs 60-69 baseline
lincom _b[3.age_cat] - _b[2.age_cat]

* 80+ vs 70-79 baseline
lincom _b[4.age_cat] - _b[3.age_cat]

* 70-79 vs 60-69 slope
lincom _b[3.age_cat#c.time] - _b[2.age_cat#c.time]

* 80+ vs 70-79 slope
lincom _b[4.age_cat#c.time] - _b[3.age_cat#c.time]

* 80+ vs 50–59 difference at Wave 5
lincom (_b[4.age_cat] + 17*_b[4.age_cat#c.time]) - ///
       (_b[1.age_cat] + 17*_b[1.age_cat#c.time])

* Female–male gap at Wave 5 (time = 17)
lincom 2.sex + 17*(2.sex#c.time)

*-------------------------------------------------------------------------------
* Plotting predicted FI by sex and age category
*------------------------------------------------------------------------------
		
* Lines without CIs
preserve
    clear
    svmat margins_sex_age_t
    
    * Rename variables to match r(table) rows
    rename margins_sex_age_t1 b
    rename margins_sex_age_t2 se
    rename margins_sex_age_t3 t
    rename margins_sex_age_t4 pvalue
    rename margins_sex_age_t5 ll
    rename margins_sex_age_t6 ul
    rename margins_sex_age_t7 df
    rename margins_sex_age_t8 crit
    rename margins_sex_age_t9 eform
    
    * Time, sex and age-category identifiers
    gen obs = _n
    gen time = .
    gen sex = .
    gen age_cat = .
    
    * 2 sexes × 4 age categories × 5 time points = 40 rows (8 sex#age_cat combinations per time point)
    forvalues i = 1/40 {
        local time_index = ceil(`i'/8)  // time point (1-5)
        local timevals "0 2 11 14 17"
        local tval : word `time_index' of `timevals'
        replace time = `tval' in `i'
        
        local within_time = mod(`i'-1, 8) + 1  // sex#age combination (1-8)
        replace sex = ceil(`within_time'/4) in `i'  // 1-4 = male, 5-8 = female
        replace age_cat = mod(`within_time'-1, 4) + 1 in `i'
    }
    
    twoway ///
        (line b time if sex==1 & age_cat==1, lcolor(green) lpattern(dash)) ///
        (line b time if sex==2 & age_cat==1, lcolor(green) lpattern(solid)) ///
        (line b time if sex==1 & age_cat==2, lcolor(blue) lpattern(dash)) ///
        (line b time if sex==2 & age_cat==2, lcolor(blue) lpattern(solid)) ///
        (line b time if sex==1 & age_cat==3, lcolor(purple) lpattern(dash)) ///
        (line b time if sex==2 & age_cat==3, lcolor(purple) lpattern(solid)) ///
        (line b time if sex==1 & age_cat==4, lcolor(red) lpattern(dash)) ///
        (line b time if sex==2 & age_cat==4, lcolor(red) lpattern(solid)), ///
        title("Predicted FI Trajectories by Sex and Age Category") ///
        xtitle("Time since baseline (years)") ytitle("Frailty Index") ///
        xlabel(0 2 11 14 17, nogrid) ylabel(, nogrid) ///
        legend(order(1 "Male, 50-59" 2 "Female, 50-59" ///
                     3 "Male, 60-69" 4 "Female, 60-69" ///
                     5 "Male, 70-79" 6 "Female, 70-79" ///
                     7 "Male, 80+" 8 "Female, 80+") ///
               cols(2) pos(6)) ///
        name(age_band_combined, replace)
graph export "$outputgraphs/fi_trajectories_age_bands_by_sex.jpg", replace
graph export "$outputgraphs/fi_trajectories_age_bands_by_sex.tif", replace width(4000)
restore

* Lines with shaded CIs
preserve
    clear
    svmat margins_sex_age_t
    
    * Rename variables to match r(table) rows
    rename margins_sex_age_t1 b
    rename margins_sex_age_t2 se
    rename margins_sex_age_t3 t
    rename margins_sex_age_t4 pvalue
    rename margins_sex_age_t5 ll
    rename margins_sex_age_t6 ul
    rename margins_sex_age_t7 df
    rename margins_sex_age_t8 crit
    rename margins_sex_age_t9 eform
    
    * Time, sex and age-category identifiers
    gen obs = _n
    gen time = .
    gen sex = .
    gen age_cat = .
    
    forvalues i = 1/40 {
        local time_index = ceil(`i'/8)  // time point (1-5)
        local timevals "0 2 11 14 17"
        local tval : word `time_index' of `timevals'
        replace time = `tval' in `i'
        
        local within_time = mod(`i'-1, 8) + 1  // sex#age combination (1-8)
        replace sex = ceil(`within_time'/4) in `i'  // 1-4 = male, 5-8 = female
        replace age_cat = mod(`within_time'-1, 4) + 1 in `i'
    }
    
    twoway ///
        (rarea ll ul time if sex==1 & age_cat==1, fcolor(green%20) lwidth(none)) ///
        (rarea ll ul time if sex==2 & age_cat==1, fcolor(green%20) lwidth(none)) ///
        (rarea ll ul time if sex==1 & age_cat==2, fcolor(blue%20) lwidth(none)) ///
        (rarea ll ul time if sex==2 & age_cat==2, fcolor(blue%20) lwidth(none)) ///
        (rarea ll ul time if sex==1 & age_cat==3, fcolor(purple%20) lwidth(none)) ///
        (rarea ll ul time if sex==2 & age_cat==3, fcolor(purple%20) lwidth(none)) ///
        (rarea ll ul time if sex==1 & age_cat==4, fcolor(red%20) lwidth(none)) ///
        (rarea ll ul time if sex==2 & age_cat==4, fcolor(red%20) lwidth(none)) ///
        (line b time if sex==1 & age_cat==1, lcolor(green) lpattern(dash)) ///
        (line b time if sex==2 & age_cat==1, lcolor(green) lpattern(solid)) ///
        (line b time if sex==1 & age_cat==2, lcolor(blue) lpattern(dash)) ///
        (line b time if sex==2 & age_cat==2, lcolor(blue) lpattern(solid)) ///
        (line b time if sex==1 & age_cat==3, lcolor(purple) lpattern(dash)) ///
        (line b time if sex==2 & age_cat==3, lcolor(purple) lpattern(solid)) ///
        (line b time if sex==1 & age_cat==4, lcolor(red) lpattern(dash)) ///
        (line b time if sex==2 & age_cat==4, lcolor(red) lpattern(solid)), ///
        title("Predicted FI Trajectories by Sex and Age Category") ///
        xtitle("Time since baseline (years)") ytitle("Frailty Index") ///
        xlabel(0 2 11 14 17, nogrid) ylabel(, nogrid) ///
        legend(order(9 "Male, 50-59" 10 "Female, 50-59" ///
                     11 "Male, 60-69" 12 "Female, 60-69" ///
                     13 "Male, 70-79" 14 "Female, 70-79" ///
                     15 "Male, 80+" 16 "Female, 80+") ///
               cols(2) pos(6)) ///
        name(age_band_combined, replace)
    graph export "$outputgraphs/fi_trajectories_age_bands_by_sex_CIover.jpg", replace
    graph export "$outputgraphs/fi_trajectories_age_bands_by_sex_CIover.tif", replace width(4000)
restore
		
*------------------------------------------------------------------------
* Predicted FI by age category, panels by sex
*------------------------------------------------------------------
preserve
    clear
    svmat margins_sex_age_t
    
    * Rename variables to match r(table) rows
    rename margins_sex_age_t1 b
    rename margins_sex_age_t2 se
    rename margins_sex_age_t3 t
    rename margins_sex_age_t4 pvalue
    rename margins_sex_age_t5 ll
    rename margins_sex_age_t6 ul
    rename margins_sex_age_t7 df
    rename margins_sex_age_t8 crit
    rename margins_sex_age_t9 eform
    
    * Time, sex and age-category identifiers
    gen obs = _n
    gen time = .
    gen sex = .
    gen age_cat = .
    
    forvalues i = 1/40 {
        local time_index = ceil(`i'/8)
        local timevals "0 2 11 14 17"
        local tval : word `time_index' of `timevals'
        replace time = `tval' in `i'
        
        local within_time = mod(`i'-1, 8) + 1
        replace sex = ceil(`within_time'/4) in `i'
        replace age_cat = mod(`within_time'-1, 4) + 1 in `i'
    }
    
    * Sex labels for panel titles
    label define sexlbl 1 "Male" 2 "Female"
    label values sex sexlbl
    
    twoway ///
        (rarea ll ul time if age_cat==1, fcolor(green%20) lwidth(none)) ///
        (rarea ll ul time if age_cat==2, fcolor(blue%20) lwidth(none)) ///
        (rarea ll ul time if age_cat==3, fcolor(purple%20) lwidth(none)) ///
        (rarea ll ul time if age_cat==4, fcolor(red%20) lwidth(none)) ///
        (line b time if age_cat==1, lcolor(green) lwidth(medium)) ///
        (line b time if age_cat==2, lcolor(blue) lwidth(medium)) ///
        (line b time if age_cat==3, lcolor(purple) lwidth(medium)) ///
        (line b time if age_cat==4, lcolor(red) lwidth(medium)), ///
        by(sex, title("Predicted FI Trajectories by Age Category") note("") legend(pos(6))) ///
        xtitle("Time since baseline (years)") ///
        ytitle("Frailty Index") ///
        xlabel(0 2 11 14 17, nogrid) ylabel(, nogrid) ///
        legend(order(5 "Age 50-59" 6 "Age 60-69" 7 "Age 70-79" 8 "Age 80+") ///
               cols(4) size(small) symxsize(5) region(lcolor(none))) ///
        name(age_band_by_sex, replace)
    graph export "$outputgraphs/fi_trajectories_agecat_by_sex_CIover.jpg", replace
    graph export "$outputgraphs/fi_trajectories_agecat_by_sex_CIover.tif", replace width(4000)
restore

*===============================================================================
* 7. Interaction terms (MI): sex × education, SRH, wealth
*===============================================================================
* Each tested in the linear-time model and again allowing for non-linear time

*1) Sex # Education (#time)

capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots: ///
	mixed wfi i.sex c.time i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.educ#c.time  i.sex#i.educ  i.sex#i.educ#c.time   ///
	 i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml pw(w1respwt) nolog stddev iterate(300) vce(robust)
	 
* Differential effect of education on FI trajectory in men vs women
test 2.sex#2.educ#c.time 2.sex#3.educ#c.time	

* Differential effect of education on FI in men vs women (baseline or follow-up)
test 2.sex#2.educ#c.time 2.sex#3.educ#c.time  2.sex#2.educ 2.sex#3.educ

* Allowing for non-linear time
capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots: ///
	mixed wfi i.sex c.time i.sex#c.time c.time#c.time c.time#c.time#c.time ///
     c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.educ#c.time  i.sex#i.educ  i.sex#i.educ#c.time   ///
	 i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
	 
test 2.sex#2.educ#c.time 2.sex#3.educ#c.time
test 2.sex#2.educ#c.time 2.sex#3.educ#c.time  2.sex#2.educ 2.sex#3.educ

*2) Sex # SRH (#time)

capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots: ///
	mixed wfi i.sex c.time i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status ///
	 i.srh i.srh#c.time i.srh#i.sex i.srh#i.sex#c.time  ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
	 
* Differential effect of SRH on FI trajectory in men vs women
testparm i.srh#i.sex#c.time 

* Differential effect of SRH on FI in men vs women (baseline or follow-up)
testparm i.srh#i.sex#c.time i.srh#i.sex

* Allowing for non-linear time
capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots: ///
	mixed wfi i.sex c.time i.sex#c.time c.time#c.time c.time#c.time#c.time ///
     c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status ///
	 i.srh i.srh#c.time i.srh#i.sex i.srh#i.sex#c.time  ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
	 
testparm i.srh#i.sex#c.time 
testparm i.srh#i.sex#c.time i.srh#i.sex

*3) Sex # Wealth (#time)

capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots: ///
	mixed wfi i.sex c.time i.sex#c.time c.age_c c.age_c#c.time ///
	c.wealth_gbp_c c.wealth_gbp_c#c.time ///
	c.wealth_gbp_c#i.sex c.wealth_gbp_c#i.sex#c.time ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* Differential effect of wealth on FI in men vs women
test 2.sex#c.wealth_gbp_c 

* Differential effect of wealth on FI in men vs women (baseline or follow-up)
test 2.sex#c.wealth_gbp_c#c.time 2.sex#c.wealth_gbp_c 
	 
* Allowing for non-linear time
capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots: ///
	mixed wfi i.sex c.time i.sex#c.time c.time#c.time c.time#c.time#c.time ///
	c.age_c c.age_c#c.time ///
	c.wealth_gbp_c c.wealth_gbp_c#c.time ///
	c.wealth_gbp_c#i.sex c.wealth_gbp_c#i.sex#c.time ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

test 2.sex#c.wealth_gbp_c 
test 2.sex#c.wealth_gbp_c#c.time 2.sex#c.wealth_gbp_c

*===============================================================================
* 7. Session end
*===============================================================================
display as text "------------------------------------------------------------"
display as text "MHAS_MI_completed"
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close