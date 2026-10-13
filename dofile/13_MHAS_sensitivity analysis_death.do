********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Sensitivity analyses for attrition and terminal decline among decedents
********************************************************************************

* ==============================================================================
* 3. Load Data
* ==============================================================================

use "$cleanmhasdata/wave1to5_long.dta", clear

* Set seed
set seed 14122025

xtset patid time

* ===========================================================================
* 4. Sensitivity Analysis - attrition
* ===========================================================================  

/* Participant status codes:
Status	Meaning
1	Alive and observed at wave 5
2	Known death
3	LFU–dead (special death codes, no wave 5 FI)
4	LFU–unknown (dropout before wave 5, not known dead)
*/

preserve 

*--------------------------------------------------------------------------------
* Complete-survivor model: excludes deaths and dropouts
*--------------------------------------------------------------------------------

mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
		i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
		i.educ i.rural i.vigact i.smoke_status i.srh ///
		i.drink_status i.drink_status#c.time i.higov ///
		i.cogstat i.cogstat#c.time if status == 1 ///
		|| patid: female time_female, noconstant covariance(unstructured) ///
		|| patid: male time_male, noconstant covariance(unstructured) ///
		ml pw(w1respwt) nolog stddev base
		
estat ic
esti store m_lfu_aliveonly

* -------------------------------------------------------------------------------
* Death-exclusion model (LFU-death): removes mortality-related attrition, keeps non-death dropout
* ------------------------------------------------------------------------------
mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
		i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
		i.educ i.rural i.vigact i.smoke_status i.srh ///
		i.drink_status i.drink_status#c.time i.higov ///
		i.cogstat i.cogstat#c.time if !(status == 2 | status == 3) ///
		|| patid: female time_female, noconstant covariance(unstructured) ///
		|| patid: male time_male, noconstant covariance(unstructured) ///
		ml pw(w1respwt) nolog stddev base

estat ic
esti store m_lfu_death

*-------------------------------------------------------------------------------
* Dropout-exclusion model (LFU-dropout): removes non-response attrition, keeps deaths
*-------------------------------------------------------------------------------

mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
		i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
		i.educ i.rural i.vigact i.smoke_status i.srh ///
		i.drink_status i.drink_status#c.time i.higov ///
		i.cogstat i.cogstat#c.time if !(status == 4) /// 
		|| patid: female time_female, noconstant covariance(unstructured) ///
		|| patid: male time_male, noconstant covariance(unstructured) ///
		ml nolog stddev pweight(w1respwt) base 
estat ic
esti store m_lfu_dropout

* ===========================================================================================
* 8 FI trajectories among participants alive at end of study vs participants who died
* ============================================================================================
* Note: descriptive comparison by vital status (not a sensitivity analysis)

*--------------------------------------------------------------------------------------------
* Model among participants alive at end of study (excluding dropouts)
*--------------------------------------------------------------------------------------------

mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
		i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
		i.educ i.rural i.vigact i.smoke_status i.srh ///
		i.drink_status i.drink_status#c.time i.higov ///
		i.cogstat i.cogstat#c.time if status == 1 ///
		|| patid: female time_female, noconstant covariance(unstructured) ///
		|| patid: male time_male, noconstant covariance(unstructured) ///
		ml pw(w1respwt) nolog stddev

* Linear combinations
lincom 2.sex
lincom 2.sex#c.time - 1.sex#c.time
lincom 2.sex + 2*(2.sex#c.time)    
lincom 2.sex + 11*(2.sex#c.time)   
lincom 2.sex + 14*(2.sex#c.time)   
lincom 2.sex + 17*(2.sex#c.time)   
lincom c.time + 2.sex#c.time

* Slope ratio: female / male
nlcom (slope_ratio: ( _b[c.time] + _b[2.sex#c.time] ) / _b[c.time])
nlcom (pct_faster: 100 * ( _b[2.sex#c.time] / _b[c.time] ))

* Relative female/male gap at t = 17 (ratio of predicted means)
nlcom (gap_ratio_t17: ( _b[2.sex] + 17*_b[2.sex#c.time] + _b[_cons] + 17*_b[c.time] ) / ///
                      ( _b[_cons] + 17*_b[c.time] ))

* -- Predicted FI by sex at each time point (alive)
margins sex, vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_sex_alive_wave1to5.dta", replace)

marginsplot, title("Predicted FI trajectories by sex among alive") ///
    ylabel(, nogrid) xlabel(,nogrid) ///
    name(margins_sex_alive, replace)
graph export "$outputgraphs/margins_sex_alive.jpg", replace

restore 

*-------------------------------------------------------------------------------------------------
* Model among participants who died (calendar time scale)
*-------------------------------------------------------------------------------------------------

/* Did not converge:

mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
		i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
		i.educ i.rural i.vigact i.smoke_status i.srh ///
		i.drink_status i.drink_status#c.time i.higov ///
		i.cogstat i.cogstat#c.time if death == 1 ///
		|| patid: female time_female, noconstant covariance(unstructured) ///
		|| patid: male time_male, noconstant covariance(unstructured) ///
		ml pw(w1respwt) nolog stddev 
*/

*-------------------------------------------------------------------------------
* Terminal decline among decedents (time-to-death scale)
* ------------------------------------------------------------------------------
* Descriptive model among participants who died during follow-up, using years until death
* (backward time). Not a sensitivity analysis; characterises terminal decline.

* Interview year = baseline year (2001) + time since baseline
local baseline_year 2001
gen intyear = `baseline_year' + time

local interview_year_var  intyear
local death_year_var      dyear
local death_ind_var       death

capture drop time_to_death
gen time_to_death = `death_year_var' - `interview_year_var' if `death_ind_var' == 1
label var time_to_death "Years until death (0=death year; higher=further from death)"

* Drop impossible values (interview after recorded death year)
drop if `death_ind_var' == 1 & time_to_death < 0

summ time_to_death if death==1, detail

* Distribution of time to death
tabulate time_to_death if death == 1, missing
bysort sex: tabulate time_to_death if death == 1, missing

histogram time_to_death if death == 1, ///
    width(1) percent ///
    xlabel(0(1)20) ///
    title("Distribution of time to death (years)") ///
    xtitle("Years until death") ytitle("Percent") ///
    name(ttd_hist, replace)

*-----------------------------------------------------------------------------------------------
* Random-effects design variables on the time-to-death scale
*-----------------------------------------------------------------------------------------------
* sex: 1 = Male, 2 = Female
capture drop female_ttd male_ttd time_female_ttd time_male_ttd
gen female_ttd = (sex == 2)
gen male_ttd   = (sex == 1)

gen time_female_ttd = female_ttd * time_to_death
gen time_male_ttd   = male_ttd   * time_to_death

label var female_ttd      "Female indicator"
label var male_ttd        "Male indicator"
label var time_female_ttd "Female x time_to_death"
label var time_male_ttd   "Male x time_to_death"

count if death==1

*-----------------------------------------------------------------------------------------------
* Mixed model among decedents (time-to-death as the time scale)
*-----------------------------------------------------------------------------------------------

mixed wfi i.sex ///
	 c.time_to_death c.time_to_death#c.time_to_death c.time_to_death#c.time_to_death#c.time_to_death ///
     i.sex#c.time_to_death c.age_c c.age_c#c.time_to_death ///
     c.wealth_gbp_c i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status i.drink_status#c.time_to_death i.higov ///
     i.cogstat i.cogstat#c.time_to_death if death == 1 ///
     || patid: female_ttd time_female_ttd, noconstant covariance(unstructured) ///
     || patid: male_ttd   time_male_ttd,   noconstant covariance(unstructured) ///
     ml pw(w1respwt) nolog stddev

estimates store m_dead_ttd

* Test sex difference in slope on the time-to-death scale
test 2.sex#c.time_to_death

*-----------------------------------------------------------------------------------------------
* Post-estimation: sex contrasts and predicted trajectories (decedents)
*-----------------------------------------------------------------------------------------------

* Sex difference in the year of death (time_to_death = 0)
lincom 2.sex

* Female vs male slope difference on the time-to-death scale
lincom 2.sex#c.time_to_death - 1.sex#c.time_to_death

* Female–male difference at selected years before death
lincom 2.sex + 2*(2.sex#c.time_to_death)
lincom 2.sex + 5*(2.sex#c.time_to_death)
lincom 2.sex + 10*(2.sex#c.time_to_death)
lincom 2.sex + 11*(2.sex#c.time_to_death)
lincom 2.sex + 15*(2.sex#c.time_to_death)
lincom 2.sex + 17*(2.sex#c.time_to_death)

* Predicted FI by sex at selected years before death (0, 2, 5, 10, 15)
margins sex, vce(unconditional) at(time_to_death=(0 2 5 10 15)) ///
    saving("$working/margins_sex_dead_time_to_death.dta", replace)

marginsplot, title("Predicted FI by sex among decedents (time-to-death scale)") ///
    ytitle("Frailty Index") xtitle("Years until death") ///
    ylabel(, nogrid) xlabel(, nogrid) ///
    name(margins_sex_dead_ttd, replace)
	
graph export "$outputgraphs/margins_sex_dead_time_to_death.jpg", replace
graph export "$outputgraphs/margins_sex_dead_time_to_death.tif", replace width(4000)
	
* Predicted FI by sex at selected years before death (0, 2, 11, 17)
margins sex, vce(unconditional) at(time_to_death=(0 2 11 17)) ///
    saving("$working/margins_sex_dead_time_to_death_bywave0-17.dta", replace)

marginsplot, title("Predicted FI by sex among decedents (time-to-death scale WAVE)") ///
    ytitle("Frailty Index") xtitle("Years until death") ///
    ylabel(, nogrid) xlabel(, nogrid) ///
    name(margins_sex_dead_ttd, replace)

graph export "$outputgraphs/margins_sex_dead_time_to_death_wave)0-17.jpg", replace
graph export "$outputgraphs/margins_sex_dead_time_to_death_wave)0-17.tif", replace width (4000)

* Predicted FI by sex at selected years before death (0, 2, 5, 11, 15, 17)
margins sex, vce(unconditional) at(time_to_death=(0 2 5 11 15 17)) ///
    saving("$working/margins_sex_dead_time_to_death_bywave0-17.dta", replace)

marginsplot, title("Predicted FI by sex among decedents (time-to-death scale WAVE)") ///
    ytitle("Frailty Index") xtitle("Years until death") ///
    ylabel(, nogrid) xlabel(, nogrid) ///
    name(margins_sex_dead_ttd, replace)

graph export "$outputgraphs/margins_sex_dead_time_to_death_wave+15)0-15&17.jpg", replace
graph export "$outputgraphs/margins_sex_dead_time_to_death_wave+15)0-15&17.tif", replace width(4000)

preserve

*-----------------------------------------------------------------------------------------------
* Predicted FI by sex with shaded CIs (years until death)
*-----------------------------------------------------------------------------------------------

* Load margins results
use "$working/margins_sex_dead_time_to_death_bywave0-17.dta", clear

rename _margin pred_fi
rename _ci_lb ci_lb
rename _ci_ub ci_ub
rename _at2 time_to_death    // time_to_death stored in _at2
rename _m1 sex               // sex stored in _m1

gen sex_label = "Male" if sex == 1
replace sex_label = "Female" if sex == 2

twoway ///
    (rarea ci_lb ci_ub time_to_death if sex == 1, color(blue%20) lwidth(none)) ///
    (rarea ci_lb ci_ub time_to_death if sex == 2, color(red%20) lwidth(none)) ///
    (line pred_fi time_to_death if sex == 1, lcolor(blue) lwidth(medium)) ///
    (line pred_fi time_to_death if sex == 2, lcolor(red) lwidth(medium)), ///
    title("Predicted FI by Sex Over Time") ///
    ytitle("Predicted Frailty Index") ///
    xtitle("Time until death (years)") ///
    ylabel(, nogrid) ///
    xlabel(0 2 5 11 15 17, nogrid) ///
    legend(order(3 "Male" 4 "Female") position(6) ring(1) cols(2) ) ///
    name(fi_sex_ttd_custom, replace)

graph export "$outputgraphs/margins_sex_dead_time_to_death_custom.jpg", replace
graph export "$outputgraphs/margins_sex_dead_time_to_death_custom.tif", replace width(4000)

restore

*-----------------------------------------------------------------------------------------------
* Predicted FI by sex with shaded CIs (reversed x-axis: 17 to 0 years before death)
*-----------------------------------------------------------------------------------------------

* Load margins results
use "$working/margins_sex_dead_time_to_death_bywave0-17.dta", clear

rename _margin pred_fi
rename _ci_lb ci_lb
rename _ci_ub ci_ub
rename _at2 time_to_death    // time_to_death stored in _at2
rename _m1 sex               // sex stored in _m1

gen sex_label = "Male" if sex == 1
replace sex_label = "Female" if sex == 2

twoway ///
    (rarea ci_lb ci_ub time_to_death if sex == 1, color(blue%20) lwidth(none)) ///
    (rarea ci_lb ci_ub time_to_death if sex == 2, color(red%20) lwidth(none)) ///
    (line pred_fi time_to_death if sex == 1, lcolor(blue) lwidth(medium)) ///
    (line pred_fi time_to_death if sex == 2, lcolor(red) lwidth(medium)), ///
    title("Predicted FI by Sex Over Time to Death") ///
    ytitle("Predicted Frailty Index") ///
    xtitle("Years before death") ///
    ylabel(, nogrid) ///
    xlabel(17 15 11 5 2 0, nogrid) ///
    xscale(reverse) ///
    legend(order(3 "Male" 4 "Female") position(6) ring(1) cols(2) ) ///
    name(fi_sex_ttd_custom, replace)

graph export "$outputgraphs/margins_sex_dead_time_to_death_17to0.jpg", replace
graph export "$outputgraphs/margins_sex_dead_time_to_death_17to0.tif", replace width(4000)

restore

* ==============================================================================
* 10 Finish
* ==============================================================================

display as text "------------------------------------------------------------"
display "==== MHAS_LMMs_Post-estimation_COMPLETED  ===="
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close