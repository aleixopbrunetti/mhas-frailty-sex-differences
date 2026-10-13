********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Sensitivity Analysis for MHAS Frailty Index Study
********************************************************************************

* ==============================================================================
* 3. Load Data
* ==============================================================================

use "$cleanmhasdata/wave1to5_long.dta", clear

xtset patid time

* ===========================================================================
* 4. Perform Sensitivity Analysis
* ===========================================================================  

* ---------------------------------------------------------------------------
* Source of model missing data (participants excluded from the final model)
* ---------------------------------------------------------------------------

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 

* -- Flag rows used in the model
gen byte inmodel = e(sample)

preserve 
bysort patid: egen byte kept_inmodel = max(inmodel)
count if kept_inmodel==0

* Keep only excluded participants
keep if kept_inmodel==0
keep patid
duplicates drop

* Variables to check
local X "wfi time age_c wealth_gbp_c educ rural vigact smoke_status srh drink_status higov cogstat sex"

* Return to full dataset
restore

* Row-level missing indicator
egen row_miss = rowmiss(wfi time age_c wealth_gbp_c educ rural vigact smoke_status srh drink_status higov cogstat sex)
gen okrow = (row_miss==0 & w1respwt>0 & !missing(w1respwt))

* Participant-level analysis
bysort patid: egen any_ok = max(okrow)
egen tagpat = tag(patid)

* Participants excluded because a variable is missing at every wave
foreach v of local X {
    bysort patid: egen allmiss_`v' = min(missing(`v'))
    quietly count if tagpat & any_ok==0 & allmiss_`v'==1
    display "Variable `v': " r(N) " patients excluded (all observations missing)"
}

* Missing patterns among excluded participants
list patid allmiss_* if tagpat & any_ok==0 in 1/1000, noobs clean

/* Not included: pwscale sensitivity analysis not appropriate (level-1 weights w_ij not available)
* ====================================================================
* 5. Model with different pwscales
* ====================================================================

// pwscale() controls how sampling weights are scaled in two-level models: size, effective, or gk
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) pwscale(size) base
esti store pwsize
estat ic

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_sensitivity_analysis.xlsx", ///
     sheet("Model pwsize") firstrow(variables) sheetreplace
restore

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) pwscale(effective) base
esti store pweffective
estat ic

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_sensitivity_analysis.xlsx", ///
     sheet("Model pweffective") firstrow(variables) sheetreplace
restore

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) pwscale(gk) base
esti store pwgk
estat ic

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_sensitivity_analysis.xlsx", ///
     sheet("Model pwsize") firstrow(variables) sheetreplace
restore

// Final model with pwscale variants
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) pwscale(size) base
estat ic
esti store m_size

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) pwscale(effective) base
estat ic
esti store m_effective

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) pwscale(gk) base
estat ic
esti store m_gk

estimates stats pwsize pweffective pwgk m_size m_effective m_gk

*/

* =======================================================================
* 6. Model with BMI participants (anthropometric sub-sample)
* =======================================================================

tab bmi, missing
sum bmi
gen bmi_c = bmi - r(mean)

* Saturated model with BMI category
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time ///
     i.bmicat##c.time if !(bmi ==.s) ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
testparm i.bmicat##c.time

* Final model with BMI category
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time ///
     i.bmicat##c.time if !(bmi ==.s) ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
estat ic
esti store m_bmicat

*-- Level-1 residuals (m_bmicat)
predict m_bmicat, rstandard            // standardized residuals
predict yhat_m_bmicat                   // fitted values

hist m_bmicat, normal title("Level-1 standardized residuals (m_bmicat)")
graph export "$outputgraphs/m_bmicat_histogram.jpg", replace
qnorm m_bmicat, title("Q–Q plot: Level-1 residuals (m_bmicat)")
graph export "$outputgraphs/m_bmicat_qnorm.jpg", replace

twoway (scatter m_bmicat yhat_m_bmicat, msymbol(o) msize(vsmall)) ///
       (lowess  m_bmicat yhat_m_bmicat), ///
       yline(0) ///
       title("Residuals vs fitted (m_bmicat)") ///
       ytitle("Std. residuals") xtitle("Fitted values")

*-- Level-2 residuals (BLUPs, m_bmicat)
predict fem_m_bmicat_intercept fem_m_bmicat_slope mal_m_bmicat_intercept mal_m_bmicat_slope, reffects

preserve
gen pickone if e(sample)=1
contract patid fem_m_bmicat_intercept fem_m_bmicat_slope mal_m_bmicat_intercept mal_m_bmicat_slope female male

* Female
qnorm fem_m_bmicat_intercept if female, title("Q–Q: random intercept (female) (m_bmicat)")
graph export "$outputgraphs/BLUP_intercept_m_bmicat_qnorm_female.jpg", replace
qnorm fem_m_bmicat_slope if female, title("Q–Q: random slope (time_female) (m_bmicat)")
graph export "$outputgraphs/BLUP_slope_m_bmicat_qnorm_female.jpg", replace

* Male
qnorm mal_m_bmicat_intercept if male, title("Q–Q: random intercept (male) (m_bmicat)")
graph export "$outputgraphs/BLUP_intercept_m_bmicat_qnorm_male.jpg", replace
qnorm mal_m_bmicat_slope if male, title("Q–Q: random slope (time_male) (m_bmicat)")
graph export "$outputgraphs/BLUP_slope_m_bmicat_qnorm_male.jpg", replace

restore 

* ===================================================================
* 7. Modeling with different missing data
* ===================================================================

*-------------------------------------------------------------------------------------------------------------
* Excluding LFU - dropouts (status 4)
*-------------------------------------------------------------------------------------------------------------

/*mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time if !(status == 4) ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

did not converge
*/
/*
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time if !(status == 4) ///
     || patid: time, covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) residuals(independent, by(sex)) base 
*/

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time if !(status == 4) /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 
estat ic
esti store m_lfu_dropout

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_sensitivity_analysis.xlsx", ///
     sheet("Model LFU_dropouts") firstrow(variables) replace
restore

* Model check 
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time if !(status == 4) /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base  
estat ic

*-- Level-1 residuals (m_lfu_dropout)
predict m_lfu_dropout, rstandard      // standardized residuals
predict yhat_m_lfu_dropout              // fitted values
predict m_lfu_dropout_re_female ///
        m_lfu_dropout_re_time_female ///
        m_lfu_dropout_re_male ///
        m_lfu_dropout_re_time_male, reffects

hist m_lfu_dropout, normal  title("Level-1 standardized residuals (m_lfu_dropout)")
qnorm m_lfu_dropout, title("Q–Q plot: Level-1 residuals (m_lfu_dropout)")

twoway ///
    (scatter m_lfu_dropout yhat_m_lfu_dropout, msymbol(o) msize(vsmall)) ///
    (lowess  m_lfu_dropout yhat_m_lfu_dropout), ///
    yline(0) ///
    title("Residuals vs fitted (m_lfu_dropout)") ///
    ytitle("Std. residuals") xtitle("Fitted values")

* Residuals by sex
twoway (scatter m_lfu_dropout yhat_m_lfu_dropout if sex==1, msymbol(o) msize(vsmall)) ///
       (lowess  m_lfu_dropout yhat_m_lfu_dropout if sex==1), ///
       yline(0) title("Residuals vs time (males, m_lfu_dropout)")

twoway (scatter m_lfu_dropout yhat_m_lfu_dropout if sex==2, msymbol(o) msize(vsmall)) ///
       (lowess  m_lfu_dropout yhat_m_lfu_dropout if sex==2), ///
       yline(0) title("Residuals vs time (females, m_lfu_dropout)")

*-- Level-2 residuals (BLUPs, m_lfu_dropout)
hist  m_lfu_dropout_re_female, normal title("Random intercepts (female, BLUPs)")
qnorm m_lfu_dropout_re_female, title("Q–Q plot: random intercepts (female, BLUPs)")
hist  m_lfu_dropout_re_time_female, normal title("Random slopes (time_female, BLUPs)")
qnorm m_lfu_dropout_re_time_female, title("Q–Q plot: random slopes (time_female, BLUPs)")
hist  m_lfu_dropout_re_male, normal title("Random intercepts (male, BLUPs)")
qnorm m_lfu_dropout_re_male, title("Q–Q plot: random intercepts (male, BLUPs)")
hist  m_lfu_dropout_re_time_male, normal title("Random slopes (time_male, BLUPs)")
qnorm m_lfu_dropout_re_time_male, title("Q–Q plot: random slopes (time_male, BLUPs)")

* ---------------------------------------------------------------------------------------------
* Excluding deaths and LFU - death (status 2 and 3)
* ---------------------------------------------------------------------------------------------

*/mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time if !(status == 2 | status == 3) ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
*/

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time if !(status == 2 | status == 3) /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 

estat ic
esti store m_lfu_death

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_sensitivity_analysis.xlsx", ///
     sheet("Model LFU_death") firstrow(variables) sheetreplace
restore

* Model check 
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time if !(status == 2 | status == 3) /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 
estat ic

*-- Level-1 residuals (m_lfu_death)
predict m_lfu_death, rstandard            // standardized residuals
predict yhat_m_lfu_death                   // fitted values

hist m_lfu_death, normal title("Level-1 standardized residuals (m_lfu_death)")
graph export "$outputgraphs/m_lfu_death_histogram.jpg", replace
qnorm m_lfu_death, title("Q–Q plot: Level-1 residuals (m_lfu_death)")
graph export "$outputgraphs/m_lfu_death_qnorm.jpg", replace

twoway (scatter m_lfu_death yhat_m_lfu_death, msymbol(o) msize(vsmall)) ///
       (lowess  m_lfu_death yhat_m_lfu_death), ///
       yline(0) ///
       title("Residuals vs fitted (m_lfu_death)") ///
       ytitle("Std. residuals") xtitle("Fitted values")

*-- Level-2 residuals (BLUPs, m_lfu_death)
predict fem_m_lfu_death_intercept fem_m_lfu_death_slope mal_m_lfu_death_intercept mal_m_lfu_death_slope, reffects

preserve
gen pickone if e(sample)=1
contract patid fem_m_lfu_death_intercept fem_m_lfu_death_slope mal_m_lfu_death_intercept mal_m_lfu_death_slope female male

* Female
qnorm fem_m_lfu_death_intercept if female, title("Q–Q: random intercept (female) (m_lfu_death)")
graph export "$outputgraphs/BLUP_intercept_m_lfu_death_qnorm_female.jpg", replace
qnorm fem_m_lfu_death_slope if female, title("Q–Q: random slope (time_female) (m_lfu_death)")
graph export "$outputgraphs/BLUP_slope_m_lfu_death_qnorm_female.jpg", replace

* Male
qnorm mal_m_lfu_death_intercept if male, title("Q–Q: random intercept (male) (m_lfu_death)")
graph export "$outputgraphs/BLUP_intercept_m_lfu_death_qnorm_male.jpg", replace
qnorm mal_m_lfu_death_slope if male, title("Q–Q: random slope (time_male) (m_lfu_death)")
graph export "$outputgraphs/BLUP_slope_m_lfu_death_qnorm_male.jpg", replace

restore

*-------------------------------------------------------------------------------------------------------------
* Alive at end of study only (status 1, including intermittent missing)
*-------------------------------------------------------------------------------------------------------------

*/mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time if status == 1 ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
*/

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time if status == 1 /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 
estat ic
esti store m_lfu_aliveonly

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_sensitivity_analysis.xlsx", ///
     sheet("Model Alive_Intermittent") firstrow(variables) sheetreplace
restore

* Model check 
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time if status == 1 /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 
estat ic

*-- Level-1 residuals (m_lfu_aliveonly)
predict m_lfu_aliveonly, rstandard            // standardized residuals
predict yhat_m_lfu_aliveonly                   // fitted values

hist m_lfu_aliveonly, normal title("Level-1 standardized residuals (m_lfu_aliveonly)")
graph export "$outputgraphs/m_lfu_aliveonly_histogram.jpg", replace
qnorm m_lfu_aliveonly, title("Q–Q plot: Level-1 residuals (m_lfu_aliveonly)")
graph export "$outputgraphs/m_lfu_aliveonly_qnorm.jpg", replace

twoway (scatter m_lfu_aliveonly yhat_m_lfu_aliveonly, msymbol(o) msize(vsmall)) ///
       (lowess  m_lfu_aliveonly yhat_m_lfu_aliveonly), ///
       yline(0) ///
       title("Residuals vs fitted (m_lfu_aliveonly)") ///
       ytitle("Std. residuals") xtitle("Fitted values")

*-- Level-2 residuals (BLUPs, m_lfu_aliveonly)
predict fem_m_lfu_aliveonly_intercept fem_m_lfu_aliveonly_slope mal_m_lfu_aliveonly_intercept mal_m_lfu_aliveonly_slope, reffects

preserve
gen pickone if e(sample)=1
contract patid fem_m_lfu_aliveonly_intercept fem_m_lfu_aliveonly_slope mal_m_lfu_aliveonly_intercept mal_m_lfu_aliveonly_slope female male

* Female
qnorm fem_m_lfu_aliveonly_intercept if female, title("Q–Q: random intercept (female) (m_lfu_aliveonly)") 
graph export "$outputgraphs/BLUP_intercept_m_lfu_aliveonly_qnorm_female.jpg", replace
qnorm fem_m_lfu_aliveonly_slope if female, title("Q–Q: random slope (time_female) (m_lfu_aliveonly)")
graph export "$outputgraphs/BLUP_slope_m_lfu_aliveonly_qnorm_female.jpg", replace

* Male
qnorm mal_m_lfu_aliveonly_intercept if male, title("Q–Q: random intercept (male) (m_lfu_aliveonly)")
graph export "$outputgraphs/BLUP_intercept_m_lfu_aliveonly_qnorm_male.jpg", replace
qnorm mal_m_lfu_aliveonly_slope if male, title("Q–Q: random slope (time_male) (m_lfu_aliveonly)")
graph export "$outputgraphs/BLUP_slope_m_lfu_aliveonly_qnorm_male.jpg", replace

restore 

* ==================================================================
* 9 Checking for gamma distribution
* ==================================================================

* Requirement: ssc install gammafit

preserve 
use "$cleanmhasdata/wave1to5_long.dta", clear

* Fit Gamma distribution to FI (excluding zeros)
gammafit wfi if wfi>0

* Gamma parameters from gammafit
local shape = e(alpha)
local scale = e(beta) 
sort wfi
gen long r = _n
gen long n = _N
gen double u = (r - 0.5) / n    // plotting positions in (0,1)
gen double gamma_q = invgammap(`shape', u) * `scale'
summ wfi gamma_q
assert r(N)>0

* Q–Q plot (observed vs theoretical Gamma quantiles)
twoway ///
  (scatter wfi gamma_q, msymbol(o) msize(small) mcolor(blue) ///
      legend(label(1 "Observed quantiles"))) ///
  (line gamma_q gamma_q, sort lcolor(red) ///
      legend(label(2 "Reference line"))), ///
  xtitle("Theoretical Gamma quantiles") ///
  ytitle("Observed Frailty Index quantiles") ///
  title("Gamma Q–Q plot for FI") ///
  name(gammaqq, replace)

* Histogram with fitted Gamma density
twoway ///
  (hist wfi, density width(0.01) start(0) ///
         color(blue) fcolor(%70) lcolor(none) legend(label(1 "Empirical density"))) ///
  (function y = ( x^(`shape'-1) * exp(-x/`scale') ) / ///
                 ( exp(lngamma(`shape')) * `scale'^`shape' ), ///
        range(0 1) n(400) lwidth(medthick) lcolor(red) ///
        legend(label(2 "Fitted Gamma"))), ///
  xscale(range(0 1)) xlab(0(0.1)1) ///
  ytitle("Density") xtitle("FI") ///
  title("Histogram with fitted Gamma density") ///
  name(fi_gamma_overlay, replace)

*/
restore

* ==================================================================
* 10 Checking for change in the estimate from SRH
* ==================================================================

* Cubic time model without self-rated health
mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
		i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
		i.educ i.rural i.vigact i.smoke_status ///
		i.drink_status i.drink_status#c.time i.higov ///
		i.cogstat i.cogstat#c.time ///
		|| patid: female time_female, noconstant covariance(unstructured) ///
		|| patid: male time_male, noconstant covariance(unstructured) ///
		ml nolog stddev pweight(w1respwt) base 
estat ic

* Cubic time model with self-rated health
mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
		i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
		i.educ i.rural i.vigact i.smoke_status ///
		i.drink_status i.drink_status#c.time i.higov i.srh ///
		i.cogstat i.cogstat#c.time ///
		|| patid: female time_female, noconstant covariance(unstructured) ///
		|| patid: male time_male, noconstant covariance(unstructured) ///
		ml nolog stddev pweight(w1respwt) base 
estat ic

*===============================================================================
* 11 Session end
*===============================================================================

display as text "------------------------------------------------------------"
display as text "MHAS_Sensitivity_Analysis_completed"
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close