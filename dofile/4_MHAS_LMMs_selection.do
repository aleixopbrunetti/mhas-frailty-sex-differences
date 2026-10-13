********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Mixed regression models analysis for MHAS Frailty Index Study
********************************************************************************

* ==============================================================================
* 3. Load Data
* ==============================================================================
use "$cleanmhasdata/wave1to5_wide_deaths_code.dta", clear

* -- Display data summary
describe
summarize   
* -- Check for missing values
misstable summarize
* -- Check for duplicates
duplicates report
* -- Check for completely missing data across all waves
tab completely_lost
list w*fi completely_lost
tab completely_lost, missing
tab completely_lost sex, col chi

*=================================================================================
*3. Convert data from wide to long format & create time variable
*=================================================================================

* Create time variables for each wave (years since baseline)
gen time1 = 0   // Wave 1: 2001 - baseline
gen time2 = 2   // Wave 2: 2003
gen time3 = 11  // Wave 3: 2012
gen time4 = 14  // Wave 4: 2015
gen time5 = 17  // Wave 5: 2018

*-- Reshape the data to long format
drop (wave)
reshape long w@fi w@fc time@, i(patid) j(wave) string

label variable time "Time since baseline (in years)"

*========================================================================
* 4 Treatment of Variables
*========================================================================

* Variable centering
sum age, meanonly
gen age_c = age - r(mean)
sum wealth_gbp, meanonly
gen wealth_gbp_c = wealth_gbp - r(mean)

* Create age categories
gen age_cat = .
replace age_cat = 1 if age >= 50 & age < 60
replace age_cat = 2 if age >= 60 & age < 70
replace age_cat = 3 if age >= 70 & age < 80
replace age_cat = 4 if age >= 80
label define age_lbl 1 "50–59" 2 "60–69" 3 "70–79" 4 "80+"
label values age_cat age_lbl

* Create sex-specific time variables for level-2 heteroskedasticity
gen female = (sex == 2)
gen male = (sex == 1)
gen time_female = female * c.time
gen time_male = male * c.time

*========================================================================
* 5 Setting up the LMM time model
*========================================================================
xtset patid time

save "$cleanmhasdata/wave1to5_long.dta", replace

*========================================================================
* 6 Check for assumptions and skewness
*========================================================================
* -- Fit mixed-effects models with random intercepts

// Requirement: ssc install parmest

* Null model and ICC
mixed wfi ||patid: ,pweight(w1respwt) mle
estat icc 

* Model 0: simple model with random intercepts
mixed wfi age_c i.sex || patid: , pweight(w1respwt) mle nolog stddev base
esti store m0
estat ic
estat icc 

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_results.xlsx", ///
     sheet("Model 0") firstrow(variables) replace
restore

* Checking for assumptions
mixed wfi age_c i.sex || patid: , pweight(w1respwt) mle nolog stddev base
esti store m0
estat ic
estat icc 

* 1) Linearity and homoscedasticity
predict residuals, residuals
predict fitted, xb

scatter residuals fitted || lowess residuals fitted, title("Residuals vs Fitted Values")
graph export "$outputgraphs/residuals_vs_fitted_wholesample.jpg", replace
scatter residuals age_c || lowess residuals age_c, title("Residuals vs Age")
graph box residuals, over(sex) title("Residuals by Sex") ///
    name(resid_sex, replace)

* 2) Normality of residuals
histogram residuals, normal title("Histogram and Normal Density of Residuals") ///
     ylabel(,nogrid) xlabel(,nogrid)
    graph export "$outputgraphs/residuals_hist_wholesample.jpg", replace
qnorm residuals, title("Normal Q-Q Plot of Residuals")
graph export "$outputgraphs/qnorm_wholesample.jpg", replace

* Standardised residuals
predict std_residuals, rstandard
hist std_residuals, normal name(residuals_hist, replace) title("Histogram Density of Residuals")
graph export "$outputgraphs/std_residuals_hist_wholesample.jpg", replace
qnorm std_residuals, name(residuals_qnorm, replace) title("Normal Q-Q Plot of Residuals")
graph export "$outputgraphs/std_qnorm_wholesample.jpg", replace

* 3) Random effects - BLUPs
predict re_patid, reffects
histogram re_patid, normal
graph export "$outputgraphs/re_patid_hist_wholesample.jpg", replace
qnorm re_patid
graph export "$outputgraphs/re_patid_qnorm_wholesample.jpg", replace
qnorm re_patid if sex == 1, title("Normal Q-Q Plot of Random Effects (Male)")
graph export "$outputgraphs/re_patid_qnorm_male.jpg", replace
qnorm re_patid if sex == 2, title("Normal Q-Q Plot of Random Effects (Female)")
graph export "$outputgraphs/re_patid_qnorm_female.jpg", replace

* Mean FI over time by sex
preserve
collapse (mean) wfi, by(sex time)
twoway (line wfi time if sex==1) (line wfi time if sex==2)
restore

* ========================================================================
* 7 Weighted Model Comparison Analysis
* ==========================================================================
// Covariance structures compared: unstructured, independent, exchangeable, identity
// Structures assuming equally spaced measurements (AR, MA, Toeplitz) not used: waves are unequally spaced
// Banded structure not used: time modelled as continuous, not per-wave random effects
// Model selection: AIC/BIC and likelihood ratio tests
// Assumption: missing data MAR
// Interaction terms selected using Wald tests (testparm)

* --------------------------------------------------------------------------------
* Model structure: Fixed effects and interactions
* --------------------------------------------------------------------------------

* -- Exploring models with different covariance structures

*- Model 1: unstructured covariance structure
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
esti store m1
estat ic
display "AIC: " %20.0f -2*e(ll) + 2*e(rank)
display "BIC: " %20.0f -2*e(ll) + e(rank)*ln(e(N))

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_results.xlsx", ///
     sheet("Model 1") firstrow(variables) sheetreplace
restore

*-- Model 2: independent covariance structure
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time ///
     || patid: time, covariance(independent) pweight(w1respwt) ml nolog stddev
esti store m2
estat ic
display "AIC: " %20.0f -2*e(ll) + 2*e(rank)
display "BIC: " %20.0f -2*e(ll) + e(rank)*ln(e(N))

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_results.xlsx", ///
     sheet("Model 2") firstrow(variables) sheetreplace
restore

*-- Model 3: exchangeable covariance structure
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time ///
     || patid: time, covariance(exchangeable) pweight(w1respwt) ml nolog stddev
esti store m3
estat ic

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_results.xlsx", ///
     sheet("Model 3") firstrow(variables) sheetreplace
restore

*-- Model 4: identity covariance structure
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time ///
     || patid: time, covariance(identity) pweight(w1respwt) ml nolog stddev
esti store m4
estat ic
// Unstructured covariance (Model 1) selected: best fit

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_results.xlsx", ///
     sheet("Model 4") firstrow(variables) sheetreplace
restore

*---------------------------------------------------------------------------------
* Selection of the interaction terms
*---------------------------------------------------------------------------------
// Interactions with time tested sequentially in Model 1; non-significant terms removed

* Sex-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm i.sex#c.time

* Age-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm c.age_c#c.time

* Wealth-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm c.wealth_gbp_c#c.time

* Education-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c i.educ##c.time ///
     i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm i.educ#c.time

* Rural-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm i.rural#c.time

* Vigorous activity-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c i.educ ///
     i.rural i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm i.vigact#c.time

* Smoking status-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm i.smoke_status#c.time

* Self-rated health-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm i.srh#c.time

* Drinking status-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm i.drink_status#c.time

* Health insurance-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm i.higov#c.time

* Cognitive status-time
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev
testparm i.cogstat#c.time

*-- Model 5: Model with significant interactions (sex, age, drinking status, cognitive status x time)
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev base 
estimates store m5
estat ic

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_results.xlsx", ///
     sheet("Model 5") firstrow(variables) sheetreplace
restore

*----------------------------------------------------------------------------------------
* Diagnostic check - Model 5
*---------------------------------------------------------------------------------------
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) pweight(w1respwt) ml nolog stddev base

* Level-1 residuals
predict e_m5, rstandard   // standardized residuals
predict yhat_m5           // fitted values

hist e_m5, normal title("Level-1 standardized residuals (m5)")
qnorm e_m5, title("Q–Q plot: Level-1 residuals (m5)")
twoway (scatter e_m5 yhat_m5, msymbol(o) msize(vsmall)) (lowess e_m5 yhat_m5), ///
    yline(0) ///
    title("Residuals vs fitted values (m5)") ///
    ytitle("Std. residuals") ///
    xtitle("Fitted values")
scatter e_m5 time, yline(0) title("Residuals vs time (m5)") ///
    ytitle("Std. residuals") xtitle("Time")

* Level-2 residuals (BLUPs)
predict re_int_m5 re_time_m5, reffects
hist re_int_m5, normal title("Random intercepts (BLUPs, m5)")
qnorm re_int_m5, title("Q–Q plot: random intercepts (m5)")
hist re_time_m5, normal title("Random slopes for time (BLUPs, m5)")
qnorm re_time_m5, title("Q–Q plot: random slopes for time (m5)")

*=======================================================================================
* 8 Model with Heteroskedasticity level 1 and 2 - Growth-curve model
*=======================================================================================

*-------------------------------------------------------------------------------------
* Heteroskedasticity level 1
*-------------------------------------------------------------------------------------

* -- Heteroskedasticity level 1 (residual variance by sex) - saturated model
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time ///
     || patid: time, covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) residuals(independent, by(sex)) base 
esti store sat_m6_hetl1
estat ic

/* Simplify model mean structure by checking interactions:
testparm i.sex#c.time
testparm c.age_c#c.time
testparm c.wealth_gbp_c#c.time
testparm i.educ#c.time
testparm i.rural#c.time
testparm i.vigact#c.time
testparm i.smoke_status#c.time
testparm i.srh#c.time
testparm i.drink_status#c.time
testparm i.higov#c.time
testparm i.cogstat#c.time
*/

* -- Heteroskedasticity level 1 - final model
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) residuals(independent, by(sex)) base 
esti store m6_hetl1
estat ic
display "AIC: " %20.0f -2*e(ll) + 2*e(rank)
display "BIC: " %20.0f -2*e(ll) + e(rank)*ln(e(N))

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_results.xlsx", ///
     sheet("Model 6_hetl1") firstrow(variables) sheetreplace
restore

*-------------------------------------------------------------------------------------
* Heteroskedasticity level 1 - model check
*-------------------------------------------------------------------------------------
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: time, covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) residuals(independent, by(sex)) base 
esti store m6_hetl1
estat ic

*-- Level-1 residuals (m6_hetl1)
predict e_m6, rstandard      // standardized residuals
predict yhat_m6              // fitted values

hist e_m6, normal  title("Level-1 standardized residuals (m6_hetl1)")
qnorm e_m6, title("Q–Q plot: Level-1 residuals (m6_hetl1)")

twoway ///
    (scatter e_m6 yhat_m6, msymbol(o) msize(vsmall)) ///
    (lowess  e_m6 yhat_m6), ///
    yline(0) ///
    title("Residuals vs fitted (m6_hetl1)") ///
    ytitle("Std. residuals") xtitle("Fitted values")

* Residuals by sex
twoway (scatter e_m6 yhat_m6 if sex==1, msymbol(o) msize(vsmall)) ///
       (lowess  e_m6 yhat_m6 if sex==1), ///
       yline(0) title("Residuals vs time (males, m6_hetl1)")

twoway (scatter e_m6 yhat_m6 if sex==2, msymbol(o) msize(vsmall)) ///
       (lowess  e_m6 yhat_m6 if sex==2), ///
       yline(0) title("Residuals vs time (females, m6_hetl1)")

*-- Level-2 residuals (BLUPs, m6_hetl1)
predict re_int_m6 re_time_m6, reffects
hist re_int_m6, normal title("Random intercepts (BLUPs, m6_hetl1)")
qnorm re_int_m6, title("Q–Q plot: random intercepts (m6_hetl1)")

hist re_time_m6, normal title("Random slopes for time (BLUPs, m6_hetl1)")
qnorm re_time_m6, title("Q–Q plot: random slopes for time (m6_hetl1)")

*-------------------------------------------------------------------------------------
* Heteroskedasticity level 2
*-------------------------------------------------------------------------------------

* -- Heteroskedasticity level 2 (random effects by sex) - saturated model
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c##c.time ///
     i.educ##c.time i.rural##c.time i.vigact##c.time i.smoke_status##c.time i.srh##c.time ///
     i.drink_status##c.time i.higov##c.time i.cogstat##c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 
esti store sat_m7_hetl2
estat ic

/* Simplify model mean structure by checking interactions:
testparm i.sex#c.time
testparm c.age_c#c.time
testparm c.wealth_gbp_c#c.time
testparm i.educ#c.time
testparm i.rural#c.time
testparm i.vigact#c.time
testparm i.smoke_status#c.
testparm i.srh#c.time
testparm i.drink_status#c.time
testparm i.higov#c.time
testparm i.cogstat#c.time
*/

* -- Heteroskedasticity level 2 - final model 
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 
esti store m7_hetl2
estat ic 
display "AIC: " %20.0f -2*e(ll) + 2*e(rank)
display "BIC: " %20.0f -2*e(ll) + e(rank)*ln(e(N))

preserve 
estat sd, post
parmest, norestore
export excel using "$outputtables/mixed_results.xlsx", ///
     sheet("Model 7_hetl2") firstrow(variables) sheetreplace
restore

*===========================================================================================
* 9 Comparing models for best fit
*===========================================================================================

estimates stats m0 sat_m6_hetl1 sat_m7_hetl2 m1 m2 m3 m4 m5 m6_hetl1 m7_hetl2

// Final model: m7_hetl2 (lowest AIC and simpler mean structure)

/* Sensitivity check - likelihood ratio tests:
lrtest m0 m1, force
lrtest m1 m5, force
lrtest m5 m6_hetl1, force
lrtest m6_hetl1 m7_hetl2, force
*/

* ===========================================================================================
* 10 Final Model Assumptions and Diagnostics: Model m7_hetl2
* ===========================================================================================

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

*-- Level-1 residuals (m7)
predict e_m7, rstandard            // standardized residuals
predict yhat_m7                    // fitted values

hist e_m7, normal ///
     ylabel(, nogrid) xlabel(, nogrid)
graph export "$outputgraphs/e_m7_histogram.jpg", replace
qnorm e_m7, title("Q–Q plot: Level-1 residuals (m7)")
graph export "$outputgraphs/e_m7_qnorm.jpg", replace

twoway (scatter e_m7 yhat_m7, msymbol(o) msize(vsmall)) ///
       (lowess  e_m7 yhat_m7, legend(label(2 "LOWESS"))), ///
       yline(0) ///
       title("Residuals vs fitted") ///
       ytitle("Std. residuals") xtitle("Fitted values") ///
       ylabel(, nogrid) ///
       xlabel(, nogrid)
graph export "$outputgraphs/e_m7_residuals_vs_fitted.jpg", replace

* Residuals by sex
twoway (scatter e_m7 yhat_m7 if sex==1, msymbol(o) msize(vsmall)) ///
       (lowess  e_m7 yhat_m7 if sex==1), ///
       yline(0) ///
       title("Residuals vs fitted — Males (m7)")
twoway (scatter e_m7 yhat_m7 if sex==2, msymbol(o) msize(vsmall)) ///
       (lowess  e_m7 yhat_m7 if sex==2), ///
       yline(0) ///
       title("Residuals vs fitted — Females (m7)")

*-- Level-2 residuals (BLUPs, m7)
predict re_femaleintercept re_femaleslope re_maleintercept re_maleslope, reffects

preserve
gen pickone if e(sample)=1
contract patid re_femaleintercept re_femaleslope re_maleintercept re_maleslope female male

* Female
hist re_femaleintercept if female, normal title("BLUPs: random intercept (female) (m7)")
graph export "$outputgraphs/BLUP_finalmodel_hist_female.jpg", replace
qnorm re_femaleintercept if female, title("Q–Q: random intercept (female) (final model-level2)")
graph export "$outputgraphs/BLUP_intercept_finalmodel_qnorm_female.jpg", replace

hist re_femaleslope if female, normal title("BLUPs: random slope time_female (m7)")
qnorm re_femaleslope if female, title("Q–Q: random slope time_female (final model-level2)")
graph export "$outputgraphs/BLUP_slope_finalmodel_qnorm_female.jpg", replace

* Male
hist re_maleintercept if male, normal title("BLUPs: random intercept (male) (m7)")
qnorm re_maleintercept if male, title("Q–Q: random intercept (male) (m7)")
graph export "$outputgraphs/BLUP_intercept_finalmodel_qnorm_male.jpg", replace
hist re_maleslope if male, normal title("BLUPs: random slope time_male (m7)")
qnorm re_maleslope if male, title("Q–Q: random slope time_male (final model-level2)")
graph export "$outputgraphs/BLUP_slope_finalmodel_qnorm_male.jpg", replace

restore

*===============================================================================
* 11 Session end
*===============================================================================

display as text "------------------------------------------------------------"
display as text "MHAS_LMMs_selection_completed"
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close