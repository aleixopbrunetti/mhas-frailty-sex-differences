********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: LMM with non-linear time and sex interaction tests (final model)
********************************************************************************

* ==============================================================================
* 3. Load Data
* ==============================================================================
use "$cleanmhasdata/wave1to5_wide_deaths_code.dta", clear

* -- Display data summary
describe
summarize   
* -- Check for duplicates
duplicates report

*=================================================================================
*3. Convert data from wide to long format & create time variable
*=================================================================================

* Time variables for each wave (years since baseline)
gen time1 = 0   // Wave 1: 2001 - baseline
gen time2 = 2   // Wave 2: 2003
gen time3 = 11  // Wave 3: 2012
gen time4 = 14  // Wave 4: 2015
gen time5 = 17  // Wave 5: 2018

*-- Reshape to long format
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

* Age categories
gen age_cat = .
replace age_cat = 1 if age >= 50 & age < 60
replace age_cat = 2 if age >= 60 & age < 70
replace age_cat = 3 if age >= 70 & age < 80
replace age_cat = 4 if age >= 80
label define age_lbl 1 "50–59" 2 "60–69" 3 "70–79" 4 "80+"
label values age_cat age_lbl

* Sex-specific time variables for level-2 heteroskedasticity
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
* 6 LMM model with interactions
*========================================================================

* Previous final model (linear time)
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* Same model with interactions written out separately
mixed wfi i.sex c.time i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
esti store m7
estat ic 
display "AIC: " %20.0f -2*e(ll) + 2*e(rank)
display "BIC: " %20.0f -2*e(ll) + e(rank)*ln(e(N))

* -----------------------------------------------------------------------------------------
* 6.1 Non-linear time: quadratic and cubic terms
* ------------------------------------------------------------------------------------------

* Quadratic and cubic time terms, allowed to differ by sex
mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time  ///
	i.sex#c.time i.sex#c.time#c.time i.sex#c.time#c.time#c.time ///
	c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* Test for sex difference in non-linearity
test 2.sex#c.time#c.time 2.sex#c.time#c.time#c.time
* Test for overall non-linearity
test c.time#c.time c.time#c.time#c.time

* Common non-linear time terms for men and women
mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
	i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
	 
test c.time#c.time c.time#c.time#c.time

* -----------------------------------------------------------------------------------------
* 6.2 Further interaction terms: sex × (education, SRH, wealth, age, rural)
* ------------------------------------------------------------------------------------------
* Each tested in the linear-time model and again allowing for non-linear time
* Trajectory test: sex × covariate × time; overall test: adds sex × covariate (baseline)

*1) Sex # Education (#time)

mixed wfi i.sex c.time i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.educ#c.time  i.sex#i.educ  i.sex#i.educ#c.time   ///
	 i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
	 
* Differential effect of education on FI trajectory in men vs women
test 2.sex#2.educ#c.time 2.sex#3.educ#c.time	

* Differential effect of education on FI in men vs women (baseline or follow-up)
test 2.sex#2.educ#c.time 2.sex#3.educ#c.time  2.sex#2.educ 2.sex#3.educ

* Allowing for non-linear time
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

mixed wfi i.sex c.time i.sex#c.time c.age_c c.age_c#c.time ///
	c.wealth_gbp_c c.wealth_gbp_c#c.time ///
	c.wealth_gbp_c#i.sex c.wealth_gbp_c#i.sex#c.time ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* Differential effect of wealth on FI trajectory in men vs women
test 2.sex#c.wealth_gbp_c#c.time

* Differential effect of wealth on FI in men vs women (baseline or follow-up)
test 2.sex#c.wealth_gbp_c#c.time 2.sex#c.wealth_gbp_c 
	 
* Allowing for non-linear time
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

test 2.sex#c.wealth_gbp_c#c.time 
test 2.sex#c.wealth_gbp_c#c.time 2.sex#c.wealth_gbp_c 

*4) Sex # Age category (#time)

* Linear-time model
mixed wfi i.sex c.time i.sex#c.time c.age_c c.age_c#c.time ///
     i.age_cat i.age_cat#c.time i.sex#i.age_cat i.sex#i.age_cat#c.time ///
     c.wealth_gbp_c i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status i.drink_status#c.time i.higov ///
     i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* Trajectory
testparm i.sex#i.age_cat#c.time

* Overall (baseline + trajectory)
testparm i.sex#i.age_cat i.sex#i.age_cat#c.time

* Allowing for non-linear time
mixed wfi i.sex c.time i.sex#c.time c.time#c.time c.time#c.time#c.time ///
     c.age_c c.age_cat#c.time ///
     i.age_cat i.age_cat#c.time i.sex#i.age_cat i.sex#i.age_cat#c.time ///
     c.wealth_gbp_c i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status i.drink_status#c.time i.higov ///
     i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

testparm i.sex#i.age_cat#c.time
testparm i.sex#i.age_cat i.sex#i.age_cat#c.time

*5) Sex # Rural (#time)

mixed wfi i.sex c.time i.sex#c.time c.age_c c.age_c#c.time ///
	i.rural i.rural#c.time i.sex#i.rural i.sex#i.rural#c.time ///
	c.wealth_gbp_c i.educ i.vigact i.smoke_status i.srh ///
	i.drink_status i.drink_status#c.time i.higov ///
	i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
	 
* Trajectory
test 2.sex#1.rural#c.time

* Overall (baseline + trajectory)
test 2.sex#1.rural#c.time 2.sex#1.rural

* Allowing for non-linear time
mixed wfi i.sex c.time i.sex#c.time c.time#c.time c.time#c.time#c.time ///
     c.age_c c.age_c#c.time ///
     i.rural i.rural#c.time i.sex#i.rural i.sex#i.rural#c.time ///
     c.wealth_gbp_c i.educ i.vigact i.smoke_status i.srh ///
     i.drink_status i.drink_status#c.time i.higov ///
     i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base
	 
test 2.sex#1.rural#c.time
test 2.sex#1.rural#c.time 2.sex#1.rura

// No additional sex interactions retained (education, SRH, wealth, age, rural)

*===============================================================================
* 7. Final model and diagnostics (m11 specification)
*===============================================================================

mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
	i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

estat ic 
display "AIC: " %20.0f -2*e(ll) + 2*e(rank)
display "BIC: " %20.0f -2*e(ll) + e(rank)*ln(e(N))

*-- Level-1 residuals (m11)
predict e_m11, rstandard            // standardized residuals
predict yhat_m11                    // fitted values

hist e_m11, normal ///
     ylabel(, nogrid) xlabel(, nogrid)
graph export "$outputgraphs/e_m11_histogram.jpg", replace
qnorm e_m11, title("Q–Q plot: Level-1 residuals (m11)")
graph export "$outputgraphs/e_m11_qnorm.jpg", replace

twoway (scatter e_m11 yhat_m11, msymbol(o) msize(vsmall)) ///
       (lowess  e_m11 yhat_m11, legend(label(2 "LOWESS"))), ///
       yline(0) ///
       title("Residuals vs fitted") ///
       ytitle("Std. residuals") xtitle("Fitted values") ///
       ylabel(, nogrid) ///
       xlabel(, nogrid)
graph export "$outputgraphs/e_m11_residuals_vs_fitted.jpg", replace

*-- Level-2 residuals (BLUPs, m11)
predict re_femaleintercept re_femaleslope re_maleintercept re_maleslope, reffects

preserve
gen pickone if e(sample)=1
contract patid re_femaleintercept re_femaleslope re_maleintercept re_maleslope female male

* Female
hist re_femaleintercept if female, normal title("BLUPs: random intercept (female) (m11)")
graph export "$outputgraphs/BLUP_finalmodel_hist_female.jpg", replace
qnorm re_femaleintercept if female, title("Q–Q: random intercept (female) (m11-level2)")
graph export "$outputgraphs/BLUP_intercept_finalmodel_qnorm_female.jpg", replace

hist re_femaleslope if female, normal title("BLUPs: random slope time_female (m11)")
qnorm re_femaleslope if female, title("Q–Q: random slope time_female (m11-level2)")
graph export "$outputgraphs/BLUP_slope_finalmodel_qnorm_female.jpg", replace

* Male
hist re_maleintercept if male, normal title("BLUPs: random intercept (male) (m11)")
qnorm re_maleintercept if male, title("Q–Q: random intercept (male) (m11)")
graph export "$outputgraphs/BLUP_intercept_finalmodel_qnorm_male.jpg", replace
hist re_maleslope if male, normal title("BLUPs: random slope time_male (m11)")
qnorm re_maleslope if male, title("Q–Q: random slope time_male (m11-level2)")
graph export "$outputgraphs/BLUP_slope_finalmodel_qnorm_male.jpg", replace

restore

*===============================================================================
* 7. Final model without SRH
*===============================================================================

mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
	i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

*===============================================================================
* 8. VIF diagnostic 
*===============================================================================

* Centre time and create polynomial terms
summarize time
gen time_c = time - r(mean)
gen time_c2 = time_c^2
gen time_c3 = time_c^3

* Correlations among centred polynomial terms
pwcorr time_c time_c2 time_c3

* mixed does not support estat vif: linear regression used as diagnostic proxy
regress wfi i.sex c.time_c c.time_c#c.time_c c.time_c#c.time_c#c.time_c ///
    i.sex#c.time_c c.age_c i.educ i.rural i.vigact i.smoke_status i.srh ///
    i.drink_status i.higov i.cogstat [pweight=w1respwt]
estat vif

*===============================================================================
* 9. Orthogonal polynomials for time
*===============================================================================

orthpoly time, generate(otime1 otime2 otime3) degree(3)

* Refit with orthogonal time basis (otime1 used for all linear-time interactions)
mixed wfi i.sex otime1 otime2 otime3 ///
    i.sex#c.otime1 c.age_c c.age_c#c.otime1 c.wealth_gbp_c ///
    i.educ i.rural i.vigact i.smoke_status i.srh ///
    i.drink_status i.drink_status#c.otime1 i.higov ///
    i.cogstat i.cogstat#c.otime1 ///
    || patid: female time_female, noconstant covariance(unstructured) ///
    || patid: male time_male, noconstant covariance(unstructured) ///
    ml nolog stddev pweight(w1respwt) base
	
estat ic 
display "AIC: " %20.0f -2*e(ll) + 2*e(rank)
display "BIC: " %20.0f -2*e(ll) + e(rank)*ln(e(N))

* Female–male difference at baseline (time = 0)
preserve 
sum otime1 if time==0, meanonly
local v0 = r(mean)
lincom 2.sex + `v0'*2.sex#c.otime1
restore 

*===============================================================================
* 10. Post-estimation
*===============================================================================

mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
	i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* Predicted FI by sex at each wave
margins sex, vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_sex_wave1to5.dta", replace)
	
marginsplot, title("Predicted FI trajectories by sex") ///
    ylabel(, nogrid) xlabel(,nogrid) ///
    name(margins_sex, replace)
graph export "$outputgraphs/margins_sex.jpg", replace

* Slope ratio: female / male
nlcom (slope_ratio: ( _b[c.time] + _b[2.sex#c.time] ) / _b[c.time])
nlcom (pct_faster: 100 * ( _b[2.sex#c.time] / _b[c.time] ))

* Relative female/male gap at t = 17 (ratio of predicted means)
nlcom (gap_ratio_t17: ( _b[2.sex] + 17*_b[2.sex#c.time] + _b[_cons] + 17*_b[c.time] ) / ///
                      ( _b[_cons] + 17*_b[c.time] ))

*===============================================================================
* Age categorical by sex
*===============================================================================
*------------------------------------------------------------------------------
* Post-estimation for age categories
* -----------------------------------------------------------------------------
 
mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
	i.sex#c.time i.age_cat i.age_cat#c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* Predicted FI trajectories by age band and sex (waves 1 to 5)
margins sex#age_cat, at(time=(0 2 11 14 17)) ///
    vce(unconditional) ///
    saving("$working/margins_sex_age_cat_wave1to5.dta", replace)

marginsplot, xdimension(time) plotdimension(age_cat sex) recast(line) noci ///
    plot1opts(lcolor(green)   lpattern(longdash) lwidth(medthin)) ///
    plot2opts(lcolor(green)   lpattern(solid)) ///
    plot3opts(lcolor(blue)   lpattern(longdash) lwidth(medthin)) ///
    plot4opts(lcolor(blue)   lpattern(solid)) ///
    plot5opts(lcolor(purple) lpattern(longdash) lwidth(medthin)) ///
    plot6opts(lcolor(purple) lpattern(solid)) ///
    plot7opts(lcolor(red) lpattern(longdash) lwidth(medthin)) ///
    plot8opts(lcolor(red) lpattern(solid)) ///
    legend(order(1 "Male, 50–59" 2 "Female, 50–59" 3 "Male, 60–69" 4 "Female, 60–69" ///
                 5 "Male, 70–79" 6 "Female, 70–79" 7 "Male, 80+" ///
                 8 "Female, 80+") ///
           cols(2) region(lstyle(none)) ///
           position(6)) /// 
    ytitle("Frailty Index") xtitle("Time since baseline (years)") ///
    ylabel(, nogrid) xlabel(, nogrid) ///
    title("Predicted FI Trajectories by Sex and Age Category") ///
    name(age_band_combined, replace)
graph export "$outputgraphs/fi_trajectories_age_bands_by_sex.jpg", replace

*-----------------------------------------------------------------------------
* Sex # age category (#time) with non-linear time
*----------------------------------------------------------------------------

mixed wfi i.sex c.time i.sex#c.time c.time#c.time c.time#c.time#c.time ///
     i.age_cat i.age_cat#c.time i.sex#i.age_cat i.sex#i.age_cat#c.time ///
	 c.wealth_gbp_c i.educ i.rural i.vigact i.smoke_status i.srh ///
	 i.drink_status i.drink_status#c.time i.higov ///
	 i.cogstat i.cogstat#c.time ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* Predicted FI trajectories by age band and sex (waves 1 to 5)
margins sex#age_cat, at(time=(0 2 11 14 17)) ///
    vce(unconditional) ///
    saving("$working/margins_sex_age_cat_non-linearity.dta", replace)

marginsplot, xdimension(time) plotdimension(age_cat sex) recast(line) noci ///
    plot1opts(lcolor(green)   lpattern(longdash) lwidth(medthin)) ///
    plot2opts(lcolor(green)   lpattern(solid)) ///
    plot3opts(lcolor(blue)   lpattern(longdash) lwidth(medthin)) ///
    plot4opts(lcolor(blue)   lpattern(solid)) ///
    plot5opts(lcolor(purple) lpattern(longdash) lwidth(medthin)) ///
    plot6opts(lcolor(purple) lpattern(solid)) ///
    plot7opts(lcolor(red) lpattern(longdash) lwidth(medthin)) ///
    plot8opts(lcolor(red) lpattern(solid)) ///
    legend(order(1 "Male, 50–59" 2 "Female, 50–59" 3 "Male, 60–69" 4 "Female, 60–69" ///
                 5 "Male, 70–79" 6 "Female, 70–79" 7 "Male, 80+" ///
                 8 "Female, 80+") ///
           cols(2) region(lstyle(none)) ///
           position(6)) /// 
    ytitle("Frailty Index") xtitle("Time since baseline (years)") ///
    ylabel(, nogrid) xlabel(, nogrid) ///
    title("Predicted FI Trajectories by Sex and Age Category interaction non-linearity") ///
    name(age_band_combined, replace)
graph export "$outputgraphs/fi_trajectories_age_bands_by_sex_non-linearity.jpg", replace

* ==============================================================================
* Supplementary Figure S1 — FI distribution by sex, overlaid by wave
* ==============================================================================

use "$cleanmhasdata/wave1to5_wide_deaths_code.dta", clear

* Keep variables needed for the figure and reshape to long
preserve

keep patid sex w1fi w2fi w3fi w4fi w5fi

reshape long w@fi, i(patid) j(wave)
rename wfi fi

* Labels for legend and panel titles
label define sexlbl 1 "Men" 2 "Women", replace
label values sex sexlbl

label define wavelbl 1 "Wave 1 (2001)" 2 "Wave 2 (2003)" ///
                     3 "Wave 3 (2012)" 4 "Wave 4 (2015)" ///
                     5 "Wave 5 (2018)", replace
label values wave wavelbl

* Overlaid kernel densities by sex, one panel per wave
twoway ///
    (kdensity fi if sex==1, lcolor(navy)      lpattern(solid) lwidth(medthick)) ///
    (kdensity fi if sex==2, lcolor(cranberry) lpattern(dash)  lwidth(medthick)) ///
    , by(wave, cols(2) xrescale note("") ///                   
         title("Frailty Index distribution by sex across waves", size(medium)) ///
         legend(pos(6) order(1 "Men" 2 "Women") rows(1) ///
                size(*1.15) symxsize(*1.7) keygap(*1.3) colgap(*1.5) ///
                region(margin(medsmall))) ///
         graphregion(color(white)) plotregion(color(white))) ///
      xtitle("Frailty index (0-1)") ytitle("Probability density") ///
      xlabel(0(.2)1) xscale(range(0 1)) ///  
	  legend(order(1 "Men" 2 "Women") rows(1)) ///
      scheme(s1color) ///
      subtitle(, fcolor(gs14) lcolor(black)) ///
      xsize(8.27) ysize(11.69)                                  
	  
graph export "$outputgraphs/figS1_fi_by_sex_overlay.pdf", replace
graph export "$outputgraphs/mfigS1_fi_by_sex_overlay.tif", replace width(4000)

restore 

*===============================================================================
* 11 Session end
*===============================================================================

display as text "------------------------------------------------------------"
display as text "MHAS_LMMs_interactions_completed"
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close