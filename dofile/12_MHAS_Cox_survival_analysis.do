********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Survival analysis - frailty index, sex and mortality (Cox models)
********************************************************************************

* ==============================================================================
* 3. Load Data
* ==============================================================================
use "$cleanmhasdata/wave1to5_long.dta", clear 

* -- Display data summary
describe
* -- Check for duplicates
duplicates report

*=================================================================================
*4. Setting the time intervals for survival analysis
*=================================================================================

* Key variables:
* - patid: person ID
* - time: years since baseline (0, 2, 11, 14, 17)
* - death: 1 = dead, 0 = alive
* - dyear, dmonth: year/month of death (month may be missing)

*------------------------
* 4.1 Study baseline year & end
*------------------------
local baseyear   = 2001
local study_end  = 18

* -----------------------
* 4.2 Build (t0, t1] rows
* -----------------------
sort patid time
by patid: gen double t0 = time
by patid: gen double t1 = time[_n+1]                 // next wave
by patid: egen double last_time = max(time)
replace t1 = last_time if missing(t1)                // close final interval

* -----------------------
* 4.3 Exit time per person
* -----------------------
gen byte died = (death==1)

* Time to death in fractional years since baseline (mid-month; mid-year if month missing)
gen double death_time = .
replace death_time = (dyear - `baseyear') + (dmonth-0.5)/12 if died==1 & dyear<. & dmonth<. 
replace death_time = (dyear - `baseyear') + 0.5             if died==1 & dyear<. & missing(dmonth)

* Check bounds
replace death_time = 0                  if death_time<0 & death_time<.
replace death_time = `study_end'        if death_time>`study_end' & death_time<.

* Exit time: death time if died, otherwise last observed wave
gen double exit = cond(died==1 & death_time<., death_time, last_time)

* Cap at study end
replace exit = `study_end' if exit>`study_end' & exit<.

* ----------------------------
* 4.4 Mark the failure interval
* ----------------------------

replace t1 = exit if t1 > exit

gen byte fail = (died==1 & exit<. & t1==exit)

* Ensure max one failure per subject
bysort patid: egen nfail = total(fail)
drop nfail

* =======================================================================
*5. stset the data
* =======================================================================
stset t1 [pweight = w1respwt], enter(t0) failure(fail) scale(1)

stsum, by(sex)
* Distribution of event times
tab _t if _d==1

* Drop intervals starting at or after exit
drop if t0 >= exit

* --- Sanity checks
count if death==1
display as result "Rows with death==1: " r(N)
count if fail==1
display as result "Failure intervals (fail==1): " r(N)
stsum

* ---------------------------------------------------------------------
* 5.1 Frailty index
* ---------------------------------------------------------------------

* Centering
sum wfi
gen wfi_c = wfi - r(mean)

sum wfi_c
* Rescale FI for interpretability
gen wfi_c10 = wfi_c*10   // HR per 0.1 increase in FI
gen wfi_c1 = wfi_c*34    // HR per additional deficit (34-item FI; 1 deficit ≈ 0.029 FI units)

* ======================================================================
* 6. Survival Analysis
* ======================================================================

*----------------------------------------------------------------------------
* Kaplan-Meier survival curves
*----------------------------------------------------------------------------
sts graph, by(sex) ///
    title("Kaplan-Meier Survival by Sex") ///
    ytitle("Survival Probability") ///
    xtitle("Time since baseline (years)") ///
    ylabel(, nogrid) xlabel(, nogrid) ///
    legend(order(1 "Male" 2 "Female") cols(2) position(6) ring(1) /// 
	region(lstyle(none))) ///
    name(km_sex, replace)
graph export "$outputgraphs/km_survival_curves.jpg", replace
graph export "$outputgraphs/km_survival_curves.tif", replace width(4000)
// Note: confidence intervals not available with pweights

* KM cumulative mortality (1 - S(t))
sts graph, failure by(sex) ///
    title("Kaplan–Meier Cumulative Mortality by Sex") ///
    ytitle("Cumulative mortality (1 - S(t))") ///
    xtitle("Time since baseline (years)") ///
	ylabel(, nogrid) xlabel(, nogrid) ///
    legend(order(1 "Male" 2 "Female") cols(2) position(6) ring(1) /// 
	region(lstyle(none))) ///
    name(km_mort_rt, replace)
graph export "$outputgraphs/km_mortality_risktable.jpg", replace
graph export "$outputgraphs/km_mortality_risktable.tif", replace width(4000)

* Survey-weighted test of survival difference by sex (Wald test from Cox model; Figure 2)
svy: stcox i.sex, nohr
test 2.sex
local lr_p = r(p)
display "Log-rank p-value = `lr_p'"

*----------------------------------------------------------------------------
* 6.2 Cox proportional hazards model
*----------------------------------------------------------------------------

* FI (original scale) with sex interaction
svy: stcox  i.sex c.wfi i.sex#c.wfi  c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            base

* FI per 0.1 with sex interaction
svy: stcox  i.sex c.wfi_c10 i.sex#c.wfi_c10  c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            base
// Sex × FI interaction not retained (not statistically significant)

* Final model: FI per 0.1, no interaction
svy: stcox  i.sex c.wfi_c10  c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            base

* ======================================================================
* 7 Checking proportional hazards assumption
* ====================================================================== 

stset t1 [pweight = w1respwt], enter(t0) failure(fail) scale(1)

stcox  i.sex c.wfi_c10 c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
       i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
       base
predict mresid, mgale

* Schoenfeld residuals test
estat phtest, detail

* Functional form of FI (martingale residuals)
scatter mresid wfi_c10, mcolor(blue) ///
    || lowess mresid wfi_c10, lcolor(red) ///
    yline(0, lpattern(dash)) ///
    ytitle("Martingale residuals") ///
    xtitle("Frailty Index (per 0.1 increase)") ///
    title("Functional form check for frailty index") ///
    legend(order(1 "Martingale residual" 2 "LOWESS smoothed curve") cols(1) position(4)) ///
    ylabel(, nogrid) xlabel(, nogrid)  
graph export "$outputgraphs/martingale_residuals_frailty.jpg", replace

* Schoenfeld residual plots
estat phtest, plot(2.sex) ylabel(, nogrid) xlabel(, nogrid)
graph export "$outputgraphs/schoenfeld_residuals_sex.jpg", replace
estat phtest, plot(c.age_c) ylabel(, nogrid) xlabel(, nogrid)
estat phtest, plot(c.wfi_c10) ylabel(, nogrid) xlabel(, nogrid)
graph export "$outputgraphs/schoenfeld_residuals_frailty.jpg", replace
estat phtest, plot(2.sex#c.wfi_c10 ) ylabel(, nogrid) xlabel(, nogrid)
graph export "$outputgraphs/schoenfeld_residuals_interaction.jpg", replace

* Log-log plot of survival
stphplot, by(sex) ylabel(, nogrid) xlabel(, nogrid)
graph export "$outputgraphs/loglog_survival.jpg", replace

* Observed (KM) vs predicted (Cox) survival
stcoxkm, by(sex) ylabel(, nogrid) xlabel(, nogrid)
graph export "$outputgraphs/cox_km_survival.jpg", replace

* Time-varying effect of FI
stcox   i.sex c.wfi_c10  c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
        i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
        tvc(wfi_c10) vce(cluster patid) base

* ----------------------------------------------------------------------
* 7.1 Cox time-varying effects
* ----------------------------------------------------------------------  

stcox   i.sex c.wfi_c10 c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
        i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
        tvc(wfi_c10) vce(cluster patid)  base

* Time-varying effects of sex, FI, age and drinking status
stcox  i.sex c.wfi_c10  c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            tvc(i.sex c.wfi_c10 c.age_c i.drink_status) ///
			vce(cluster patid) base	

* Adjusting for non-linear time
svy: stcox  i.sex c.wfi_c10 c.time c.time#c.time c.time#c.time#c.time c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            base

* ======================================================================
* 8 Post-estimation: final Cox model
* =================================================================
svy: stcox  i.sex c.wfi_c10 c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            base	
			
* Adjusted cumulative mortality by sex
stcurve, failure at1(sex=1) at2(sex=2) ///
    title("Adjusted Cumulative Mortality by Sex (Cox model)") ///
    ytitle("Adjusted cumulative mortality") ///
    xtitle("Time since baseline (years)") ///
	ylabel(, nogrid) xlabel(, nogrid) ///
    legend(order(1 "Male" 2 "Female") cols(2) position(6) ring(1) /// 
	region(lstyle(none))) ///
    name(adj_cox_mort, replace)
graph export "$outputgraphs/adj_cox_mort.jpg", replace
graph export "$outputgraphs/adj_cox_mort.tif", replace width(4000)

* Relative hazard by sex
margins sex, expression(exp(predict(xb)))

margins r.sex, expression(exp(predict(xb)))
nlcom pct_lower: 100*(1 - exp(_b[2.sex]))

margins, dydx(wfi_c10)

nlcom HR_overall: exp(_b[wfi_c10])

* HR per 0.1 FI
lincom wfi_c10, eform

* Adjusted survival curves by sex
stcurve, survival at1(sex=1) at2(sex=2) ///
     title("Adjusted Survival Curves by Sex (Cox, weighted)") ///
     ytitle("Survival probability") xtitle("Time since baseline (years)") ///
     ylabel(0(0.25)1, nogrid) xlabel(, nogrid) ///
     legend(order(1 "Male" 2 "Female"))
graph export "$outputgraphs/cox_survival_curves.jpg", replace 
graph export "$outputgraphs/cox_survival_curves.tif", replace width(4000)

* ======================================================================
* 8.1 Adjusted cumulative mortality by sex with 95% CI
* =================================================================

preserve 

* 1. Fit with female (sex=2) as reference
stcox ib2.sex c.wfi_c10 c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
      i.smoke_status i.srh i.drink_status i.higov i.cogstat, vce(robust)

* 2. Baseline cumulative hazard and male coefficient
predict double H0_f, basechazard
local b_m  = _b[1.sex]
local se_m = _se[1.sex]

* 3. Male cumulative mortality and CI
gen double f_male_v2    = 1 - exp(-H0_f * exp(`b_m'))
gen double f_male_lo    = 1 - exp(-H0_f * exp(`b_m' - 1.96*`se_m'))
gen double f_male_hi    = 1 - exp(-H0_f * exp(`b_m' + 1.96*`se_m'))

* 4. Fit with male as reference
stcox i.sex c.wfi_c10 c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
      i.smoke_status i.srh i.drink_status i.higov i.cogstat, vce(robust)

predict double H0_m, basechazard
local b_f  = _b[2.sex]
local se_f = _se[2.sex]

* 5. Female cumulative mortality and CI
gen double f_female    = 1 - exp(-H0_m * exp(`b_f'))
gen double f_female_lo = 1 - exp(-H0_m * exp(`b_f' - 1.96*`se_f'))
gen double f_female_hi = 1 - exp(-H0_m * exp(`b_f' + 1.96*`se_f'))

* Male point estimate
gen double f_male = 1 - exp(-H0_f * exp(`b_m'))

* 6. One row per time point
keep _t f_male f_male_lo f_male_hi f_female f_female_lo f_female_hi
bysort _t: keep if _n == _N
drop if missing(_t)

* 7. Plot with CIs
twoway ///
    (rarea f_male_lo    f_male_hi    _t, color(navy%20)     fintensity(50)) ///
    (rarea f_female_lo  f_female_hi  _t, color(cranberry%20) fintensity(50)) ///
    (line f_male    _t, lcolor(navy)      lwidth(medthick)) ///
    (line f_female  _t, lcolor(cranberry) lwidth(medthick)) ///
    , ylabel(, nogrid) xlabel(, nogrid) ///
    title("Adjusted cumulative mortality by sex") ///
    ytitle("Adjusted cumulative mortality") ///
    xtitle("Time since baseline (years)") ///
    legend(order(3 "Male" 1 "Male 95% CI" 4 "Female" 2 "Female 95% CI") ///
           cols(2) position(6) ring(1) region(lstyle(none))) ///
    name(adj_cox_mort_ci, replace)

graph export "$outputgraphs/adj_cox_mort_ci.jpg", replace
graph export "$outputgraphs/adj_cox_mort_ci.tif", replace width(4000)

restore

* ==============================================================================
* 9 Finish
* ==============================================================================

display as text "------------------------------------------------------------"
display "==== Cox survival analysis_ COMPLETED SUCCESSFULLY ===="
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close