********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Survival analysis - frailty index and mortality (Cox models)
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
* ----------------------
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

*----------------------------------------------------------------------------
* Cox proportional hazards model
*----------------------------------------------------------------------------

* -- Continuous FI (per 0.1)
svy: stcox   i.sex##c.wfi_c10  c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
              base

stcurve, survival at1(sex=1) at2(sex=2) ///
     title("Adjusted Survival Curves by Sex (Cox, weighted)") ///
     ytitle("Survival probability") xtitle("Time since baseline (years)") ///
     ylabel(0(0.25)1, nogrid) xlabel(, nogrid) ///
     legend(order(1 "Male" 2 "Female"))
graph export "$outputgraphs/cox_survival_curves.jpg", replace

* -- Frailty category (frail vs not frail)
svy: stcox i.sex##i.wfc c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
             base

* -- Continuous FI (per additional deficit)
svy: stcox   i.sex##c.wfi_c1  c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
              base

* ----------------------------------------------------------------------
* Checking proportional hazards assumption
* ----------------------------------------------------------------------    

stset t1 [pweight = w1respwt], enter(t0) failure(fail) scale(1)

stcox  i.sex##c.wfi_c10  c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
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
stcox   i.sex##c.wfi_c10 c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            tvc(wfi_c10) base

stcox   i.sex##c.wfi_c10 c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            tvc(wfi_c10) texp(_t) base

* Flexible parametric (Royston-Parmar) models
ssc install stpm2, replace
ssc install rcsgen, replace

stpm2   i.sex##c.wfi_c10 c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            scale(hazard) df(3) eform

* With time-varying effect of FI
stpm2   i.sex##c.wfi_c10 c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            scale(hazard) df(3) tvc(wfi_c10) dftvc(2) eform

predict hr_wfi, hrnumerator(wfi_c10 0.1) hrdenominator(wfi_c10 0) timevar(_t)
twoway line hr_wfi _t, ytitle("HR per 0.1 FI") xtitle("Years")

* ----------------------------------------------------------------------
* Post-estimation: hazard ratios by sex
* ----------------------------------------------------------------------    
svy: stcox   i.sex##c.wfi_c10  c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
              base

margins sex, expression(exp(predict(xb)))

margins r.sex, expression(exp(predict(xb)))
nlcom pct_lower: 100*(1 - exp(_b[2.sex]))

margins, dydx(wfi_c10)

nlcom HR_overall: exp(_b[wfi_c10])

* Men (reference): HR per 0.1 FI
lincom wfi_c10, eform

* Women: HR per 0.1 FI = exp( b[wfi_c10] + b[2.sex#c.wfi_c10] )
lincom wfi_c10 + 2.sex#c.wfi_c10, eform

* ======================================================================
* 8 Cox adjusted survival models   
* =================================================================

stcox   i.sex c.wfi_c10  c.age_c c.wealth_gbp_c i.educ i.rural i.vigact ///
            i.smoke_status i.srh i.drink_status i.higov i.cogstat, ///
            tvc(c.wfi_c10) texp(time) ///
            base

* ==============================================================================
* 9 Finish
* ==============================================================================

display as text "------------------------------------------------------------"
display "==== Cox survival analysis_ COMPLETED SUCCESSFULLY ===="
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close