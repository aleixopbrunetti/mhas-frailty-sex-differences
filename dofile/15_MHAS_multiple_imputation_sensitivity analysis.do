********************************************************************************
* Project: MHAS Frailty Index Study
* Author:  Aleixo Brunetti
* Purpose: MI mixed models - orthogonal polynomial time and model without SRH
********************************************************************************

*===============================================================================
* 3 Refit MI mixed model with orthogonal polynomials in time 
*===============================================================================

use "$cleanmhasdata/wave1to5_long_MI.dta", clear

*-------------------------------------------------------------------------------
* 3.1 Orthogonal polynomial terms for time
*-------------------------------------------------------------------------------

* orthpoly built on the full long MI dataset (time is not imputed)
orthpoly time, degree(3) generate(otime1 otime2 otime3)

* Random slopes on raw linear time (as in primary model)
capture drop time_female time_male
gen double time_female = female * time
gen double time_male   = male   * time

mi update

* Checks: orthogonality and consistency across imputations
mi xeq 1: correlate otime1 otime2 otime3    
mi xeq 1: regress otime3 otime1 otime2      
mi xeq 1:  summarize otime1 otime2 otime3
mi xeq 40: summarize otime1 otime2 otime3

save "$cleanmhasdata/wave1to5_long_MI.dta", replace

*-------------------------------------------------------------------------------
* 3.3 Mixed model on the imputed datasets (orthogonal time)
*-------------------------------------------------------------------------------
capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots ///
    saving(miests_orth, replace): ///
    mixed wfi i.sex c.otime1 c.otime2 c.otime3 ///
        i.sex#c.otime1 c.age_c c.age_c#c.otime1 c.wealth_gbp_c ///
        i.educ i.rural i.vigact i.smoke_status i.srh ///
        i.drink_status i.drink_status#c.otime1 i.higov ///
        i.cogstat i.cogstat#c.otime1 ///
        || patid: female time_female, noconstant covariance(unstructured) ///
        || patid: male   time_male,   noconstant covariance(unstructured) ///
        ml pw(w1respwt) nolog stddev iterate(300) vce(robust)

estimates store mi_model_orth
estat ic 
display "AIC: " %20.0f -2*e(ll) + 2*e(rank)
display "BIC: " %20.0f -2*e(ll) + e(rank)*ln(e(N))

* Female–male FI difference at baseline (time = 0)
preserve 
sum otime1 if time==0, meanonly
local v0 = r(mean)
lincom 2.sex + `v0'*2.sex#c.otime1
restore 

* Female–male FI difference at each wave
levelsof time, local(times)
foreach t of local times {
    quietly sum otime1 if time==`t', meanonly
    di _n "=== Female-male FI gap at time `t' (otime1=" r(mean) ") ==="
    lincom 2.sex + `=r(mean)'*2.sex#c.otime1
}

*===============================================================================
* 4 Refit MI mixed model without SRH
*===============================================================================

use "$cleanmhasdata/wave1to5_long_MI.dta", clear

capture noisily mi estimate, post mcerror vartable esampvaryok errorok dots saving(miests, replace): ///
		mixed wfi i.sex c.time c.time#c.time c.time#c.time#c.time ///
		i.sex#c.time c.age_c c.age_c#c.time c.wealth_gbp_c ///
		i.educ i.rural i.vigact i.smoke_status ///
		i.drink_status i.drink_status#c.time i.higov ///
		i.cogstat i.cogstat#c.time ///
		|| patid: female time_female, noconstant covariance(unstructured) ///
		|| patid: male time_male, noconstant covariance(unstructured) ///
		ml pw(w1respwt) nolog stddev iterate(300) vce(robust)
		
estimates store mi_model_nosrh

*===============================================================================
* 7. Session end
*===============================================================================
display as text "------------------------------------------------------------"
display as text "MHAS_MI_completed"
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close