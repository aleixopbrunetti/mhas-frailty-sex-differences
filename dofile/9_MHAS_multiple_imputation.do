********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Multiple Imputation for MHAS Frailty Index Study
********************************************************************************

* ==============================================================================
* 3. Load Data
* ==============================================================================

use "$cleanmhasdata/wave1to5_long.dta", clear

set seed 20250920

* ==============================================================================
* 4 Multiple Imputation for missing data
* ==============================================================================
*---------------------------
* 1) Set MI structure
*---------------------------
mi set wide
mi svyset [pweight=w1respwt]
mi register regular patid sex time female male ///
     time_female time_male w1respwt age_c rural
mi register imputed wealth_gbp educ vigact ///
    smoke_status srh drink_status higov cogstat wfi
mi register passive  ///
    wealth_gbp_c

*------------------------------------------------------------
* 2) Multiple imputation by chained equations (MICE)
*---------------------------------------------------------------
// Truncated regression for FI (bounded 0–1); PMM for wealth
// mlogit for nominal; ologit for ordered; logit for binary variables

* Number of imputations and burn-in
local Mtarget 20
local Burnin 30

mi impute chained /// 
    (truncreg, ll(0) ul(1)) wfi ///
    (pmm, knn(10))           wealth_gbp /// 
    (mlogit, augment)       smoke_status cogstat educ ///
    (ologit, augment)       srh ///
    (logit, augment)        vigact drink_status higov ///
    = i.sex c.time c.age_c i.rural [pweight=w1respwt], ///
    add(`Mtarget') force burnin(`Burnin') savetrace("$working/mitrace_MI.dta", replace) ///
    dots noisily

* Passive variable: wealth centred on observed mean (same constant in all imputations)
mi xeq 0: quietly summarize wealth_gbp
scalar WEALTH_MEAN = r(mean)
mi passive: replace wealth_gbp_c = wealth_gbp - WEALTH_MEAN

*----------------------------------------------------------
* 3) Chain diagnostics
*----------------------------------------------------------

* -------------------------------
* FCS trace: mean of wfi by imputation, with observed mean as reference
* -------------------------------
mi xeq 0: summarize wfi
scalar wfi_obs = r(mean)

preserve
    * Load trace file
    use "$working/mitrace_MI.dta", clear
    isid m iter, sort
    xtset m iter

    * Trace of mean(wfi) for each imputation
    twoway line wfi_mean iter, ///
        by(m, cols(3) note("") legend(off) title("Trace: mean(wfi) across iterations")) ///
        ylabel(, nogrid) xlabel(, nogrid) ///
        ytitle("mean(wfi)") xtitle("Iteration") ///
        yline(`=wfi_obs', lpattern(dash)) ///
        name(tr_wfi_by, replace)
    graph export "$outputgraphs/mi_trace_wfi_by.jpg", replace
restore 

/* Not used: autocorrelation plots are not standard for MICE diagnostics
* ------------------------------------------
* 4) ACF trace
* -----------------------------------------

preserve
use "$working/mitrace_MI.dta", clear
isid m iter, sort
xtset m iter
levelsof m, local(mlist)

foreach mm of local mlist {
    keep if m==`mm'
    tsset iter
    quietly su iter, meanonly
    local T = r(max) - r(min) + 1
    local L = max(1, `T' - 1)

    ac wfi_mean, lags(`L') ///
        xlabel(, nogrid) ylabel(, nogrid) ///
        title("ACF: mean(wfi) — m=`mm'") ///
        name(ac_wfi_`mm', replace)
    graph export "$outputgraphs/mi_trace_wfi_acf_m`mm'.jpg", replace

    restore, preserve 
}

restore
*/

* --------------------------
* Overlay trace plots (all imputations)
*---------------------------

* Mean
preserve 
    use "$working/mitrace_MI.dta", clear
    isid m iter, sort
    xtset m iter
    levelsof m, local(mlist)

    local layers
    foreach mm of local mlist {
        local layers `"`layers' (line wfi_mean iter if m==`mm', sort)"'
    }
    twoway `layers', ///
        title("Trace: mean(wfi), all imputations overlaid") ///
        ytitle("mean(wfi)") xtitle("Iteration") legend(off) ///
		xlabel(,nogrid) ylabel(,nogrid) ///
        name(tr_wfi_overlay, replace)
    graph export "$outputgraphs/mi_trace_wfi_mean_overlay.jpg", replace
restore 

* Standard deviation
preserve
    use "$working/mitrace_MI.dta", clear
    isid m iter, sort
    levelsof m, local(mlist)

    local layers_sd
    foreach mm of local mlist {
        local layers_sd `"`layers_sd' (line wfi_sd iter if m==`mm', sort)"'
    }
    twoway `layers_sd', ///
        title("Trace: sd(wfi), all imputations overlaid") ///
        ytitle("sd(wfi)") xtitle("Iteration") legend(off) ///
		xlabel(,nogrid) ylabel(,nogrid) ///
        name(tr_wfi_sd_overlay, replace)
    graph export "$outputgraphs/mi_trace_wfi_sd_overlay.jpg", replace
restore

*---------------------------------------------------
* 5) Diagnostics with midiagplots (observed vs imputed)
*---------------------------------------------------

* Kernel density overlay for FI
mi describe
local M = r(M) 
midiagplots wfi, m(1/`M') combine ///
    xlabel(, nogrid) ylabel(, nogrid)
graph export "$outputgraphs/mi_midiagplots_wfi_kdens.jpg", replace

* Cumulative distribution
midiagplots wfi, combine plottype(cumul) ///
	 xlabel(, nogrid) ylabel (, nogrid)
graph export "$outputgraphs/mi_midiagplots_wfi_cumul.jpg", replace

* Two-sample Kolmogorov–Smirnov test (observed vs imputed distributions)
midiagplots wfi, ksmirnov
midiagplots wfi, m(3) ksmirnov 

mi describe
local M = r(M) 
midiagplots wfi, m(1/`M') ksmirnov combine plottype(cumul) ///
	 xlabel(, nogrid) ylabel (, nogrid)

* Kernel density: observed vs all imputations
mi describe
local M = r(M)

preserve
    keep wfi _*_wfi w1respwt
    local layers "(kdensity wfi, lcolor(black) lwidth(thin))"
    forvalues m = 1/`M' {
        local layers `"`layers' (kdensity _`m'_wfi, lwidth(thin))"'
    }

    twoway `layers', ///
        title("Observed vs ALL imputations — wfi (kernel density)") ///
        ytitle("Density") xtitle("wfi") ///
        legend(off) note("Observed = black; Imputations m=1..`M'") ///
		xlabel(,nogrid) ylabel(,nogrid) ///
        name(wfi_kdens_all, replace)
    graph export "$outputgraphs/mi_overlay_wfi_kdensity_all.jpg", replace
restore

*===============================================================================
* 5) Final mixed model on imputed data
*===============================================================================

capture noisily mi estimate, mcerror vartable esampvaryok errorok dots: ///
    mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
          i.educ i.rural i.vigact i.smoke_status i.srh ///
          i.drink_status##c.time i.higov i.cogstat##c.time ///
          || patid: female time_female, noconstant covariance(unstructured) ///
          || patid: male   time_male, noconstant covariance(unstructured) ///
          ml pw(w1respwt) nolog stddev difficult iterate(200)

estimates store mi_model

*------------------------------------------------------------------------------
* Complete-case analysis for comparison
*------------------------------------------------------------------------------
preserve

mi xeq 0: mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
    i.educ i.rural i.vigact i.smoke_status i.srh ///
    i.drink_status##c.time i.higov i.cogstat##c.time ///
    || patid: female time_female, noconstant covariance(unstructured) ///
    || patid: male   time_male,   noconstant covariance(unstructured) ///
    ml pw(w1respwt) nolog stddev difficult iterate(200) vce(robust)
estimates store cc_model
estat ic 

restore 

* Comparison table: complete case vs multiple imputation
esttab cc_model mi_model, ///
    ci b(%9.4f) ci(%9.4f) keep(`keep') nomtitles label nonotes compress ///
    coeflabels( ///
        2.sex                      "Female (vs male)" ///
        time                       "Time (years)" ///
        2.sex#c.time               "Female × Time" ///
        age_c                      "Age (centered)" ///
        c.age_c#c.time             "Age × Time" ///
        wealth_gbp_c               "Wealth (centered)" ///
        2.educ                     "Education: Upper sec." ///
        3.educ                     "Education: Tertiary" ///
        1.rural                    "Rural (yes)" ///
        1.vigact                   "Vigorous activity (yes)" ///
        1.smoke_status             "Former smoker" ///
        2.smoke_status             "Current smoker" ///
        2.srh                      "SRH: Very good" ///
        3.srh                      "SRH: Good" ///
        4.srh                      "SRH: Fair" ///
        5.srh                      "SRH: Poor" ///
        1.drink_status             "Current drinker" ///
        1.drink_status#c.time      "Current drinker × Time" ///
        1.higov                    "Gov't insurance (yes)" ///
        1.cogstat                  "Instr. impairment" ///
        2.cogstat                  "CIND" ///
        3.cogstat                  "Dementia" ///
        1.cogstat#c.time           "Instr. impairment × Time" ///
        2.cogstat#c.time           "CIND × Time" ///
        3.cogstat#c.time           "Dementia × Time" ///
        _cons                      "Constant" ///
    ) ///
    title("Mixed-effects regression of FI: Complete case vs Multiple imputation") ///
    stats(N ll df_m, fmt(%9.0g %9.0g %9.0g) labels("N" "LogLik" "df(model)"))

*===============================================================================
* 6) Session end
*===============================================================================
display as text "------------------------------------------------------------"
display as text "MHAS_MI_completed"
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close