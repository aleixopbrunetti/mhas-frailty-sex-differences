********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: LMM post-estimation (contrasts, margins, trajectories, subgroups)
********************************************************************************

* ==============================================================================
* 3. Load Data
* ==============================================================================

use "$cleanmhasdata/wave1to5_long.dta", clear

xtset patid time

* ==============================================================================
* 4. LMM Postestimation
* ==============================================================================

* Final model (m7_hetl2) with continuous age
mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 

* -----------------------------------------------------------
* Individual fitted trajectories and smoothed average by sex
* -----------------------------------------------------------
predict fi_indiv, fitted

* Males
twoway ///
    (line fi_indiv time if sex==1, lcolor(blue%1) lwidth(vthin) connect(ascending)) ///
    (lowess fi_indiv time if sex==1, lcolor(black) lwidth(midthick) lpattern(solid)), ///
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
    (line fi_indiv time if sex==2, lcolor(red%1) lwidth(vthin) connect(ascending)) ///
    (lowess fi_indiv time if sex==2, lcolor(black) lwidth(midthick) lpattern(solid)), ///
    title("Individual FI Trajectories (Females)") ///
    ytitle("Frailty Index") ///
    xtitle("Years since baseline") ///
    ylabel(0(0.2) 0.8, nogrid) ///
    xlabel(0 2 11 14 17, nogrid) ///
    legend(off) ///
    name(fi_trajectories_female, replace)
graph export "$outputgraphs/fi_trajectories_female.jpg", replace

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

* -----------------------------------------------------------------------------------
* FI trajectories of 60 random individuals per panel
* ----------------------------------------------------------------------------------

sort patid time
local perpanel = 60     
capture drop panel
gen byte panel = .

* -- Sample MALES -> panels 1 & 2
preserve
    keep if male==1
    bysort patid: keep if _n==1               // one row per id
    set seed 62026
    gen ur = runiform()
    sort ur
    count
    local N = r(N)
    local K = cond(`N' < 2*`perpanel', `N', 2*`perpanel')
    keep in 1/`K'
    gen byte panel_m = cond(_n<=`perpanel', 1, 2)
    keep patid panel_m
    tempfile maleids
    save `maleids'
restore
merge m:1 patid using `maleids', nogen
replace panel = panel_m if panel_m < .
drop panel_m

* -- Sample FEMALES -> panels 3 & 4 
preserve
    keep if female==1
    bysort patid: keep if _n==1
    set seed 62026
    gen ur = runiform()
    sort ur
    count
    local N = r(N)
    local K = cond(`N' < 2*`perpanel', `N', 2*`perpanel')
    keep in 1/`K'
    gen byte panel_f = cond(_n<=`perpanel', 3, 4)
    keep patid panel_f
    tempfile femids
    save `femids'
restore
merge m:1 patid using `femids', nogen
replace panel = panel_f if panel_f < .
drop panel_f
keep if panel < .

* Four subplots
xtset patid time

xtline fi_indiv if panel==1, overlay legend(off) ///
    tline(0, lpattern(shortdash) lcolor(black)) ///
    xlabel(0 2 11 14 17, nogrid) ///
    ylabel(0(0.2) 0.8, nogrid) ///
    ytitle("Frailty Index") xtitle("Years since baseline") ///
    title("Male") name(g1, replace)

xtline fi_indiv if panel==2, overlay legend(off) ///
    tline(0, lpattern(shortdash) lcolor(black)) ///
    xlabel(0 2 11 14 17, nogrid) ///
    ylabel(0(0.2) 0.8, nogrid) ///
    ytitle("Frailty Index") xtitle("Years since baseline") ///
    title("Male") name(g2, replace)

xtline fi_indiv if panel==3, overlay legend(off) ///
    tline(0, lpattern(shortdash) lcolor(black)) ///
    xlabel(0 2 11 14 17, nogrid) ///
    ylabel(0(0.2) 0.8, nogrid) ///
    ytitle("Frailty Index") xtitle("Years since baseline") ///
    title("Female") name(g3, replace)

xtline fi_indiv if panel==4, overlay legend(off) ///
    tline(0, lpattern(shortdash) lcolor(black)) ///
    xlabel(0 2 11 14 17, nogrid) ///
    ylabel(0(0.2) 0.8, nogrid) ///
    ytitle("Frailty Index") xtitle("Years since baseline") ///
    title("Female") name(g4, replace)

* Combine panels
graph combine g1 g2 g3 g4, cols(2) imargin(vsmall) ///
    title("Fitted FI trajectories by sex") ///
    name(spaghetti4, replace)
capture graph export "$outputgraphs/fi_spaghetti_4panels.png", replace

* ===================================================================================
* 4 Margins data extraction  
* ===================================================================================
use "$cleanmhasdata/wave1to5_long.dta", clear

xtset patid time

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 

* -- Predicted FI by sex at each time point
margins sex, vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_sex_wave1to5.dta", replace)
marginsplot, title("Predicted FI trajectories by sex") ///
    ylabel(, nogrid) xlabel(,nogrid) ///
    name(margins_sex, replace)
graph export "$outputgraphs/margins_sex.jpg", replace

* -- Female–male difference at each time point
margins, dydx(sex) at(time=(0 2 11 14 17)) vce(unconditional) ///
    saving("$working/margins_dydx_sex_wave1to5.dta", replace)
marginsplot, title("Marginal effect of time on FI by sex") name(margins_time_sex, replace)

* -- Average over time by sex
margins sex, vce(unconditional) ///
    saving("$working/margins_sex_average.dta", replace)
marginsplot, title("Predicted FI trajectories by sex (Average)") ///
    ylabel(, nogrid) xlabel(,nogrid) ///
    name(margins_sex_average, replace)

margins, dydx(sex) vce(unconditional) ///
    saving("$working/margins_dydx_sex_average.dta", replace)

/* Not used for the analysis, but kept for reference:
* -- Average marginal effects of covariates
margins sex#educ sex#rural sex#vigact ///
    sex#smoke_status sex#srh sex#drink_status sex#higov sex#cogstat, ///
    vce(unconditional) ///
    saving("$working/margins_overall_average.dta", replace)

margins, dydx(sex) over(i.educ) vce(unconditional) ///
    saving("$working/margins_dydx_educ.dta", replace)
margins, dydx(sex) over(i.rural) vce(unconditional) ///
    saving("$working/margins_dydx_rural.dta", replace)
margins, dydx(sex) over(i.vigact) vce(unconditional) ///
    saving("$working/margins_dydx_vigact.dta", replace)
margins, dydx(sex) over(i.smoke_status) vce(unconditional) ///
    saving("$working/margins_dydx_smoke_status.dta", replace)
margins, dydx(sex) over(i.srh) vce(unconditional) ///
    saving("$working/margins_dydx_srh.dta", replace)
margins, dydx(sex) over(i.drink_status) vce(unconditional) ///
    saving("$working/margins_dydx_drink_status.dta", replace)
margins, dydx(sex) over(i.higov) vce(unconditional) ///
    saving("$working/margins_dydx_higov.dta", replace)
margins, dydx(sex) over(i.cogstat) vce(unconditional) ///
    saving("$working/margins_dydx_cogstat.dta", replace)

* -- at all waves marginal effects of covariates:
margins sex#educ sex#rural sex#vigact sex#smoke_status ///
    sex#srh sex#drink_status sex#higov sex#cogstat, ///
    at(time=(0 2 11 14 17)) vce(unconditional) ///
    saving("$working/margins_wave1to5.dta", replace)

margins, dydx(sex) over(i.educ) vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_dydx_educ_wave1to5.dta", replace)
margins, dydx(sex) over(i.rural) vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_dydx_rural_wave1to5.dta", replace)
margins, dydx(sex) over(i.vigact) vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_dydx_vigact_wave1to5.dta", replace)
margins, dydx(sex) over(i.smoke_status) vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_dydx_smoke_status_wave1to5.dta", replace)
margins, dydx(sex) over(i.srh) vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_dydx_srh_wave1to5.dta", replace)
margins, dydx(sex) over(i.drink_status) at(time=(0 2 11 14 17)) vce(unconditional) ///
    saving("$working/margins_dydx_drink_status_wave1to5.dta", replace)
margins, dydx(sex) over(i.higov) vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_dydx_higov_wave1to5.dta", replace)
margins, dydx(sex) over(i.cogstat) at(time=(0 2 11 14 17)) vce(unconditional) ///
    saving("$working/margins_dydx_cogstat_wave1to5.dta", replace)


* -- Margins main effect average and wave1to5 female vs males
preserve
    use "$working/margins_sex_wave1to5.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("FI sex time") firstrow(variables) replace

    use "$working/margins_dydx_sex_wave1to5.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("Dydx FI sexdiff 1to5") firstrow(variables) sheetreplace

    use "$working/margins_sex_average.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("FI sex average") firstrow(variables) sheetreplace

    use "$working/margins_dydx_sex_average.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("Dydx FI sexdiff average") firstrow(variables) sheetreplace

restore 

* remove main effect male vs female average and wave 1to5
capture erase "$working/margins_sex_wave1to5.dta"
capture erase "$working/margins_dydx_sex_wave1to5.dta"
capture erase "$working/margins_sex_average.dta"
capture erase "$working/margins_dydx_sex_average.dta"

*-- at average savings
preserve
    use "$working/margins_overall_average.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("fi_sex_overall") firstrow(variables) replace

    use "$working/margins_dydx_educ.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("Dydx FI sexdiff educ") firstrow(variables) sheetreplace

    use "$working/margins_dydx_rural.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("Dydx FI sexdiff rural") firstrow(variables) sheetreplace

    use "$working/margins_dydx_vigact.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("Dydx FI sexdiff vigact") firstrow(variables) sheetreplace

    use "$working/margins_dydx_smoke_status.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("Dydx FI sexdiff smoke_status") firstrow(variables) sheetreplace

    use "$working/margins_dydx_srh.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("Dydx FI sexdiff srh") firstrow(variables) sheetreplace

    use "$working/margins_dydx_drink_status.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("Dydx FI sexdiff drink_status") firstrow(variables) sheetreplace

    use "$working/margins_dydx_higov.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("Dydx FI sexdiff higov") firstrow(variables) sheetreplace

    use "$working/margins_dydx_cogstat.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("Dydx FI sexdiff cogstat") firstrow(variables) sheetreplace
restore

*/
* Remove temporary files (average)
capture erase "$working/margins_overall_average.dta"
capture erase "$working/margins_dydx_educ.dta"
capture erase "$working/margins_dydx_rural.dta"
capture erase "$working/margins_dydx_vigact.dta"
capture erase "$working/margins_dydx_smoke_status.dta"
capture erase "$working/margins_dydx_srh.dta"
capture erase "$working/margins_dydx_drink_status.dta"
capture erase "$working/margins_dydx_higov.dta"
capture erase "$working/margins_dydx_cogstat.dta"

/*
* -- at all waves savings

preserve
    use "$working/margins_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("fi_sex_wave1to5") firstrow(variables) replace

    use "$working/margins_dydx_educ_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("Dydx FI educ 1to5") firstrow(variables) sheetreplace
    
    use "$working/margins_dydx_rural_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("Dydx FI rural 1to5") firstrow(variables) sheetreplace

    use "$working/margins_dydx_vigact_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("Dydx FI vigact 1to5") firstrow(variables) sheetreplace

    use "$working/margins_dydx_smoke_status_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("Dydx FI smoke_status 1to5") firstrow(variables) sheetreplace

    use "$working/margins_dydx_srh_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("Dydx FI srh 1to5") firstrow(variables) sheetreplace

    use "$working/margins_dydx_drink_status_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("Dydx FI drink_status 1to5") firstrow(variables) sheetreplace

    use "$working/margins_dydx_higov_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("Dydx FI higov 1to5") firstrow(variables) sheetreplace

    use "$working/margins_dydx_cogstat_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("Dydx FI cogstat 1to5") firstrow(variables) sheetreplace
restore
*/

* Remove temporary files (wave 1 to 5)
capture erase "$working/margins_wave1to5.dta"
capture erase "$working/margins_dydx_educ_wave1to5.dta"
capture erase "$working/margins_dydx_rural_wave1to5.dta"
capture erase "$working/margins_dydx_vigact_wave1to5.dta"
capture erase "$working/margins_dydx_smoke_status_wave1to5.dta"
capture erase "$working/margins_dydx_srh_wave1to5.dta"
capture erase "$working/margins_dydx_drink_status_wave1to5.dta"
capture erase "$working/margins_dydx_higov_wave1to5.dta"
capture erase "$working/margins_dydx_cogstat_wave1to5.dta"

* =====================================================================================
* 5. Marginal effects by wealth quintile
* ======================================================================================
mixed wfi i.sex##c.time c.age_c##c.time i.wealth_q ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* -- Average
margins sex#wealth_q, vce(unconditional) ///
    saving("$working/margins_sex_wealth_q.dta", replace)

margins, dydx(sex) over(wealth_q) vce(unconditional) ///
    saving("$working/margins_dydx_wealth_q.dta", replace)

* -- Wave 1 to 5
margins sex#wealth_q, vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_sex_wealth_q_wave1to5.dta", replace)

marginsplot, title("Predicted FI trajectories by sex and wealth quintile") ///
    ylabel(, nogrid) xlabel(,nogrid) ///
    name(margins_sex_wealth_q_time, replace)
graph export "$outputgraphs/margins_sex_wealth_q_time.jpg", replace

margins, dydx(sex) over(i.wealth_q) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_dydx_sex_wealth_q_wave1to5.dta", replace)

preserve
    use "$working/margins_sex_wealth_q.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("sexdiff_wealth_q") firstrow(variables) sheetreplace
    
    use "$working/margins_dydx_wealth_q.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("sexdiff_dydx_wealth_q") firstrow(variables) sheetreplace

    use "$working/margins_sex_wealth_q_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("sexdiff_wealth_q_wave1to5") firstrow(variables) sheetreplace

    use "$working/margins_dydx_sex_wealth_q_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("sexdiff_dydx_wealth_q_wave1to5") firstrow(variables) sheetreplace
restore

capture erase "$working/margins_sex_wealth_q.dta"
capture erase "$working/margins_dydx_wealth_q.dta"
capture erase "$working/margins_sex_wealth_q_wave1to5.dta"
capture erase "$working/margins_dydx_sex_wealth_q_wave1to5.dta"

* =======================================================================================
* 6 Marginal effects of the model with age categories
* =======================================================================================

mixed wfi i.sex##c.time i.age_cat##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
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

margins, dydx(sex) over(age_cat) at(time=(0 2 11 14 17)) ///
    vce(unconditional) ///
    saving("$working/margins_dydx_agecat.dta", replace)

* Age categories - average over time
margins sex#age_cat, vce(unconditional) ///
    saving("$working/margins_sex_age_cat.dta", replace)

margins, dydx(sex) over(age_cat)  ///
    vce(unconditional) ///
    saving("$working/margins_dydx_agecat.dta", replace)

preserve
    use "$working/margins_sex_age_cat_wave1to5.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("sexdiff_sex_agecat") firstrow(variables) sheetreplace
    
    use "$working/margins_sex_age_cat.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("sexdiff_age_cat") firstrow(variables) sheetreplace

    use "$working/margins_dydx_agecat.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("sexdiff_dydx_age_cat") firstrow(variables) sheetreplace
restore

capture erase "$working/margins_sex_age_cat_wave1to5.dta"
capture erase "$working/margins_dydx_agecat_wave1to5.dta"
capture erase "$working/margins_sex_age_cat.dta"
capture erase "$working/margins_dydx_agecat.dta"
capture erase "$working/margins_dydx_sex_alive_average.dta"
capture erase "$working/margins_dydx_sex_average.dta"
capture erase "$working/margins_dydx_sex_wave1to5.dta"
capture erase "$working/margins_sex_average.dta"
capture erase "$working/margins_sex_wave1to5.dta"

* ===========================================================================================
* 7 Trajectory model with BMI subsample
* ===========================================================================================

mixed wfi i.sex##c.time i.age_cat##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time i.bmicat##c.time if !(bmi ==.s) ///
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* Average marginal effects
margins sex#bmicat, ///
    vce(unconditional) ///
    saving("$working/margins_sex_bmicat.dta", replace)

margins, dydx(sex) over(bmicat) ///
    vce(unconditional) ///
    saving("$working/margins_dydx_sex_bmicat.dta", replace)

* Wave 1 to 5
margins sex#bmicat, at(time=(0 2 11 14 17)) vce(unconditional) ///
    saving("$working/margins_sex_bmicat_wave1to5.dta", replace)

marginsplot, xdimension(time) plotdimension(bmicat sex) recast(line) noci ///
    plot1opts(lcolor(green)   lpattern(longdash) lwidth(medthin)) ///
    plot2opts(lcolor(green)   lpattern(solid)) ///
    plot3opts(lcolor(blue)    lpattern(longdash) lwidth(medthin)) ///
    plot4opts(lcolor(blue)    lpattern(solid)) ///
    plot5opts(lcolor(purple)  lpattern(longdash) lwidth(medthin)) ///
    plot6opts(lcolor(purple)  lpattern(solid)) ///
    plot7opts(lcolor(orange)  lpattern(longdash) lwidth(medthin)) ///
    plot8opts(lcolor(orange)  lpattern(solid)) ///
    plot9opts(lcolor(red)     lpattern(longdash) lwidth(medthin)) ///
    plot10opts(lcolor(red)    lpattern(solid)) ///
    plot11opts(lcolor(brown)  lpattern(longdash) lwidth(medthin)) ///
    plot12opts(lcolor(brown)  lpattern(solid)) ///
    legend(order(1 "Male, Underweight" 2 "Female, Underweight" ///
                 3 "Male, Normal" 4 "Female, Normal" ///
                 5 "Male, Pre-obesity" 6 "Female, Pre-obesity" ///
                 7 "Male, Obesity I" 8 "Female, Obesity I" ///
                 9 "Male, Obesity II" 10 "Female, Obesity II" ///
                 11 "Male, Obesity III" 12 "Female, Obesity III") ///
           cols(2) region(lstyle(none)) ///
           position(6)) ///
    ytitle("Frailty Index") xtitle("Time since baseline (years)") ///
    ylabel(, nogrid) xlabel(, nogrid) ///
    title("Predicted FI Trajectories by Sex and BMI Category") ///
    name(bmi_band_combined, replace)
graph export "$outputgraphs/fi_trajectories_BMI_bands_by_sex.jpg", replace

margins, dydx(sex) over(bmicat) at(time=(0 2 11 14 17)) ///
    vce(unconditional) ///
    saving("$working/margins_dydx_sex_bmicat_time.dta", replace)

preserve
    use "$working/margins_sex_bmicat.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("sexdiff_bmicat") firstrow(variables) sheetreplace

    use "$working/margins_dydx_sex_bmicat.dta", clear
    export excel using "$outputtables/margins_average_bysex.xlsx", ///
        sheet("sexdiff_dydx_bmicat") firstrow(variables) sheetreplace

    use "$working/margins_sex_bmicat_wave1to5.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("sexdiff_sex_wave1to5") firstrow(variables) sheetreplace

    use "$working/margins_dydx_sex_bmicat_time.dta", clear
    export excel using "$outputtables/margins_wave1to5_bysex.xlsx", ///
        sheet("sexdiff_dydx_bmicat") firstrow(variables) sheetreplace
restore

* Remove temporary files
capture erase "$working/margins_sex_bmicat.dta"
capture erase "$working/margins_dydx_sex_bmicat.dta"
capture erase "$working/margins_sex_bmicat_wave1to5.dta"
capture erase "$working/margins_dydx_sex_bmicat_time.dta"

* ===========================================================================================
* 8 Trajectory model among participants who died vs alive until the end of the study
* ============================================================================================
* Note: descriptive comparison of FI trajectories by vital status (not a sensitivity analysis)

*--------------------------------------------------------------------------------------------
* Model among participants alive until the end of the study (excluding dropouts)
*--------------------------------------------------------------------------------------------

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time if status == 1 /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 

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

* -- Female–male difference at each time point
margins, dydx(sex) at(time=(0 2 11 14 17)) vce(unconditional) ///
    saving("$working/margins_dydx_sex_alive_wave1to5.dta", replace)
marginsplot, title("Marginal effect of time on FI by sex among alive") ///
    name(margins_time_sex_alive, replace)

* -- Average over time by sex
margins sex, vce(unconditional) ///
    saving("$working/margins_sex_alive_average.dta", replace)
marginsplot, title("Predicted FI trajectories by sex among alive (Average)") ///
    ylabel(, nogrid) xlabel(,nogrid) ///
    name(margins_sex_alive_average, replace)

margins, dydx(sex) vce(unconditional) ///
    saving("$working/margins_dydx_sex_alive_average.dta", replace)

* -- Export (alive)
preserve
    use "$working/margins_sex_alive_wave1to5.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("FI_sex_alive") firstrow(variables) sheetreplace

    use "$working/margins_dydx_sex_alive_wave1to5.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("Dydx_FI_1to5_alive") firstrow(variables) sheetreplace

    use "$working/margins_sex_alive_average.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("FI_sex_alive") firstrow(variables) sheetreplace

    use "$working/margins_dydx_sex_alive_average.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("Dydx_FI_average_alive") firstrow(variables) sheetreplace

restore

* Remove temporary files
capture erase "$working/margins_sex_alive_wave1to5.dta"
capture erase "$working/margins_dydx_sex_alive_wave1to5.dta"
capture erase "$working/margins_sex_alive_average.dta"
capture erase "$working/margins_dydx_sex_alive_average.dta"

*-------------------------------------------------------------------------------------------------
* Model among participants who died
*-------------------------------------------------------------------------------------------------

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time if death == 1 /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base 

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

* -- Predicted FI by sex at each time point (died)
margins sex, vce(unconditional) at(time=(0 2 11 14 17)) ///
    saving("$working/margins_sex_dead_wave1to5.dta", replace)

marginsplot, title("Predicted FI trajectories by sex among dead") ///
    ylabel(, nogrid) xlabel(,nogrid) ///
    name(margins_sex_dead, replace)
graph export "$outputgraphs/margins_sex_dead.jpg", replace

* -- Female–male difference at each time point
margins, dydx(sex) at(time=(0 2 11 14 17)) vce(unconditional) ///
    saving("$working/margins_dydx_sex_dead_wave1to5.dta", replace)
marginsplot, title("Marginal effect of time on FI by sex among dead") ///
    name(margins_time_sex_dead, replace)

* -- Average over time by sex
margins sex, vce(unconditional) ///
    saving("$working/margins_sex_dead_average.dta", replace)
marginsplot, title("Predicted FI trajectories by sex among dead (Average)") ///
    ylabel(, nogrid) xlabel(,nogrid) ///
    name(margins_sex_dead_average, replace)

margins, dydx(sex) vce(unconditional) ///
    saving("$working/margins_dydx_sex_dead_average.dta", replace)

* -- Export (died)
preserve
    use "$working/margins_sex_dead_wave1to5.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("FI_sex_time_dead") firstrow(variables) sheetreplace

    use "$working/margins_dydx_sex_dead_wave1to5.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("Dydx_FI_1to5_dead") firstrow(variables) sheetreplace

    use "$working/margins_sex_dead_average.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("FI_sex_dead") firstrow(variables) sheetreplace

    use "$working/margins_dydx_sex_dead_average.dta", clear
    export excel using "$outputtables/margins_bysex.xlsx", ///
        sheet("Dydx_FI_sexdiff_dead") firstrow(variables) sheetreplace

restore

* Remove temporary files
capture erase "$working/margins_sex_dead_wave1to5.dta"
capture erase "$working/margins_dydx_sex_dead_wave1to5.dta"
capture erase "$working/margins_sex_dead_average.dta"
capture erase "$working/margins_dydx_sex_dead_average.dta"

* ==============================================================================
* 9 Distribution Check and correlation matrix
* ==============================================================================
use "$cleanmhasdata/wave1to5_long.dta", clear

xtset patid time

mixed wfi i.sex##c.time c.age_c##c.time c.wealth_gbp_c ///
     i.educ i.rural i.vigact i.smoke_status i.srh ///
     i.drink_status##c.time i.higov i.cogstat##c.time /// 
     || patid: female time_female, noconstant covariance(unstructured) ///
     || patid: male time_male, noconstant covariance(unstructured) ///
     ml nolog stddev pweight(w1respwt) base

* -------------------------------------------------------------------------------
* Correlation Matrix
* -------------------------------------------------------------------------------

* Correlation of parameter estimates (multicollinearity check)
estat vce, correlation

*-- Pairwise correlations of FI across waves
preserve

use "$cleanmhasdata/wave1to5_wide_deaths_code.dta", clear
pwcorr w*fi, obs
pwcorr w*fi if sex == 1, obs // Males
pwcorr w*fi if sex == 2, obs // Females

restore

* ------------------------------------------------------------------
* Checking for gamma distribution
* ------------------------------------------------------------------

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

restore

* ==============================================================================
* 10 Finish
* ==============================================================================

display as text "------------------------------------------------------------"
display "==== MHAS_LMMs_Post-estimation_COMPLETED  ===="
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close