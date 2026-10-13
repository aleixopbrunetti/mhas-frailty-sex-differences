********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Missing Data Analysis for MHAS Frailty Index
********************************************************************************

* ==============================================================================
* 3. Check the missing data of those who have ONLY baseline data
* ==============================================================================

* Load the cleaned MHAS data
use "$cleanmhasdata/wave1to5_wide.dta", clear

tabstat w1fi w2fi w3fi w4fi w5fi, stats(mean sd p25 p50 p75 min max count) longstub columns(statistics)
misstable sum w*fi

* -- General overview of FI missing data across all waves
misstable patterns w*fi
misstable sum w*fi if dyear== .x
misstable sum w*fi if !missing(dyear)

* ---------------------------------------------------------------
* 3.1 Weighted Summary statistics for 1to5 waves
* ---------------------------------------------------------------

svyset [pweight=w1respwt]

foreach wave of varlist w1fi w2fi w3fi w4fi w5fi {
    display as text "=== `wave' ==="
    svy: mean `wave'
    estat sd
}

* Set the survey design
svyset [pweight=w1respwt]

* List of FI variables
local fi_vars w1fi w2fi w3fi w4fi w5fi
local fi_vars w1fi w2fi w3fi w4fi w5fi

* Weighted FI by sex and test of sex difference
foreach var of local fi_vars {

    display as text "===== `var' in MEN (sex = 1) ====="
    svy, subpop(if sex == 1): mean `var'
    estat sd

    display as text "===== `var' in WOMEN (sex = 2) ====="
    svy, subpop(if sex == 2): mean `var'
    estat sd

    display as text "===== P-VALUE for sex difference in `var' ====="
    svy: regress `var' i.sex
    test 1.sex
}

*---------------------------------------------------------------
*3.1 Summary statistics for 1to5 waves
*---------------------------------------------------------------

preserve
tempname results
tempfile postresults

postfile `results' str20 statistic str25 wave1 str25 wave2 str25 wave3 str25 wave4 str25 wave5 using `postresults', replace

local waves "w1fi w2fi w3fi w4fi w5fi"

* Calculate statistics for all waves
forvalues i = 1/5 {
    local w : word `i' of `waves'
    quietly summarize `w'
    local mean`i' = string(r(mean), "%5.2f")
    local sd`i' = string(r(sd), "%5.2f")
    local min`i' = string(r(min), "%5.2f")
    local max`i' = string(r(max), "%5.2f")
    
    quietly count if missing(`w')
    local nmiss`i' = r(N)
    quietly count
    local total`i' = r(N)
    local pctmiss`i' = string(100 * `nmiss`i'' / `total`i'', "%5.1f")

    quietly centile `w', centile(25 50 75)
    local p25_`i' = string(r(c_1), "%5.2f")
    local p50_`i' = string(r(c_2), "%5.2f")
    local p75_`i' = string(r(c_3), "%5.2f")

    * Format 
    local mean_sd`i' `"`mean`i'' (`sd`i'')"'
    local median_iqr`i' `"`p50_`i'' (`p25_`i'',`p75_`i'')"'
    local range`i' "`min`i'', `max`i''"
    local missing`i' `"`nmiss`i'' (`pctmiss`i'')"'
}

* Post rows 
post `results' ("Mean (SD)") (`"`mean_sd1'"') (`"`mean_sd2'"') (`"`mean_sd3'"') (`"`mean_sd4'"') (`"`mean_sd5'"')
post `results' ("Median (IQR)") (`"`median_iqr1'"') (`"`median_iqr2'"') (`"`median_iqr3'"') (`"`median_iqr4'"') (`"`median_iqr5'"')
post `results' ("Range") (`"`range1'"') (`"`range2'"') (`"`range3'"') (`"`range4'"') (`"`range5'"')
post `results' ("Missing, n (%)") (`"`missing1'"') (`"`missing2'"') (`"`missing3'"') (`"`missing4'"') (`"`missing5'"')

postclose `results'

use `postresults', clear

forvalues i = 1/4 {
    local stat_val = statistic[`i']
    local w1_val = wave1[`i']
    local w2_val = wave2[`i']
    local w3_val = wave3[`i']
    local w4_val = wave4[`i']
    local w5_val = wave5[`i']
    
    display as text %16s "`stat_val'" %16s "`w1_val'" %16s "`w2_val'" %16s "`w3_val'" %16s "`w4_val'" %16s "`w5_val'"
}

display as text "{hline 96}"
display as text ""

* Dataset for export
gen row_order = _n
order statistic wave1 wave2 wave3 wave4 wave5 row_order

* Add identifier for overall results
gen sex_group = "Overall"
order sex_group statistic wave1 wave2 wave3 wave4 wave5 row_order

* Store overall results for later combination
tempfile overall_results
save `overall_results'

* Export overall table
export delimited "$outputtables/fi_missing_summary.csv", replace

restore

* -----------------------------------------------------------------------------
* 3.2 Summary table for MALES ONLY (sex == 1) waves 1 to 5
* -----------------------------------------------------------------------------

preserve
tempname results_male
tempfile postresults_male

postfile `results_male' str20 statistic str25 wave1 str25 wave2 str25 wave3 str25 wave4 str25 wave5 using `postresults_male', replace

local waves "w1fi w2fi w3fi w4fi w5fi"

* Calculate statistics for all waves for MALES only
forvalues i = 1/5 {
    local w : word `i' of `waves'
    
    quietly summarize `w' if sex == 1
    local mean`i' = string(r(mean), "%5.2f")
    local sd`i' = string(r(sd), "%5.2f")
    local min`i' = string(r(min), "%5.2f")
    local max`i' = string(r(max), "%5.2f")
    
    quietly count if missing(`w') & sex == 1
    local nmiss`i' = r(N)
    quietly count if sex == 1
    local total`i' = r(N)
    local pctmiss`i' = string(100 * `nmiss`i'' / `total`i'', "%5.1f")

    quietly centile `w' if sex == 1, centile(25 50 75)
    local p25_`i' = string(r(c_1), "%5.2f")
    local p50_`i' = string(r(c_2), "%5.2f")
    local p75_`i' = string(r(c_3), "%5.2f")

    * Format 
    local mean_sd`i' `"`mean`i'' (`sd`i'')"'
    local median_iqr`i' `"`p50_`i'' (`p25_`i'',`p75_`i'')"'
    local range`i' "`min`i'', `max`i''"
    local missing`i' `"`nmiss`i'' (`pctmiss`i'')"'
}

* Post rows for MALES
post `results_male' ("Mean (SD)") (`"`mean_sd1'"') (`"`mean_sd2'"') (`"`mean_sd3'"') (`"`mean_sd4'"') (`"`mean_sd5'"')
post `results_male' ("Median (IQR)") (`"`median_iqr1'"') (`"`median_iqr2'"') (`"`median_iqr3'"') (`"`median_iqr4'"') (`"`median_iqr5'"')
post `results_male' ("Range") (`"`range1'"') (`"`range2'"') (`"`range3'"') (`"`range4'"') (`"`range5'"')
post `results_male' ("Missing, n (%)") (`"`missing1'"') (`"`missing2'"') (`"`missing3'"') (`"`missing4'"') (`"`missing5'"')

postclose `results_male'

use `postresults_male', clear

* Table for Males
display as text %16s "Statistic" %16s "Wave 1" %16s "Wave 2" %16s "Wave 3" %16s "Wave 4" %16s "Wave 5"
display as text "{hline 96}"

forvalues i = 1/4 {
    local stat_val = statistic[`i']
    local w1_val = wave1[`i']
    local w2_val = wave2[`i']
    local w3_val = wave3[`i']
    local w4_val = wave4[`i']
    local w5_val = wave5[`i']
    
    display as text %16s "`stat_val'" %16s "`w1_val'" %16s "`w2_val'" %16s "`w3_val'" %16s "`w4_val'" %16s "`w5_val'"
}

display as text "{hline 96}"
display as text ""

* Identifier for males
gen sex_group = "Males"
gen row_order = _n
order sex_group statistic wave1 wave2 wave3 wave4 wave5 row_order

* Save males results temporarily
tempfile males_results
save `males_results'
restore

* -----------------------------------------------------------------------------------------------
* 3.3 Summary table for FEMALES ONLY (sex == 2) waves 1 to 5
* -----------------------------------------------------------------------------------------------

preserve
tempname results_female
tempfile postresults_female

postfile `results_female' str20 statistic str25 wave1 str25 wave2 str25 wave3 str25 wave4 str25 wave5 using `postresults_female', replace

local waves "w1fi w2fi w3fi w4fi w5fi"

* Calculate statistics for all waves for FEMALES only
forvalues i = 1/5 {
    local w : word `i' of `waves'
    
    quietly summarize `w' if sex == 2
    local mean`i' = string(r(mean), "%5.2f")
    local sd`i' = string(r(sd), "%5.2f")
    local min`i' = string(r(min), "%5.2f")
    local max`i' = string(r(max), "%5.2f")
    
    quietly count if missing(`w') & sex == 2
    local nmiss`i' = r(N)
    quietly count if sex == 2
    local total`i' = r(N)
    local pctmiss`i' = string(100 * `nmiss`i'' / `total`i'', "%5.1f")

    quietly centile `w' if sex == 2, centile(25 50 75)
    local p25_`i' = string(r(c_1), "%5.2f")
    local p50_`i' = string(r(c_2), "%5.2f")
    local p75_`i' = string(r(c_3), "%5.2f")

    * Format 
    local mean_sd`i' `"`mean`i'' (`sd`i'')"'
    local median_iqr`i' `"`p50_`i'' (`p25_`i'',`p75_`i'')"'
    local range`i' "`min`i'', `max`i''"
    local missing`i' `"`nmiss`i'' (`pctmiss`i'')"'
}

* Post rows 
post `results_female' ("Mean (SD)") (`"`mean_sd1'"') (`"`mean_sd2'"') (`"`mean_sd3'"') (`"`mean_sd4'"') (`"`mean_sd5'"')
post `results_female' ("Median (IQR)") (`"`median_iqr1'"') (`"`median_iqr2'"') (`"`median_iqr3'"') (`"`median_iqr4'"') (`"`median_iqr5'"')
post `results_female' ("Range") (`"`range1'"') (`"`range2'"') (`"`range3'"') (`"`range4'"') (`"`range5'"')
post `results_female' ("Missing, n (%)") (`"`missing1'"') (`"`missing2'"') (`"`missing3'"') (`"`missing4'"') (`"`missing5'"')

postclose `results_female'

use `postresults_female', clear

* Table for Females
display as text %16s "Statistic" %16s "Wave 1" %16s "Wave 2" %16s "Wave 3" %16s "Wave 4" %16s "Wave 5"
display as text "{hline 96}"

forvalues i = 1/4 {
    local stat_val = statistic[`i']
    local w1_val = wave1[`i']
    local w2_val = wave2[`i']
    local w3_val = wave3[`i']
    local w4_val = wave4[`i']
    local w5_val = wave5[`i']
    
    display as text %16s "`stat_val'" %16s "`w1_val'" %16s "`w2_val'" %16s "`w3_val'" %16s "`w4_val'" %16s "`w5_val'"
}

display as text "{hline 96}"
display as text ""

* Identifier for females and combine all three groups
gen sex_group = "Females"
gen row_order = _n
order sex_group statistic wave1 wave2 wave3 wave4 wave5 row_order
tempfile females_results
save `females_results'
use `overall_results', clear
append using `males_results'
append using `females_results'

* Export combined table (overall, males, females)
export delimited "$outputtables/fi_missing_summary.csv", replace

restore

* ------------------------------------------------------------------------
* 3.4 Analysis of Frailty Categories (w*fc) waves 1 to 5
*-------------------------------------------------------------------------

preserve
tempname results_cat
tempfile postresults_cat

postfile `results_cat' str25 category str25 wave1 str25 wave2 str25 wave3 str25 wave4 str25 wave5 using `postresults_cat', replace

local waves "w1fc w2fc w3fc w4fc w5fc"

* Calculate categorical distributions for all waves
forvalues i = 1/5 {
    local w : word `i' of `waves'
    display as text "Debug for `w':"
    tab `w', missing
    
    * Total sample size (all observations)
    quietly count
    local total_all`i' = r(N)
    
    * Not Frail (0) vs Frail (1)
    quietly count if `w' == 0
    local not_frail`i' = r(N)
    local pct_not_frail`i' = 100 * `not_frail`i'' / `total_all`i''
    quietly count if `w' == 1
    local frail`i' = r(N)
    local pct_frail`i' = 100 * `frail`i'' / `total_all`i''
    
    * Missing (as percentage of total sample)
    quietly count if missing(`w')
    local missing`i' = r(N)
    local pct_missing`i' = 100 * `missing`i'' / `total_all`i''
    display as text "  Not frail: `not_frail`i'' (`pct_not_frail`i')"
    display as text "  Frail: `frail`i'' (`pct_frail`i')"
    display as text "  Missing: `missing`i'' (`pct_missing`i')"
    
    * Verification: percentages sum to 100%
    local total_pct = `pct_not_frail`i'' + `pct_frail`i'' + `pct_missing`i''
    display as text "  Total percentage: `total_pct'% (should be 100%)"
    
    * Format strings
    local not_frail_str`i' "`not_frail`i'' (`=string(`pct_not_frail`i'', "%5.1f")')"
    local frail_str`i' "`frail`i'' (`=string(`pct_frail`i'', "%5.1f")')"
    local missing_str`i' "`missing`i'' (`=string(`pct_missing`i'', "%5.1f")')"
}

* Post categorical results
post `results_cat' ("Not Frail, n (%)") (`"`not_frail_str1'"') (`"`not_frail_str2'"') (`"`not_frail_str3'"') (`"`not_frail_str4'"') (`"`not_frail_str5'"')
post `results_cat' ("Frail, n (%)") (`"`frail_str1'"') (`"`frail_str2'"') (`"`frail_str3'"') (`"`frail_str4'"') (`"`frail_str5'"')
post `results_cat' ("Missing, n (%)") (`"`missing_str1'"') (`"`missing_str2'"') (`"`missing_str3'"') (`"`missing_str4'"') (`"`missing_str5'"')

postclose `results_cat'

use `postresults_cat', clear

forvalues i = 1/3 {
    local cat_val = category[`i']
    local w1_val = wave1[`i']
    local w2_val = wave2[`i']
    local w3_val = wave3[`i']
    local w4_val = wave4[`i']
    local w5_val = wave5[`i']
    
    display as text %20s "`cat_val'" %16s "`w1_val'" %16s "`w2_val'" %16s "`w3_val'" %16s "`w4_val'" %16s "`w5_val'"
}

display as text "{hline 100}"
display as text ""

* Add identifier and export
gen analysis_type = "Categorical_Overall"
gen row_order = _n + 12  // after the 12 continuous rows (3 groups × 4 statistics)
order analysis_type category wave1 wave2 wave3 wave4 wave5 row_order

* Rename category to statistic to match existing structure
rename category statistic
gen sex_group = "Categorical"
order sex_group statistic wave1 wave2 wave3 wave4 wave5 row_order

tempfile categorical_results
save `categorical_results'

* Append to existing combined results
import delimited "$outputtables/fi_missing_summary.csv", clear
append using `categorical_results'

export delimited "$outputtables/fi_missing_summary.csv", replace

display as text "Categorical analysis appended to: $outputtables/fi_missing_summary.csv"

restore

* ------------------------------------------------------------------------
* 3.5 Categorical Analysis - Frailty Categories (w*fc) MALES ONLY
*-------------------------------------------------------------------------

preserve
tempname results_cat_male
tempfile postresults_cat_male

postfile `results_cat_male' str25 category str25 wave1 str25 wave2 str25 wave3 str25 wave4 str25 wave5 using `postresults_cat_male', replace

local waves "w1fc w2fc w3fc w4fc w5fc"

* Calculate categorical distributions for all waves for MALES only
forvalues i = 1/5 {
    local w : word `i' of `waves'
    display as text "Debug for `w' (Males only):"
    tab `w' if sex == 1, missing
    quietly count if sex == 1
    local total_all`i' = r(N)
    
    * Not Frail (0) vs Frail (1)
    quietly count if `w' == 0 & sex == 1
    local not_frail`i' = r(N)
    local pct_not_frail`i' = 100 * `not_frail`i'' / `total_all`i''
    quietly count if `w' == 1 & sex == 1
    local frail`i' = r(N)
    local pct_frail`i' = 100 * `frail`i'' / `total_all`i''
    
    * Missing (as percentage of total male sample)
    quietly count if missing(`w') & sex == 1
    local missing`i' = r(N)
    local pct_missing`i' = 100 * `missing`i'' / `total_all`i''
    display as text "  Not frail (Males): `not_frail`i'' (`pct_not_frail`i')"
    display as text "  Frail (Males): `frail`i'' (`pct_frail`i')"
    display as text "  Missing (Males): `missing`i'' (`pct_missing`i')"
    local total_pct = `pct_not_frail`i'' + `pct_frail`i'' + `pct_missing`i''
    display as text "  Total percentage (Males): `total_pct'% (should be 100%)"
    local not_frail_str`i' "`not_frail`i'' (`=string(`pct_not_frail`i'', "%5.1f")')"
    local frail_str`i' "`frail`i'' (`=string(`pct_frail`i'', "%5.1f")')"
    local missing_str`i' "`missing`i'' (`=string(`pct_missing`i'', "%5.1f")')"
}

* Post categorical results for males 
post `results_cat_male' ("Not Frail, n (%)") (`"`not_frail_str1'"') (`"`not_frail_str2'"') (`"`not_frail_str3'"') (`"`not_frail_str4'"') (`"`not_frail_str5'"')
post `results_cat_male' ("Frail, n (%)") (`"`frail_str1'"') (`"`frail_str2'"') (`"`frail_str3'"') (`"`frail_str4'"') (`"`frail_str5'"')
post `results_cat_male' ("Missing, n (%)") (`"`missing_str1'"') (`"`missing_str2'"') (`"`missing_str3'"') (`"`missing_str4'"') (`"`missing_str5'"')

postclose `results_cat_male'

use `postresults_cat_male', clear

* Display categorical results
forvalues i = 1/3 {
    local cat_val = category[`i']
    local w1_val = wave1[`i']
    local w2_val = wave2[`i']
    local w3_val = wave3[`i']
    local w4_val = wave4[`i']
    local w5_val = wave5[`i']
    
    display as text %20s "`cat_val'" %16s "`w1_val'" %16s "`w2_val'" %16s "`w3_val'" %16s "`w4_val'" %16s "`w5_val'"
}

display as text "{hline 100}"
display as text ""

* Add identifier and export for males
gen analysis_type = "Categorical_Males"
gen row_order = _n + 15  // after 15 rows (12 continuous + 3 overall categorical)
order analysis_type category wave1 wave2 wave3 wave4 wave5 row_order
rename category statistic
gen sex_group = "Categorical_Males"
order sex_group statistic wave1 wave2 wave3 wave4 wave5 row_order

tempfile categorical_males_results
save `categorical_males_results'
import delimited "$outputtables/fi_missing_summary.csv", clear
append using `categorical_males_results'

export delimited "$outputtables/fi_missing_summary.csv", replace

display as text "Males categorical analysis appended to: $outputtables/fi_missing_summary.csv"

restore

* ------------------------------------------------------------------------
* 3.6 Categorical Analysis - Frailty Categories (w*fc) FEMALES ONLY
*-------------------------------------------------------------------------

preserve
tempname results_cat_female
tempfile postresults_cat_female

postfile `results_cat_female' str25 category str25 wave1 str25 wave2 str25 wave3 str25 wave4 str25 wave5 using `postresults_cat_female', replace

local waves "w1fc w2fc w3fc w4fc w5fc"

* Calculate categorical distributions for all waves for FEMALES only
forvalues i = 1/5 {
    local w : word `i' of `waves'
    display as text "Debug for `w' (Females only):"
    tab `w' if sex == 2, missing
    
    quietly count if sex == 2
    local total_all`i' = r(N)
    quietly count if `w' == 0 & sex == 2
    local not_frail`i' = r(N)
    local pct_not_frail`i' = 100 * `not_frail`i'' / `total_all`i''
    quietly count if `w' == 1 & sex == 2
    local frail`i' = r(N)
    local pct_frail`i' = 100 * `frail`i'' / `total_all`i''
    
    * Missing (as percentage of total female sample)
    quietly count if missing(`w') & sex == 2
    local missing`i' = r(N)
    local pct_missing`i' = 100 * `missing`i'' / `total_all`i''
    display as text "  Not frail (Females): `not_frail`i'' (`pct_not_frail`i')"
    display as text "  Frail (Females): `frail`i'' (`pct_frail`i')"
    display as text "  Missing (Females): `missing`i'' (`pct_missing`i')"
    local total_pct = `pct_not_frail`i'' + `pct_frail`i'' + `pct_missing`i''
    display as text "  Total percentage (Females): `total_pct'% (should be 100%)"
    local not_frail_str`i' "`not_frail`i'' (`=string(`pct_not_frail`i'', "%5.1f")')"
    local frail_str`i' "`frail`i'' (`=string(`pct_frail`i'', "%5.1f")')"
    local missing_str`i' "`missing`i'' (`=string(`pct_missing`i'', "%5.1f")')"
}

* Post categorical results for females 
post `results_cat_female' ("Not Frail, n (%)") (`"`not_frail_str1'"') (`"`not_frail_str2'"') (`"`not_frail_str3'"') (`"`not_frail_str4'"') (`"`not_frail_str5'"')
post `results_cat_female' ("Frail, n (%)") (`"`frail_str1'"') (`"`frail_str2'"') (`"`frail_str3'"') (`"`frail_str4'"') (`"`frail_str5'"')
post `results_cat_female' ("Missing, n (%)") (`"`missing_str1'"') (`"`missing_str2'"') (`"`missing_str3'"') (`"`missing_str4'"') (`"`missing_str5'"')

postclose `results_cat_female'

use `postresults_cat_female', clear

* Display categorical results for females
forvalues i = 1/3 {
    local cat_val = category[`i']
    local w1_val = wave1[`i']
    local w2_val = wave2[`i']
    local w3_val = wave3[`i']
    local w4_val = wave4[`i']
    local w5_val = wave5[`i']
    
    display as text %20s "`cat_val'" %16s "`w1_val'" %16s "`w2_val'" %16s "`w3_val'" %16s "`w4_val'" %16s "`w5_val'"
}

display as text "{hline 100}"
display as text ""

* Add identifier and export for females
gen analysis_type = "Categorical_Females"
gen row_order = _n + 18  // after 18 rows (12 continuous + 6 categorical)
order analysis_type category wave1 wave2 wave3 wave4 wave5 row_order
rename category statistic
gen sex_group = "Categorical_Females"
order sex_group statistic wave1 wave2 wave3 wave4 wave5 row_order

tempfile categorical_females_results
save `categorical_females_results'
import delimited "$outputtables/fi_missing_summary.csv", clear
append using `categorical_females_results'

export delimited "$outputtables/fi_missing_summary.csv", replace

display as text "Females categorical analysis appended to: $outputtables/fi_missing_summary.csv"

restore

* ------------------------------------------------------------------------
* 3.7 Histograms of FI by wave (baseline cohort)
*-------------------------------------------------------------------------

* Histogram for Wave 1 Frailty Index
histogram w1fi, frequency ///
    title("Wave 1", size(medium)) ///
    xtitle("Frailty Index") ///
    ytitle("Frequency") ///
    bin(30) ///
    fcolor(eltblue) ///
    lcolor(black) ///
    lwidth(thin) ///
    graphregion(color(white)) ///
    plotregion(color(white)) ///
    scheme(s1mono)

* Histograms for w1fi to w5fi
forvalues w = 1/5 {
    histogram w`w'fi, frequency ///
        title("Wave `w'", size(medium)) ///
        xtitle("Frailty Index") ///
        ytitle("Frequency") ///
        bin(30) ///
        fcolor(eltblue) ///
        lcolor(black) ///
        lwidth(thin) ///
        graphregion(color(white)) ///
        plotregion(color(white)) ///
        scheme(s1mono)
    graph export "$outputgraphs/hist_baseline_w`w'fi.pdf", replace
}

* Histograms for w1fi to w5fi by sex
forvalues w = 1/5 {
    histogram w`w'fi, by(sex, note("")) frequency ///
        title("Wave `w'", size(medium)) ///
        xtitle("Frailty Index") ///
        ytitle("Frequency") ///
        bin(30) ///
        fcolor(eltblue) ///
        lcolor(black) ///
        lwidth(thin) ///
        graphregion(color(white)) ///
        plotregion(color(white)) ///
        scheme(s1mono)
    graph export "$outputgraphs/hist_baseline_w`w'fi_by_sex.pdf", replace
}

* ==============================================================================
* 4. Analysing if missing Data due to death or dropout   
* ==============================================================================

* ------------------------------------------------------------------------------
* 4.1 Pattern of the extended codes for death year
* ------------------------------------------------------------------------------

*-- Death-year special missing codes:
*   .d = don't know (participant died, year unknown)
*   .r = refused
*   .s = skipped
*   .x = alive

tab dyear, missing

list patid dyear w1fi w2fi w3fi w4fi w5fi if dyear == .d
list patid dyear w1fi w2fi w3fi w4fi w5fi if dyear == .r
list patid dyear w1fi w2fi w3fi w4fi w5fi if dyear == .s

* ------------------------------------------------------------------------------
* 4.2  Completely missing data across 5 waves
* ------------------------------------------------------------------------------

tab dyear if missing(w1fi), missing

* Exploring .s & missing(w1fi)
list patid dyear w1fi w2fi w3fi w4fi w5fi if dyear ==.s & missing(w1fi)

* Exploring alive (.x) & missing(w1fi)
list patid dyear w1fi w2fi w3fi w4fi w5fi if dyear ==.x & missing(w1fi)

* -- Completely missing across 5 waves (removed in data cleaning; expected none)
tab completely_lost, missing
list w*fi completely_lost

* ==============================================================================
* 5. Analysing Dropout
* ==============================================================================

* Definition: participants not known to have died who stop participating permanently after a given wave

* -------
* wave 1
* -------

* -- Dropout after Wave 1
gen dropout_after_w1 = missing(dyear) & !missing(w1fi) & (missing(w2fi) & missing(w3fi) & missing(w4fi) & missing(w5fi))
tab dropout_after_w1 dyear, missing

* Consistency check
list dropout_after_w1 dyear w*fi if dropout_after_w1 == 1, sum mean

* -------
* wave 2
* -------

* -- Dropout after Wave 2
gen dropout_after_w2 = missing(dyear) & !missing(w2fi) & (missing(w3fi) & missing(w4fi) & missing(w5fi))
tab dropout_after_w2 dyear, missing

* Consistency check
list dropout_after_w2 dyear w*fi if dropout_after_w2 == 1

* -------
* wave 3
* -------

* -- Dropout after Wave 3
gen dropout_after_w3 = missing(dyear) & !missing(w3fi) & (missing(w4fi) & missing(w5fi))
tab dropout_after_w3 dyear, missing

* Consistency check
list dropout_after_w3 dyear w*fi if dropout_after_w3 ==1

* -------
* wave 4
* -------

* -- Dropout after Wave 4
gen dropout_after_w4 = missing(dyear) & !missing(w4fi) & missing(w5fi)
tab  dropout_after_w4 dyear, missing

* Consistency check
list dropout_after_w4 dyear w*fi if dropout_after_w4 == 1, sum mean

* -------
* wave 5
* -------

* -- Dropout after Wave 5
* No dropout can be detected after wave 5 (end of study); included for code symmetry (expected zero)
gen dropout_after_w5 = missing(dyear) & dropout_after_w4 == 1 & !missing(w5fi) 
tab dropout_after_w5 dyear, missing

* Consistency check
list dropout_after_w5 dyear w*fi if dropout_after_w5 == 1, mean sum

*==================================================================================
* 6. Analysing Death Events
*================================================================================== 
* -------
* wave 1
* -------

* -- Death between wave 1 and wave 2 (2001–2002)
* Check inconsistencies between year of death and FI data
list patid w*fi dyear if (dyear <2003) 
list patid w*fi dyear if dyear <2003 & !missing(w2fi) 
list patid w1fi w2fi dyear if (w2fi >0 & dyear ==2001) & !missing(w2fi)

gen death_event_w1 = (dyear < 2003) 
tab death_event_w1 dyear, missing

* Consistency check
list death_event_w1 dyear w*fi if death_event_w1 == 1, sum mean

* -------
* wave 2
* -------

* -- Death between wave 2 and wave 3 (2003–2011)
list patid w*fi dyear if (dyear >=2003 & dyear <2012) & (!missing(w3fi) | !missing(w4fi, w5fi))

tab dyear if missing(w2fi), missing

list dyear w*fi if dyear == .r

gen death_event_w2 = (dyear >= 2003 & dyear < 2012) 
tab death_event_w2 dyear, missing

* Consistency check
list death_event_w2 dyear w*fi if death_event_w2 ==1 , sum mean

* -------
* Wave 3
* -------

* -- Death between wave 3 and wave 4 (2012–2014)
list patid w*fi if (dyear >= 2012 & dyear < 2015) & (!missing(w4fi) | !missing(w5fi))

tab dyear if missing(w3fi), missing
gen death_event_w3 = (dyear >= 2012 & dyear <2015)
tab death_event_w3 dyear, missing

* Consistency check
list death_event_w3 dyear w*fi if death_event_w3 == 1, sum mean

* -------
* wave 4
* -------

* -- Death between wave 4 and wave 5 (2015–2017)
list patid w*fi if (dyear >=2015 & dyear <2018) & !missing(w5fi)
tab dyear if missing(w4fi), missing
gen death_event_w4 = (dyear>=2015 & dyear < 2018) 

* Consistency check
list death_event_w4 dyear w*fi if death_event_w4 == 1, sum mean

* -------
* wave 5
* -------

* -- Death in 2018 (wave 5)
list patid w*fi dyear if (dyear == 2018) & !missing(w5fi) 

gen death_event_w5 = (dyear ==2018)
tab death_event_w5 dyear, missing

* Consistency check
list death_event_w5 dyear w*fi if death_event_w5 == 1, mean sum

* -------------------------------------------------------------------------------
* 5.1 Summary of Dropout and Death by Wave  
* -------------------------------------------------------------------------------
* -- Deaths by inter-wave period
gen death_period = .
replace death_period = 1 if dyear >= 2001 & dyear < 2003
replace death_period = 2 if dyear >= 2003 & dyear < 2012
replace death_period = 3 if dyear >= 2012 & dyear < 2015
replace death_period = 4 if dyear >= 2015 & dyear < 2018
replace death_period = 5 if dyear >= 2018 & dyear <= 2020

label define death_period 1 "2001-2002 (W1)" 2 "2003-2011 (W2-W3)" ///
                         3 "2012-2014 (W3-W4)" 4 "2015-2017 (W4-W5)" ///
                         5 "2018-2020 (Post-W5)"
label values death_period death_period

tab death_period, missing
display as text "Deaths by wave periods shown above"

* Graph of deaths by wave period
preserve

contract death_period if !missing(death_period)

* Bar graph
twoway (bar _freq death_period, barwidth(0.5) lcolor(black) lwidth(medium)) ///
       , title("Deaths by Wave Periods") ///
         xtitle("Wave Period") ytitle("Number of Deaths") ///
         xlabel(1 "W1-W2" 2 "W2-W3" 3 "W3-W4" 4 "W4-W5" 5 "Post-W5", angle(45)) ///
         legend(off) ///
         graphregion(color(white)) bgcolor(white) ///
         note("Based on death year and wave periods")

graph export "$outputgraphs/deaths_by_wave_period.png", replace

* Line graph (alternative)
twoway (line _freq death_period, lcolor(black) lwidth(medium) msymbol(circle) mcolor(black)) ///
    , title("Deaths by Wave Periods") ///
      xtitle("Wave Period") ytitle("Number of Deaths") ///
      xlabel(1 "W1-W2" 2 "W2-W3" 3 "W3-W4" 4 "W4-W5" 5 "Post-W5", angle(45)) ///
      legend(off) ///
      graphregion(color(white)) bgcolor(white) ///
      note("Based on death year and wave periods")

restore

* ------------------------------------------------------------------------------
* 5.2 Number and percentage missing due to death and dropout by wave
* ------------------------------------------------------------------------------

preserve

tempname results
tempfile postresults

postfile `results' str3 wave int n_death int n_dropout ///
                   float pct_death float pct_dropout using `postresults', replace

foreach w in w1 w2 w3 w4 w5 {
    quietly count if dropout_after_`w' == 1
    local n_dropout = r(N)

    quietly count if inlist(dropout_after_`w', 0, 1)
    local n_atrisk = r(N)

    quietly count if death_event_`w' == 1
    local n_death = r(N)

    quietly count if inlist(death_event_`w', 0, 1)
    local n_death_atrisk = r(N)

    local pct_dropout = 100 * `n_dropout' / `n_atrisk'
    local pct_death   = 100 * `n_death' / `n_death_atrisk'

    post `results' ("`w'") (`n_death') (`n_dropout') ///
         (`pct_death') (`pct_dropout')
}

postclose `results'

use `postresults', clear
gen n_atrisk = n_death + n_dropout
gen wave_num = _n

twoway (line pct_death wave_num, sort lpattern(dash) lwidth(medthick)) ///
       (line pct_dropout wave_num, sort lpattern(solid) lwidth(medthick)), ///
       title("Dropout and Death by Wave") ///
       ytitle("Percent") xtitle("Wave") ///
       xlabel(1 "w1" 2 "w2" 3 "w3" 4 "w4" 5 "w5") ///
       legend(order(1 "Death" 2 "Dropout"))

use `postresults', clear
list, sep(0)

restore

* ===============================================================================
* 6.0 Summary of Dropout and Death per Wave by sex
* ===============================================================================

*-- Participant status at end of follow-up

gen status = .
* (1) Alive: alive (.x) or refused (.r) with data at wave 5
replace status = 1 if inlist(dyear, .x, .r) & !missing(w5fi)

* (2) Died: known year of death
replace status = 2 if !missing(dyear)

* (3) Lost to follow-up, known dead (LFU-dead): death year unknown (.d) or skipped (.s), no wave 5 data
replace status = 3 if inlist(dyear, .d, .s) & missing(w5fi)

* (4) Lost to follow-up, vital status unknown (dropout): alive/refused with no wave 5 data
replace status = 4 if inlist(dyear, .x, .r) & missing(w5fi)

label define status 1 "Alive" 2 "Deaths" 3 "LFU-dead" 4 "LFU-unknown (dropout)"
label values status status

* Verify status assignment
tab status, missing
tab status dyear, missing

* ------------------------------------------------------------------------------
* 6.1 Create summary table by sex
* ------------------------------------------------------------------------------

preserve

tempname results
tempfile postresults

postfile `results' str20 category str20 total str20 females str20 males using `postresults', replace

quietly count
local total_n = r(N)

quietly count if sex == 2  // Females
local female_n = r(N)

quietly count if sex == 1  // Males  
local male_n = r(N)

* Header row
post `results' ("N") ("`total_n'") ("`female_n'") ("`male_n'")

* Each status category
foreach stat_num in 1 2 3 4 {
    local stat_label : label status `stat_num'
    
    * Total
    quietly count if status == `stat_num'
    local total_count = r(N)
    local total_pct = 100 * `total_count' / `total_n'
    local total_str "`total_count' (" + string(`total_pct', "%5.2f") + ")"
    
    * Females
    quietly count if status == `stat_num' & sex == 2
    local female_count = r(N)
    local female_pct = 100 * `female_count' / `female_n'
    local female_str "`female_count' (" + string(`female_pct', "%5.2f") + ")"
    
    * Males
    quietly count if status == `stat_num' & sex == 1
    local male_count = r(N)
    local male_pct = 100 * `male_count' / `male_n'
    local male_str "`male_count' (" + string(`male_pct', "%5.2f") + ")"
    
    post `results' ("`stat_label', N (percent)") ("`total_str'") ("`female_str'") ("`male_str'")
}

* Combined LFU (categories 3 and 4)
quietly count if inlist(status, 3, 4)
local lfu_total = r(N)
local lfu_total_pct = 100 * `lfu_total' / `total_n'
local lfu_total_str "`lfu_total' (" + string(`lfu_total_pct', "%5.2f") + ")"

quietly count if inlist(status, 3, 4) & sex == 2
local lfu_female = r(N)
local lfu_female_pct = 100 * `lfu_female' / `female_n'
local lfu_female_str "`lfu_female' (" + string(`lfu_female_pct', "%5.2f") + ")"

quietly count if inlist(status, 3, 4) & sex == 1
local lfu_male = r(N)
local lfu_male_pct = 100 * `lfu_male' / `male_n'
local lfu_male_str "`lfu_male' (" + string(`lfu_male_pct', "%5.2f") + ")"

post `results' ("Lost to follow-up (LFU), N (percent)") ("`lfu_total_str'") ("`lfu_female_str'") ("`lfu_male_str'")
postclose `results'

use `postresults', clear

* Display table
forvalues i = 1/`=_N' {
    local cat = category[`i']
    local tot = total[`i']
    local fem = females[`i']
    local mal = males[`i']
    
    display as text %35s "`cat'" %20s "`tot'" %20s "`fem'" %20s "`mal'"
}

display as text "{hline 95}"

* Export to CSV
export delimited "$outputtables/participant_status_summary.csv", replace
gen row_order = _n
order category total females males row_order

restore

*--------------------------------------------------------------------------------
* Proportion of deaths in total sample
*--------------------------------------------------------------------------------
gen death = .
replace death = 1 if inlist(status, 2, 3)
replace death = 0 if inlist(status, 1, 4)
label variable death "Death Status (1 = Death, 0 = Alive or Other)"

* Proportion of deaths by sex (weighted and unweighted)
svy: tab death sex, col percent 
tab death sex, col 

* Save dataset with death/status coding (used in subsequent analyses)
save "$cleanmhasdata/wave1to5_wide_deaths_code.dta", replace

* ==============================================================================
* 8. Missing Data Analysis - MCAR, MAR, MNAR
* ==============================================================================

* ------------------------------------------------------------------------------
* 8.1 Count participants with complete FI data in all 5 waves
* ------------------------------------------------------------------------------
misstable summarize w*fi
count if !missing(w1fi, w2fi, w3fi, w4fi, w5fi)
local complete_fi = r(N)

* Total participants
count
local total = r(N)

local pct = 100 * `complete_fi' / `total'

display as text "Participants with complete FI data: `complete_fi' out of `total' (`pct'%)"

*--- Missing-data patterns (excluding deaths)
* Pattern string: 1 = observed, 0 = missing (waves 1–5)
gen pattern = ""
replace pattern = pattern + cond(missing(w1fi), "0", "1")
replace pattern = pattern + cond(missing(w2fi), "0", "1") 
replace pattern = pattern + cond(missing(w3fi), "0", "1")
replace pattern = pattern + cond(missing(w4fi), "0", "1")
replace pattern = pattern + cond(missing(w5fi), "0", "1")

* Total participants not known to have died (excluding status 2 and 3)
count if !inlist(status, 2, 3)
local total_alive = r(N)

* Complete cases
count if pattern == "11111" & !inlist(status, 2, 3)
local complete = r(N)

* Monotone missing (dropout)
count if (pattern == "10000" | pattern == "11000" | pattern == "11100" | pattern == "11110") & !inlist(status, 2, 3)
local monotone = r(N)

* Intermittent missing (observed after a missing wave)
count if pattern != "11111" & pattern != "10000" & pattern != "11000" & pattern != "11100" & pattern != "11110" & pattern != "00000" & !inlist(status, 2, 3)
local intermittent = r(N)

* All missing
count if pattern == "00000" & !inlist(status, 2, 3)
local all_missing = r(N)

* Deaths
count if inlist(status, 2, 3)
local deaths = r(N)

* Percentages (denominator: participants not known to have died)
local pct_complete = 100 * `complete' / `total_alive'
local pct_monotone = 100 * `monotone' / `total_alive' 
local pct_intermittent = 100 * `intermittent' / `total_alive'
local pct_all_missing = 100 * `all_missing' / `total_alive'

tab pattern if pattern != "11111" & pattern != "10000" & pattern != "11000" & pattern != "11100" & pattern != "11110" & pattern != "00000" & !inlist(status, 2, 3), sort

* ==============================================================================
* 8.2 Intermittent missingness vs dropout: predictors of missingness
* ==============================================================================

preserve

use "$cleanmhasdata/wave1to5_wide_deaths_code.dta", clear

* Reshape to person–wave
capture drop wave
reshape long w@fi, i(patid) j(wave)    
recode wave (1=0)(2=2)(3=11)(4=14)(5=17), gen(time)
label var time "Years since baseline (0,2,11,14,17)"
xtset patid wave
svyset patid [pweight=w1respwt]

* Age centered
sum age, meanonly
gen age_c = age - r(mean)

* Person-level tags
bys patid: gen byte _tag = _n==1
quietly count if _tag
local N0 = r(N)                             
bys patid: egen byte ever_dead = max(inlist(status,2,3))   // ever died (Deaths or LFU-dead)
quietly count if _tag & ever_dead==0
local N_alive = r(N)                             
tempvar alive_end                               // alive at end of study (may be missing w5fi)
gen byte `alive_end' = inlist(dyear,.x,.r)
tempvar later_obs
gen byte `later_obs' = (!missing(F.wfi) | !missing(F2.wfi) | !missing(F3.wfi) | !missing(F4.wfi)) if wave<=4

* Death before next interview 
tempvar y_next death_next any_future
gen int  `y_next' = 2001 + F.time if wave<5
gen byte `death_next' = 0
replace `death_next' = 1 if wave<5 & !missing(dyear) & dyear <= `y_next'
capture confirm variable dmonth
if _rc==0 {
    replace `death_next' = 1 if wave<5 & inlist(dyear,.d,.r) & !missing(dmonth)   // special code + month => death
}
gen byte `any_future' = (!missing(F.wfi)|!missing(F2.wfi)|!missing(F3.wfi)|!missing(F4.wfi)) if wave<=4
replace `death_next' = 1 if wave<5 & dyear==.s & `any_future'==0                 // .s treated as death if never returns

*-- A) Intermittent missingness (missing now, observed later)
tempvar inter_t
gen byte `inter_t' = missing(wfi) & `later_obs' if wave<=4 & ever_dead==0
bys patid: egen byte intermittent = max(`inter_t')

* Per-wave logistic model, clustered by person
logit `inter_t' L.wfi i.wave i.sex c.age_c i.educ i.rural i.srh c.wealth_gbp if wave<=4 & ever_dead==0 [pweight=w1respwt], vce(cluster patid) or base
testparm L1.wfi i.wave i.sex c.age_c i.educ i.rural i.srh c.wealth_gbp

*-- B) Dropout (monotone): observed at w, missing at w+1 and never returns
tempvar drop_t
gen byte `drop_t' = missing(F.wfi) & !missing(wfi)  & !(!missing(F2.wfi) | !missing(F3.wfi) | !missing(F4.wfi)) if wave<5
replace `drop_t' = . if `death_next'                 // exclude death-related missingness 
replace `drop_t' = . if ever_dead==1                 // exclude anyone who ever dies

* Wave 4: refusal at wave 5 treated as non-response, not dropout
replace `drop_t' = 0 if wave==4 & `alive_end'==1 & inlist(dyear,.r)

bys patid: egen byte dropout = max(`drop_t')

* Per-wave logistic model, clustered by person
logit `drop_t' L.wfi i.wave i.sex c.age_c i.educ i.rural i.srh c.wealth_gbp if wave<5 & !missing(`drop_t') [pweight=w1respwt], vce(cluster patid) or base
testparm L1.wfi i.wave i.sex c.age_c i.educ i.rural i.srh c.wealth_gbp

* Sanity check by wave
tab wave `inter_t', row m
tab wave `drop_t', row m

restore

* ==============================================================================
* 11. Session end
* ==============================================================================
display as text "------------------------------------------------------------"
display "==== Missing data analysis_ COMPLETED SUCCESSFULLY ===="
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close