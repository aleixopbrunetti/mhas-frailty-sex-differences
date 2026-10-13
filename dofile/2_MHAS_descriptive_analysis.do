********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Descriptive analysis (Table 1)
********************************************************************************

* ==============================================================================
* 3. Descriptive analysis
* ==============================================================================

* -- Load baseline data
use "$cleanmhasdata/wave1to5_wide.dta", clear

* ------------------------------------------------------------------------------
* 3.1 Demographic Characteristics
* ------------------------------------------------------------------------------

* age
sum age, d

* sex
tab sex, missing

* education
tab educ, missing

* rural
tab rural, missing  

* ------------------------------------------------------------------------------
* 3.2 Health Characteristics
* ------------------------------------------------------------------------------
* frailty index
sum w1fi, d
tab w1fc, missing
tab w1fc sex, chi2

* BMI characteristics
sum bmi, d  
tab bmicat sex, missing
tab bmicat sex, chi2

* adl functional limitations
tab adl, missing
replace adl = . if adl == .d
replace adl = . if adl == .m
replace adl = . if adl == .r
tab adl, missing
tab adl if sex ==1, missing
tab adl if sex ==2, missing
tab adl sex, chi2

* Percentage with adl >= 3 (whole sample, males, females)
tab adl if adl>=3
display as text "Percentage with adl >= 3:"
count if adl >= 3 & !missing(adl)
local n_above3 = r(N)
count
local total = r(N)
display as text "  " %6.2f (100 * `n_above3' / `total') "%"

tab adl if adl>=3 & sex==1
display as text "Percentage with adl >= 3:"
count if adl >= 3 & !missing(adl) & sex == 1
local n_above3 = r(N)
count
local total = r(N)
display as text "  " %6.2f (100 * `n_above3' / `total') "%"

tab adl if adl>=3 & sex==2
display as text "Percentage with adl >= 3:"
count if adl >= 3 & !missing(adl) & sex == 2
local n_above3 = r(N)
count
local total = r(N)
display as text "  " %6.2f (100 * `n_above3' / `total') "%"

* Self-reported health
tab srh, missing    
tab srh sex, chi2

* diabetes
tab diabe, missing

* cognitive status
tab cogstat, missing  

* adl functional limitations
tab adl, missing

* ------------------------------------------------------------------------------
* 3.3 Health Behaviors and Lifestyle
* ------------------------------------------------------------------------------

* physical activity
tab vigact, missing

* smoking status
tab smoke_status, missing  
tab smokev, missing
tab smoke, missing 

* alcohol consumption
tab drink_status, missing
tab drink, missing
sum drinkdw, d
summarize drinkdw if drink == 1, detail

* health insurance
tab higov, missing

* ------------------------------------------------------------------------------
* 3.4 Wealth Characteristics
* ------------------------------------------------------------------------------
sum wealth_gbp, d
tab wealth_cat, missing   

* ===============================================================================
* 4. Extracting summary statistics for continuous variables
* ==============================================================================

* List of variables to summarize
local vars_demogr w1fi age bmi wealth_gbp

* Create matrices to store results
local wordcount: word count `vars_demogr'

* Matrices for whole sample, males, and females
matrix results_whole = J(`wordcount',10,.)
matrix results_male = J(`wordcount', 10, .)
matrix results_female = J(`wordcount', 10, .)
matrix ttest_results = J(`wordcount', 1, .)

matrix colnames results_whole = "Mean" "SD" "Median" "Min" "Max" "Skewness" "Kurtosis" "N" "Missing" "%_Missing"
matrix colnames results_male = "Mean" "SD" "Median" "Min" "Max" "Skewness" "Kurtosis" "N" "Missing" "%_Missing"
matrix colnames results_female = "Mean" "SD" "Median" "Min" "Max" "Skewness" "Kurtosis" "N" "Missing" "%_Missing"
matrix colnames ttest_results = "p_value"

matrix rownames results_whole = `vars_demogr'
matrix rownames results_male = `vars_demogr'
matrix rownames results_female = `vars_demogr'
matrix rownames ttest_results = `vars_demogr'

* Initialize row counter
local row = 1

* Loop over demographic variables
foreach var of local vars_demogr {
    * Whole sample statistics
    quietly sum `var', detail
    matrix results_whole[`row', 1] = r(mean)
    matrix results_whole[`row', 2] = r(sd)
    matrix results_whole[`row', 3] = r(p50)
    matrix results_whole[`row', 4] = r(min)
    matrix results_whole[`row', 5] = r(max)
    matrix results_whole[`row', 6] = r(skewness)
    matrix results_whole[`row', 7] = r(kurtosis)
    quietly count
    local total = r(N)
    matrix results_whole[`row', 8] = `total'   // N: total sample size including missing
    quietly count if missing(`var')
    matrix results_whole[`row', 9] = r(N)
    local miss = r(N)
    matrix results_whole[`row', 10] = 100 * `miss' / `total'

    * Male statistics (sex == 1)
    quietly count if sex == 1
    local total_male = r(N)
    quietly sum `var' if sex == 1, detail
    matrix results_male[`row', 1] = r(mean)
    matrix results_male[`row', 2] = r(sd)
    matrix results_male[`row', 3] = r(p50)
    matrix results_male[`row', 4] = r(min)
    matrix results_male[`row', 5] = r(max)
    matrix results_male[`row', 6] = r(skewness)
    matrix results_male[`row', 7] = r(kurtosis)
    matrix results_male[`row', 8] = `total_male'   // N: total males including missing
    quietly count if missing(`var') & sex == 1
    matrix results_male[`row', 9] = r(N)
    local miss_male = r(N)
    matrix results_male[`row', 10] = cond(`total_male' > 0, 100 * `miss_male' / `total_male', .)
    
    * Female statistics (sex == 2)
    quietly count if sex == 2
    local total_female = r(N)
    quietly sum `var' if sex == 2, detail
    matrix results_female[`row', 1] = r(mean)
    matrix results_female[`row', 2] = r(sd)
    matrix results_female[`row', 3] = r(p50)
    matrix results_female[`row', 4] = r(min)
    matrix results_female[`row', 5] = r(max)
    matrix results_female[`row', 6] = r(skewness)
    matrix results_female[`row', 7] = r(kurtosis)
    matrix results_female[`row', 8] = `total_female'   // N: total females including missing
    quietly count if missing(`var') & sex == 2
    matrix results_female[`row', 9] = r(N)
    local miss_female = r(N)
    matrix results_female[`row', 10] = cond(`total_female' > 0, 100 * `miss_female' / `total_female', .)
    
    * T-test between sexes
    quietly ttest `var', by(sex)
    matrix ttest_results[`row', 1] = r(p)
    
    local ++row
}

* Display results
display as text "=== DESCRIPTIVE STATISTICS_01 ==="
display as text "=== Results Whole Sample ==="
matrix list results_whole, format(%9.3f)                                                

display as text ""
display as text "=== Results for males (sex = 1) ==="
matrix list results_male, format(%9.3f)

display as text ""
display as text "=== Results for females (sex = 2) ==="
matrix list results_female, format(%9.3f)

display as text ""
display as text "=== T-test (Male vs Female) ==="
matrix list ttest_results, format(%12.10f)

* Display detailed p-values
display as text ""
display as text "=== T-test P-values (Detailed) ==="
local row = 1
foreach var of local vars_demogr {
    scalar pval_scalar = ttest_results[`row', 1]
    local pval = scalar(pval_scalar)
    display as text "`var': p = " %10.6f `pval'
    local ++row
}

* Export results to Excel
local filename "MHAS_descriptive_analysis_01"
putexcel set "$outputtables/`filename'.xlsx", sheet("Whole Sample") replace
putexcel A1 = matrix(results_whole), names

putexcel set "$outputtables/`filename'.xlsx", sheet("Male") modify
putexcel A1 = matrix(results_male), names

putexcel set "$outputtables/`filename'.xlsx", sheet("Female") modify
putexcel A1 = matrix(results_female), names

putexcel set "$outputtables/`filename'.xlsx", sheet("Ttest Pvalues") modify
putexcel A1 = matrix(ttest_results), names

*===============================================================================
* 5. Categorical Analysis extraction
*===============================================================================

* Categorical variables to summarize:
* adl          - Activities of Daily Living
* educ         - Education level
* vigact       - Vigorous physical activity
* rural        - Rural residence indicator
* smoke_status - Smoking status
* smokev       - Ever smoked
* smoke        - Current smoker
* srh          - Self-reported health
* drink        - Alcohol consumption
* drinkdw      - Drinks per week
* drink_status - Drinking status
* higov        - Health insurance (government)
* diabe        - Diabetes status
* bmicat       - BMI category
* w1fc         - Frailty category (wave 1)
* cogstat      - Cognitive status
* wealth_cat   - Wealth category

* -------------------------------------------------------------------------
* 5.1 CATEGORICAL VARIABLES ANALYSIS
* -------------------------------------------------------------------------
* Summarize: frequencies, percentages, N, and Chi2 by sex, showing subcategories    

* Define categorical variables list
local vars_cat adl educ vigact rural smoke_status smokev smoke srh drink drinkdw ///
              drink_status higov diabe bmicat w1fc cogstat wealth_cat

* Create temporary files to store results
tempname catresults
tempfile categorical_results
postfile `catresults' str20 variable str20 category str10 sample int frequency float percentage int total_n float chi2 float p_value str20 freq_pct using `categorical_results', replace

* Loop through each categorical variable
foreach var of local vars_cat {
    
    display as text ""
    display as text "--- ANALYZING VARIABLE: `var' ---"
    
    * Get variable label
    local varlabel : variable label `var'
    if "`varlabel'" == "" local varlabel "`var'"
    display as text "Variable label: `varlabel'"

    * Whole sample analysis
    display as text ""
    display as text "WHOLE SAMPLE:"
    tab `var', missing
    
    * Store whole sample results
    quietly levelsof `var', local(levels) missing
    quietly count
    local total_all = r(N)
    foreach level of local levels {
        if "`level'" == "." {
            local level_label "Missing"
            quietly count if missing(`var')
            local freq = r(N)
            local pct = cond(`total_all' > 0, 100 * `freq' / `total_all', 0)
            local freq_pct = "`freq'(" + string(`pct', "%6.2f") + ")"
            post `catresults' ("`var'") ("`level_label'") ("Whole") (`freq') (`pct') (`total_all') (.) (.) ("`freq_pct'")
        }
        else {
            local level_label : label (`var') `level'
            if "`level_label'" == "" local level_label "`level'"
            quietly count if `var' == `level'
            local freq = r(N)
            local pct = cond(`total_all' > 0, 100 * `freq' / `total_all', 0)
            local freq_pct = "`freq'(" + string(`pct', "%6.2f") + ")"
            post `catresults' ("`var'") ("`level_label'") ("Whole") (`freq') (`pct') (`total_all') (.) (.) ("`freq_pct'")
        }
        display as text "  `level_label': " %6.0f `freq' " (" %6.2f `pct' "%)"
    }
    
    * Male analysis (sex == 1)
    display as text ""
    display as text "MALES (sex = 1):"
    tab `var' if sex == 1, missing
    
    quietly count if sex == 1
    local total_all_male = r(N)
    quietly levelsof `var' if sex == 1, local(levels_male) missing
    foreach level of local levels_male {
        if "`level'" == "." {
            local level_label "Missing"
            quietly count if missing(`var') & sex == 1
            local freq = r(N)
            local pct = cond(`total_all_male' > 0, 100 * `freq' / `total_all_male', 0)
            local freq_pct = "`freq'(" + string(`pct', "%6.2f") + ")"
            post `catresults' ("`var'") ("`level_label'") ("Male") (`freq') (`pct') (`total_all_male') (.) (.) ("`freq_pct'")
        }
        else {
            local level_label : label (`var') `level'
            if "`level_label'" == "" local level_label "`level'"
            quietly count if `var' == `level' & sex == 1
            local freq = r(N)
            local pct = cond(`total_all_male' > 0, 100 * `freq' / `total_all_male', 0)
            local freq_pct = "`freq'(" + string(`pct', "%6.2f") + ")"
            post `catresults' ("`var'") ("`level_label'") ("Male") (`freq') (`pct') (`total_all_male') (.) (.) ("`freq_pct'")
        }
        display as text "  `level_label': " `freq' " (" %4.1f `pct' "%)"
    }
    
    * Female analysis (sex == 2)
    display as text ""
    display as text "FEMALES (sex = 2):"
    tab `var' if sex == 2, missing
    
    quietly count if sex == 2
    local total_all_female = r(N)
    quietly levelsof `var' if sex == 2, local(levels_female) missing
    foreach level of local levels_female {
        if "`level'" == "." {
            local level_label "Missing"
            quietly count if missing(`var') & sex == 2
            local freq = r(N)
            local pct = cond(`total_all_female' > 0, 100 * `freq' / `total_all_female', 0)
            local freq_pct = "`freq'(" + string(`pct', "%6.2f") + ")"
            post `catresults' ("`var'") ("`level_label'") ("Female") (`freq') (`pct') (`total_all_female') (.) (.) ("`freq_pct'")
        }
        else {
            local level_label : label (`var') `level'
            if "`level_label'" == "" local level_label "`level'"
            quietly count if `var' == `level' & sex == 2
            local freq = r(N)
            local pct = cond(`total_all_female' > 0, 100 * `freq' / `total_all_female', 0)
            local freq_pct = "`freq'(" + string(`pct', "%6.2f") + ")"
            post `catresults' ("`var'") ("`level_label'") ("Female") (`freq') (`pct') (`total_all_female') (.) (.) ("`freq_pct'")
        }
        display as text "  `level_label': " `freq' " (" %4.1f `pct' "%)"
    }
    
    * Chi2 Test
    display as text ""
    capture quietly tab `var' sex, chi2
    local chi2_val = .
    local p_val = .
    if _rc == 0 {
        local chi2_val = r(chi2)
        local p_val = r(p)
        display as text "Chi-squared test (variable vs sex): chi2 = " %6.3f `chi2_val' ", p = " %6.4f `p_val'
    }
    else {
        display as text "Chi-squared test could not be performed for `var'"
    }
    
    * Post chi2 and p-value for the variable (once per variable, category = "Chi2")
    post `catresults' ("`var'") ("Chi2") ("Test") (.) (.) (.) (`chi2_val') (`p_val') ("")
    
    display as text "----------------------------------------"
}

postclose `catresults'

* Load and display summary results
preserve
use `categorical_results', clear

* Format results
format percentage %6.1f
format chi2 %8.3f
format p_value %8.4f

* Save the results
save "$outputtables/categorical_analysis_detailed.dta", replace
export delimited "$outputtables/categorical_analysis_detailed.csv", replace

* Display summary
display as text ""
display as text "=== SUMMARY OF ALL CATEGORICAL VARIABLES ==="
list, sepby(variable) abbreviate(15)

restore

* ------------------------------------------------------------------------------
* 5.2 Creating tables
* ------------------------------------------------------------------------------

display as text ""
display as text "Creating publication-ready summary tables..."

* Summary table for each variable showing all categories by sex
preserve

* Excel workbook with separate sheets for each variable
local excel_file "$outputtables/categorical_variables_by_sex.xlsx"

use `categorical_results', clear

* Get list of unique variables
quietly levelsof variable, local(variables)

* Create summary table for each variable
local sheet_counter = 0
foreach var of local variables {
    
    local ++sheet_counter
    
    * Keep only data for current variable
    tempfile temp_var_data
    keep if variable == "`var'"
    save `temp_var_data'
    
    * Reshape to have sample types as columns
    capture drop if category == "Chi2"  // Remove Chi2 rows before reshape
    
    capture reshape wide frequency percentage total_n chi2 p_value freq_pct, i(category) j(sample) string
    
    if _rc == 0 {
        * Clean up variable names
        capture rename frequencyWhole freq_whole
        capture rename percentageWhole pct_whole  
        capture rename total_nWhole n_whole
        capture rename freq_pctWhole freq_pct_whole
        capture rename frequencyMale freq_male
        capture rename percentageMale pct_male
        capture rename total_nMale n_male
        capture rename freq_pctMale freq_pct_male
        capture rename frequencyFemale freq_female
        capture rename percentageFemale pct_female
        capture rename total_nFemale n_female
        capture rename freq_pctFemale freq_pct_female
        
        * Export to Excel
        if `sheet_counter' == 1 {
            export excel using "`excel_file'", sheet("`var'") firstrow(variables) replace
        }
        else {
            export excel using "`excel_file'", sheet("`var'") firstrow(variables) sheetmodify
        }
        
        display as text "Exported `var' to Excel sheet"
    }
    else {
        display as text "Warning: Could not reshape data for variable `var'"
    }
    
    * Restore data for next variable
    use `temp_var_data', clear
}

restore

display as text "Categorical analysis complete!"
display as text "Results saved to:"
display as text "  - $outputtables/categorical_analysis_detailed.csv"
display as text "  - $outputtables/categorical_variables_by_sex.xlsx"

* ==============================================================================
* 6. Session End
* ==============================================================================

display as text ""
display as text "=========================================================================="
display as text "MHAS_Descriptive_Analysis _COMPLETED"
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "=========================================================================="

log close