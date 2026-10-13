********************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: Descriptive Analysis for MHAS Frailty Index Study
********************************************************************************

* ==============================================================================
* 2. Weighted Descriptive analysis
* ==============================================================================

* -- Load baseline data
use "$cleanmhasdata/wave1to5_wide.dta", clear
 

*===============================================================================
* 3. Declaring the weight variable
*===============================================================================

* aim: to adjust for unequal probabilities of selection. In this case, age < 50 years old. 

svyset [pweight=w1respwt]

sum w1respwt, d

*---------------------------------------------------------------
* 3.1 Weighted Continuous Variables Analysis
*---------------------------------------------------------------
svy: mean age // Sd error
estat sd // SD 
// https://stats.oarc.ucla.edu/stata/seminars/survey-data-analysis-in-stata-17/

* Weighted mean difference in age by sex
svy: mean age, over(sex)
svy: regress age i.sex


svy: tab educ sex, row percent missing
svy: mean bmi, over(sex)
svy: regress bmi i.sex

svy: mean w1fi

svy: regress w1fi i.sex

local vars_demogr w1fi age bmi wealth_gbp
foreach var of local vars_demogr {
    svy: mean `var'
    svy: mean `var', over(sex)
    svy: regress `var' i.sex
    test 1.sex
}


*---------------------------------------------------------------
* 3.2 Weighted Categorical Variables Analysis
*---------------------------------------------------------------
* Categorical variables

local vars_cat adl educ vigact rural smoke_status smokev smoke srh drink drinkdw ///
              drink_status higov diabe bmicat w1fc cogstat wealth_cat
foreach var of local vars_cat {
    svy: tabulate `var' sex, row percent missing 
    local N = r(N)
    display as text "Total N for `var': " `N'
}


local vars_cat adl educ vigact rural smoke_status smokev smoke srh drink drinkdw ///
              drink_status higov diabe bmicat w1fc cogstat wealth_cat
foreach var of local vars_cat {
    count if !missing(sex, `var')
    local N = r(N)
    di "==== `var' by sex ===="
    di "Sample size: `N'"
    svy: tabulate `var' sex, row percent missing
}

svy: tab wealth_cat if sex==1, percent missing
svy: tab wealth_cat if sex==2, percent missing


*===============================================================================
* 4. Extracting Results
*===============================================================================

*----------------------------------------------------------------------
* 4.1 Extracting Weighted Continuous Variables Analysis
*----------------------------------------------------------------------

* Ensure survey design is properly set
svyset [pweight=w1respwt]

* List of variables to summarize
local vars_demogr w1fi age bmi wealth_gbp

* Create matrices to store weighted results
local wordcount: word count `vars_demogr'

* Matrices for weighted whole sample, males, and females
matrix weighted_results_whole = J(`wordcount',11,.)
matrix weighted_results_male = J(`wordcount', 11, .)
matrix weighted_results_female = J(`wordcount', 11, .)
matrix weighted_ttest_results = J(`wordcount', 3, .)

matrix colnames weighted_results_whole = "Mean" "SE" "SD" "Median" "Min" "Max" "Skewness" "Kurtosis" "N" "Missing" "%_Missing"
matrix colnames weighted_results_male = "Mean" "SE" "SD" "Median" "Min" "Max" "Skewness" "Kurtosis" "N" "Missing" "%_Missing"
matrix colnames weighted_results_female = "Mean" "SE" "SD" "Median" "Min" "Max" "Skewness" "Kurtosis" "N" "Missing" "%_Missing"
matrix colnames weighted_ttest_results = "F_stat" "p_value" "df"

matrix rownames weighted_results_whole = `vars_demogr'
matrix rownames weighted_results_male = `vars_demogr'
matrix rownames weighted_results_female = `vars_demogr'
matrix rownames weighted_ttest_results = `vars_demogr'

* Initialize row counter
local row = 1

* Loop over demographic variables for weighted analysis
foreach var of local vars_demogr {
    
    display as text ""
    display as text "--- Weighted Analysis for `var' ---"
    
    * Check if variable has any non-missing values
    quietly count if !missing(`var') & !missing(w1respwt) & w1respwt > 0
    if r(N) == 0 {
        display as error "No valid observations for `var' with positive weights"
        local ++row
        continue
    }
    
    * Whole sample weighted statistics
    capture svy: mean `var'
    if _rc == 0 {
        matrix temp_mean = e(b)
        matrix temp_V = e(V)
        local weighted_mean = temp_mean[1,1]
        local weighted_se = sqrt(temp_V[1,1])
        
        * Store weighted mean and SE
        matrix weighted_results_whole[`row', 1] = `weighted_mean'
        matrix weighted_results_whole[`row', 2] = `weighted_se'
        
        * Weighted SD
        quietly estat sd
        matrix weighted_results_whole[`row', 3] = r(sd)
    }
    else {
        display as error "Error in weighted mean calculation for `var'"
        matrix weighted_results_whole[`row', 1] = .
        matrix weighted_results_whole[`row', 2] = .
        matrix weighted_results_whole[`row', 3] = .
    }
    
    * For median, min, max, skewness, kurtosis - use unweighted as survey commands don't provide these
    quietly sum `var', detail
    matrix weighted_results_whole[`row', 4] = r(p50)
    matrix weighted_results_whole[`row', 5] = r(min)
    matrix weighted_results_whole[`row', 6] = r(max)
    matrix weighted_results_whole[`row', 7] = r(skewness)
    matrix weighted_results_whole[`row', 8] = r(kurtosis)
    
    * Sample size and missing (unweighted counts)
    quietly count
    local total = r(N)
    matrix weighted_results_whole[`row', 9] = `total'
    quietly count if missing(`var')
    matrix weighted_results_whole[`row', 10] = r(N)
    local miss = r(N)
    matrix weighted_results_whole[`row', 11] = 100 * `miss' / `total'
    
    * Male weighted statistics (sex == 1)
    capture svy, subpop(if sex == 1): mean `var'
    if _rc == 0 & e(N_sub) > 0 {
        matrix temp_mean_m = e(b)
        matrix temp_V_m = e(V)
        local weighted_mean_m = temp_mean_m[1,1]
        local weighted_se_m = sqrt(temp_V_m[1,1])
        matrix weighted_results_male[`row', 1] = `weighted_mean_m'
        matrix weighted_results_male[`row', 2] = `weighted_se_m'
        
        * Weighted SD for males
        quietly estat sd
        matrix weighted_results_male[`row', 3] = r(sd)
    }
    else {
        matrix weighted_results_male[`row', 1] = .
        matrix weighted_results_male[`row', 2] = .
        matrix weighted_results_male[`row', 3] = .
    }
    
    * Male unweighted descriptives for other statistics
    quietly count if sex == 1
    local total_male = r(N)
    quietly sum `var' if sex == 1, detail
    matrix weighted_results_male[`row', 4] = r(p50)
    matrix weighted_results_male[`row', 5] = r(min)
    matrix weighted_results_male[`row', 6] = r(max)
    matrix weighted_results_male[`row', 7] = r(skewness)
    matrix weighted_results_male[`row', 8] = r(kurtosis)
    matrix weighted_results_male[`row', 9] = `total_male'
    quietly count if missing(`var') & sex == 1
    matrix weighted_results_male[`row', 10] = r(N)
    local miss_male = r(N)
    matrix weighted_results_male[`row', 11] = cond(`total_male' > 0, 100 * `miss_male' / `total_male', .)
    
    * Female weighted statistics (sex == 2)
    capture svy, subpop(if sex == 2): mean `var'
    if _rc == 0 & e(N_sub) > 0 {
        matrix temp_mean_f = e(b)
        matrix temp_V_f = e(V)
        local weighted_mean_f = temp_mean_f[1,1]
        local weighted_se_f = sqrt(temp_V_f[1,1])
        matrix weighted_results_female[`row', 1] = `weighted_mean_f'
        matrix weighted_results_female[`row', 2] = `weighted_se_f'
        
        * Weighted SD for females
        quietly estat sd
        matrix weighted_results_female[`row', 3] = r(sd)
    }
    else {
        matrix weighted_results_female[`row', 1] = .
        matrix weighted_results_female[`row', 2] = .
        matrix weighted_results_female[`row', 3] = .
    }
    
    * Female unweighted descriptives for other statistics
    quietly count if sex == 2
    local total_female = r(N)
    quietly sum `var' if sex == 2, detail
    matrix weighted_results_female[`row', 4] = r(p50)
    matrix weighted_results_female[`row', 5] = r(min)
    matrix weighted_results_female[`row', 6] = r(max)
    matrix weighted_results_female[`row', 7] = r(skewness)
    matrix weighted_results_female[`row', 8] = r(kurtosis)
    matrix weighted_results_female[`row', 9] = `total_female'
    quietly count if missing(`var') & sex == 2
    matrix weighted_results_female[`row', 10] = r(N)
    local miss_female = r(N)
    matrix weighted_results_female[`row', 11] = cond(`total_female' > 0, 100 * `miss_female' / `total_female', .)
    
    * Weighted test of difference between sexes
    capture svy: regress `var' i.sex
    if _rc == 0 {
        capture test 2.sex
        if _rc == 0 {
            matrix weighted_ttest_results[`row', 1] = r(F)
            matrix weighted_ttest_results[`row', 2] = r(p)
            matrix weighted_ttest_results[`row', 3] = r(df_r)
            
            display as text "Weighted F-test for `var': F(" r(df) "," r(df_r) ") = " %6.3f r(F) ", p = " %8.6f r(p)
        }
        else {
            display as error "Error in F-test for `var'"
            matrix weighted_ttest_results[`row', 1] = .
            matrix weighted_ttest_results[`row', 2] = .
            matrix weighted_ttest_results[`row', 3] = .
        }
    }
    else {
        display as error "Error in weighted regression for `var'"
        matrix weighted_ttest_results[`row', 1] = .
        matrix weighted_ttest_results[`row', 2] = .
        matrix weighted_ttest_results[`row', 3] = .
    }
    
    local ++row
}

* Display weighted results
display as text ""
display as text "=== WEIGHTED DESCRIPTIVE STATISTICS ==="
display as text "=== Weighted Results Whole Sample ==="
matrix list weighted_results_whole, format(%9.3f)

display as text ""
display as text "=== Weighted Results for Males (sex = 1) ==="
matrix list weighted_results_male, format(%9.3f)

display as text ""
display as text "=== Weighted Results for Females (sex = 2) ==="
matrix list weighted_results_female, format(%9.3f)

display as text ""
display as text "=== Weighted F-tests (Male vs Female) ==="
matrix list weighted_ttest_results, format(%12.6f)

* Display significant differences with error checking
display as text ""
display as text "=== Weighted F-test P-values (Detailed) ==="
local row = 1
foreach var of local vars_demogr {
    scalar pval_scalar = weighted_ttest_results[`row', 2]
    scalar fstat_scalar = weighted_ttest_results[`row', 1]
    
    if !missing(pval_scalar) & !missing(fstat_scalar) {
        local pval = scalar(pval_scalar)
        local fstat = scalar(fstat_scalar)
        local significance = cond(`pval' < 0.001, "***", cond(`pval' < 0.01, "**", cond(`pval' < 0.05, "*", "")))
        display as text "`var': F = " %8.3f `fstat' ", p = " %8.6f `pval' " `significance'"
    }
    else {
        display as text "`var': F-test failed or missing data"
    }
    local ++row
}

*-----------------------------------------------------------------------
* 4.2 Diagnostic checks for survey weights
*-----------------------------------------------------------------------
svyset [pweight=w1respwt]

display as text ""
display as text "=== SURVEY WEIGHT DIAGNOSTICS ==="

* Check weight distribution
summarize w1respwt, detail
display as text "Weight summary statistics:"
display as text "  Mean: " %8.3f r(mean)
display as text "  Min: " %8.3f r(min) 
display as text "  Max: " %8.3f r(max)
display as text "  Zero weights: " 
count if w1respwt == 0
display r(N) " observations"

* Check for missing weights
count if missing(w1respwt)
display as text "Missing weights: " r(N) " observations"

* Effective sample size
quietly svy: mean w1fi
display as text "Effective sample size: " e(N_pop)
display as text "Actual sample size: " e(N)


* Export weighted results to Excel
local filename_weighted "MHAS_descriptive_analysis_weighted"
putexcel set "$outputtables/`filename_weighted'.xlsx", sheet("Weighted_Whole_Sample") replace
putexcel A1 = matrix(weighted_results_whole), names

putexcel set "$outputtables/`filename_weighted'.xlsx", sheet("Weighted_Male") modify
putexcel A1 = matrix(weighted_results_male), names

putexcel set "$outputtables/`filename_weighted'.xlsx", sheet("Weighted_Female") modify
putexcel A1 = matrix(weighted_results_female), names

putexcel set "$outputtables/`filename_weighted'.xlsx", sheet("Weighted_Ftests") modify
putexcel A1 = matrix(weighted_ttest_results), names

*----------------------------------------------------------------------
* 4.1 Extracting Weighted Categorical Variables Analysis
*----------------------------------------------------------------------

*----------------------------------------------------------------------
* 4.1 Weighted Categorical Variables Analysis - Whole Sample
*----------------------------------------------------------------------

* Ensure survey design is properly set
svyset [pweight=w1respwt]

* Define categorical variables list
local vars_cat adl educ vigact rural smoke_status smokev smoke srh drink drinkdw ///
              drink_status higov diabe bmicat w1fc cogstat wealth_cat

* Create temporary file to store results
tempname catresults
tempfile categorical_results
postfile `catresults' str20 variable str20 category ///
    int sample_n double weighted_pct double missing_pct str30 n_wpct ///
    using `categorical_results', replace

* Loop through each categorical variable
foreach var of local vars_cat {
    
    display as text ""
    display as text "=== ANALYZING VARIABLE: `var' ==="
    
    * Get variable label
    local varlabel : variable label `var'
    if "`varlabel'" == "" local varlabel "`var'"
    display as text "Variable label: `varlabel'"
    
    * Basic diagnostics
    quietly count
    local total_obs = r(N)
    quietly count if !missing(`var')
    local valid_obs = r(N)
    quietly count if missing(`var')
    local missing_obs = r(N)
    local missing_pct = 100 * `missing_obs' / `total_obs'
    
    display as text "Total observations: " `total_obs'
    display as text "Valid observations: " `valid_obs'
    display as text "Missing observations: " `missing_obs' " (" %4.1f `missing_pct' "%)"
    
    * Get all levels (excluding missing for weighted calculations)
    quietly levelsof `var' if !missing(`var'), local(levels)
    
    if "`levels'" != "" {
        display as text "Categories found: `levels'"
        
        * Calculate weighted percentages for each category
        foreach level of local levels {
            
            * Get category label
            local level_label : label (`var') `level'
            if "`level_label'" == "" local level_label "`level'"
            
            * Get sample size (unweighted N)
            quietly count if `var' == `level'
            local sample_n = r(N)
            
            if `sample_n' > 0 {
                * Create indicator variable for this category (including all observations)
                quietly gen temp_indicator = (`var' == `level')
                
                * Calculate weighted proportion (of total sample, not just valid)
                capture svy: mean temp_indicator
                if _rc == 0 {
                    matrix temp_prop = e(b)
                    local weighted_prop = temp_prop[1,1]
                    local weighted_pct = `weighted_prop' * 100
                    
                    * Format combined N(%) with commas
                    local n_wpct = string(`sample_n', "%12.0gc") + " (" + string(`weighted_pct', "%5.2f") + ")"
                    
                    * Post results
                    post `catresults' ("`var'") ("`level_label'") ///
                        (`sample_n') (`weighted_pct') (`missing_pct') ("`n_wpct'")
                    
                    display as text "  `level_label': N=" string(`sample_n', "%12.0gc") " (Weighted %: " %5.2f `weighted_pct' ")"
                }
                else {
                    display as error "  Failed to calculate weighted percentage for `level_label'"
                }
                
                quietly drop temp_indicator
            }
            else {
                display as text "  `level_label': N=0 (no observations)"
            }
        }
        
        * Add missing values entry if any
        if `missing_obs' > 0 {
            * Calculate weighted percentage for missing values
            quietly gen temp_missing = missing(`var')
            capture svy: mean temp_missing
            if _rc == 0 {
                matrix temp_prop_miss = e(b)
                local weighted_prop_miss = temp_prop_miss[1,1]
                local weighted_pct_miss = `weighted_prop_miss' * 100
                
                local missing_formatted = string(`missing_obs', "%12.0gc") + " (" + string(`weighted_pct_miss', "%5.2f") + ")"
                post `catresults' ("`var'") ("Missing") ///
                    (`missing_obs') (`weighted_pct_miss') (.) ("`missing_formatted'")
                display as text "  Missing: N=" string(`missing_obs', "%12.0gc") " (Weighted %: " %4.2f `weighted_pct_miss' ")"
            }
            else {
                * Fallback to unweighted percentage if weighted fails
                local missing_formatted = string(`missing_obs', "%12.0gc") + " (" + string(`missing_pct', "%5.2f") + ")"
                post `catresults' ("`var'") ("Missing") ///
                    (`missing_obs') (`missing_pct') (.) ("`missing_formatted'")
                display as text "  Missing: N=" string(`missing_obs', "%12.0gc") " (Unweighted %: " %4.2f `missing_pct' ")"
            }
            quietly drop temp_missing
        }
    }
    else {
        display as error "No valid categories found for `var'"
        * Still post missing data if present
        if `missing_obs' > 0 {
            quietly gen temp_missing = missing(`var')
            capture svy: mean temp_missing
            if _rc == 0 {
                matrix temp_prop_miss = e(b)
                local weighted_prop_miss = temp_prop_miss[1,1]
                local weighted_pct_miss = `weighted_prop_miss' * 100
                
                local missing_formatted = string(`missing_obs', "%12.0gc") + " (" + string(`weighted_pct_miss', "%5.2f") + ")"
                post `catresults' ("`var'") ("Missing") ///
                    (`missing_obs') (`weighted_pct_miss') (.) ("`missing_formatted'")
            }
            else {
                local missing_formatted = string(`missing_obs', "%12.0gc") + " (" + string(`missing_pct', "%5.2f") + ")"
                post `catresults' ("`var'") ("Missing") ///
                    (`missing_obs') (`missing_pct') (.) ("`missing_formatted'")
            }
            quietly drop temp_missing
        }
    }
    
    display as text "----------------------------------------"
}

postclose `catresults'

*----------------------------------------------------------------------
* 4.2 Weighted Categorical Variables Analysis - Males Only (sex = 1)
*----------------------------------------------------------------------

* Create temporary file to store results for males
tempname catresults_male
tempfile categorical_results_male
postfile `catresults_male' str20 variable str20 category ///
    int sample_n double weighted_pct double missing_pct str30 n_wpct ///
    using `categorical_results_male', replace

* Loop through each categorical variable for males only
foreach var of local vars_cat {
    
    display as text ""
    display as text "=== ANALYZING VARIABLE: `var' (MALES ONLY) ==="
    
    * Get variable label
    local varlabel : variable label `var'
    if "`varlabel'" == "" local varlabel "`var'"
    display as text "Variable label: `varlabel'"
    
    * Basic diagnostics for males only
    quietly count if sex == 1
    local total_obs_male = r(N)
    quietly count if !missing(`var') & sex == 1
    local valid_obs_male = r(N)
    quietly count if missing(`var') & sex == 1
    local missing_obs_male = r(N)
    local missing_pct_male = 100 * `missing_obs_male' / `total_obs_male'
    
    display as text "Total male observations: " `total_obs_male'
    display as text "Valid male observations: " `valid_obs_male'
    display as text "Missing male observations: " `missing_obs_male' " (" %4.1f `missing_pct_male' "%)"
    
    * Get all levels for males (excluding missing for weighted calculations)
    quietly levelsof `var' if !missing(`var') & sex == 1, local(levels_male)
    
    if "`levels_male'" != "" {
        display as text "Categories found in males: `levels_male'"
        
        * Calculate weighted percentages for each category in males
        foreach level of local levels_male {
            
            * Get category label
            local level_label : label (`var') `level'
            if "`level_label'" == "" local level_label "`level'"
            
            * Get sample size (unweighted N) for males
            quietly count if `var' == `level' & sex == 1
            local sample_n_male = r(N)
            
            if `sample_n_male' > 0 {
                * Create indicator variable for this category in males
                quietly gen temp_indicator_male = (`var' == `level') if sex == 1
                
                * Calculate weighted proportion using subpopulation estimation
                capture svy, subpop(if sex == 1): mean temp_indicator_male
                if _rc == 0 & e(N_sub) > 0 {
                    matrix temp_prop_male = e(b)
                    local weighted_prop_male = temp_prop_male[1,1]
                    local weighted_pct_male = `weighted_prop_male' * 100
                    
                    * Format combined N(%) with commas
                    local n_wpct_male = string(`sample_n_male', "%12.0gc") + " (" + string(`weighted_pct_male', "%5.2f") + ")"
                    
                    * Post results
                    post `catresults_male' ("`var'") ("`level_label'") ///
                        (`sample_n_male') (`weighted_pct_male') (`missing_pct_male') ("`n_wpct_male'")
                    
                    display as text "  `level_label': N=" string(`sample_n_male', "%12.0gc") " (Weighted %: " %5.2f `weighted_pct_male' ")"
                }
                else {
                    display as error "  Failed to calculate weighted percentage for `level_label' in males"
                }
                
                quietly drop temp_indicator_male
            }
            else {
                display as text "  `level_label': N=0 (no male observations)"
            }
        }
        
        * Add missing values entry for males if any
        if `missing_obs_male' > 0 {
            * Calculate weighted percentage for missing values in males
            quietly gen temp_missing_male = missing(`var') if sex == 1
            capture svy, subpop(if sex == 1): mean temp_missing_male
            if _rc == 0 & e(N_sub) > 0 {
                matrix temp_prop_miss_male = e(b)
                local weighted_prop_miss_male = temp_prop_miss_male[1,1]
                local weighted_pct_miss_male = `weighted_prop_miss_male' * 100
                
                local missing_formatted_male = string(`missing_obs_male', "%12.0gc") + " (" + string(`weighted_pct_miss_male', "%5.2f") + ")"
                post `catresults_male' ("`var'") ("Missing") ///
                    (`missing_obs_male') (`weighted_pct_miss_male') (.) ("`missing_formatted_male'")
                display as text "  Missing: N=" string(`missing_obs_male', "%12.0gc") " (Weighted %: " %4.2f `weighted_pct_miss_male' ")"
            }
            else {
                * Fallback to unweighted percentage if weighted fails
                local missing_formatted_male = string(`missing_obs_male', "%12.0gc") + " (" + string(`missing_pct_male', "%5.2f") + ")"
                post `catresults_male' ("`var'") ("Missing") ///
                    (`missing_obs_male') (`missing_pct_male') (.) ("`missing_formatted_male'")
                display as text "  Missing: N=" string(`missing_obs_male', "%12.0gc") " (Unweighted %: " %4.2f `missing_pct_male' ")"
            }
            quietly drop temp_missing_male
        }
    }
    else {
        display as error "No valid categories found for `var' in males"
        * Still post missing data if present
        if `missing_obs_male' > 0 {
            quietly gen temp_missing_male = missing(`var') if sex == 1
            capture svy, subpop(if sex == 1): mean temp_missing_male
            if _rc == 0 & e(N_sub) > 0 {
                matrix temp_prop_miss_male = e(b)
                local weighted_prop_miss_male = temp_prop_miss_male[1,1]
                local weighted_pct_miss_male = `weighted_prop_miss_male' * 100
                
                local missing_formatted_male = string(`missing_obs_male', "%12.0gc") + " (" + string(`weighted_pct_miss_male', "%5.2f") + ")"
                post `catresults_male' ("`var'") ("Missing") ///
                    (`missing_obs_male') (`weighted_pct_miss_male') (.) ("`missing_formatted_male'")
            }
            else {
                local missing_formatted_male = string(`missing_obs_male', "%12.0gc") + " (" + string(`missing_pct_male', "%5.2f") + ")"
                post `catresults_male' ("`var'") ("Missing") ///
                    (`missing_obs_male') (`missing_pct_male') (.) ("`missing_formatted_male'")
            }
            quietly drop temp_missing_male
        }
    }
    
    display as text "----------------------------------------"
}

postclose `catresults_male'

*----------------------------------------------------------------------
* 4.3 Weighted Categorical Variables Analysis - Females Only (sex = 2)
*----------------------------------------------------------------------

* Create temporary file to store results for females
tempname catresults_female
tempfile categorical_results_female
postfile `catresults_female' str20 variable str20 category ///
    int sample_n double weighted_pct double missing_pct str30 n_wpct ///
    using `categorical_results_female', replace

* Loop through each categorical variable for females only
foreach var of local vars_cat {
    
    display as text ""
    display as text "=== ANALYZING VARIABLE: `var' (FEMALES ONLY) ==="
    
    * Get variable label
    local varlabel : variable label `var'
    if "`varlabel'" == "" local varlabel "`var'"
    display as text "Variable label: `varlabel'"
    
    * Basic diagnostics for females only
    quietly count if sex == 2
    local total_obs_female = r(N)
    quietly count if !missing(`var') & sex == 2
    local valid_obs_female = r(N)
    quietly count if missing(`var') & sex == 2
    local missing_obs_female = r(N)
    local missing_pct_female = 100 * `missing_obs_female' / `total_obs_female'
    
    display as text "Total female observations: " `total_obs_female'
    display as text "Valid female observations: " `valid_obs_female'
    display as text "Missing female observations: " `missing_obs_female' " (" %4.1f `missing_pct_female' "%)"
    
    * Get all levels for females (excluding missing for weighted calculations)
    quietly levelsof `var' if !missing(`var') & sex == 2, local(levels_female)
    
    if "`levels_female'" != "" {
        display as text "Categories found in females: `levels_female'"
        
        * Calculate weighted percentages for each category in females
        foreach level of local levels_female {
            
            * Get category label
            local level_label : label (`var') `level'
            if "`level_label'" == "" local level_label "`level'"
            
            * Get sample size (unweighted N) for females
            quietly count if `var' == `level' & sex == 2
            local sample_n_female = r(N)
            
            if `sample_n_female' > 0 {
                * Create indicator variable for this category in females
                quietly gen temp_indicator_female = (`var' == `level') if sex == 2
                
                * Calculate weighted proportion using subpopulation estimation
                capture svy, subpop(if sex == 2): mean temp_indicator_female
                if _rc == 0 & e(N_sub) > 0 {
                    matrix temp_prop_female = e(b)
                    local weighted_prop_female = temp_prop_female[1,1]
                    local weighted_pct_female = `weighted_prop_female' * 100
                    
                    * Format combined N(%) with commas
                    local n_wpct_female = string(`sample_n_female', "%12.0gc") + " (" + string(`weighted_pct_female', "%5.2f") + ")"
                    
                    * Post results
                    post `catresults_female' ("`var'") ("`level_label'") ///
                        (`sample_n_female') (`weighted_pct_female') (`missing_pct_female') ("`n_wpct_female'")
                    
                    display as text "  `level_label': N=" string(`sample_n_female', "%12.0gc") " (Weighted %: " %5.2f `weighted_pct_female' ")"
                }
                else {
                    display as error "  Failed to calculate weighted percentage for `level_label' in females"
                }
                
                quietly drop temp_indicator_female
            }
            else {
                display as text "  `level_label': N=0 (no female observations)"
            }
        }
        
        * Add missing values entry for females if any
        if `missing_obs_female' > 0 {
            * Calculate weighted percentage for missing values in females
            quietly gen temp_missing_female = missing(`var') if sex == 2
            capture svy, subpop(if sex == 2): mean temp_missing_female
            if _rc == 0 & e(N_sub) > 0 {
                matrix temp_prop_miss_female = e(b)
                local weighted_prop_miss_female = temp_prop_miss_female[1,1]
                local weighted_pct_miss_female = `weighted_prop_miss_female' * 100
                
                local missing_formatted_female = string(`missing_obs_female', "%12.0gc") + " (" + string(`weighted_pct_miss_female', "%5.2f") + ")"
                post `catresults_female' ("`var'") ("Missing") ///
                    (`missing_obs_female') (`weighted_pct_miss_female') (.) ("`missing_formatted_female'")
                display as text "  Missing: N=" string(`missing_obs_female', "%12.0gc") " (Weighted %: " %4.2f `weighted_pct_miss_female' ")"
            }
            else {
                * Fallback to unweighted percentage if weighted fails
                local missing_formatted_female = string(`missing_obs_female', "%12.0gc") + " (" + string(`missing_pct_female', "%5.2f") + ")"
                post `catresults_female' ("`var'") ("Missing") ///
                    (`missing_obs_female') (`missing_pct_female') (.) ("`missing_formatted_female'")
                display as text "  Missing: N=" string(`missing_obs_female', "%12.0gc") " (Unweighted %: " %4.2f `missing_pct_female' ")"
            }
            quietly drop temp_missing_female
        }
    }
    else {
        display as error "No valid categories found for `var' in females"
        * Still post missing data if present
        if `missing_obs_female' > 0 {
            quietly gen temp_missing_female = missing(`var') if sex == 2
            capture svy, subpop(if sex == 2): mean temp_missing_female
            if _rc == 0 & e(N_sub) > 0 {
                matrix temp_prop_miss_female = e(b)
                local weighted_prop_miss_female = temp_prop_miss_female[1,1]
                local weighted_pct_miss_female = `weighted_prop_miss_female' * 100
                
                local missing_formatted_female = string(`missing_obs_female', "%12.0gc") + " (" + string(`weighted_pct_miss_female', "%5.2f") + ")"
                post `catresults_female' ("`var'") ("Missing") ///
                    (`missing_obs_female') (`weighted_pct_miss_female') (.) ("`missing_formatted_female'")
            }
            else {
                local missing_formatted_female = string(`missing_obs_female', "%12.0gc") + " (" + string(`missing_pct_female', "%5.2f") + ")"
                post `catresults_female' ("`var'") ("Missing") ///
                    (`missing_obs_female') (`missing_pct_female') (.) ("`missing_formatted_female'")
            }
            quietly drop temp_missing_female
        }
    }
    
    display as text "----------------------------------------"
}

postclose `catresults_female'

*----------------------------------------------------------------------
* Display and Export Results for All Three Groups
*----------------------------------------------------------------------

* Export Whole Sample Results
preserve
use `categorical_results', clear
format weighted_pct %6.2f
format missing_pct %6.2f

display as text ""
display as text "=== WEIGHTED CATEGORICAL ANALYSIS SUMMARY - WHOLE SAMPLE ==="
list variable category sample_n weighted_pct missing_pct, sepby(variable) abbreviate(15)

save "$outputtables/weighted_categorical_whole_sample.dta", replace
export delimited "$outputtables/weighted_categorical_whole_sample.csv", replace
export excel using "$outputtables/weighted_categorical_whole_sample.xlsx", ///
    sheet("Whole_Sample_Summary") firstrow(variables) replace
restore

* Export Male Results
preserve
use `categorical_results_male', clear
format weighted_pct %6.2f
format missing_pct %6.2f

display as text ""
display as text "=== WEIGHTED CATEGORICAL ANALYSIS SUMMARY - MALES ONLY ==="
list variable category sample_n weighted_pct missing_pct, sepby(variable) abbreviate(15)

save "$outputtables/weighted_categorical_males_only.dta", replace
export delimited "$outputtables/weighted_categorical_males_only.csv", replace
export excel using "$outputtables/weighted_categorical_males_only.xlsx", ///
    sheet("Males_Only_Summary") firstrow(variables) replace
restore

* Export Female Results
preserve
use `categorical_results_female', clear
format weighted_pct %6.2f
format missing_pct %6.2f

display as text ""
display as text "=== WEIGHTED CATEGORICAL ANALYSIS SUMMARY - FEMALES ONLY ==="
list variable category sample_n weighted_pct missing_pct, sepby(variable) abbreviate(15)

save "$outputtables/weighted_categorical_females_only.dta", replace
export delimited "$outputtables/weighted_categorical_females_only.csv", replace
export excel using "$outputtables/weighted_categorical_females_only.xlsx", ///
    sheet("Females_Only_Summary") firstrow(variables) replace
restore

display as text ""
display as text "ALL WEIGHTED CATEGORICAL ANALYSES COMPLETE!"
display as text ""
display as text "=== SUMMARY OF EXPORTED FILES ==="
display as text "Whole Sample:"
display as text "  - $outputtables/weighted_categorical_whole_sample.xlsx"
display as text ""
display as text "Males Only:"
display as text "  - $outputtables/weighted_categorical_males_only.xlsx"
display as text ""
display as text "Females Only:"
display as text "  - $outputtables/weighted_categorical_females_only.xlsx"

*===============================================================================
* 5. Session end
*===============================================================================

display as text "------------------------------------------------------------"
display as text "MHAs weighted descriptives analysis completed"
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close