****************************************************************************
* Project: MHAS frailty trajectories
* Author:  Aleixo Brunetti
* Purpose: MHAS data cleaning
****************************************************************************
* ==============================================================================
* 3. Data Cleaning & Exploring Wave 1 - Baseline
* ==============================================================================
set scheme stcolor

* -- Load the data
use "mhas_wave1_2001.dta", clear

codebook, compact
codebook, problems

* -- Decision to delete the variables that are not needed for the analysis:

* Education variables 
codebook raedyrs raeducl, compact
list patid raedyrs raeducl in 1/20
drop raedyrs
// Harmonised education (raeducl) retained

* ADL variables:
codebook r1adla r1walkr r1mobila r1lowermoba r1uppermoba r1adlfive r1adlwa r1adltot6, compact
list patid r1adla r1walkr r1mobila r1lowermoba r1uppermoba r1adlfive r1adlwa r1adltot6 in 1/20
drop r1adla r1walkr r1mobila r1lowermoba r1uppermoba r1adlfive r1adlwa
// r1adltot6 retained: RwADLTOT6 = sum(RwBATHA, RwDRESSA, RwEATA, RwBEDA, RwWALKRA, RwTOILTA)

* Wealth variables:
codebook h1atotb h1atotf
list patid h1atotb h1atotf in 1/20
drop h1atotf
// Total wealth (h1atotb) retained

* Alcohol consumption variables:
codebook r1drink r1drinkd r1drinkb r1cage
list patid r1drink r1drinkd r1drinkb r1cage in 1/20
list patid r1drink r1drinkd r1drinkb r1cage if r1cage == 4
// Current drinking status (r1drink) and days per week (r1drinkd) retained
drop r1drinkb r1cage

* Health care coverage variables:
codebook r1higov r1covr_m
list patid r1higov r1covr_m 
tab r1higov r1covr_m
// r1covr_m adds no information beyond r1higov; r1higov retained
drop r1covr_m

* Smoking variables:
list r1smokev r1smoken r1quitsmok r1strtsmok in 1/20
drop r1strtsmok r1quitsmok
// r1smokev and r1smoken retained to derive never/former/current smoking

* Participants weight calibaration variable:
list patid r1wtresp in 1/10 
codebook r1wtresp
tab r1wtresp if r1wtresp ==0
label variable r1wtresp "Wave 1 R Personal-level Calibration Weight (0=Not age-eligible)"
// MHAS calibration weight = 0 when respondent is not age-eligible (born 1951 or later at wave 1)

* -- Rename variables to make them more consistent and descriptive

* Participants weight calibaration variable:
rename r1wtresp w1respwt
label variable w1respwt "Wave 1 R Personal-level Calibration Weight (0=Not age-eligible)"

* Demographics
rename r1agey age
label variable age "age (years) at ivw wave 1"

rename rabyear byear
label variable byear "R Birth Year (years)"

rename rabmonth bmonth
label variable bmonth "R Birth Month"

rename radyear dyear
label variable dyear "R death year (Next-of-Kin)"

rename radmonth dmonth
label variable dmonth "R Death Month (Next-of-Kin)"

rename ragender sex
label variable sex "W1 R sex 1=Male & 2=Female"
label define sex 1 "1.Male" 2 "2.Female", replace
label values sex sex

rename raeducl educ
label variable educ "W1R Hamonized Education (1 < up sec,  2 = Up sec and voc & 3 = tertiary)"

rename h1rural rural
label variable rural "W1 R live in 0 = urban & 1= rural"

* Body Measurements
rename r1mbmi bmi
label variable bmi "W1 R Body Mass Index (BMI) at wave 1"

rename r1mbmicat bmicat
label variable bmicat "W1 R Body Mass Index (BMI) Categorization at wave 1 (1=Underweight, 2=Normal, 3=Pre-obesity, 4=Obesity class 1, 5=Obesity class 2, 6=Obesity class 3)"

* Activity of Daily Living (ADL)
rename r1adltot6 adl
label variable adl "W1 R ADL Total Score 6 (0-6)"

* Health behaviour & Lifestyle
rename r1vigact vigact
label variable vigact "W1 R Vigorous Activity (1=Yes & 0=No)"

rename r1drink drink
label variable drink "W1 R Alcohol Consumption Status (1=Yes & 0=No)"

rename r1drinkd drinkdw
label variable drinkdw "W1 R Number of Days per Week (0-7)" 

rename r1smokev smokev
label variable smokev "W1 R Ever Smoked (1=Yes & 0=No)"

rename r1smoken smoke
label variable smoke "W1 R Smoking Status Now (1= Yes & 0=No)"

* Health care coverage
rename r1higov higov
label variable higov "W1 R Health Insurance Coverage (1=Yes & 0=No)"

* Wealth
rename h1atotb wealth
label variable wealth "W1 R Total Wealth (in 1000s of pesos)"

* Health
rename r1shlt srh
label variable srh "W1 R self-report of health (1= Ex,  2=VG, 3=G, 4=F & 5=P)"

rename r1diabe diabe
label variable diabe "W1 R Diabetes 1=Yes & 0=No"

rename cognitive_status_01 cogstat
label variable cogstat "Cognitive Status (0 = Normal, 1 = Normal with Instrumental Impairment, 2 = CIND, 3 = Dementia)"

* Frailty Index_ Baseline
rename r1frailty_10 w1fi
label variable w1fi "W1 R Frailty Index"

rename r1frailty_10_cat w1fc
label variable w1fc "W1 R Frailty Category 0=Not Frail & 1=Frail"

misstable sum
tab dyear, missing
tab dmonth, missing
tab bmonth, missing
tab educ, missing
tab w1fi, missing 
tab w1fc, missing

misstable nested
list w1fi w1fc if missing(w1fi) | missing(w1fc)

*-------------------------------------------------------------------------------
* 3.1 Inspect baseline labels
*-------------------------------------------------------------------------------
label list
// Recode 1/2 coded variables to 0/1

replace vigact = 0 if vigact == 1
replace vigact = 1 if vigact == 2
label define vigact 0 "0.No" 1 "1.Yes"
label values vigact vigact

replace rural = 0 if rural == 1
replace rural = 1 if rural == 2
label define rural 0 "0.No" 1 "1.Yes"
label values rural rural

replace smoke = 0 if smoke == 1
replace smoke = 1 if smoke == 2
label define smoke 0 "0.No" 1 "1.Yes"
label values smoke smoke

replace smokev = 0 if smokev == 1
replace smokev = 1 if smokev == 2
label define smokev 0 "0.No" 1 "1.Yes"
label values smokev smokev

replace drink = 0 if drink == 1
replace drink = 1 if drink == 2
label define drink 0 "0.No" 1 "1.Yes"
label values drink drink

replace higov = 0 if higov == 1
replace higov = 1 if higov == 2
label define higov 0 "0.No" 1 "1.Yes"
label values higov higov

replace diabe = 0 if diabe == 1
replace diabe = 1 if diabe == 2
label define diabe 0 "0.No" 1 "1.Yes"
label values diabe diabe

replace cogstat = 0 if cogstat == 1
replace cogstat = 1 if cogstat == 2
replace cogstat = 2 if cogstat == 3
replace cogstat = 3 if cogstat == 4
label define cogstat 0 "0.Normal" 1 "1.Normal with Instrumental Impairment" 2 "2.CIND" 3 "3.Dementia"
label values cogstat cogstat

replace w1fc = 0 if w1fc == 1
replace w1fc = 1 if w1fc == 2
label define w1fc 0 "0.Not frail" 1 "1.Frail"
label values w1fc w1fc

label list

inspect 

*===============================================================================
* 11 Check for errors and making changes in the BMI, Wealth and weight calibration variables
*===============================================================================
*--- Check BMI outliers
list patid bmi bmic if bmi > 70 & !missing(bmi)

list patid bmi if bmi <14 & !missing(bmi)

* drop implausible BMI values
drop if (bmi < 14 | bmi > 70) & !missing(bmi)   // implausible BMI (Freedman et al.)

* recode BMI values as missing
tab bmi, missing
replace bmi = . if bmi ==.m
replace bmi = . if bmi ==.n
replace bmi = . if bmi ==.r
replace bmi = . if bmi ==.x
tab bmi, missing
// .s retained: participants not part of the anthropometric sub-sample

sum bmi
qnorm bmi
histogram bmi, normal

*-- Check Wealth outliers
sum wealth, detail
local rate 25.17
// Exchange rate on 10 June 2025: 1 GBP = 25.17 MXN
gen wealth_gbp = wealth / `rate'
label variable wealth_gbp "W1 R Total Wealth (in 1000s of GBP)"
sum wealth_gbp, detail
qnorm wealth_gbp
histogram wealth_gbp, normal
list patid wealth_gbp if wealth_gbp > 1036000 & !missing(wealth_gbp)

* -- Weight calibration variable
list patid w1respwt in 1/10
codebook w1respwt
tab w1respwt if w1respwt ==0
sum age if age < 50
sum age w1fc w1respwt if age < 50 & w1respwt ==0
list patid age w1fc w1respwt if age < 50 & w1respwt ==0 in 1/1000
count if age >= 50 & w1respwt ==0 
drop if w1respwt == 0
// Participants not age-eligible for population representation (weight = 0) excluded

*==============================================================================
* 4. Creating variables 
*==============================================================================

* -- Wealth categories:
* Wealth quintiles
tab wealth_gbp, missing
xtile wealth_q = wealth_gbp, nq(5)
label define wealth_lbl 1 "Q1.Poorest" 2 "Q2.Poorer" 3 "Q3.Middle" 4 "Q4.Richer" 5 "Q5.Richest"
label values wealth_q wealth_lbl
label variable wealth_q "Wealth quintiles (1=Q1.Poorest, 2=Q2.Poorer, 3=Q3.Middle, 4=Q4.Richer, 5=Q5.Richest)"
* Wealth categories
gen wealth_cat = .
replace wealth_cat = 1 if wealth_q == 1 | wealth_q == 2   // Q1 or Q2 → Low
replace wealth_cat = 2 if wealth_q == 3                   // Q3 → Middle
replace wealth_cat = 3 if wealth_q == 4 | wealth_q == 5   // Q4 or Q5 → High
replace wealth_cat = . if missing(wealth_q) // if wealth_q is missing
label define wealthcat_lbl 1 "1.Low wealth (Q1&Q2)" 2 "2.Middle wealth (Q3)" 3 "3.High wealth (Q4&Q5)"
label values wealth_cat wealthcat_lbl
label variable wealth_cat "Wealth category (1=Low wealth, 2=Middle wealth, 3=High wealth)"

list wealth wealth_cat wealth_q in 1/30
order wealth_cat wealth_q wealth_gbp, after (wealth)

* -- Smoking status
tab smoke, missing
tab smokev, missing
list patid smoke smokev if smoke ==0 & smokev==0 in 1/20
count if smoke ==0 & smokev==0

list patid smoke smokev if smoke ==0 & smokev>0 in 1/20
count if smoke ==0 & smokev>0

* -- Create a new variable for smoking status
gen smoke_status = .
replace smoke_status = 0 if smoke == 0 & (smokev == 0 | missing(smokev)) // 0.Never
replace smoke_status = 1 if (smoke == 0 & smokev > 0 & !missing(smoke) & !missing(smokev)) // 1.Former
replace smoke_status = 2 if smoke == 1 & !missing(smokev) // 2.Current smokers
replace smoke_status = . if missing(smoke) & missing(smokev) // if both are missing              

* -- Label the new variable
label define smokel_status 0 "0.Never" 1 "1.Former" 2 "2.Current Smokers"
label variable smoke_status "smoke status (0=Never, 1=Former, 2=Current Smokers)"
label values smoke_status smokel_status 

tab smoke_status, missing

* -- Drinking status

tab drink, missing
tab drinkdw, missing
list patid drink drinkdw if drink ==0 & drinkdw==0 in 1/20
list patid drink drinkdw if drink ==0 & drinkdw>0
list drink drinkdw if drink == 1 & drinkdw==0 in 1/20
count if drink == 1 & drinkdw==0
// Assumption: current drinkers drinking < 1 day per week classified as current drinkers

* -- Create a new variable for drink status
gen drink_status = .
replace drink_status = 0 if drink == 0 & (drinkdw == 0 | missing(drinkdw)) | (drink == 0 & drinkdw > 0 & !missing(drink) & !missing(drinkdw)) // 0.Never
replace drink_status = 1 if drink == 1 & !missing(drinkdw) // 1.Current drinkers
replace drink_status = . if missing(drink) & missing(drinkdw) // if both are missing
tab drink_status, missing

* -- Label the new variable 
label define drinkl_status 0 "0.Never" 1 "1.Current Drinkers"
label variable drink_status "drink status (0=Never, 1=Current Drinkers)"
label values drink_status drinkl_status

tab drink_status, missing

*-------------------------------------------------------------------------------
*4.1 Change special extended missing values to regular missing values
*-------------------------------------------------------------------------------
* Activities of Daily Living
replace adl = . if adl == .d
replace adl = . if adl == .m
replace adl = . if adl == .r

* ==============================================================================
* 5. Save Cleaned Data and Results Wave 1 -Baseline
* ==============================================================================
order w1fi w1fc, last

save "$cleanmhasdata/baseline_only.dta", replace

* ==============================================================================
* 6. Preparing & Exploring Waves 2-5
* ==============================================================================

* ------ 
* Wave 1 (cleaned) 
* ------

use "$cleanmhasdata/baseline_only.dta", clear
codebook patid
sum w1fi, d
gen wave = 1
tab w1fi, missing
tab w1fc, missing
* Males
sum w1fi if sex ==1, d
tab w1fc if sex ==1, missing
*Females
sum w1fi if sex ==2, d
tab w1fc if sex ==2, missing

* histogram for Wave 1

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

graph export "$outputgraphs/hist_w1fi.pdf", replace

save "$working/wave1_temp.dta", replace

* histogram for Wave 1 by sex
histogram w1fi, by(sex, note("")) frequency ///
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
graph export "$outputgraphs/hist_w1fi_by_sex.pdf", replace

* ------
* Wave 2 
* ------
use "$rawmhasdata/mhas_wave2_2003.dta", clear
gen wave = 2
rename r2frailty_10 w2fi
label variable w2fi "W2 R Frailty Index"
rename r2frailty_10_cat w2fc
label variable w2fc "W2 R Frailty Category 0=Not Frail & 1=Frail"

* -- missingness wave 2
codebook patid
misstable sum
misstable nested
sum w2fi, d
tab w2fi, missing
tab w2fc, missing 

* histogram for Wave 2

histogram w2fi, frequency ///
    title("Wave 2", size(medium)) ///
    xtitle("Frailty Index") ///
    ytitle("Frequency") ///
    bin(30) ///
    fcolor(eltblue) ///
    lcolor(black) ///
    lwidth(thin) ///
    graphregion(color(white)) ///
    plotregion(color(white)) ///
    scheme(s1mono)
    
graph export "$outputgraphs/hist_w2fi.pdf", replace

save "$working/wave2_temp.dta", replace

* ------
* Wave 3 
* ------

use "$rawmhasdata/mhas_wave3_2012.dta", clear
gen wave = 3
rename r3frailty_10 w3fi
label variable w3fi "W3 R Frailty Index"
rename r3frailty_10_cat w3fc
label variable w3fc "W3 R Frailty Category 0=Not Frail & 1=Frail"

* -- missingness Wave 3
codebook patid
misstable sum
misstable nested
sum w3fi, d
tab w3fi, missing
tab w3fc, missing

* histogram for Wave 3

histogram w3fi, frequency ///
    title("Wave 3", size(medium)) ///
    xtitle("Frailty Index") ///
    ytitle("Frequency") ///
    bin(30) ///
    fcolor(eltblue) ///
    lcolor(black) ///
    lwidth(thin) ///
    graphregion(color(white)) ///
    plotregion(color(white)) ///
    scheme(s1mono)

graph export "$outputgraphs/hist_w3fi.pdf", replace

save "$working/wave3_temp.dta", replace

* ------
* Wave 4 
* ------

use "$rawmhasdata/mhas_wave4_2015.dta", clear
gen wave = 4
rename r4frailty_10 w4fi
label variable w4fi "W4 R Frailty Index"
rename r4frailty_10_cat w4fc
label variable w4fc "W4 R Frailty Category 0=Not Frail & 1=Frail"

* -- missingness Wave 4
codebook patid
misstable sum
misstable nested
sum w4fi, d
tab w4fi, missing
tab w4fc, missing 

* histogram for Wave 4
histogram w4fi, frequency ///
    title("Wave 4", size(medium)) ///
    xtitle("Frailty Index") ///
    ytitle("Frequency") ///
    bin(30) ///
    fcolor(eltblue) ///
    lcolor(black) ///
    lwidth(thin) ///
    graphregion(color(white)) ///
    plotregion(color(white)) ///
    scheme(s1mono)
    
graph export "$outputgraphs/hist_w4fi.pdf", replace

save "$working/wave4_temp.dta", replace

* ------
* Wave 5
* ------

use "$rawmhasdata/mhas_wave5_2018.dta", clear
gen wave = 5
rename r5frailty_10 w5fi
label variable w5fi "W5 R Frailty Index"
rename r5frailty_10_cat w5fc
label variable w5fc "W5 R Frailty Category 0=Not Frail & 1=Frail"

* -- Missingness Wave 5
codebook patid
misstable sum
misstable nested
sum w5fi, detail
tab w5fi, missing
tab w5fc, missing 

* histogram for Wave 5

histogram w5fi, frequency ///
    title("Wave 5", size(medium)) ///
    xtitle("Frailty Index") ///
    ytitle("Frequency") ///
    bin(30) ///
    fcolor(eltblue) ///
    lcolor(black) ///
    lwidth(thin) ///
    graphregion(color(white)) ///
    plotregion(color(white)) ///
    scheme(s1mono)
    
graph export "$outputgraphs/hist_w5fi.pdf", replace

save "$working/wave5_temp.dta", replace

* ==============================================================================
* 7. MERGING MHAS Waves 1–5 wide formnat
* ==============================================================================

* Start with wave 1
use "$working/wave1_temp.dta", clear
duplicates report patid
duplicates report patid w1fi

* Merge with wave 2
merge 1:1 patid using "$working/wave2_temp.dta"

tab _merge
list patid w1fi w2fi _merge in 1/10 if _merge != 3
count if _merge == 1  // Only in wave 1
count if _merge == 2  // Only in wave 2  
count if _merge == 3  // Matched in both waves

* Check for unexpected patterns
list patid w1fi w2fi if _merge == 2 & !missing(w1fi)  // Should be empty
assert _merge != 2 | missing(w1fi)  // Participants only in wave 2 should not have w1fi

duplicates report patid w1fi w2fi
drop _merge

* Merge with wave 3
merge 1:1 patid using "$working/wave3_temp.dta"

tab _merge
count if _merge == 1  // In previous waves but not wave 3
count if _merge == 2  // Only in wave 3 (new entrants)
count if _merge == 3  // Matched

* Participants with w3fi but missing earlier waves
list patid w1fi w2fi w3fi if _merge == 2 & !missing(w3fi) in 1/10

duplicates report patid w1fi w2fi w3fi
drop _merge

* Merge with wave 4
merge 1:1 patid using "$working/wave4_temp.dta"

duplicates report patid w1fi w2fi w3fi w4fi
drop _merge

* Merge with wave 5
merge 1:1 patid using "$working/wave5_temp.dta"
drop _merge

* -------------------------------------------------------------------------------
* 7.1 Check the final merged dataset (error checking)
* -------------------------------------------------------------------------------

duplicates report patid
duplicates report patid w1fi w2fi w3fi w4fi w5fi 
duplicates tag patid w1fi w2fi w3fi w4fi w5fi, gen(_dup)
capture duplicates list patid if _dup > 0
drop _dup

order wave, before(w1fi)
summarize w*fi
list patid wave w1fi w2fi w3fi w4fi w5fi in 1/20

sort wave
label variable wave "First wave in which person entered (wave 1 to 5)"
list patid wave w1fi w2fi w3fi w4fi w5fi in 1/20 

tab wave
tab wave, missing
assert !missing(wave) 

* data integrity across waves
misstable patterns w*fi  // Missing patterns
misstable sum w*fi      // Missing counts per wave

* Checking for frailty index ranges
display as text "---- Frailty Index Range Validation ----"
forvalues w = 1/5 {
    quietly summarize w`w'fi
    if r(min) < 0 | r(max) > 1 {
        display as error "ERROR: w`w'fi has values outside [0,1] range"
        count if w`w'fi < 0 | w`w'fi > 1
        display "Number of invalid values: " r(N)
    }
    else {
        display as text "✓ w`w'fi values in valid range [0,1]"
    }
}

*===============================================================================
* 8. Saving baseline data
*===============================================================================

drop if wave >1
// Analytic cohort restricted to participants who entered at wave 1

*===============================================================================
* 9. Checking after deletion of wave > 1 participants errors
*===============================================================================

order cogstat w1fi w1fc w2fi w2fc w3fi w3fc w4fi w4fc w5fi w5fc wave, after(bmicat)

tab wave, missing
codebook patid, compact
count if missing(patid)
duplicates report patid
duplicates tag patid w1fi w2fi w3fi w4fi w5fi, gen(_dup)
capture duplicates list patid if _dup > 0
drop _dup

count if missing(w1fi)
tab wave 
codebook, compact
codebook, problems
list patid wave w*fi in 1/20

inspect 

*-------------------------------------------------------------------------------
* 9.1 Label inconsistent values
*-------------------------------------------------------------------------------

label list

replace w2fc = 0 if w2fc == 1
replace w2fc = 1 if w2fc == 2
label define w2fc 0 "0.Not frail" 1 "1.Frail"
label values w2fc w2fc

replace w3fc = 0 if w3fc == 1
replace w3fc = 1 if w3fc == 2
label define w3fc 0 "0.Not frail" 1 "1.Frail"
label values w3fc w3fc

replace w4fc = 0 if w4fc == 1
replace w4fc = 1 if w4fc == 2
label define w4fc 0 "0.Not frail" 1 "1.Frail"
label values w4fc w4fc

replace w5fc = 0 if w5fc == 1
replace w5fc = 1 if w5fc == 2
label define w5fc 0 "0.Not frail" 1 "1.Frail"
label values w5fc w5fc

label list

* ------------------------------------------------------------------------------
* 9.2  Completely missing data across 5 waves
* ------------------------------------------------------------------------------

tab dyear if missing(w1fi), missing

* Exploring the .s & missing(w1fi)
list patid dyear w1fi w2fi w3fi w4fi w5fi if dyear ==.s & missing(w1fi)

* Exploring among alive (.x) and missing (w1fi)
list patid dyear w1fi w2fi w3fi w4fi w5fi if dyear ==.x & missing(w1fi)

* -- Explore completely missing across 5 waves

gen completely_lost = missing(w1fi) + missing(w2fi) + missing(w3fi) + missing(w4fi) + missing(w5fi)
tab completely_lost
list w*fi completely_lost
tab completely_lost, missing

drop if completely_lost == 5
// Participants with no frailty index measurement in any wave excluded

tab completely_lost, missing

*==============================================================================
* 10 Save the cleaned wide dataset
*==============================================================================
* Save the wide dataset with all waves

order drink_status, after(drinkdw)
order smoke_status, after(smokev)

save "$cleanmhasdata/wave1to5_wide.dta", replace

*===============================================================================
* 11 Session end
*===============================================================================

display as text "------------------------------------------------------------"
display as text "MHAS_data_cleaning completed"
display as text "Log closed on `c(current_date)' at `c(current_time)'"
display as text "------------------------------------------------------------"
log close
