# Sex Differences in Frailty Trajectories and Terminal Decline in Older Adults in Mexico
 
Stata code for analysing sex differences in frailty trajectories over 17 years of follow-up, and in the years before death, among older adults in the Mexican Health and Aging Study (MHAS).
 
## Background
 
Women accumulate more frailty than men yet live longer (the **health-survival paradox**), and frailty accelerates sharply in the years before death (**terminal decline**). Most evidence on these sex differences comes from high-income countries. This project examines whether they hold in Mexico, part of Latin America, one of the fastest-ageing regions in the world.
 
## Aims
 
1. Investigate sex differences in frailty trajectories over 17 years of follow-up among older adults in Mexico.
2. Investigate sex differences in frailty trajectories in the years preceding death.
## Data
 
- **Source:** [Mexican Health and Aging Study (MHAS)](https://www.mhasweb.org/)
- **Waves:** 2001, 2003, 2012, 2015, 2018
- **Frailty measure:** Frailty Index (FI) based on deficit accumulation
> **Note:** Raw data are not included in this repository. Register and download the data from the official MHAS website.
 
## Repository structure
 
All analysis scripts are Stata do-files in [`dofile/`](dofile/), numbered in the order they should be run.
 
| # | Script | Purpose |
|---|--------|---------|
| 1 | `1_MHAS_data_cleaning.do` | Data cleaning and Frailty Index construction |
| 2 | `2_MHAS_descriptive_analysis.do` | Descriptive statistics |
| 3 | `3_MHAS_weighted_analysis.do` | Survey-weighted analysis |
| 4 | `4_MHAS_LMMs_selection.do` | Linear mixed model selection |
| 5 | `5_MHAS_LMMs_Post-estimation.do` | LMM post-estimation |
| 6 | `6_MHAS_Cox_survival_analysis.do` | Cox survival analysis |
| 7 | `7_MHAS_missing_data_analysis.do` | Missing data analysis |
| 8 | `8_MHAS_General_sensitivity_analysis.do` | General sensitivity analyses |
| 9 | `9_MHAS_multiple_imputation.do` | Multiple imputation |
| 10 | `10_MHAS_LMMs_interaction.do` | LMMs with sex interaction terms |
| 11 | `11_MHAS_multiple_imputation_interaction.do` | Multiple imputation for interaction models |
| 12 | `12_MHAS_Cox_survival_analysis.do` | Cox survival analysis (updated) |
| 13 | `13_MHAS_sensitivity analysis_death.do` | Terminal decline (time-to-death) analysis |
| 14 | `14_MHAS_MI_time-to-death.do` | Multiple imputation for time-to-death models |
| 15 | `15_MHAS_multiple_imputation_sensitivity analysis.do` | Multiple imputation sensitivity analyses |
 
## How to use
 
1. Clone this repository:
```bash
   git clone https://github.com/aleixopbrunetti/mhas-frailty-sex-differences.git
```
2. Download the MHAS data.
3. Set your working directory and define the folder paths used by the do-files. Run this once at the start of your Stata session, replacing the paths with your own:
```stata
   cd "path/to/raw_mhas_data"
 
   global rawmhasdata   "path/to/raw_mhas_data"
   global cleanmhasdata "path/to/clean_data"
   global working       "path/to/working"
   global outputtables  "path/to/outputs/tables"
   global outputgraphs  "path/to/outputs/graphs"
```
4. Install the required user-written package:
```stata
   ssc install parmest   
   ssc install estout
   ssc install gammafit
   net from https://www.andrew.cmu.edu/user/bjones/traj
   net install traj
   net install st0263, from(http://www.stata-journal.com/software/sj12-3)
   net install mimrgns, from(http://fmwww.bc.edu/RePEc/bocode/m)
```
5. Run the do-files in numerical order.

## Author

Aleixo P Brunetti<br>
PhD Researcher<br>
London School of Hygiene & Tropical Medicine<br>
Faculty of Epidemiology and Population Health<br>
Department of Non-Communicable Disease Epidemiology
