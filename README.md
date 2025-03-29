# The determinants of the diagnostic accuracy of chest Xray for silicosis: A systematic review. meta-analysis and modelling study

## Overview
This contains the code and data for the meta-analysis and modelling of CXR diagnostic accuracy for silicosis. It relates to an upcoming research article

## Data loading and running scripts
The CT.xlsx, HRCT.xlsx and Autopsy.xlsx files are to be used with the MA_forest_plots.R script. This allows you to reproduce the meta-analysis and forest plots for the main meta-analysis.

The CT_all.xlsx, HRCT_all.xlsx and Autopsy.xlsx files are to be used with the MA_forest_all.R script. This allows you to reproduce the sensitivity analyses for the main meta-analysis. 

The MA_forest_diff_cutoffs.R script requires no data inputs. This allows you to run repeated meta-analyses at increasing cut-offs of the reference test. 

The props.xlsx, pred_mine_tab.csv, pred_mine.csv, pred_non_mine_tab.csv and pred_nonmine.csv are to be used with the mreg_missed.R script. This is the longest script, and performs the meta-regression of CXR severity and CXR sensitivity, and then model the impact of both relative and fixed sensitivity scenarios on the relationship between cumulative silica exposure and silica cases. 

The new_missed_calcs.R requires not inputs. It allows you to calculate the number of missed cases and number needed to treat. 

The extraction_data.csv file includes all of the original data extracted from the studies, including quality assessment. 




