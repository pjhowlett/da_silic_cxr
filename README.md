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

## Copyright

Copyright (c) 2022–2026, Imperial College London.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

Neither the name of Imperial College London nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.




