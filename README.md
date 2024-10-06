# Estimation and prediction in $t$ linear mixed models with non-ignorable dropout

## Supplementary Material for "Estimation and prediction in $t$ linear mixed models with non-ignorable dropout" by Yu-Chen Yang, Wan-Lun Wang, Luis M. Castro, and Tsung-I Lin

### Author responsible for the code
For questions, comments or remarks about the code please contact responsible author, Yu-Chen Yang (d107053002@smail.nchu.edu.tw)

### Configurations
The code was written/evaluated in R with the following software versions:
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.utf8  LC_CTYPE=Chinese (Traditional)_Taiwan.utf8    LC_MONETARY=Chinese (Traditional)_Taiwan.utf8
[4] LC_NUMERIC=C                                  LC_TIME=Chinese (Traditional)_Taiwan.utf8    

time zone: Asia/Taipei
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] glue_1.7.0    pROC_1.18.5   gridExtra_2.3 ggtext_0.1.2  rlang_1.1.4   cowplot_1.1.3 mvtnorm_1.3-1 ggplot2_3.5.1

loaded via a namespace (and not attached):
 [1] crayon_1.5.3      vctrs_0.6.5       nlme_3.1-164      cli_3.6.3         generics_0.1.3    labeling_0.4.3    colorspace_2.1-1 
 [8] plyr_1.8.9        gridtext_0.1.5    scales_1.3.0      fansi_1.0.6       munsell_0.5.1     tibble_3.2.1      lifecycle_1.0.4  
[15] compiler_4.4.1    dplyr_1.1.4       Rcpp_1.0.13       pkgconfig_2.0.3   rstudioapi_0.16.0 lattice_0.22-6    farver_2.1.2     
[22] R6_2.5.1          tidyselect_1.2.1  utf8_1.2.4        pillar_1.9.0      magrittr_2.0.3    tools_4.4.1       withr_3.0.1      
[29] gtable_0.3.5      xml2_1.3.6 

### Descriptions of the codes 
Please copy the files to the "current working directory" of the R package.
The 'getwd()' function shall determine an absolute pathname of the "current working directory".

Before running all of the codes, one needs to install the following R packages:

    install.packages("mvtnorm")  Version: 1.3-1
    install.packages("ggplot2") Version: 3.5.1
    install.packages("pROC") Version: 1.18.5
    install.packages("gridExtra") Version: 2.3
    install.packages("ggtext") Version: 0.1.2
    install.packages("rlang") Version: 1.1.4
    install.packages("cowplot") Version: 1.1.3
    install.packages("glue") Version: 1.7.0
    install.packages("nlme") Version: 3.1-164
    

R codes for the implementation of our methodology are provided.

#### Subfolder: ./Function ####
'./Function/'
       contains 
       
	(1) the subfolder 'simulation', in this folder contains '.R' file code for simulation studies that can fit the LME and tLME models with non-ignorable dropout consisting of different missing mechanisms;
	(2) the subfolder 'fixed_alpha', in this folder contains '.R' file code for sensitivity analysis;
	(3) 'tLMMmissingSEM.R' main script for fitting the tLME model with non-ignorable dropout;
	(4) 'LMMmissingSEM.R' main script for fitting the LME model with with non-ignorable dropout;
	(5) 'computer.pvi.R' calculate the probability of dropout for each subject at each time point;
	(6) 'multiplot.R'  main script for combining multiple plots by ggplot2 package;
	(7) 'functionR' that can assist for fitting the LME and tLME models with non-ignorable dropout.

'./Function/simulation'
	subfolder collects functions for maximum likelihood (ML) estimation for LME and tLME models with non-ignorable dropout in the simulation study, including
	
	(1) 'tLMMmissingSEM.R' main script for fitting the tLME model with non-ignorable dropout;
	(2) 'LMMmissingSEM.R' main script for fitting the LME model with with non-ignorable dropout;
	(3) 'functionR' that can assist for fitting the LME and tLME models with non-ignorable dropout

'./Function/fixed_alpha'
	subfolder collects functions for maximum likelihood (ML) estimation for LME and tLME models with non-ignorable dropout in the sensitivity analysis, including
	
	(1) 'tLMMmissingSEM.R' main script for fitting the tLME model with non-ignorable dropout;
	(2) 'LMMmissingSEM.R' main script for fitting the LME model with with non-ignorable dropout;
	(3) 'functionR' that can assist for fitting the LME and tLME models with non-ignorable dropout


#### Subfolder: ./Code ####
'./Code/'
       contains 
       
	(1) 'Fig1.R' main script for reproducing Figure 1;
	(2) 'Fig2.R' main script for reproducing Figure 2; 
	(3) 'Fig3.R' main script for reproducing Figure 3; 
	(4) 'Fig4.R' main script for reproducing Figure 4; 
	(5) 'Fig5.R' main script for reproducing Figure 5; 
	(6) 'Tab1.R' main script for Table 1;
	(7) 'Tab2.R' main script for Table 2;
	(8) 'Tab3.R' main script for Table 3;
	(9) 'Tab4.R' main script for Table 4;
	(10) 'Tab5.R' main script for Table 5;

###### Note for Section 5 - Illustrative examples - simulated data:
1. The 'simSEM25.R', 'simSEM50.R', and 'simSEM75.R' codes require a significant amount of time to run the MCECM procedure for fitting the LME and tLME models with non-ignorable dropout.
2. Every 100 replications result is recorded in the subfolders 'SS-simulationSEM-t25', 'SS-simulationSEM-t50', and 'SS-simulationSEM-t75' from './Data/Simulation/' respectively.
3. The R codes 'Fig1.R', 'Fig2.R', 'Tab1', 'Tab2', and 'Tab3' can be used to obtain the final results immediately.
4. To reproduce the results presented in Figures 1 and 2 and Tables 1 to 3: (i) Load the files from the subfolders 'SS-simulationSEM-t25', 'SS-simulationSEM-t50', and 'SS-simulationSEM-t75' located in './Data/Simulation/'. (ii) Then, run the scripts 'Fig1.R', 'Fig2.R', 'Tab1', 'Tab2', and 'Tab3' from the subfolder './Code/'.



       (11) The main script 'simSEM25.R' is used to fit LME and tLME models under three missing data mechanisms to simulated datasets with 25% dropout rate; and
       (12) The main script 'simSEM50.R' is used to fit LME and tLME models under three missing data mechanisms to simulated datasets with 50% dropout rate; and
       (13) The main script 'simSEM75.R' is used to fit LME and tLME models under three missing data mechanisms to simulated datasets with 75% dropout rate.



###### Note for Section 6 - Illustrative examples - AIDS clinical trial:
1. The 'fit_aids.R' code takes a significant amount of time to run the MCECM procedure for fitting the LME and tLME models with non-ignorable dropout.
2. Intermediate results are recorded in 'fit_result.RData' located in './Data' so that the R codes 'Fig3.R', 'Fig4.R', 'Fig5.R', 'Tab4.R', and 'Tab5.R' can be used to obtain the final results immediately.
3. Steps to reproduce the results:
   - **For Figure 3:** (i) Load the 'aids.RData' file from './Data/source'. (ii) Run the 'Fig3.R' script from the './Code/' subfolder.
   - **For Figure 4:** (i) Load the 'fit_result.RData' file from './Data/'. (ii) Run the 'Fig4.R' script from the './Code/' subfolder.
   - **For Figure 5:** (i) Load the 'fit_result.RData' file from './Data/'. (ii) Load the 'fix_alpha' files: 'fix_alpha1.RData' through 'fix_alpha22.RData' from the './Data/fix_alpha/' directory. (iii) Run the 'Fig4.R' script from the './Code/' subfolder.
   - **For Table 4 & 5:** (i) Load the 'fit_result.RData' file from './Data/'. (ii) Run the 'Tab4.R' and 'Tab5.R' scripts from the './Code/' subfolder.



       (14) 'fit_aids.R' main script for fitting the LME and tLME models under the three mechanisms. It includes four structures for within-patient autocorrelation.


#### Subfolder: ./Data ####
'./Data/'
       contains 
       
       (1) 'fit_result.RData' collecting the fitting results of the LME and tLME models with the three missing mechanisms and four within-patient autocorrelation structures to the AIDS data;
       (2) The subfolder 'Simulation', which contains the three simulation studies used in Section 5, collecting the 100 replications result for each components;
       (3) The subfolder “fixed_alpha”, which contains the results of the sensitivity analysis, collecting fixed alpha 2 value within the range of -16, -8, -2, -1, -0.5, -0.1, -0.05, -0.01, -0.001, 0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 8, 16.
       (4) The subfolder 'source', which contains the AIDS dataset used in Section 6.

'./Data/source'
	subfolder contains
	
    	(1) 'aids.RData' is the AIDS dataset used in Section 6.

'./Data/Simulation'
	subfolder contains
	
	(1) The subfolder “SS-simulationSEM-t25” contains the fitting results of the LME and tLME models for data with a 25% dropout rate, where three missing data mechanisms are based on 100 repetitions under various sample sizes.
	(2)  The subfolder “SS-simulationSEM-t50” contains the fitting results of the LME and tLME models for data with a 50% dropout rate, where three missing data mechanisms are based on 100 repetitions under various sample sizes.
	(3)  The subfolder “SS-simulationSEM-t75” contains the fitting results of the LME and tLME models for data with a 75% dropout rate, where three missing data mechanisms are based on 100 repetitions under various sample sizes.

'./Data/Simulation/SS-simulationSEM-t25'
	subfolder contains
	
	(1) 'SIM1', 'SIM2', 'SIM3', 'SIM4', and 'SIM5' the fitting results of the LME and tLME models for data with a 25% dropout rate, where three missing data mechanisms are based on 100 repetitions under various sample sizes.

'./Data/Simulation/SS-simulationSEM-t50'
	subfolder contains
	
	(1) 'SIM1', 'SIM2', 'SIM3', 'SIM4', and 'SIM5' the fitting results of the LME and tLME models for data with a 50% dropout rate, where three missing data mechanisms are based on 100 repetitions under various sample sizes.

 './Data/Simulation/SS-simulationSEM-t75'
	subfolder contains
	
	(1) 'SIM1', 'SIM2', 'SIM3', 'SIM4', and 'SIM5' the fitting results of the LME and tLME models for data with a 75% dropout rate, where three missing data mechanisms are based on 100 repetitions under various sample sizes.

 './Data/fixed_alpha'
 	subfolder contains
	
	(1) 'fix_alpha1.RData', 'fix_alpha2.RData', 'fix_alpha3.RData', 'fix_alpha4.RData', 'fix_alpha5.RData', 'fix_alpha6.RData', 'fix_alpha7.RData', 'fix_alpha8.RData', 'fix_alpha9.RData', 'fix_alpha10.RData', 'fix_alpha13.RData', 'fix_alpha14.RData', 'fix_alpha15.RData', 'fix_alpha16.RData', 'fix_alpha17.RData', 'fix_alpha18.RData', 'fix_alpha19.RData', 'fix_alpha20.RData', 'fix_alpha21.RData', 'fix_alpha22.RData' the fitting results of the tLME model under the MNAR mechanism, with the fixed alpha2 parameter values ranging from -16 to 16.
#### Subfolder: ./Results ####
'./Results/'
       contains 
       
       (1) 'sim-trajectory.eps' shows the trajectory plot of responses for one simulated case of size N = 100 (left panel), box plots for responses for subjects in the two prespecified groups (right top panel), and bar chart of frequency distribution of numbers of observed responses (right bottom panel).
       (2) 'sim-fig2.eps' shows the MSE scores for the estimated parameters under fitted tLME model with the MNAR mechanism for three dropout rates across various sample sizes.
       (3) 'AIDS_Figure1.eps' displays the trajectory plots of the square roots of CD4 observations for 467 patients (top panel) and bar chart of frequency distribution of numbers of observed the square roots of CD4 (bottom panel) for the two drug treatments.
       (4) 'AIDS_Figure2.eps' displays ROC curves obtained from the best fitted LME (dashed lines) and tLME (sold lines) models under the three mechanisms.
       (5) 'AIDS_Figure3.eps' displays estimated fixed effects parameters (solid line) and 95% confidence intervals (dashed lines) from the fitted tLME model with MNAR.
       (6) 'Table1.csv' lists the simulation results for assessing the asymptotic standard errors (IM SE) and empirical standard deviations (MC Sd) of parameters estimates under fitted tLME model with the MNAR mechanism across various sample sizes.
       (7) 'Table2.csv' reports the average AIC and BIC scores together with frequencies (in parentheses) supported by the two criteria for the LME and tLME models with various sample sizes and dropout rates. The result for the best performance per row is highlighted in bold.
       (8) 'Table3.csv' reports the comparison of predictive accuracies of missing responses in terms of MSPE between the LME and tLME models. The frequencies (Freq) favored by the model with a lower MSPE value are also recorded. The result for the best performance per row is highlighted in bold.
       (9) 'Table4.csv' lists the fitting results for the 24 candidate models, including the numbers of unknown parameters, maximized log-likelihood values together with AIC and BIC scores for determining the preferred model.
       (10) 'Table5.csv' presents a comparison of the ML estimates of parameters together with their standard errors (SE) obtained from the tLME model with AR(1) errors under the three mechanisms.

# Additional Remark 
One can directly run each "source(.)" described in 'master.r' file in the seperate R session to obtain the results.
