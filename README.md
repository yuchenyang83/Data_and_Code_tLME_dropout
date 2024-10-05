# tLME with with non-ignorable dropout

## Supplementary Material for "Estimation and prediction in $t$ linear mixed models with non-ignorable dropout" by Yu-Chen Yang and Tsung-I Lin

### Author responsible for the code
For questions, comments or remarks about the code please contact responsible author, 

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
    

R codes for the implementation of our methodology are provided.

#### Subfolder: ./Function ####
'./Function/'
       contains 
       
       (1) the subfolder 'simulation', in this folder contains '.R' file code for simulation studies that can fit the 2-, 3- and 4-component FM-NLME and EFM-NLME models with three different mixtures consisting of poorly separated (PS) and well separated (WS) components and components that are not affected by explanatory covariates (NIC).
       (2) 'Efmnlmm.fn.R' main script for fitting the EFM-NLME model; and
       (3) 'fmnlmm.fn.R' main script for fitting the FM-NLME model

'./Function/simulation'
	subfolder collects functions for maximum likelihood (ML) estimation for FM-NLME and EFM-NLME models, including
	
    	(1) 'Efmnlmm.fn.R' main script for fitting the EFM-tLME model
        (2) 'fmnlmm.fn.R' main script for fitting the FM-tLME model
	(3) 'multiplot.R'  main script for combining multiple plots by ggplot2 package.
	

#### Subfolder: ./Code ####
'./Code/'
       contains 
       
       (1) 'Fig1.R' main script for reproducing Figure 1;
       (2) 'Fig2.R' main script for reproducing Figure 2; 
       (3) 'Fig3.R' main script for reproducing Figure 3; 
       (4) 'Fig4.R' main script for reproducing Figure 4; 
       (5) 'Fig5.R' main script for reproducing Figure 5; 
       (6) 'Fig6.R' main script for reproducing Figure 6; 
       (7) 'Fig7.R' main script for reproducing Figure 7; 
       (8) 'Fig8.R' main script for reproducing Figure 8; 
       (9) 'Tab1.R' main script for Table 1;
       (10) 'Tab2.R' main script for Table 2;
       (11) 'Tab3.R' main script for Table 3;
       (12) 'Tab4.R' main script for Table 4;
       (13) 'Tab5.R' main script for Table 5;
       (14) 'Tab6.R' main script for Table 6;

###### Note for Section 4 - Illustrative examples - simulated data:
Because the 'simulation.case1.R', 'simulation.case2.R' and 'simulation.case3.R' codes take a huge amount of time to run the ECM procedure for fitting the EFM-NLME and FM-NLME models, we record every 100 replications result in the subfolder ''SS-simulation1'', ''SS-simulation2'', and ''SS-simulation3'' from ./Data/Simulation/' respectively so that one can use the R codes 'Fig1.R', 'Fig2.R', 'Fig3.R', 'Fig4.R', 'Tab1', 'Tab2', and 'Tab3' to obtain the final results immediately.
To reproduce the results presented in Figures 1 to 4 and Tables 1 to 3, just  load the subfolder 'SS-simulation1', 'SS-simulation2', and 'SS-simulation3' files in the './Data/Simulation/' and then run the scripts 'Fig1.R', 'Fig2.R', 'Fig3.R', 'Fig4.R', 'Tab1', 'Tab2', and 'Tab3' in the subfolder './Code/';

       (15) ''simulation.case1.R' main script for fitting the FM-NLME and EFM-NLME models with g = 2 to 4 to simulated datasets with PS component; and
       (16) ''simulation.case2.R' main script for fitting the FM-NLME and EFM-NLME models with g = 2 to 4 to simulated datasets with WS component; and
       (17) ''simulation.case3.R' main script for fitting the FM-NLME and EFM-NLME models with g = 2 to 4 to simulated datasets with NIC component.
       

###### Note for Section 5 - Illustrative examples - ACTG 315 data:
Because the 'actg315_code.R' code takes a huge amount of time to run the ECM procedure for fitting the FM-NLME and EFM-NLME models, we record these intermediate results in 'ACTG315result.RData' so that one can use the R codes 'Fig5.R', 'Fig6.R', 'Fig7.R', 'Fig8.R', 'Tab4.R', 'Tab5.R', and 'Tab6.R' to obtain the final results immediately.
To reproduce the results presented in Figures 5, just load 'actg315.txt' file in the './Data/source' and then run the script 'Fig5.R' in the subfolder './Code/';
To reproduce the results presented in Figures 6, 7, and 8, just load 'actg315_code.RData' file in the './Data/' and then run the scripts 'Fig6.R', 'Fig7.R', and 'Fig8.R' in the subfolder './Code/';  
To reproduce the results presented in Tables 4, 5, and 6, just load 'actg315_code.RData' file in the './Data/' and then run the scripts 'Tab4.R', 'Tab5.R', and 'Tab6.R' in the subfolder './Code/'

       (18) 'actg315_code.R' main script for fitting the 1-, 2- and 3-component FM-NLME and the 2- and 3-component EFM-NLME models. It includes four scenarios to incorporate fixed and random effects, along with three structures for within-patient autocorrelation.


#### Subfolder: ./Data ####
'./Data/'
       contains 
       
       (1) 'ACTG315result.RData' collecting the fitting results of the the 1-, 2- and 3-component FM-NLME models and EFM-NLME models with the four scenarios and three within-patient autocorrelation structures to the ACTG 315 data.
       (2) The subfolder 'Simulation', which contains the three simulation studies used in Section 4, collecting the 100 replications result for each components;
       (3) The subfolder 'source', which contains the ACTG 315 dataset used in Section 5.

'./Data/source'
	subfolder contains
	
    	(1) 'actg315.txt' is the ACTG 315 dataset used in Section 5.

'./Data/Simulation'
	subfolder contains
	
	(1) The subfolder 'SS-simulation1' collects the fitting results for the FM-NLME and EFM-NLME models with g = 2 to 4 based on 100 replications under various sample sizes and PS component.
	(2)  The subfolder 'SS-simulation1' collects the fitting results for the FM-NLME and EFM-NLME models with g = 2 to 4 based on 100 replications under various sample sizes and WS component.
	(3)  The subfolder 'SS-simulation1' collects the fitting results for the FM-NLME and EFM-NLME models with g = 2 to 4 based on 100 replications under various sample sizes and NIC component.

'./Data/Simulation/SS-simulation1'
	subfolder contains
	
	(1) SIM2', 'SIM3', 'SIM4', 'SIM5', and 'SIM6' collect the fitting results from 100 replications obtained by fitting the FM-NLME and EFM-NLME models with g = 2 to 4 under various sample sizes and PS component.

'./Data/Simulation/SS-simulation2'
	subfolder contains
	
	(1) SIM2', 'SIM3', 'SIM4', 'SIM5', and 'SIM6' collect the fitting results from 100 replications obtained by fitting the FM-NLME and EFM-NLME models with g = 2 to 4 under various sample sizes and PS component.

'./Data/Simulation/SS-simulation3'
	subfolder contains
 
	(1) SIM2', 'SIM3', 'SIM4', 'SIM5', and 'SIM6' collect the fitting results from 100 replications obtained by fitting the FM-NLME and EFM-NLME models with g = 2 to 4 under various sample sizes and NIC component.

#### Subfolder: ./Results ####
'./Results/'
       contains 
       
       (1) 'Figure1.eps' shows the trajectory plots for one simulated dataset of size n = 90 generated from EFM-NLME model with true parameters specified in cases (i)–(iii) corresponding to the PS, WS and NIC components.
       (2) 'Figure2.eps' shows the box plots of CCR and ARI values obtained from fitting the FM-NLME and EFM-NLME models to 100 simulated datasets generated under three considered cases across various sample sizes.
       (3) 'Figure3.eps' displays the empirical mean absolute relative bias (MARB) for the estimated parameters obtained from fitting the 3-component FM-NLME and EFM-NLME models to 100 datasets simulated under three considered cases with various sample sizes n = 90, 150, 300, 450 and 600.
       (4) 'Figure4.eps' displays the mean squared errors (MSE) for the estimated parameters obtained from fitting the 3-component FM-NLME and EFM-NLME models to 100 datasets simulated under three considered cases with various sample sizes n = 90, 150, 300, 450 and 600.
       (5) 'Figure5.eps' displays the trajectory plot of log10(RNA) observations for 48 patients (left panel), bar chart of frequency distribution of numbers of observations (right top panel), and violin plots for log10(RNA) levels for patients in the two pre-specified groups (right bottom panel).
       (6) 'Figure6.eps' displays estimated mean curves of log10(RNA) across days (left panel) and violin plots of fitted responses (right panel) obtained form the best fittedFM-NLMEandEFM-NLMEmodels and their predicted clustering labels.
       (7) 'Figure7.eps' displays plots of CCR scores for the 24 candidate models with g = 2 (left panel) and ROC curves obtained from the best fitted FM-NLME and EFM-NLME models (right panel).
       (8) 'Figure8.eps' presents the temporal trajectory of log10(RNA) for six randomly selected patients along with their fitted values under the best fitted EFM-LME, EFM-tLME, and EFM-NLME models.
       (9) 'Table1.csv' lists the average BIC scores for the two fittedmodels together with the frequencies preferred by the criterion under every scenario considered.
       (10) 'Table2.csv' reports the average BIC values of both models with g = 2–4 components fitted to the data with PS, WS and NIC components and different sample sizes, and the frequencies of each model being selected.
       (11) 'Table3.csv' reports the MC Sd of each parameter along with the average values of IM SE obtained from the fitted 3-component EFM-NLME model over 100 trials with PS samples.
       (12) 'Table4.csv' lists the fitting results for the 60 candidate models, including the numbers of unknown parameters, maximized log-likelihood values together with BIC scores for determining the preferred model.
       (13) 'Table5.csv' presents a comparison of the ML estimates of parameters together with their standard errors (SE) obtained from the best EFM-NLME model and FM-NLME analogue.
       (14) 'Table6.csv' lists the mean values of MSD, MAD and MARD over 48 patients together with their standard deviations (Std) for comparing the fitting accuracies of the EFM-LME, EFM-tLME, and EFM-NLME models.

# Additional Remark 
One can directly run each "source(.)" described in 'master.r' file in the seperate R session to obtain the results.
