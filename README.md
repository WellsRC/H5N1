# Spatial risk assessment of H5N1 outbreaks among poultry farms, dairy farms, and humans in the United States
Chad R. Wells <sup>1</sup>,✝ , Michael Cairo <sup>1,2</sup>, Sarah Galvani-Townsend <sup>1</sup>, Seyed M. Moghadas <sup>3</sup>, Meagan C. Fitzpatrick <sup>4</sup>, Abhishek Pandey <sup>1</sup>,✝

<sup>1</sup> Center for Infectious Disease Modeling and Analysis, Yale School of Public Health, New Haven, CT, USA <br />
<sup>2</sup> Northeastern University, Boston, MA, USA <br />
<sup>3</sup> Agent-Based Modelling Laboratory, York University, Toronto, Ontario, Canada <br />
<sup>4</sup> Center for Vaccine Development and Global Health, University of Maryland School of Medicine, Baltimore, MD, USA <br />

✝Corresponding authors: abhishek.pandey@yale.edu and chad.richard.wells@gmail.com <br />
## OS System requirements
The codes developed here are tested on Windows operating system (Windows 10 Home: 64-bit). Optimization was also conducted on a 48 core nodes with 384GB RAM, 2TB SSD, dual Xeon Platinum 8268 2.9GHz CPUs.

## Installation guide
### MATLAB
The code was run with Both MATLAB R2022a and R2024a. Installation instruction for MATLAB can be found at https://www.mathworks.com/help/install/install-products.html. Typical install time for MATLAB on a "normal" desktop is around 30-40 minutes. 

## Demo
Fit_Dairy_Measure runs the optimization for all the models computing the risk for spillover events associated with dairy farms. Aftr the fitting is complete, Compute_Risk_Dairy.m computes the risk and averag risk across the different models to be plotted by Figure_Risk_Map.m (i.e., Figure_Risk_Map('Dairy'));

## Instructions for use
To generate the Figures and output of the calculations, select a script from Figures and Tables section to run in MATLAB and enter the name in the command line. All mat file are availble to generate figures and conduct the calculations.

## Data
Contains the raw data files that were downlaoded and the scripts used to scrape and aggregate it into a single structure at the county-level for use in the analysis
### Age_Sex
ACSST5Y2010.S0101-Data.csv - County-level age and sex demographics obtained from the US Cenesus Bureau for 2010 <br />
ACSST5Y2020.S0101-Data.csv - County-level age and sex demographics obtained from the US Cenesus Bureau for 2020 <br />
ACSST5Y2022.S0101-Data.csv - County-level age and sex demographics obtained from the US Cenesus Bureau for 2022 <br />
### Bird_Stopover
Stopover_Data.m - Scrapes the stopover data from the tif files <br /> 
County_Level_Stopover.mat - The county level data obtained from Stopover_Data.m <br />
spring_stopover_2500_v9_265_class.tif - Tif file for spring (DOI:10.6084/m9.figshare.24438280.v2)
fall_stopover_2500_v9_265_class.tif - Tif file for fall (DOI:10.6084/m9.figshare.24438280.v2) <br />
### COVID-19
OxCGRTUS_timeseries_all.xlsx - Stringency index for each state <br />
AH_County_of_Residence_COVID-19_Deaths_Counts__2020_Provisional.csv - COVID-19 provisional mortality data at the county-level <br />
### Dairy
Infer_County_Milk_Cow_Inventory.m - Script used to infer the dairy cattle inventory for counties in which values were supressed to maintain farm anonymity <br />
County_Milk_Cow_Inventory_2022.csv - The raw dairy cattle inventory data for counties <br />
County_Operations_with_Inventory_Cattle_Milk_2022.csv - The raw data for dairy cattle operations with inventory stratefied by numebr of head for counties <br />
Inferred_County_Milk_Cow_Inventory_2022.csv - The inferred and raw data for dairy cattle inventory <br />
State_Milk_Cow_Inventory_2022.csv - State-level dairy cattle inventory
### Dairy_Network
Average_In_Out_County.m - Computes the average in and out flow of each county
Dairy_County_Network.mat - The average daiy network
FromToKernelGenhwgaallX.txt- The output of the network from the Xth ensemble of network construction (https://mountainscholar.org/items/75fd3841-11c3-4cec-91b5-45d36b8c6b5a)
### Education
ACSST5Y2010.S1501-Data.csv - County-level education demographics obtained from the US Cenesus Bureau for 2010 <br />
ACSST5Y2020.S1501-Data.csv - County-level education demographics obtained from the US Cenesus Bureau for 2020 <br />
ACSST5Y2022.S1501-Data.csv - County-level education demographics obtained from the US Cenesus Bureau for 2022 <br />
### Gini_Index
ACSDT5Y2010.B19083-Data.csv - County-level GINI index obtained from the US Cenesus Bureau for 2010 <br />
ACSDT5Y2020.B19083-Data.csv - County-level GINI index obtained from the US Cenesus Bureau for 2020 <br />
ACSDT5Y2022.B19083-Data.csv - County-level GINI index obtained from the US Cenesus Bureau for 2022 <br />
### H1N1
h1n1-070609.xlsx - H1N1 county-level incidence data
### H5N1_Outbreaks
Event_4451.pdf - Raw data specifying county-level information for H5N1 outbreaks among  dairy cattle  <br />
HPAI Detections in Wild Birds.csv - H5N1 surviellane among wild birds  <br />
H5N1_State_Human.csv - H5N1 spillover events linked to dairy and poultry  <br />
HPAI_Poulty_Farms.xlsx - H5N1 outbreaks among poultry farms  <br />
HPAI_Dairy_Farms.xlsx - H5N1 outbreaks among dairy farms  <br />
### Household_Income
ACSST5Y2010.S1901-Data.csv - County-level household income demographics obtained from the US Cenesus Bureau for 2010 <br />
ACSST5Y2020.S1901-Data.csv - County-level household income demographics obtained from the US Cenesus Bureau for 2020 <br />
ACSST5Y2022.S1901-Data.csv - County-level household income demographics obtained from the US Cenesus Bureau for 2022 <br />
### Influenza_Testing
Scrape_State_Testing.m - Script used to scrape state-level testing data for influenza <br />
State_Level_Influenza_Testing.mat - The state-level testing data <br />
WHO_NREVSS_Clinical_Labs.csv - State-level testing among clinical labs <br />
WHO_NREVSS_Public_Health_Labs.csv - State-level testing among public health labs <br />
### Poultry
Infer_County_Chicken_Inventory.m -  Script used to infer the poultry inventory for counties in which values were supressed to maintain farm anonymity <br />
Inferred_Rooster_Inventory_County_2022.csv -  The inferred and raw data for rooster inventory <br />
Inferred_Pullet_Inventory_County_2022.csv -  The inferred and raw data for pullet inventory <br />
Inferred_Layers_Inventory_County_2022.csv -  The inferred and raw data for layer inventory <br />
Inferred_Broiler_Inventory_County_2022.csv -  The inferred and raw data for broiler inventory <br />
County_Turkey_Farms_with_Inventory_2022.sv - County-level raw data for the numebr of turkey operations with invenotry <br />
County_Operations_with_Inventory_Poultry_2022.sv - County-level raw data for the numebr of poultry operations with invenotry <br />
County_Operations_with_Inventory_Layer_Chickens_2022.sv - County-level raw data for the numebr of layer operations with invenotry stratiied by the number of head<br />
County_Level_Pullet_Operations_with_Inventory_2022.sv - County-level raw data for the numebr of pullet operations with invenotry <br />
County_Level_Layer_Operations_with_Inventory_2022.sv - County-level raw data for the numebr of layer operations with invenotry <br />
County_Level_Broiler_Operations_with_Inventory_2022.sv - County-level raw data for the numebr of broiler operations with invenotry <br />
County_Chicken_Roosters_Inventory_2022.csv - The raw data for rooster inventory <br />
County_Chicken_Pullets_Inventory_2022.csv - The raw data for pullet inventory <br />
County_Chicken_Layers_Inventory_2022.csv - The raw data for layer inventory <br />
County_Chicken_Broiler_Inventory_2022.csv - The raw data for broiler inventory <br />
### Rural_Population
Rural_Population_County_2010.xls - Rural and urban demographics for 2010 from the US Census Bureau <br /> 
Rural_Population_County_2020.xlsx - Rural and urban demographics for 2020 from the US Census Bureau <br />
### Slaughterhouse
County_FIPS_Slaughterhouse.m - Script used to pull the location of the slaughterhouses <br />
geoCode.m - Looks up the geographical coordinates of an address (https://www.mathworks.com/matlabcentral/fileexchange/37860-geocode) <br />
ZIP-COUNTY-FIPS_2017-06.csv - Used in the mapping of zip code to county fips <br />
slaughter_productionData.xlsx - Location/address of slaughterhouses <br />
Location_Slaughterhouse.mat - The county fips of the slaughter houses obtained from County_FIPS_Slaughterhouse.m <br />
### Swine
Infer_Hog_County.m - Script used to infer the hog inventory for counties in which values were supressed to maintain farm anonymity <br />
Inferred_Hogs_Inventory_County_2022.csv - The inferred and raw data for hog inventory <br />
County_Level_Hog_Inventory_2022.csv - The raw hog inventory data for counties <br />
County_Level_Operations_Inventory_2022.csv - The raw data for hog operations with inventory stratefied by numebr of head for counties <br />
## Shapefile
Contains the files for the shapefile used in the analysis and the generation of the figures
## Fitting
Risk_Assesment_Farms.m - Logistic risk function used in the models
### Poultry
Fit_Poultry_Measure.m - Fits all of the poultry models <br />
Optimize_Poultry_Farm_Risk.m - Runs the optimization for a specified poultry model <br />
Objective_Function_Poultry_Farm.m - Objective function to be optimzied for the poultry model <br />
Poultry_Covariates.m - Covariates for the poultry models <br />
### Dairy
Fit_Dairy_Measure.m - Fits all of the dairy models <br />
Optimize_Dairy_Farm_Risk.m - Runs the optimization for a specified dairy model <br />
Objective_Function_Dairy_Farm.m - Objective function to be optimzied for the dairy model <br />
Dairy_Covariates.m - Covariates for the dairy models <br />
### Humans
Fit_Population_Measure.m - Fits all of the human population models <br />
Population_Covariates.m - Populatino covariates for 2010 and 2020 <br />
Population_Covariates_H5N1.m - Populatino covariates for 2022 <br />
Optimize_Population_Farm_Risk.m - Runs the optimization for a specified human population model <br /> 
Objective_Function_Population.m - Objective function to be optimzied for the human population model <br />
## Risk computation
Compute_Risk_Swine.m - Computes the risk for swine based on the inference from poultry and dairy cattle <br />
Swine_Covariates.m - Covariates for the swine model <br />
Compute_Underreporting.m - Computes the probability of underreporting at the county-level <br />
Compute_Surveillance.m - Computes the likelihood of a H5N1 case going undetected within a county  <br /> 
Compute_Risk_Poultry.m - Computes the risks for counties with poultry operations <br />
Compute_Risk_Dairy.m - Computes the risks for counties with dairy operations <br />
Compute_Risk_Population.m - Computes the risks for counties for H5N1 outreaks among humans based on inference using H1N1 and COVID-19 data <br />
Compute_Overall_Spillover_Risk.m - Computes the risk of spillover from either poultry or dairy cattle <br />
Compute_Risk_Dairy_Reduced_Connectivity.m - Computes the risks for counties with dairy operations for a specified level of reduction in dairy cattle transportation <br />
## Figures and Tables
Table_Risk_Model.m - Outputs the AIC model weighted output for the specified model <br /> 
Figure_Change_Susceptibility_Risk_Dairy_Farm.m - Plots the the relative reduction dairy cattle transportation <br /> 
Figure_Exacerbate_Spillover.m - Plots the risk assocaited with swine and the number of slaughterhouses in a county <br /> 
Figure_Risk_Map.m - Plots the risks for poultry, dairy, and the scenario for reduction in dairy cattle transportation <br /> 
Figure_Spillover_Risk_Map.m - Plots the results for spillover risk and subsequent localcised transmission <br /> 
Figure_Surveillance_Map.m - Plots the likelihood of a H5N1 case going undetected  <br />
Figure_Susceptible_Risk_Map.m - Plots the susceptiblity of countrys for H5N1 outbreaks among humans  <br />
Figure_Underreporting_Map.m - Plots the results for underreporting among dairy cattle operations   <br />
Supplementary_Figure_Spillover_Risk_Map.m - Plots the results for spillvoer from a specified farm type <br />

 
