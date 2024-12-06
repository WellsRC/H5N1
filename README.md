# Spatial risk assessment of H5N1 outbreaks among poultry farms, dairy farms, and humans in the United States: a modelling study
Chad R. Wells <sup>1</sup>,✝ , Michael Cairo <sup>1,2</sup>, Sarah Galvani-Townsend <sup>1</sup>, Zvonimir Poljak <sup>3,4</sup>, Seyed M. Moghadas <sup>5</sup>, Meagan C. Fitzpatrick <sup>6</sup>, Abhishek Pandey <sup>1</sup>,✝

<sup>1</sup> Center for Infectious Disease Modeling and Analysis, Yale School of Public Health, New Haven, CT, USA <br />
<sup>2</sup> Northeastern University, Boston, MA, USA <br />
<sup>3</sup> Department of Population Medicine, Ontario Veterinary College, University of Guelph, 50 Stone Road East, Guelph, Ontario, Canada <br />
<sup>4</sup> Centre for Public Health and Zoonoses, Ontario Veterinary College, University of Guelph, 50 Stone Road East, Guelph, Ontario, Canada <br />
<sup>5</sup> Agent-Based Modelling Laboratory, York University, Toronto, Ontario, Canada <br />
<sup>6</sup> Center for Vaccine Development and Global Health, University of Maryland School of Medicine, Baltimore, MD, USA <br />

✝Corresponding authors: abhishek.pandey@yale.edu and chad.richard.wells@gmail.com <br />

## Data
Contains the raw data files that were downlaoded and the scripts used to scrape and aggregate it into a single structure at the county-level for use in the analysis
### Swine
Infer_Hog_County.m - Script used to infer the hog inventory for counties in which values were supressed to maintain farm anonymity <br />
Inferred_Hogs_Inventory_County_2022.csv - The inferred and raw data for hog inventory <br />
County_Level_Hog_Inventory_2022.csv - The raw hog inventory data for counties <br />
County_Level_Operations_Inventory_2022.csv - The raw data for hog operations with inventory stratefied by numebr of head for counties <br />
### Slaughterhouse
County_FIPS_Slaughterhouse.m - Script used to pull the location of the slaughterhouses <br />
geoCode.m - Looks up the geographical coordinates of an address (https://www.mathworks.com/matlabcentral/fileexchange/37860-geocode) <br />
ZIP-COUNTY-FIPS_2017-06.csv - Used in the mapping of zip code to county fips <br />
slaughter_productionData.xlsx - Location/address of slaughterhouses <br />
Location_Slaughterhouse.mat - The county fips of the slaughter houses obtained from County_FIPS_Slaughterhouse.m <br />
### Rural_Population
### Influenza_Testing
### Household_Income
### H5N1_Outbreaks
Event_4451.pdf - Raw data specifying county-level information for H5N1 outbreaks among  dairy cattle
HPAI Detections in Wild Birds.csv - H5N1 surviellane among wild birds
H5N1_State_Human.csv - H5N1 spillover events linked to dairy and poultry
HPAI_Poulty_Farms.xlsx - H5N1 outbreaks among poultry farms
HPAI_Dairy_Farms.xlsx - H5N1 outbreaks among dairy farms
### H1N1
h1n1-070609.xlsx - H1N1 county-level incidence data
### Gini_Index
ACSDT5Y2010.B19083-Data.csv - County-level GINI index obtained from the US Cenesus Bureau for 2010 <br />
ACSDT5Y2020.B19083-Data.csv - County-level GINI index obtained from the US Cenesus Bureau for 2020 <br />
ACSDT5Y2022.B19083-Data.csv - County-level GINI index obtained from the US Cenesus Bureau for 2022 <br />
### Education
ACSST5Y2010.S1501-Data.csv - County-level education demographics obtained from the US Cenesus Bureau for 2010 <br />
ACSST5Y2020.S1501-Data.csv - County-level education demographics obtained from the US Cenesus Bureau for 2020 <br />
ACSST5Y2022.S1501-Data.csv - County-level education demographics obtained from the US Cenesus Bureau for 2022 <br />
### Dairy_Network
Average_In_Out_County.m - Computes the average in and out flow of each county
Dairy_County_Network.mat - The average daiy network
FromToKernelGenhwgaallX.txt- The output of the network from the Xth ensemble of network construction (https://mountainscholar.org/items/75fd3841-11c3-4cec-91b5-45d36b8c6b5a)
### Dairy
Infer_County_Milk_Cow_Inventory.m - Script used to infer the dairy cattle inventory for counties in which values were supressed to maintain farm anonymity <br />
County_Milk_Cow_Inventory_2022.csv - The raw dairy cattle inventory data for counties <br />
County_Operations_with_Inventory_Cattle_Milk_2022.csv - The raw data for dairy cattle operations with inventory stratefied by numebr of head for counties <br />
Inferred_County_Milk_Cow_Inventory_2022.csv - The inferred and raw data for dairy cattle inventory <br />
State_Milk_Cow_Inventory_2022.csv - State-level dairy cattle inventory
### COVID-19
OxCGRTUS_timeseries_all.xlsx - Stringency index for each state <br />
AH_County_of_Residence_COVID-19_Deaths_Counts__2020_Provisional.csv - COVID-19 provisional mortality data at the county-level <br />
### Bird_Stopover
Stopover_Data.m - Scrapes the stopover data from the tif files <br /> 
County_Level_Stopover.mat - The county level data obtained from Stopover_Data.m <br />
spring_stopover_2500_v9_265_class.tif - Tif file for spring (DOI:10.6084/m9.figshare.24438280.v2)
fall_stopover_2500_v9_265_class.tif - Tif file for fall (DOI:10.6084/m9.figshare.24438280.v2) <br />
### Age_Sex
ACSST5Y2010.S0101-Data.csv - County-level age and sex demographics obtained from the US Cenesus Bureau for 2010 <br />
ACSST5Y2020.S0101-Data.csv - County-level age and sex demographics obtained from the US Cenesus Bureau for 2020 <br />
ACSST5Y2022.S0101-Data.csv - County-level age and sex demographics obtained from the US Cenesus Bureau for 2022 <br />
## Shapefile
Contains the files for the shapefile used in the analysis and the generation of the figures
## Fitting

## Computation

## Figures

 
