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
### H1N1
### Gini_Index
### Flyway
### Education
### Dairy_Network
### Dairy
### COVID-19
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

 
