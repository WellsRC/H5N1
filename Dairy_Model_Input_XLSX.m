clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');

state_remove=strcmp(US_County.STATE_NAME,"Alaska");

US_County=US_County(~state_remove,:);

Number_of_Farms=US_County.TOTAL_DAIRY_OPERATIONS;
H5N1_Farms = US_County.DAIRY_HPAI_OUTBREAK_KNOWN;
H5N1_Farms_State=US_County.DAIRY_HPAI_OUTBREAK_UNKNOWN;
H5N1_Dairy_to_Human_State=US_County.SPILLOVER_DAIRY;
Head_of_Cattle=US_County.CATTLE_INVENTORY;
Operations_1_to_9=US_County.("CATTLE INVENTORY OF MILK COWS: (1 TO 9 HEAD)");
Operations_10_to_19=US_County.("CATTLE INVENTORY OF MILK COWS: (10 TO 19 HEAD)");
Operations_20_to_49=US_County.("CATTLE INVENTORY OF MILK COWS: (20 TO 49 HEAD)");
Operations_50_to_99=US_County.("CATTLE INVENTORY OF MILK COWS: (50 TO 99 HEAD)");
Operations_100_to_199=US_County.("CATTLE INVENTORY OF MILK COWS: (100 TO 199 HEAD)");
Operations_200_to_499=US_County.("CATTLE INVENTORY OF MILK COWS: (200 TO 499 HEAD)");
Operations_500_plus=US_County.("CATTLE INVENTORY OF MILK COWS: (500 OR MORE HEAD)");
Poultry_operations=US_County.POULTRY_OPR_w_INVENTORY;
Stopover_intensity=US_County.LIGHT_INT;
Mallard_population=US_County.MALLARD;
Canada_Goose_population=US_County.CANADA_GOOSE;
AGW_Teal_population=US_County.AGW_TEAL;
N_Pintail_population=US_County.NORTH_PINTAIL;
GEOID=US_County.GEOID;
County=US_County.NAME;
State=US_County.STATE_NAME;

T=table(State,County,GEOID,Number_of_Farms,Head_of_Cattle,Operations_1_to_9,Operations_10_to_19,Operations_20_to_49,Operations_50_to_99,Operations_100_to_199,Operations_200_to_499,Operations_500_plus,Poultry_operations,Stopover_intensity,Mallard_population,Canada_Goose_population,AGW_Teal_population,N_Pintail_population,H5N1_Farms,H5N1_Farms_State,H5N1_Dairy_to_Human_State);

writetable(T,'H5N1_Dairy_Risk_Model.xlsx','Sheet','Model_Inputs');