clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');

Flyway=cell(height(US_County),1);
Flyway(US_County.ATLANTIC_FLYWAY==1)={'Atlantic'};
Flyway(US_County.PACIFIC_FLYWAY==1)={'Pacific'};
Flyway(US_County.CENTRAL_FLYWAY==1)={'Central'};
Flyway(US_County.MISSISSIPPI_FLYWAY==1)={'Mississippi'};

Number_of_Farms=US_County.POULTRY_OPR_w_INVENTORY;
Number_of_Pullet_Farms=US_County.PULLET_OPR_w_INVENTORY;
Number_of_Layer_Farms=US_County.LAYER_OPR_w_INVENTORY;
Number_of_Turkey_Farms=US_County.TURKEY_OPR_w_INVENTORY;
Number_of_Broiler_Farms=US_County.BROILER_OPR_w_INVENTORY;
H5N1_Farms = US_County.POULTRY_HPAI_OUTBREAK;
H5N1_Poultry_to_Human_State=US_County.SPILLOVER_POULTRY;

Turkey_Inventory=log10(1+US_County.TURKEY_INVENTORY);
Broiler_Inventory=log10(1+US_County.BROILER_INVENTORY);
Pullet_Inventory=log10(1+US_County.PULLET_INVENTORY);
Layer_Inventory=log10(1+US_County.LAYER_INVENTORY);

Temperature=US_County.TEMP-mean(US_County.TEMP);
Stopover_intensity=US_County.LIGHT_INT;
Mallard_population=log10(US_County.MALLARD+1);
Canada_Goose_population=log10(US_County.CANADA_GOOSE+1);
AGW_Teal_population=log10(US_County.AGW_TEAL+1);
N_Pintail_population=log10(US_County.NORTH_PINTAIL+1);
GEOID=US_County.GEOID;
County=US_County.NAME;
State=US_County.STATE_NAME;

T=table(State,County,GEOID,Flyway,H5N1_Farms,H5N1_Poultry_to_Human_State,Number_of_Pullet_Farms,Pullet_Inventory,Number_of_Layer_Farms,Layer_Inventory,Number_of_Turkey_Farms,Turkey_Inventory,Number_of_Broiler_Farms,Broiler_Inventory,Stopover_intensity,Mallard_population,Canada_Goose_population,AGW_Teal_population,N_Pintail_population);

writetable(T,'H5N1_Poultry_Risk_Model.xlsx','Sheet','Model_Inputs');