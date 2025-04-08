clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');

state_remove=strcmp(US_County.STATE_NAME,"Alaska");

US_County=US_County(~state_remove,:);

Number_of_Farms=US_County.POULTRY_OPR_w_INVENTORY;
Number_of_Pullet_Farms=US_County.PULLET_OPR_w_INVENTORY;
Number_of_Layer_Farms=US_County.LAYER_OPR_w_INVENTORY;
Number_of_Turkey_Farms=US_County.TURKEY_OPR_w_INVENTORY;
Number_of_Broiler_Farms=US_County.BROILER_OPR_w_INVENTORY;
H5N1_Farms = US_County.POULTRY_HPAI_OUTBREAK;
H5N1_Poultry_to_Human_State=US_County.SPILLOVER_POULTRY;

Total_Poultry_Inventory=US_County.BROILER_INVENTORY+US_County.ROOSTER_INVENTORY+US_County.PULLET_INVENTORY+US_County.LAYER_INVENTORY;
Broiler_Inventory=US_County.BROILER_INVENTORY;
Pullet_Inventory=US_County.PULLET_INVENTORY;
Layer_Inventory=US_County.LAYER_INVENTORY;

Stopover_intensity=US_County.LIGHT_INT;
Mallard_population=US_County.MALLARD;
Canada_Goose_population=US_County.CANADA_GOOSE;
AGW_Teal_population=US_County.AGW_TEAL;
N_Pintail_population=US_County.NORTH_PINTAIL;
GEOID=US_County.GEOID;
County=US_County.NAME;
State=US_County.STATE_NAME;

T=table(State,County,GEOID,Number_of_Farms,Number_of_Pullet_Farms,Number_of_Layer_Farms,Number_of_Turkey_Farms,Number_of_Broiler_Farms,Total_Poultry_Inventory,Broiler_Inventory,Pullet_Inventory,Layer_Inventory,Stopover_intensity,Mallard_population,Canada_Goose_population,AGW_Teal_population,N_Pintail_population,H5N1_Farms,H5N1_Poultry_to_Human_State);

writetable(T,'H5N1_Poultry_Risk_Model.xlsx','Sheet','Model_Inputs');