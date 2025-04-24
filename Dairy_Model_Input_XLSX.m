clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');

Flyway=cell(height(US_County),1);
Flyway(US_County.ATLANTIC_FLYWAY==1)={'Atlantic'};
Flyway(US_County.PACIFIC_FLYWAY==1)={'Pacific'};
Flyway(US_County.CENTRAL_FLYWAY==1)={'Central'};
Flyway(US_County.MISSISSIPPI_FLYWAY==1)={'Mississippi'};

Number_of_Farms=US_County.TOTAL_DAIRY_OPERATIONS;
H5N1_Farms = US_County.DAIRY_HPAI_OUTBREAK_KNOWN;
H5N1_Farms_State=US_County.DAIRY_HPAI_OUTBREAK_UNKNOWN+US_County.DAIRY_HPAI_OUTBREAK_KNOWN;
H5N1_Dairy_to_Human_State=US_County.SPILLOVER_DAIRY;

US=unique(US_County.STATE_NAME);

for ss=1:length(US)
    tf=strcmp(US_County.STATE_NAME,US{ss});
    H5N1_Dairy_to_Human_State(tf)=sum(H5N1_Dairy_to_Human_State(tf)).*ones(size(H5N1_Dairy_to_Human_State(tf)));
    H5N1_Farms_State(tf)=sum(H5N1_Farms_State(tf)).*ones(size(H5N1_Farms_State(tf)));
end



Head_of_Cattle=log10(US_County.CATTLE_INVENTORY);
Head_of_Cattle(isinf(Head_of_Cattle))=NaN;
Operations_1_to_9=US_County.("CATTLE INVENTORY OF MILK COWS: (1 TO 9 HEAD)");
Operations_10_to_19=US_County.("CATTLE INVENTORY OF MILK COWS: (10 TO 19 HEAD)");
Operations_20_to_49=US_County.("CATTLE INVENTORY OF MILK COWS: (20 TO 49 HEAD)");
Operations_50_to_99=US_County.("CATTLE INVENTORY OF MILK COWS: (50 TO 99 HEAD)");
Operations_100_to_199=US_County.("CATTLE INVENTORY OF MILK COWS: (100 TO 199 HEAD)");
Operations_200_to_499=US_County.("CATTLE INVENTORY OF MILK COWS: (200 TO 499 HEAD)");
Operations_500_plus=US_County.("CATTLE INVENTORY OF MILK COWS: (500 OR MORE HEAD)");
Dairy_Network=US_County.CONNECT_DAIRY;
       Dairy_Network=Dairy_Network-diag(diag(Dairy_Network));
Connectivity=Dairy_Network;
Stopover_intensity=US_County.LIGHT_INT;
Mallard_population=log10(US_County.MALLARD+1);
Canada_Goose_population=log10(US_County.CANADA_GOOSE+1);
AGW_Teal_population=log10(US_County.AGW_TEAL+1);
N_Pintail_population=log10(US_County.NORTH_PINTAIL+1);
Temperature=US_County.TEMP-mean(US_County.TEMP);
GEOID=US_County.GEOID;
County=US_County.NAME;
State=US_County.STATE_NAME;

T=table(State,County,GEOID,Flyway,H5N1_Farms,H5N1_Farms_State,H5N1_Dairy_to_Human_State,Stopover_intensity,Mallard_population,Canada_Goose_population,AGW_Teal_population,N_Pintail_population,Temperature,Head_of_Cattle,Number_of_Farms,Operations_1_to_9,Operations_10_to_19,Operations_20_to_49,Operations_50_to_99,Operations_100_to_199,Operations_200_to_499,Operations_500_plus,Connectivity);

writetable(T,'H5N1_Dairy_Risk_Model.xlsx','Sheet','Model_Inputs');