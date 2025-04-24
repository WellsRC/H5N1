clear;
clc;
load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);

load('Poultry_Risk_AIC.mat');

GEOID=US_County.GEOID;
County=US_County.NAME;
State=US_County.STATE_NAME;

Outbreak=mle_outbreak_poultry_farm_County;
Outbreak_95UR=outbreak_poultry_farm_County_95;

Outbreak_Risk=mle_outbreak_risk_poultry_farm_County;
Outbreak_Risk_95UR=outbreak_risk_poultry_farm_County_95;

Spillover=mle_spillover_poultry_farm_County;
Spillover_95UR=spillover_poultry_farm_County_95;

Spillover_Risk=mle_spillover_risk_poultry_farm_County;
Spillover_Risk_95UR=spillover_risk_poultry_farm_County_95;

Onward_Transmission=mle_onward_transmission_poultry_farm_County;
Onward_Transmission(isnan(Outbreak))=NaN;
Onward_Transmission_95UR=onward_transmission_poultry_farm_County_95;
Onward_Transmission_95UR(isnan(Outbreak),:)=NaN;

Poisson_Mean=mle_potential_outbreak_poultry_farm_County;
Poisson_Mean_95UR=potential_outbreak_poultry_farm_County_95;

T=table(State,County,GEOID,Outbreak,Outbreak_95UR,Outbreak_Risk,Outbreak_Risk_95UR,Spillover,Spillover_95UR,Spillover_Risk,Spillover_Risk_95UR,Onward_Transmission,Onward_Transmission_95UR,Poisson_Mean,Poisson_Mean_95UR);

writetable(T,'H5N1_Poultry_Risk_Model.xlsx','Sheet','County_Estimates');


State=State_Name;
Outbreak=mle_outbreak_poultry_farm_State;
Outbreak_95UR=outbreak_poultry_farm_State_95;

Outbreak_Risk=mle_outbreak_risk_poultry_farm_State;
Outbreak_Risk_95UR=outbreak_risk_poultry_farm_State_95;

Spillover=mle_spillover_poultry_farm_State;
Spillover_95UR=spillover_poultry_farm_State_95;

Spillover_Risk=mle_spillover_risk_poultry_farm_State;
Spillover_Risk_95UR=spillover_risk_poultry_farm_State_95;

Onward_Transmission=mle_onward_transmission_poultry_farm_State;
Onward_Transmission_95UR=onward_transmission_poultry_farm_State_95;

T=table(State,Outbreak,Outbreak_95UR,Outbreak_Risk,Outbreak_Risk_95UR,Spillover,Spillover_95UR,Spillover_Risk,Spillover_Risk_95UR,Onward_Transmission,Onward_Transmission_95UR);

writetable(T,'H5N1_Poultry_Risk_Model.xlsx','Sheet','State_Estimates');