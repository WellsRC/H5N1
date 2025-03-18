clear;


load([pwd '/Data/Data_US_County.mat'],'US_County');
US_County=US_County(:,[5,4,3]);
Filler=table({'N/A'},{'N/A'},{'Model_Weight'});
Filler.Properties.VariableNames=US_County.Properties.VariableNames;
US_County_temp=[Filler;US_County];

load('Average_Risk_Poultry.mat','w_AIC','State_Name','outbreak_poultry_farm_County','outbreak_risk_poultry_farm_County','outbreak_layer_farm_County','outbreak_pullet_farm_County','outbreak_broiler_farm_County','outbreak_turkey_farm_County','outbreak_risk_layer_farm_County','outbreak_risk_pullet_farm_County','outbreak_risk_broiler_farm_County','outbreak_risk_turkey_farm_County','spillover_poultry_farm_County','spillover_risk_poultry_farm_County','outbreak_poultry_farm_State','outbreak_risk_poultry_farm_State','spillover_poultry_farm_State','spillover_risk_poultry_farm_State','outbreak_layer_farm_State','outbreak_pullet_farm_State','outbreak_broiler_farm_State','outbreak_turkey_farm_State','outbreak_risk_layer_farm_State','outbreak_risk_pullet_farm_State','outbreak_risk_broiler_farm_State','outbreak_risk_turkey_farm_State');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
US_County=US_County_temp;
temp_v=[1; outbreak_poultry_farm_County*w_AIC];
US_County.Model_Average=temp_v;

temp_v=[w_AIC'; outbreak_poultry_farm_County];
temp_Table=array2table(temp_v);

for mm=1:length(w_AIC)
    temp_Table.Properties.VariableNames{mm}=['Model_' num2str(mm)];
end

US_County=[US_County temp_Table];

writetable(US_County,'H5N1_Poultry_Risk_Model.xlsx','Sheet','Outbreak');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outbreak Risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
US_County=US_County_temp;
temp_v=[1; outbreak_risk_poultry_farm_County*w_AIC];
US_County.Model_Average=temp_v;

temp_v=[w_AIC'; outbreak_risk_poultry_farm_County];
temp_Table=array2table(temp_v);

for mm=1:length(w_AIC)
    temp_Table.Properties.VariableNames{mm}=['Model_' num2str(mm)];
end

US_County=[US_County temp_Table];

writetable(US_County,'H5N1_Poultry_Risk_Model.xlsx','Sheet','Outbreak_Risk');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spillover 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
US_County=US_County_temp;
temp_v=[1; spillover_poultry_farm_County*w_AIC];
US_County.Model_Average=temp_v;

temp_v=[w_AIC'; spillover_poultry_farm_County];
temp_Table=array2table(temp_v);

for mm=1:length(w_AIC)
    temp_Table.Properties.VariableNames{mm}=['Model_' num2str(mm)];
end

US_County=[US_County temp_Table];

writetable(US_County,'H5N1_Poultry_Risk_Model.xlsx','Sheet','Spillover');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spillover Risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
US_County=US_County_temp;
temp_v=[1; spillover_risk_poultry_farm_County*w_AIC];
US_County.Model_Average=temp_v;

temp_v=[w_AIC'; spillover_risk_poultry_farm_County];
temp_Table=array2table(temp_v);

for mm=1:length(w_AIC)
    temp_Table.Properties.VariableNames{mm}=['Model_' num2str(mm)];
end

US_County=[US_County temp_Table];

writetable(US_County,'H5N1_Poultry_Risk_Model.xlsx','Sheet','Spillover_Risk');