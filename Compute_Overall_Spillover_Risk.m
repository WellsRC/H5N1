clear;
clc;
load('Average_Risk_Dairy.mat','w_AIC','spillover_risk_dairy_farm_County','avg_overall_risk_dairy_farm_County','spillover_risk_dairy_farm_State');
d_AIC=w_AIC;
load('Average_Risk_Poultry.mat','w_AIC','spillover_risk_poultry_farm_County','avg_overall_risk_poultry_farm_County','spillover_risk_poultry_farm_State','State_Name');
p_AIC=w_AIC;

clear w_AIC;

[d_AIC,p_AIC]=meshgrid(d_AIC,p_AIC);

w_AIC=d_AIC.*p_AIC;
w_AIC=w_AIC./sum(w_AIC(:));

avg_spillover_risk_total_County=zeros(size(spillover_risk_poultry_farm_County,1),1);
avg_spillover_risk_total_State=zeros(size(spillover_risk_poultry_farm_State,1),1);

for dd=1:size(w_AIC,1)
    for pp=1:size(w_AIC,2)
        avg_spillover_risk_total_County=avg_spillover_risk_total_County+w_AIC(pp,dd).*(1-(1-spillover_risk_dairy_farm_County(:,dd)).*(1-spillover_risk_poultry_farm_County(:,pp)));
        avg_spillover_risk_total_State=avg_spillover_risk_total_State+w_AIC(pp,dd).*(1-(1-spillover_risk_dairy_farm_State(:,dd)).*(1-spillover_risk_poultry_farm_State(:,pp)));
    end
end


avg_spillover_risk_total_County(isnan(avg_overall_risk_dairy_farm_County) & isnan(avg_overall_risk_poultry_farm_County)) = NaN; % i.e. model suggests that since no farms/data cannot infer
save('Total_Spillover_Risk_County.mat','avg_spillover_risk_total_County','avg_spillover_risk_total_State','State_Name');


