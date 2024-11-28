clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% COVID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
load('Average_Risk_Population_COVID.mat','avg_susceptible_risk_population_County_COVID');

load('Average_Risk_Dairy.mat','w_AIC','spillover_risk_dairy_farm_County','avg_overall_risk_dairy_farm_County','spillover_risk_dairy_farm_State');
d_AIC=w_AIC;
load('Average_Risk_Poultry.mat','w_AIC','spillover_risk_poultry_farm_County','avg_overall_risk_poultry_farm_County','spillover_risk_poultry_farm_State','State_Name');
p_AIC=w_AIC;

load([pwd '/Data/Data_US_County.mat'],'US_County');
clear w_AIC;

[d_AIC,p_AIC]=meshgrid(d_AIC,p_AIC);

w_AIC=d_AIC.*p_AIC;
w_AIC=w_AIC./sum(w_AIC(:));

avg_spillover_risk_total_County=zeros(size(spillover_risk_poultry_farm_County,1),1);
avg_spillover_risk_total_State=zeros(size(spillover_risk_poultry_farm_State,1),1);

avg_localized_transmission_risk_total_State=zeros(size(spillover_risk_poultry_farm_State,1),1);


for dd=1:size(w_AIC,1)
    for pp=1:size(w_AIC,2)
        avg_spillover_risk_total_County=avg_spillover_risk_total_County+w_AIC(pp,dd).*(1-(1-spillover_risk_dairy_farm_County(:,dd)).*(1-spillover_risk_poultry_farm_County(:,pp)));
        avg_spillover_risk_total_State=avg_spillover_risk_total_State+w_AIC(pp,dd).*(1-(1-spillover_risk_dairy_farm_State(:,dd)).*(1-spillover_risk_poultry_farm_State(:,pp)));        
    end
end

avg_localized_transmission_risk_total_County=avg_spillover_risk_total_County.*avg_susceptible_risk_population_County_COVID;

for ss=1:length(State_Name)
    t_state=strcmp(State_Name{ss},US_County.STATE_NAME);

    c_r=avg_localized_transmission_risk_total_County(t_state);
    w_c=US_County.POPULATION_SIZE_2022(t_state);
    t_inc=w_c>0 & ~isnan(c_r);
    c_r=c_r(t_inc);
    w_c=w_c(t_inc);
    w_c=w_c./sum(w_c);
    if(~isempty(c_r))
        avg_localized_transmission_risk_total_State(ss)=1-exp(sum(w_c.*log(1-c_r)));
    end
end


avg_spillover_risk_total_County(isnan(avg_overall_risk_dairy_farm_County) & isnan(avg_overall_risk_poultry_farm_County)) = NaN; % i.e. model suggests that since no farms/data cannot infer
avg_localized_transmission_risk_total_County(isnan(avg_overall_risk_dairy_farm_County) & isnan(avg_overall_risk_poultry_farm_County)) = NaN;
save('Total_Spillover_Risk_County_COVID.mat','avg_spillover_risk_total_County','avg_spillover_risk_total_State','State_Name','avg_localized_transmission_risk_total_State','avg_localized_transmission_risk_total_County');

clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% H1N1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
load('Average_Risk_Population_H1N1.mat','avg_susceptible_risk_population_County_H1N1');

load('Average_Risk_Dairy.mat','w_AIC','spillover_risk_dairy_farm_County','avg_overall_risk_dairy_farm_County','spillover_risk_dairy_farm_State');
d_AIC=w_AIC;
load('Average_Risk_Poultry.mat','w_AIC','spillover_risk_poultry_farm_County','avg_overall_risk_poultry_farm_County','spillover_risk_poultry_farm_State','State_Name');
p_AIC=w_AIC;

load([pwd '/Data/Data_US_County.mat'],'US_County');
clear w_AIC;

[d_AIC,p_AIC]=meshgrid(d_AIC,p_AIC);

w_AIC=d_AIC.*p_AIC;
w_AIC=w_AIC./sum(w_AIC(:));

avg_spillover_risk_total_County=zeros(size(spillover_risk_poultry_farm_County,1),1);
avg_spillover_risk_total_State=zeros(size(spillover_risk_poultry_farm_State,1),1);

avg_localized_transmission_risk_total_State=zeros(size(spillover_risk_poultry_farm_State,1),1);


for dd=1:size(w_AIC,1)
    for pp=1:size(w_AIC,2)
        avg_spillover_risk_total_County=avg_spillover_risk_total_County+w_AIC(pp,dd).*(1-(1-spillover_risk_dairy_farm_County(:,dd)).*(1-spillover_risk_poultry_farm_County(:,pp)));
        avg_spillover_risk_total_State=avg_spillover_risk_total_State+w_AIC(pp,dd).*(1-(1-spillover_risk_dairy_farm_State(:,dd)).*(1-spillover_risk_poultry_farm_State(:,pp)));        
    end
end

avg_localized_transmission_risk_total_County=avg_spillover_risk_total_County.*avg_susceptible_risk_population_County_H1N1;

for ss=1:length(State_Name)
    t_state=strcmp(State_Name{ss},US_County.STATE_NAME);

    c_r=avg_localized_transmission_risk_total_County(t_state);
    w_c=US_County.POPULATION_SIZE_2022(t_state);
    t_inc=w_c>0 & ~isnan(c_r);
    c_r=c_r(t_inc);
    w_c=w_c(t_inc);
    w_c=w_c./sum(w_c);
    if(~isempty(c_r))
        avg_localized_transmission_risk_total_State(ss)=1-exp(sum(w_c.*log(1-c_r)));
    end
end

avg_spillover_risk_total_County(isnan(avg_overall_risk_dairy_farm_County) & isnan(avg_overall_risk_poultry_farm_County)) = NaN; % i.e. model suggests that since no farms/data cannot infer
avg_localized_transmission_risk_total_County(isnan(avg_overall_risk_dairy_farm_County) & isnan(avg_overall_risk_poultry_farm_County)) = NaN;
save('Total_Spillover_Risk_County_H1N1.mat','avg_spillover_risk_total_County','avg_spillover_risk_total_State','State_Name','avg_localized_transmission_risk_total_State','avg_localized_transmission_risk_total_County');



