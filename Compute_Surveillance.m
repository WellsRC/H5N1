clear;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Average_Risk_Population_H1N1.mat','avg_susceptible_risk_population_County_H1N1','avg_susceptible_risk_population_State_H1N1');
load('Total_Spillover_Risk_County_H1N1.mat','avg_spillover_risk_total_County','avg_spillover_risk_total_State','State_Name');
load([pwd '/Data/Influenza_Testing/State_Level_Influenza_Testing.mat'])
load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Surveillance=NaN.*zeros(size(State_Name));

County_Surveillance=(avg_susceptible_risk_population_County_H1N1.*avg_spillover_risk_total_County).^(US_County_Influenza_Test.TEST_PER_CAPITA);


for ss=1:length(State_Name)
    t_state=strcmp(State_Name{ss},US_County.STATE_NAME);

    c_r=County_Surveillance(t_state);
    w_c=US_County.POPULATION_SIZE_2022(t_state);
    t_inc=w_c>0 & ~isnan(c_r);
    c_r=c_r(t_inc);
    w_c=w_c(t_inc);
    w_c=w_c./sum(w_c);
    if(~isempty(c_r))
        State_Surveillance(ss)=1-exp(sum(w_c.*log(1-c_r)));
    end
end

save('Surveillance_H5N1_Population_H1N1.mat','County_Surveillance','State_Surveillance','State_Name');