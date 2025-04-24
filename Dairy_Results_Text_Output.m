clear;
clc;

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);

load('Dairy_Risk_AIC.mat');

US_County=US_County(~no_farms,[4 5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range of county outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb=zeros(3,1);
ub=zeros(3,1);
md=zeros(3,1);

lb(1)=min(mle_outbreak_dairy_farm_County(~no_farms));
ub(1)=max(mle_outbreak_dairy_farm_County(~no_farms));
md(1)=median(mle_outbreak_dairy_farm_County(~no_farms));

lb(2:3)=min_outbreaks_95;
ub(2:3)=max_outbreaks_95;
md(2:3)=median_outbreaks_95;

fprintf(['Minimum number of outbreaks across all counties: ' num2str(lb(1),'%4.3f') ' (95%% UR: ' num2str(lb(2),'%4.3f') char(8211) num2str(lb(3),'%4.3f') ') \n']);
fprintf(['Maximum number of outbreaks across all counties: ' num2str(ub(1),'%3.1f') ' (95%% UR: ' num2str(ub(2),'%3.1f') char(8211) num2str(ub(3),'%3.1f') ') \n']);
fprintf(['Median number of outbreaks across all counties: ' num2str(md(1),'%3.2f') ' (95%% UR: ' num2str(md(2),'%3.2f') char(8211) num2str(md(3),'%3.2f') ') \n \n']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Risk of outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
risk_t=0.5;

cc=zeros(3,1);

cc(1)=sum(mle_outbreak_risk_dairy_farm_County>risk_t);

cc(2:3)=County_outbreak_risk_over_50_95;

fprintf(['Number of counties with an outbreak risk over ' num2str(risk_t) ': ' num2str(cc(1),'%3.2f') ' (95%% UR: ' num2str(cc(2),'%3.2f') char(8211) num2str(cc(3),'%3.2f') ') \n \n']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% low risk of outbreaks but potential for large number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
risk_t=0.25;

cc=zeros(3,1);

cc(1)=prctile(mle_potntial_outbreak_dairy_farm_County(~no_farms),95);

cc(2:3)=potential_outbreak_95_prctl_95;

fprintf(['95th prentile of the potential numebr of outbreaks: ' num2str(cc(1),'%3.1f') ' (95%% UR: ' num2str(cc(2),'%3.1f') char(8211) num2str(cc(3),'%3.1f') ') \n']);

temp=mle_potntial_outbreak_dairy_farm_County>prctile(mle_potntial_outbreak_dairy_farm_County(~no_farms),95);
cc(1)=mean(mle_outbreak_risk_dairy_farm_County(temp)<risk_t);

cc(2:3)=per_risk_under_25_95;



fprintf(['Of the counties above the 95th prentile of the potential numebr of outbreaks, the percentage that have an outbreak risk below ' num2str(risk_t) ': ' num2str(100.*cc(1),'%3.2f') '%% (95%% UR: ' num2str(100.*cc(2),'%3.2f') '%%' char(8211) num2str(100.*cc(3),'%3.2f') '%%) \n \n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,IndxS]=sort(mle_outbreak_dairy_farm_State,'descend');

for ii=1:3
    cc(1)=mle_outbreak_dairy_farm_State(IndxS(ii));

    cc(2:3)=outbreak_dairy_farm_State_95(IndxS(ii),:);
    fprintf(['Expectd number of outbreaks in ' State_Name{IndxS(ii)} ':' num2str(cc(1),'%3.1f') ' (95%% UR: ' num2str(cc(2),'%3.1f') char(8211) num2str(cc(3),'%3.1f') ') \n']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spill over
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
c=zeros(3,1);

cc(1)=100.*mle_par_spillover;

cc(2:3)=100.*par_spillover_95;

fprintf(['Spillover events per 100 outbreaks: ' num2str(cc(1),'%3.2f') ' (95%% UR: ' num2str(cc(2),'%3.2f') char(8211) num2str(cc(3),'%3.2f') ') \n \n']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State spillover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,IndxS]=sort(mle_spillover_risk_dairy_farm_State,'descend');
for ii=1:3
    fprintf(['State Rank ' num2str(ii) ' for spillover-risk: '  State_Name{IndxS(ii)} ' (' num2str(mle_spillover_risk_dairy_farm_State(IndxS(ii)),'%3.2f') ': ' num2str(spillover_risk_dairy_farm_State_95(IndxS(ii),1),'%3.2f') char(8211) num2str(spillover_risk_dairy_farm_State_95(IndxS(ii),2),'%3.2f') ') \n']);
end
 fprintf('\n');

mle_spillover_risk_dairy_farm_County=mle_spillover_risk_dairy_farm_County(~no_farms);
spillover_risk_dairy_farm_County_95=spillover_risk_dairy_farm_County_95(~no_farms,:);
[~,IndxS]=sort(mle_spillover_risk_dairy_farm_County,'descend');
for ii=1:5
    fprintf(['County Rank ' num2str(ii) ' for spillover-risk: ' US_County.NAME{IndxS(ii)} ', ' US_County.STATE_NAME{IndxS(ii)}  ' (' num2str(mle_spillover_risk_dairy_farm_County(IndxS(ii)),'%3.2f') ':' num2str(spillover_risk_dairy_farm_County_95(IndxS(ii),1),'%3.2f') char(8211) num2str(spillover_risk_dairy_farm_County_95(IndxS(ii),2),'%3.2f') ') \n']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State zero outbreak likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=zeros(3,1);

cc(1)=mle_state_31_zero;

cc(2:3)=state_31_zero_95;

fprintf(['Likelihood of zero-outbreaks among the 31 states not reporting any outbraks: ' num2str(cc(1),'%3.2e') ' (95%% UR: ' num2str(cc(2),'%3.2e') char(8211) num2str(cc(3),'%3.2e') ') \n \n']);


[~,~,~,~,~,~,Affected_State_Farms,~,~,~,~,~,~]= Dairy_Covariates({},{},{});
State_Name_0=State_Name(Affected_State_Farms==0);
[mle_state_zero_0,IndxS]=sort(mle_state_zero(Affected_State_Farms==0),'descend');

state_zero_0=state_zero_95(Affected_State_Farms==0,:);
state_zero_0=state_zero_0(IndxS,:);
State_Name_0=State_Name_0(IndxS);

fprintf(['State with highest likelihood of no outbreaks: '  State_Name_0{1} ' (' num2str(mle_state_zero_0(1),'%3.2f') ': ' num2str(state_zero_0(1,1),'%3.2f') char(8211) num2str(state_zero_0(1,2),'%3.2f') ') \n']);
fprintf(['State with lowest likelihood of no outbreaks: '  State_Name_0{end} ' (' num2str(mle_state_zero_0(end),'%4.3f') ': ' num2str(state_zero_0(end,1),'%4.3f') char(8211) num2str(state_zero_0(end,2),'%4.3f') ') \n']);


