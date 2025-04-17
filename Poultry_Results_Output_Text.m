clear;
clc;
rng(12500410)
load([pwd '/Data/Data_US_County.mat'],'US_County');
load('Average_Risk_Poultry.mat');

US_County=US_County(~no_farms,[4 5]);

avg_outbreak=outbreak_poultry_farm_County*w_AIC;
avg_outbreak_risk=outbreak_risk_poultry_farm_County*w_AIC;
avg_pot_outbreak=potential_outbreak_poultry_farm_County*w_AIC;
avg_spillover_per_outbreak=par_spillover*w_AIC;
avg_spillover_risk_County=spillover_risk_poultry_farm_County*w_AIC;

avg_outbreak_state=outbreak_poultry_farm_State*w_AIC;
avg_spillover_risk_State=spillover_risk_poultry_farm_State*w_AIC;

avg_spillover_risk_County(~no_farms);
avg_outbreak=avg_outbreak(~no_farms);
avg_outbreak_risk=avg_outbreak_risk(~no_farms);
avg_pot_outbreak=avg_pot_outbreak(~no_farms);
wc=cumsum(w_AIC);


r=rand(10^3,1);
f_indx=zeros(10^3,1);
for ii=1:length(r)
   f_indx(ii) = find(r(ii)<=wc,1,"first");
end

samp_outbreak=outbreak_poultry_farm_County(~no_farms,f_indx);
samp_outbreak_risk=outbreak_risk_poultry_farm_County(~no_farms,f_indx);
samp_pot_outbreak=potential_outbreak_poultry_farm_County(~no_farms,f_indx);
samp_outbreak_state=outbreak_poultry_farm_State(:,f_indx);
samp_spillover_per_outbreak=par_spillover(f_indx);
samp_spillover_risk_County=spillover_risk_poultry_farm_County(~no_farms,f_indx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range of county outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb=zeros(3,1);
ub=zeros(3,1);
md=zeros(3,1);

lb(1)=min(avg_outbreak);
ub(1)=max(avg_outbreak);
md(1)=median(avg_outbreak);

lb(2:3)=prctile(min(samp_outbreak,[],1),[2.5 97.5]);
ub(2:3)=prctile(max(samp_outbreak,[],1),[2.5 97.5]);
md(2:3)=prctile(median(samp_outbreak,1),[2.5 97.5]);

fprintf(['Minimum number of outbreaks across all counties: ' num2str(lb(1),'%4.3f') ' (95%% UR: ' num2str(lb(2),'%4.3f') char(8211) num2str(lb(3),'%4.3f') ') \n']);
fprintf(['Maximum number of outbreaks across all counties: ' num2str(ub(1),'%3.2f') ' (95%% UR: ' num2str(ub(2),'%3.2f') char(8211) num2str(ub(3),'%3.2f') ') \n']);
fprintf(['Median number of outbreaks across all counties: ' num2str(md(1),'%3.2f') ' (95%% UR: ' num2str(md(2),'%3.2f') char(8211) num2str(md(3),'%3.2f') ') \n \n']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Risk of outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
risk_t=0.5;

cc=zeros(3,1);

cc(1)=sum(avg_outbreak_risk>risk_t);

cc(2:3)=prctile(sum(samp_outbreak_risk>risk_t,1),[2.5 97.5]);

fprintf(['Number of counties with an outbreak risk over ' num2str(risk_t) ': ' num2str(cc(1),'%3.2f') ' (95%% UR: ' num2str(cc(2),'%3.2f') char(8211) num2str(cc(3),'%3.2f') ') \n \n']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% low risk of outbreaks but potential for large number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
risk_t=0.1;

cc=zeros(3,1);

cc(1)=prctile(avg_pot_outbreak,95);

cc(2:3)=prctile(prctile(samp_pot_outbreak,95),[2.5 97.5]);

fprintf(['95th prentile of the potential numebr of outbreaks: ' num2str(cc(1),'%3.1f') ' (95%% UR: ' num2str(cc(2),'%3.1f') char(8211) num2str(cc(3),'%3.1f') ') \n']);


cc(1)=mean(avg_outbreak_risk(avg_pot_outbreak>prctile(avg_pot_outbreak,95))<risk_t);

X_temp=zeros(size(samp_pot_outbreak,2),1);

for ii=1:length(X_temp)
    X_temp(ii)=mean(samp_outbreak_risk(samp_pot_outbreak(:,ii)>prctile(samp_pot_outbreak(:,ii),95),ii)<risk_t);
end

cc(2:3)=prctile(X_temp,[2.5 97.5]);



fprintf(['Of the counties above the 95th prentile of the potential numebr of outbreaks, the percentage that have an outbreak risk below ' num2str(risk_t) ': ' num2str(100.*cc(1),'%3.2f') '%% (95%% UR: ' num2str(100.*cc(2),'%3.2f') '%%' char(8211) num2str(100.*cc(3),'%3.2f') '%%) \n \n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,IndxS]=sort(avg_outbreak_state,'descend');

for ii=1:3
    cc(1)=avg_outbreak_state(IndxS(ii));

    cc(2:3)=prctile(samp_outbreak_state(IndxS(ii),:),[2.5 97.5]);
    fprintf(['Expectd number of outbreaks in ' State_Name{IndxS(ii)} ':' num2str(cc(1),'%3.1f') ' (95%% UR: ' num2str(cc(2),'%3.1f') char(8211) num2str(cc(3),'%3.1f') ') \n']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spill over
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
c=zeros(3,1);

cc(1)=100.*avg_spillover_per_outbreak;

cc(2:3)=100.*prctile(samp_spillover_per_outbreak,[2.5 97.5]);

fprintf(['Spillover events per 100 outbreaks: ' num2str(cc(1),'%3.2f') ' (95%% UR: ' num2str(cc(2),'%3.2f') char(8211) num2str(cc(3),'%3.2f') ') \n \n']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State spillover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,IndxS]=sort(avg_spillover_risk_State,'descend');
for ii=1:3
    fprintf(['State Rank ' num2str(ii) ' for spillover-risk: '  State_Name{IndxS(ii)} ' (' num2str(avg_spillover_risk_State(IndxS(ii)),'%3.2f') ') \n']);
end
 fprintf('\n');
[~,IndxS]=sort(avg_spillover_risk_County,'descend');

for ii=1:5
    fprintf(['County Rank ' num2str(ii) ' for spillover-risk: ' US_County.NAME{IndxS(ii)} ', ' US_County.STATE_NAME{IndxS(ii)}  ' (' num2str(avg_spillover_risk_County(IndxS(ii)),'%3.2f') ') \n']);
end

