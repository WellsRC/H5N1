clear;
clc;
rng(12500410)
load('Average_Risk_Poultry.mat');

avg_outbreak=outbreak_poultry_farm_County*w_AIC;
avg_outbreak_risk=outbreak_risk_poultry_farm_County*w_AIC;
avg_pot_outbreak=potential_outbreak_poultry_farm_County*w_AIC;
avg_outbreak_state=outbreak_poultry_farm_State*w_AIC;
avg_spillover_per_outbreak=par_spillover*w_AIC;

temp_spillover=sum(outbreak_poultry_farm_County.*repmat(par_spillover,size(outbreak_poultry_farm_County,1),1),1);

avg_spillover=temp_spillover*w_AIC;


avg_outbreak=avg_outbreak(~no_farms);
avg_outbreak_risk=avg_outbreak_risk(~no_farms);
avg_pot_outbreak=avg_pot_outbreak(~no_farms);
wc=cumsum(w_AIC);


r=rand(10^4,1);
f_indx=zeros(10^4,1);
for ii=1:length(r)
   f_indx(ii) = find(r(ii)<=wc,1,"first");
end

samp_outbreak=outbreak_poultry_farm_County(~no_farms,f_indx);
samp_outbreak_risk=outbreak_risk_poultry_farm_County(~no_farms,f_indx);
samp_pot_outbreak=potential_outbreak_poultry_farm_County(~no_farms,f_indx);
samp_outbreak_state=outbreak_poultry_farm_State(:,f_indx);
samp_spillover_per_outbreak=par_spillover(f_indx);
samp_spillover=temp_spillover(f_indx);

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

fprintf(['Minimum number of outbreaks across all counties: ' num2str(lb(1),'%3.2f') ' (95%% UR: ' num2str(lb(2),'%3.2f') char(8211) num2str(lb(3),'%3.2f') ') \n']);
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

X_temp=zeros(size(samp_pot_outbreak,1),1);

for ii=1:length(X_temp)
    X_temp(ii)=mean(samp_outbreak_risk(samp_pot_outbreak(:,ii)>prctile(samp_pot_outbreak(:,ii),95),ii)<risk_t);
end

cc(2:3)=prctile(X_temp,[2.5 97.5]);



fprintf(['Of the counties above the 95th prentile of the potential numebr of outbreaks, the percentage that have an outbreak risk below ' num2str(risk_t) ': ' num2str(100.*cc(1),'%3.2f') '%% (95%% UR: ' num2str(100.*cc(2),'%3.1f') '%%' char(8211) num2str(100.*cc(3),'%3.1f') '%%) \n \n']);

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
% State outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c=zeros(3,1);

cc(1)=avg_spillover_per_outbreak;

cc(2:3)=prctile(samp_spillover_per_outbreak,[2.5 97.5]);

fprintf(['Spillover events per outbreak: ' num2str(cc(1),'%4.3f') ' (95%% UR: ' num2str(cc(2),'%4.3f') char(8211) num2str(cc(3),'%4.3f') ') \n \n']);

cc(1)=avg_spillover;

cc(2:3)=prctile(samp_spillover,[2.5 97.5]);

fprintf(['Spillover events natioanlly: ' num2str(cc(1),'%4.3f') ' (95%% UR: ' num2str(cc(2),'%4.3f') char(8211) num2str(cc(3),'%4.3f') ') \n \n']);

