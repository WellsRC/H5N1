clear;
clc;
rng(12500410)
load('Average_Risk_Poultry.mat');

avg_outbreak=outbreak_poultry_farm_County*w_AIC;
avg_outbreak_risk=outbreak_risk_poultry_farm_County*w_AIC;
avg_pot_outbreak=potential_outbreak_poultry_farm_County*w_AIC;
wc=cumsum(w_AIC);


r=rand(10^4,1);
f_indx=zeros(10^4,1);
for ii=1:length(r)
   f_indx(ii) = find(r(ii)<=wc,1,"first");
end

samp_outbreak=outbreak_poultry_farm_County(:,f_indx);
samp_outbreak_risk=outbreak_risk_poultry_farm_County(:,f_indx);
samp_pot_outbreak=potential_outbreak_poultry_farm_County(:,f_indx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range of county outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb=zeros(3,1);
ub=zeros(3,1);
md=zeros(3,1);

lb(1)=min(avg_outbreak);
ub(1)=max(avg_outbreak);
md(1)=median(avg_outbreak);

lb(2:3)=prctile(min(samp_outbreak,[],2),[2.5 97.5]);
ub(2:3)=prctile(max(samp_outbreak,[],2),[2.5 97.5]);
md(2:3)=prctile(median(samp_outbreak,2),[2.5 97.5]);

fprintf(['Minimum number of outbreaks across all counties: ' num2str(lb(1),'%3.2f') ' (95% UR: ' num2str(lb(2),'%3.2f') char(8211) num2str(lb(3),'%3.2f') ') \n']);
fprintf(['Maximum number of outbreaks across all counties: ' num2str(ub(1),'%3.2f') ' (95% UR: ' num2str(ub(2),'%3.2f') char(8211) num2str(ub(3),'%3.2f') ') \n']);
fprintf(['Median number of outbreaks across all counties: ' num2str(md(1),'%3.2f') ' (95% UR: ' num2str(md(2),'%3.2f') char(8211) num2str(md(3),'%3.2f') ') \n \n']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Risk of outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
risk_t=0.5;

cc=zeros(3,1);

cc(1)=sum(avg_outbreak_risk>risk_t);

cc(2:3)=prctile(sum(samp_outbreak_risk>risk_t,2),[2.5 97.5]);

fprintf(['Number of counties with an outbrak risk over ' num2str(risk_t) ': ' num2str(cc(1),'%3.2f') ' (95% UR: ' num2str(cc(2),'%3.2f') char(8211) num2str(cc(3),'%3.2f') ') \n \n']);
