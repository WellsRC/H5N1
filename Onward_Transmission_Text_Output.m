clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);

load('Onward_Transmission.mat')

US_County=US_County(~isnan(mle_County_onward),[4 5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range of county outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb=zeros(3,1);
ub=zeros(3,1);
md=zeros(3,1);

lb(1)=min(mle_County_onward(~isnan(mle_County_onward)));
ub(1)=max(mle_County_onward(~isnan(mle_County_onward)));
md(1)=median(mle_County_onward(~isnan(mle_County_onward)));

lb(2:3)=min_County_onward_95;
ub(2:3)=max_County_onward_95;
md(2:3)=median_County_onward_95;

fprintf(['Minimum number of onward transmission events  across all counties: ' num2str(lb(1),'%3.2e') ' (95%% UR: ' num2str(lb(2),'%3.2e') char(8211) num2str(lb(3),'%3.2e') ') \n']);
fprintf(['Maximum number of onward transmission events  across all counties: ' num2str(ub(1),'%4.3f') ' (95%% UR: ' num2str(ub(2),'%4.3f') char(8211) num2str(ub(3),'%4.3f') ') \n']);
fprintf(['Median number of onward transmission events  across all counties: ' num2str(md(1),'%3.2f') ' (95%% UR: ' num2str(md(2),'%3.2f') char(8211) num2str(md(3),'%3.2f') ') \n \n']);


