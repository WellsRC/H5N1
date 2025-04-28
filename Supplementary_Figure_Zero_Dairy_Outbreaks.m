clear;
clc;
close all;

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);

load('Dairy_Risk_AIC.mat');

US_County=US_County(~no_farms,[4 5]);

[~,~,~,~,~,~,Affected_State_Farms,~,~,~,~,~,~]= Dairy_Covariates({},{},{});
State_Name_0=State_Name(Affected_State_Farms==0);
[mle_state_zero_0,IndxS]=sort(mle_state_zero(Affected_State_Farms==0),'descend');

state_zero_0=state_zero_95(Affected_State_Farms==0,:);
state_zero_0=state_zero_0(IndxS,:);
State_Name_0=State_Name_0(IndxS);

figure('units','normalized','outerposition',[0.2 0.3 0.6 0.7])
subplot('Position',[0.08 0.2 0.9 0.75])

bar(State_Name_0,mle_state_zero_0,'k','LineStyle','none');
set(gca,'LineWidth',2,'TickDir','out','YTick',[0:0.1:1],'Fontsize',14)
ylabel('Likelihood of having zero outbreaks among dairy herds')
box off

print(gcf,['Likelihood_Zero_Dairy_Outbreaks.png'],'-dpng','-r300');