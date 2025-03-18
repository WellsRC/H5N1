clear;
clc;
close all;

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Names=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Names,"Alaska") | strcmp(State_Names,"District of Columbia");
State_Names=State_Names(~state_remove);

clearvars US_County

N_Samp=10^4;

[X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect,logic_par]= Dairy_Covariates({},{},{});
Outbreak_State=zeros(size(State_Spillover_Events));
for ss=1:length(Outbreak_State)
    Outbreak_State(ss)=state_weight_matrix(ss,:)*Affected_County_Farms;
end
Outbreak_State=Outbreak_State+Affected_State_Farms;


% [avg_outbreak_farm_County,avg_outbreak_risk_farm_County,avg_spillover_risk_farm_County,avg_outbreak_risk_farm_State,avg_spillover_risk_farm_State,Spillover_Post_State,Outbreak_Post_State]=Temporary_Results_Remove_After();
% 
% save('TEMP_FIT_REMOVE.mat','Spillover_Post_State','Outbreak_Post_State');
load('TEMP_FIT_REMOVE.mat','Spillover_Post_State');

[State_Spillover_Events,R_Indx]=sort(State_Spillover_Events,'descend');
Spillover_Post_State=Spillover_Post_State(R_Indx,:);
State_Names=State_Names(R_Indx);

nr=[5 5 5 1];
for pp=1:4
    figure('units','normalized','outerposition',[0.2 0.06 0.8 1]);
    for ii=1:nr(pp)
        for jj=1:3
            subplot('Position',[0.06+0.33.*(jj-1) 0.87-0.1975.*(ii-1) 0.27 0.1]);
            bar([0:25],Spillover_Post_State(jj+3.*(ii-1)+15.*(pp-1),:),'k');
            hold on
            mx=max(Spillover_Post_State(jj+3.*(ii-1)+15.*(pp-1),:));
            plot([State_Spillover_Events(jj+3.*(ii-1)+15.*(pp-1))  State_Spillover_Events(jj+3.*(ii-1)+15.*(pp-1))],[0 1.5],'r','LineWidth',2);
            tempx=find(cumsum(Spillover_Post_State(jj+3.*(ii-1)+15.*(pp-1),:))>=0.995,1,'first')-1;
            if(isempty(tempx))
                tempx=400;
            end
            mxl=max([tempx 5+State_Spillover_Events(jj+3.*(ii-1)+15.*(pp-1))]);
            xlim([-0.5 mxl+0.5])
            ylim([0 (ceil(10.*mx))./10]);
            dx=[0.01 0.025 0.05 0.075 0.1 0.15 0.2];
            dx=dx(find(abs((ceil(10.*mx))./10./5-dx)==min(abs((ceil(10.*mx))./10./5-dx))));
            set(gca,'LineWidth',2,'TickDir','out','Fontsize',16,'YTick',[0:dx:1]);
            xlabel('Number of outbreaks','FontSize',18)
            ylabel('Density','FontSize',18)
            box off;
            title(State_Names(jj+3.*(ii-1)+15.*(pp-1)))
        end
    end
end