clear;
clc;
close all;

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Names=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Names,"Alaska") | strcmp(State_Names,"District of Columbia");
State_Names=State_Names(~state_remove);

clearvars US_County

[X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_par] = Poultry_Covariates({},{},{});

Affected_County_Farms(isnan(Affected_County_Farms))=0;


load('Average_Risk_Poultry.mat','post_spillover_poultry_farm_State','w_AIC');


[State_Spillover_Events,R_Indx]=sort(State_Spillover_Events,'descend');
post_spillover_poultry_farm_State=post_spillover_poultry_farm_State(R_Indx,:,:);
State_Names=State_Names(R_Indx);

nr=[5 5 5 1];
for pp=1:4
    figure('units','normalized','outerposition',[0.2 0.06 0.8 1]);
    for ii=1:nr(pp)
        for jj=1:3
            Outbreak_Post_State=squeeze(post_spillover_poultry_farm_State(jj+3.*(ii-1)+15.*(pp-1),:,:))*w_AIC;
            subplot('Position',[0.06+0.33.*(jj-1) 0.87-0.1975.*(ii-1) 0.27 0.1]);
            bar([0:100],Outbreak_Post_State,'k');
            hold on
            mx=max(Outbreak_Post_State);
            plot([State_Spillover_Events(jj+3.*(ii-1)+15.*(pp-1))  State_Spillover_Events(jj+3.*(ii-1)+15.*(pp-1))],[0 1.5],'r','LineWidth',2);
            tempx=find(cumsum(Outbreak_Post_State)>=0.995,1,'first')-1;
            if(isempty(tempx))
                tempx=min(400,State_Spillover_Events(jj+3.*(ii-1)+15.*(pp-1))+5);
            end
            mxl=max([tempx 5+State_Spillover_Events(jj+3.*(ii-1)+15.*(pp-1))]);
            xlim([-0.5 mxl+0.5])
            ylim([0 mx.*1.1]);
            dx=[0.01 0.025 0.05 0.075 0.1 0.15 0.2];
            dx=dx(find(abs((ceil(10.*mx))./10./5-dx)==min(abs((ceil(10.*mx))./10./5-dx))));
            set(gca,'LineWidth',2,'TickDir','out','Fontsize',16,'YTick',[0:dx:1]);
            xlabel('Number of spillover events','FontSize',18)
            ylabel('Density','FontSize',18)
            box off;
            title(State_Names(jj+3.*(ii-1)+15.*(pp-1)))
        end
    end
end