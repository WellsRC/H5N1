clear;
clc;
close all;

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Names=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Names,"Alaska") | strcmp(State_Names,"District of Columbia");
State_Names=State_Names(~state_remove);

clearvars US_County

[F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_par] = Poultry_Covariates({},{});

Affected_County_Farms(isnan(Affected_County_Farms))=0;
Outbreak_State=zeros(size(State_Spillover_Events));
for ss=1:length(Outbreak_State)
    Outbreak_State(ss)=state_weight_matrix(ss,:)*Affected_County_Farms;
end


load('Average_Risk_Poultry.mat','post_outbreak_poultry_farm_State','w_AIC');


[Outbreak_State,R_Indx]=sort(Outbreak_State,'descend');
post_outbreak_poultry_farm_State=post_outbreak_poultry_farm_State(R_Indx,:,:);
State_Names=State_Names(R_Indx);


nr=[4 4 4 4];
for pp=1:4
    figure('units','normalized','outerposition',[0.2 0.06 0.8 1]);
    for ii=1:nr(pp)
        for jj=1:3
            Outbreak_Post_State=squeeze(post_outbreak_poultry_farm_State(jj+3.*(ii-1)+12.*(pp-1),:,w_AIC==max(w_AIC)));
            subplot('Position',[0.025+0.3275.*(jj-1) 0.80-0.24.*(ii-1) 0.30 0.17]);
            bar([0:length(Outbreak_Post_State)-1],Outbreak_Post_State,'FaceColor',hex2rgb('#011A27'));
            hold on
            mx=max(Outbreak_Post_State);
            plot([Outbreak_State(jj+3.*(ii-1)+12.*(pp-1))  Outbreak_State(jj+3.*(ii-1)+12.*(pp-1))],[0 1.5],'color',hex2rgb('#F0810F'),'LineWidth',2);
            tempx=find(cumsum(Outbreak_Post_State)>=0.995,1,'first')-1;
            f_75=find(cumsum(Outbreak_Post_State)>=0.75,1,"first");
            mxl=max([f_75 5+Outbreak_State(jj+3.*(ii-1)+12.*(pp-1))]);
            
            mxl=max(mxl,10);
            dxtick=1;
            if(mxl>10)
                mxl=ceil(mxl./10).*10;
                if(mxl<=50)
                    dxtick=5;
                elseif(mxl<150)
                    dxtick=10;
                elseif(mxl<=300)
                    dxtick=25;
                else
                    dxtick=75;
                end
            end
            xlim([-0.5 mxl+0.5])
            ylim([0 mx.*1.1]);
            dx=[0.01 0.025 0.05 0.075 0.1 0.15 0.2];
            dx=dx(find(abs((ceil(10.*mx))./10./5-dx)==min(abs((ceil(10.*mx))./10./5-dx))));
            set(gca,'LineWidth',2,'TickDir','out','Fontsize',16,'YTick',[],'XTick',[0:dxtick:mxl]);
            if(ii==nr(pp))
                xlabel('Number of outbreaks','FontSize',18)
            end
            if(jj==1)
                ylabel('Density','FontSize',18)
            end
            box off;
            xtickangle(0);
            title(State_Names(jj+3.*(ii-1)+12.*(pp-1)))
        end
    end
    print(gcf,['Model_Fit_Poultry_State_' num2str(pp) '.png'],'-dpng','-r300');
end