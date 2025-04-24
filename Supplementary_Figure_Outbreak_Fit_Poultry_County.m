clear;
clc;
close all;

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Names=US_County.STATE_NAME;
state_remove=strcmp(State_Names,"Alaska") | strcmp(State_Names,"District of Columbia");
State_Names=State_Names(~state_remove);

UState_Names=unique(US_County.STATE_NAME);
state_remove=strcmp(UState_Names,"Alaska") | strcmp(UState_Names,"District of Columbia");
UState_Names=UState_Names(~state_remove);
clearvars US_County

[F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_par,logic_temp] = Poultry_Covariates({},{});

Affected_County_Farms(isnan(Affected_County_Farms))=0;
Outbreak_State=zeros(size(State_Spillover_Events));
for ss=1:length(Outbreak_State)
    Outbreak_State(ss)=state_weight_matrix(ss,:)*Affected_County_Farms;
end


load('Average_Risk_Poultry.mat','post_outbreak_poultry_farm_County_CI','w_AIC');

post_outbreak_poultry_farm_County_CI=post_outbreak_poultry_farm_County_CI(:,:,w_AIC==max(w_AIC));


[Outbreak_State,R_Indx]=sort(Outbreak_State,'descend');
UState_Names=UState_Names(R_Indx);


nr=[4 4 4 4];
for pp=1:4
    figure('units','normalized','outerposition',[0.2 0.06 0.8 1]);
    for ii=1:nr(pp)
        for jj=1:3
            f_state=strcmp(UState_Names{jj+3.*(ii-1)+12.*(pp-1)},State_Names);
            post_outbreak_County=post_outbreak_poultry_farm_County_CI(f_state,:);
            AfC=Affected_County_Farms(f_state);
            
            [AfC,SAFc]=sort(AfC,'descend');
            post_outbreak_County=post_outbreak_County(SAFc,:);
            subplot('Position',[0.065+0.32.*(jj-1) 0.77-0.24.*(ii-1) 0.29 0.2]);
            for cc=1:size(post_outbreak_County,1)
                patch(cc+[-0.45 -0.45 0.45 0.45],[post_outbreak_County(cc,1) post_outbreak_County(cc,5) post_outbreak_County(cc,5) post_outbreak_County(cc,1)],hex2rgb('#011A27'),'FaceAlpha',0.33,'linestyle','none'); hold on;                
                plot(cc+[-0.45 0.45],(squeeze(post_outbreak_County(cc,3))).*ones(1,2),'-','Color',hex2rgb('#011A27'),'LineWidth',2);
            end
            scatter(1:size(post_outbreak_County,1),AfC,20,hex2rgb('#F0810F'),'filled');
            xlim([0.55 size(post_outbreak_County,1)+0.45]);

            mxl=ceil(max(max(AfC),max(post_outbreak_County(:))).*1.01);
            dxtick=1;
            if(mxl>=10)
                mxl=ceil(mxl./10).*10;
                if(mxl<20)
                    dxtick=2;
                elseif(mxl<=45)
                    dxtick=5;
                elseif(mxl<50)
                    dxtick=10;
                elseif(mxl<=150)
                    dxtick=25;
                else
                    dxtick=75;
                end
            end
            ylim([0 mxl])
            if(ii==nr(pp))
                xlabel('State counties','FontSize',18)
            end
            if(jj==1)
                ylabel({'Number of','outbreaks'},'FontSize',18)
            end
            set(gca,'LineWidth',2,'TickDir','out','Fontsize',16,'XTick',[],'YTick',[0:dxtick:mxl]);
            box off;
            title(UState_Names(jj+3.*(ii-1)+12.*(pp-1)))
        end
    end
    print(gcf,['Model_Fit_Poultry_County_' num2str(pp) '.png'],'-dpng','-r300');
end