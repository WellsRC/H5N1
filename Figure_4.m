% close all;
% states = shaperead('usastatelo', 'UseGeoCoords', true);
% S=shaperead([pwd '/Shapefile/cb_2021_us_county_500k.shp'],'UseGeoCoords',true);
% 
% US_County=readgeotable([pwd '/Shapefile/cb_2021_us_county_500k.shp']);
% US_County=US_County(:,[2 3 6 7 10 9 12 13]);
% Var_Name=US_County.Properties.VariableNames;
% Var_Name{end-1}='AREA_LAND';
% Var_Name{end}='AREA_WATER';
% US_County.Properties.VariableNames=Var_Name;
% 
% [US_County,Indx]=sortrows(US_County,[1 2]);
% 
% S=S(Indx);
% 
% County_remove=strcmp("AK",US_County.STUSPS) | strcmp("HI",US_County.STUSPS) | strcmp("AS",US_County.STUSPS) | strcmp("GU",US_County.STUSPS) | strcmp("MP",US_County.STUSPS) | strcmp("PR",US_County.STUSPS) | strcmp("VI",US_County.STUSPS);
% S=S(~County_remove,:);
% 
% NS=length(S);
% 
% load([pwd '/Data/Data_US_County.mat'],'US_County');
% 
% State_Name=unique(US_County.STATE_NAME);
% state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
% State_Name=State_Name(~state_remove);
% 
% 
% load('Average_Risk_Poultry.mat','post_spillover_poultry_farm_County','post_spillover_poultry_farm_State','no_farms','w_AIC');
% spillover_Poultry_County=squeeze(post_spillover_poultry_farm_County(:,:,w_AIC==max(w_AIC)));
% spillover_Poultry_State=squeeze(post_spillover_poultry_farm_State(:,:,w_AIC==max(w_AIC)));
% Poultry_no_farms=no_farms;
% 
% load('Average_Risk_Dairy.mat','post_spillover_dairy_farm_County','post_spillover_dairy_farm_State','no_farms','w_AIC');
% spillover_dairy_County=squeeze(post_spillover_dairy_farm_County(:,:,w_AIC==max(w_AIC)));
% spillover_dairy_State=squeeze(post_spillover_dairy_farm_State(:,:,w_AIC==max(w_AIC)));
% Dairy_no_farms=no_farms;
% 
% no_farms=Poultry_no_farms & Dairy_no_farms;
% 
% avg_onward_transmission_County=zeros(size(Dairy_no_farms));
% avg_onward_transmission_Dairy_County=zeros(size(Dairy_no_farms));
% avg_onward_transmission_State=zeros(size(spillover_Poultry_State,1),1);
% avg_onward_transmission_Dairy_State=zeros(size(spillover_Poultry_State,1),1);
% 
% 
% k_onward_transmission=2.69;
% R0=0.05;
% p_no_onward_transmission=nbinpdf(0,k_onward_transmission,k_onward_transmission./(k_onward_transmission+R0));
% 
% NK=zeros(size(Dairy_no_farms));
% NKS=zeros(size(spillover_Poultry_State,1),1);
% for pp=0:size(spillover_Poultry_County,2)-1
%     for dd=0:size(spillover_dairy_County,2)-1
%         NK=NK+(spillover_Poultry_County(:,pp+1).*spillover_dairy_County(:,dd+1));
%         NKS=NKS+(spillover_Poultry_State(:,pp+1).*spillover_dairy_State(:,dd+1));
%         avg_onward_transmission_County=avg_onward_transmission_County+(spillover_Poultry_County(:,pp+1).*spillover_dairy_County(:,dd+1)).*(1-p_no_onward_transmission.^(dd+pp));
%         avg_onward_transmission_Dairy_County=avg_onward_transmission_Dairy_County+(spillover_Poultry_County(:,pp+1).*spillover_dairy_County(:,dd+1)).*(1-p_no_onward_transmission.^(dd));
%         avg_onward_transmission_State=avg_onward_transmission_State+(spillover_Poultry_State(:,pp+1).*spillover_dairy_State(:,dd+1)).*(1-p_no_onward_transmission.^(dd+pp));
%         avg_onward_transmission_Dairy_State=avg_onward_transmission_Dairy_State+(spillover_Poultry_State(:,pp+1).*spillover_dairy_State(:,dd+1)).*(1-p_no_onward_transmission.^(dd));
%     end
% end
% avg_onward_transmission_County=avg_onward_transmission_County./NK;
% avg_onward_transmission_Dairy_County=avg_onward_transmission_Dairy_County./NK;
% avg_onward_transmission_State=avg_onward_transmission_State./NKS;
% avg_onward_transmission_Dairy_State=avg_onward_transmission_Dairy_State./NKS;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % Figure 4
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% figure('units','normalized','outerposition',[0.1 0.15 0.8 0.7]);
%  ax1=usamap('conus');
% 
% framem off; gridm off; mlabel off; plabel off;
% ax1.Position=[-0.3,0.4,0.6,0.6];
% 
% states = shaperead('usastatelo', 'UseGeoCoords', true);
% geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;
% 
% subplot('Position',[0.61,0.035,0.02,0.94]);
% 
% Title_Name={'Likelihood of onward transmission'};
% 
% risk_measure=log10(avg_onward_transmission_County);
% C_Risk=[hex2rgb('#fff7f3');
% hex2rgb('#fde0dd');
% hex2rgb('#fcc5c0');
% hex2rgb('#fa9fb5');
% hex2rgb('#f768a1');
% hex2rgb('#dd3497');
% hex2rgb('#ae017e');
% hex2rgb('#7a0177');
% hex2rgb('#49006a');];
% 
% y_indx=linspace(-5,0,size(C_Risk,1));
% x_risk=y_indx;
% 
%  x_indx=y_indx;
%  c_indx=linspace(y_indx(1),y_indx(end),1001);
% dx=c_indx(2)-c_indx(1);
% xlim([0 1]);
% ylim([y_indx(1) y_indx(end)+dx])
% ymin=0.75;
% dy=2/(1+sqrt(5));
% for ii=1:length(c_indx)
%     patch([0 0 dy dy],c_indx(ii)+[dx -dx -dx dx],interp1(x_risk,C_Risk,c_indx(ii)),'LineStyle','none');
% end
% 
% patch([0 0 dy dy], [y_indx(1) y_indx(end)+dx y_indx(end)+dx y_indx(1)],'k','FaceAlpha',0,'LineWidth',2);
% 
% yl_indx=[y_indx(1):y_indx(end)];
% for yy=1:length(yl_indx)
%     text(ymin,yl_indx(yy),['10^{' num2str(yl_indx(yy)) '}'],'Fontsize',11);            
% end
% 
% text(ymin+1.25,-2.5,Title_Name,'HorizontalAlignment','center','Fontsize',14,'Rotation',270);
% 
% axis off;  
% 
% NS=length(S);
% CC_Risk=interp1(x_risk,C_Risk,risk_measure);
% 
% CC_Risk(no_farms,:)=repmat([0.5 0.5 0.5],sum(no_farms),1);
% 
% CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
% geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
% geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;
% 
% subplot('Position',[0.75,0.085,0.24,0.905]);
% 
% 
% [~,indxs]=sort(avg_onward_transmission_State,'ascend');
% barh(State_Name(indxs),avg_onward_transmission_State(indxs),'k')
% box off;
% set(gca,'LineWidth',2,'tickdir','out','Fontsize',11,'XTick',[0:0.1:0.5],'XMinorTick','on');
% 
% xlabel({'Likelihood of onward transmission'},'Fontsize',14);
% text(-0.35,1.02,'B','Fontsize',28,'Units','normalized','VerticalAlignment','top','HorizontalAlignment','center');
% 
% ax1.Position=[-0.29,-0.12,1.2,1.2];
% text(-3.075,1.02,'A','Fontsize',28,'Units','normalized','VerticalAlignment','top','HorizontalAlignment','center');
% 
% print(gcf,['Figure_4.png'],'-dpng','-r300');
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % % Supplmental Figure
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
% 
% figure('units','normalized','outerposition',[0.1 0.15 0.8 0.7]);
%  ax1=usamap('conus');
% 
% framem off; gridm off; mlabel off; plabel off;
% ax1.Position=[-0.3,0.4,0.6,0.6];
% 
% states = shaperead('usastatelo', 'UseGeoCoords', true);
% geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;
% 
% subplot('Position',[0.61,0.035,0.02,0.94]);
% 
% Title_Name={'Relative contribution of dairy to onward transmission'};
% 
% risk_measure=(avg_onward_transmission_Dairy_County./avg_onward_transmission_County);
% C_Risk=[hex2rgb('#ffffff');
% hex2rgb('#f0f0f0');
% hex2rgb('#d9d9d9');
% hex2rgb('#bdbdbd');
% hex2rgb('#969696');
% hex2rgb('#737373');
% hex2rgb('#525252');
% hex2rgb('#252525');
% hex2rgb('#000000');];
% 
% y_indx=linspace(0,1,size(C_Risk,1));
% x_risk=y_indx;
% 
%  x_indx=y_indx;
%  c_indx=linspace(y_indx(1),y_indx(end),1001);
% dx=c_indx(2)-c_indx(1);
% xlim([0 1]);
% ylim([y_indx(1) y_indx(end)+dx])
% ymin=0.75;
% dy=2/(1+sqrt(5));
% for ii=1:length(c_indx)
%     patch([0 0 dy dy],c_indx(ii)+[dx -dx -dx dx],interp1(x_risk,C_Risk,c_indx(ii)),'LineStyle','none');
% end
% 
% patch([0 0 dy dy], [y_indx(1) y_indx(end)+dx y_indx(end)+dx y_indx(1)],'k','FaceAlpha',0,'LineWidth',2);
% 
% yl_indx=[0:0.1:1];
% for yy=1:length(yl_indx)
%     text(ymin,yl_indx(yy),[num2str(yl_indx(yy))],'Fontsize',11);            
% end
% 
% text(ymin+1.25,0.5,Title_Name,'HorizontalAlignment','center','Fontsize',14,'Rotation',270);
% 
% axis off;  
% 
% NS=length(S);
% CC_Risk=interp1(x_risk,C_Risk,risk_measure);
% 
% CC_Risk(no_farms,:)=repmat([0.5 0.5 0.5],sum(no_farms),1);
% 
% CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
% geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
% geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;
% 
% subplot('Position',[0.75,0.085,0.24,0.905]);
% 
% 
% [~,indxs]=sort(avg_onward_transmission_State,'ascend');
% b=barh(State_Name(indxs),[avg_onward_transmission_Dairy_State(indxs) avg_onward_transmission_State(indxs)-avg_onward_transmission_Dairy_State(indxs)],'stacked','LineStyle','none');
% b(1).FaceColor=[0 0 0];
% b(2).FaceColor=[0.7 0.7 0.7];

for ii=1:length(indxs)
    p=avg_onward_transmission_Dairy_State(indxs(ii))./avg_onward_transmission_State(indxs(ii));
    text(avg_onward_transmission_State(indxs(ii))+0.01,ii,[num2str(100.*p,'%3.1f') '%'],'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',10);
end
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',11,'XTick',[0:0.1:0.5],'XMinorTick','on');
legend({'Dairy','Poultry'},'Location','southeast','Fontsize',11);
legend boxoff;
xlabel({'Likelihood of onward transmission'},'Fontsize',14);
text(-0.35,1.02,'B','Fontsize',28,'Units','normalized','VerticalAlignment','top','HorizontalAlignment','center');

ax1.Position=[-0.29,-0.12,1.2,1.2];
text(-3.075,1.02,'A','Fontsize',28,'Units','normalized','VerticalAlignment','top','HorizontalAlignment','center');

print(gcf,['Supplemental_Figure_Dairy_Proportion.png'],'-dpng','-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Text output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fprintf(['Minimum likelihood of onward transmission: ' num2str(min(avg_onward_transmission_County(~no_farms)),'%3.2e') '\n'])
fprintf(['Maximum likelihood of onward transmission: ' num2str(max(avg_onward_transmission_County),'%4.3f') '\n'])
fprintf(['Median likelihood of onward transmission: ' num2str(median(avg_onward_transmission_County(~no_farms)),'%3.2e') '\n \n'])

US_County=US_County(~no_farms,[4 5]);

findx=find(avg_onward_transmission_County>0.01);

for jj=1:length(findx)
    fprintf(['County with onward transmission risk over 0.01: ' US_County.NAME{findx(jj)} ', ' US_County.STATE_NAME{findx(jj)}  ' (' num2str(avg_onward_transmission_County(findx(jj)),'%3.2f') ') \n']);
end

fprintf('\n');

fprintf(['State with the maximum likelihood: ' State_Name{indxs(end)} ' (' num2str(avg_onward_transmission_State(indxs(end)),'%4.3f') ') \n']);
fprintf(['State with the minimum likelihood: ' State_Name{indxs(1)} ' (' num2str(avg_onward_transmission_State(indxs(1)),'%4.3f') ') \n \n']);


fprintf(['Numebr of states with the risk ovr 0.1: ' num2str(sum(avg_onward_transmission_State>0.1),'%3.2f') ' \n']);
fprintf(['Numebr of states with the risk ovr 0.05: ' num2str(sum(avg_onward_transmission_State>0.05),'%3.2f') ' \n']);


p_dairy=avg_onward_transmission_Dairy_State./avg_onward_transmission_State;


fprintf(['Minumum dairy proportion: ' num2str(100.*min(p_dairy),'%3.2') '%% \n'])
fprintf(['Maximum dairy proportion: ' num2str(100.*max(p_dairy),'%3.2') '%% \n'])
fprintf(['Median dairy proportion: ' num2str(100.*median(p_dairy),'%3.2') '%% \n'])
