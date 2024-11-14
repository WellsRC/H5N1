function Figure_Change_Susceptibility_Risk_Dairy_Farm()

load('Average_Risk_Dairy.mat','avg_overall_risk_dairy_farm_State','avg_overall_risk_dairy_farm_County');
baseline_avg_overall_risk_dairy_farm_State=avg_overall_risk_dairy_farm_State;
baseline_avg_overall_risk_dairy_farm_County=avg_overall_risk_dairy_farm_County;

red_con=1;
load(['Average_Risk_Dairy_Reduction_Connectivity=' num2str(red_con.*100) '.mat'],'avg_overall_risk_dairy_farm_State','avg_overall_risk_dairy_farm_County','State_Name');
red_connectivity_avg_overall_risk_dairy_farm_County=avg_overall_risk_dairy_farm_County;

delta_avg_overall_risk_dairy_farm_County=1-red_connectivity_avg_overall_risk_dairy_farm_County./baseline_avg_overall_risk_dairy_farm_County;


red_con_v=[0.1:0.1:1];
red_connectivity_avg_overall_risk_dairy_farm_State=zeros(length(State_Name),length(red_con_v));
for ii=1:length(red_con_v)
    load(['Average_Risk_Dairy_Reduction_Connectivity=' num2str(red_con_v(ii).*100) '.mat'],'avg_overall_risk_dairy_farm_State','State_Name');
    red_connectivity_avg_overall_risk_dairy_farm_State(:,ii)=1-avg_overall_risk_dairy_farm_State./baseline_avg_overall_risk_dairy_farm_State;
end


% Plot_Variable
close all;
states = shaperead('usastatelo', 'UseGeoCoords', true);
S=shaperead([pwd '/Shapefile/cb_2021_us_county_500k.shp'],'UseGeoCoords',true);

US_County=readgeotable([pwd '/Shapefile/cb_2021_us_county_500k.shp']);
US_County=US_County(:,[2 3 6 7 10 9 12 13]);
Var_Name=US_County.Properties.VariableNames;
Var_Name{end-1}='AREA_LAND';
Var_Name{end}='AREA_WATER';
US_County.Properties.VariableNames=Var_Name;

[US_County,Indx]=sortrows(US_County,[1 2]);

S=S(Indx);

County_remove=strcmp("HI",US_County.STUSPS) | strcmp("AS",US_County.STUSPS) | strcmp("GU",US_County.STUSPS) | strcmp("MP",US_County.STUSPS) | strcmp("PR",US_County.STUSPS) | strcmp("VI",US_County.STUSPS);
S=S(~County_remove,:);

NS=length(S);


C_Risk=[hex2rgb('#ffffff');
        hex2rgb('#f7fcfd');
hex2rgb('#e5f5f9');
hex2rgb('#ccece6');
hex2rgb('#99d8c9');
hex2rgb('#66c2a4');
hex2rgb('#41ae76');
hex2rgb('#238b45');
hex2rgb('#006d2c');
hex2rgb('#00441b');];

figure('units','normalized','outerposition',[0 0.075 1 0.5]);
 ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,0.4,0.6,0.6];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;
       

risk_measure=delta_avg_overall_risk_dairy_farm_County;
t_upper=risk_measure>prctile(risk_measure,95);
v_upper=prctile(risk_measure,95);
t_lower=risk_measure<prctile(risk_measure,5);
v_lower=prctile(risk_measure,5);
risk_measure(t_lower)=v_lower;
risk_measure(t_upper)=v_upper;
risk_measure=log(risk_measure);
risk_measure=(risk_measure-min(risk_measure))./(max(risk_measure)-min(risk_measure));
x_risk=linspace(0,1,size(C_Risk,1));
c_indx=[1:size(C_Risk,1)];
    

subplot('Position',[0.38,0.05,0.01,0.93]);

xlim([0 1]);
ylim([0 max(c_indx)]);    
ymin=1;
dy=2/(1+sqrt(5));
for ii=1:length(c_indx)
    patch([0 0 dy dy],c_indx(ii)-[1 0 0 1],C_Risk(ii,:),'LineStyle','none');
end
        
text(ymin,0,'Lowest','Fontsize',16);
text(ymin,max(c_indx),'Highest','Fontsize',16);
text(ymin+2,max(c_indx)./2,{'Reduction in overall risk'},'HorizontalAlignment','center','Fontsize',18,'Rotation',270);

axis off;  

text(-37,0.98,'A','FontSize',32,'Units','normalized');
CC_Risk=ones(length(S),3);
        
for ii=1:length(S)
    if(~isnan(risk_measure(ii)))
        CC_Risk(ii,:)=interp1(x_risk,C_Risk,risk_measure(ii));
    else
        CC_Risk(ii,:)=[0.7 0.7 0.7];
    end
end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot uptake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});

geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% State-level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
state_nan=~isnan(mean(red_connectivity_avg_overall_risk_dairy_farm_State,2));
red_connectivity_avg_overall_risk_dairy_farm_State=red_connectivity_avg_overall_risk_dairy_farm_State(state_nan,:);
State_Name=State_Name(state_nan);

State_Risk=zeros(length(State_Name),1);
for ss=1:length(State_Name)
    State_Risk(ss)=red_connectivity_avg_overall_risk_dairy_farm_State(ss,:)*red_con_v(:);
end

[~,Indx_Risk]=sort(State_Risk,'descend');
State_Name=State_Name(Indx_Risk);
red_connectivity_avg_overall_risk_dairy_farm_State=red_connectivity_avg_overall_risk_dairy_farm_State(Indx_Risk,:);

red_connectivity_avg_overall_risk_dairy_farm_State=log(red_connectivity_avg_overall_risk_dairy_farm_State);

subplot('Position',[0.52,0.4,0.47,0.57]);

for uu=1:length(red_con_v)
    for ss=1:length(State_Name)
        patch(ss+[-0.5 -0.5 0.5 0.5],100.*(red_con_v(uu)+[-0.05 0.05 0.05 -0.05]),interp1(linspace(min(red_connectivity_avg_overall_risk_dairy_farm_State(:)),max(red_connectivity_avg_overall_risk_dairy_farm_State(:)),size(C_Risk,1)),C_Risk,red_connectivity_avg_overall_risk_dairy_farm_State(ss,uu)),'LineStyle','none');
    end
end

xlim([0.5 length(State_Name)+0.5])
ylim([5 105])
set(gca,'XTick',[1:length(State_Name)],'XTickLabel',State_Name,'YTick',100.*red_con_v,'Fontsize',14,'LineWidth',2,'TickDir','out', 'YDir','reverse');
ytickformat('percentage');

ylabel({'Reduction in','interstate connectivity'},'Fontsize',18)
xlabel('State','Fontsize',18)

text(-0.15,0.97,'B','FontSize',32,'Units','normalized');

ax1.Position=[-0.445,-0.2,1.3,1.3];


print(gcf,['Supplementary_Interstate_Transporation_Dairy.png'],'-dpng','-r300');
end