function  Figure_Susceptible_Risk_Map(Var_Plot)
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

[~,Indx]=sortrows(US_County,[1 2]);

S=S(Indx);

County_remove=strcmp("HI",US_County.STUSPS) | strcmp("AS",US_County.STUSPS) | strcmp("GU",US_County.STUSPS) | strcmp("MP",US_County.STUSPS) | strcmp("PR",US_County.STUSPS) | strcmp("VI",US_County.STUSPS);
S=S(~County_remove,:);


NS=length(S);


State_Name=unique({S.STUSPS});
avg_susceptible_risk_total_State=NaN.*zeros(size(State_Name));

load([pwd '/Data/Data_US_County.mat'],'US_County');

if(strcmp(Var_Plot,'COVID'))
    load('Average_Risk_Population_COVID.mat','avg_susceptible_risk_population_County_COVID','avg_susceptible_risk_population_State_COVID','State_Name');

    avg_susceptible_risk_total_County=avg_susceptible_risk_population_County_COVID;
    avg_susceptible_risk_total_State=avg_susceptible_risk_population_State_COVID;

elseif(strcmp(Var_Plot,'H1N1'))
    load('Average_Risk_Population_COVID.mat','avg_susceptible_risk_population_County_H1N1','avg_susceptible_risk_population_State_H1N1','State_Name');
    avg_susceptible_risk_total_County=avg_susceptible_risk_population_County_H1N1;
    avg_susceptible_risk_total_State=avg_susceptible_risk_population_State_H1N1;
end



state_nan=~isnan(avg_susceptible_risk_total_State);
avg_susceptible_risk_total_State=avg_susceptible_risk_total_State(state_nan);
State_Name=State_Name(state_nan);
[avg_susceptible_risk_total_State,Indx_Risk]=sort(avg_susceptible_risk_total_State,'descend');
State_Name=State_Name(Indx_Risk);

figure('units','normalized','outerposition',[0 0.075 1 0.5]);
 ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,0.4,0.6,0.6];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;
       
C_Risk=[hex2rgb('#ffffff');
        hex2rgb('#fff7ec');
        hex2rgb('#fee8c8');
        hex2rgb('#fdd49e');
        hex2rgb('#fdbb84');
        hex2rgb('#fc8d59');
        hex2rgb('#ef6548');
        hex2rgb('#d7301f');
        hex2rgb('#b30000');
        hex2rgb('#7f0000');];

risk_measure=avg_susceptible_risk_total_County;
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
text(ymin+2,max(c_indx)./2,{'Risk of spillover','within dairy farming counties'},'HorizontalAlignment','center','Fontsize',18,'Rotation',270);

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
avg_susceptible_risk_total_State=log(avg_susceptible_risk_total_State);
risk_measure=(avg_susceptible_risk_total_State-min(avg_susceptible_risk_total_State))./(max(avg_susceptible_risk_total_State)-min(avg_susceptible_risk_total_State));

CC_Risk=ones(length(State_Name),3);

for ii=1:length(State_Name)
    if(~isnan(risk_measure(ii)))
        CC_Risk(ii,:)=interp1(x_risk,C_Risk,risk_measure(ii));
    else
        CC_Risk(ii,:)=[0.7 0.7 0.7];
    end
end

subplot('Position',[0.52,0.45,0.47,0.52]);

b=bar([1:length(State_Name)],risk_measure);
b.FaceColor = 'flat';
for ss=1:length(State_Name)
    b.CData(ss,:) =CC_Risk(ss,:);
end
box off;
set(gca,'LineWidth',2,'tickdir','out','YTick',[0 1],'YTickLabel',{'Lowest','Highest'},'Fontsize',16,'XTick',[1:length(State_Name)],'XTickLabel',State_Name);
xlabel('State','FontSize',18);
ylabel({'State spillover risk'},'FontSize',18);
xlim([0.5 length(State_Name)+0.5]);
ylim([0 1.01])

text(-0.175,0.97,'B','FontSize',32,'Units','normalized');

ax1.Position=[-0.445,-0.2,1.3,1.3];


print(gcf,['Supplementary_Figure_Susceptible_' Var_Plot '.png'],'-dpng','-r300');
end

