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

County_remove=strcmp("AK",US_County.STUSPS) | strcmp("HI",US_County.STUSPS) | strcmp("AS",US_County.STUSPS) | strcmp("GU",US_County.STUSPS) | strcmp("MP",US_County.STUSPS) | strcmp("PR",US_County.STUSPS) | strcmp("VI",US_County.STUSPS);
S=S(~County_remove,:);

NS=length(S);

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);


load('Average_Risk_Poultry.mat','onward_transmission_poultry_farm_County','no_farms','w_AIC');

avg_onward_transmission_County=onward_transmission_poultry_farm_County(:,w_AIC==max(w_AIC));

figure('units','normalized','outerposition',[0.25 0.25 0.4 0.5]);
 ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,0.4,0.6,0.6];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;

subplot('Position',[0.85,0.035,0.03,0.94]);

Title_Name={'Likelihood of onward transmission','after spillover from poultry farm'};

risk_measure=log10(avg_onward_transmission_County);
C_Risk=[hex2rgb('#fff7f3');
hex2rgb('#fde0dd');
hex2rgb('#fcc5c0');
hex2rgb('#fa9fb5');
hex2rgb('#f768a1');
hex2rgb('#dd3497');
hex2rgb('#ae017e');
hex2rgb('#7a0177');
hex2rgb('#49006a');];

y_indx=linspace(-5,0,size(C_Risk,1));
x_risk=y_indx;

 x_indx=y_indx;
 c_indx=linspace(y_indx(1),y_indx(end),1001);
dx=c_indx(2)-c_indx(1);
xlim([0 1]);
ylim([y_indx(1) y_indx(end)+dx])
ymin=0.75;
dy=2/(1+sqrt(5));
for ii=1:length(c_indx)
    patch([0 0 dy dy],c_indx(ii)+[dx -dx -dx dx],interp1(x_risk,C_Risk,c_indx(ii)),'LineStyle','none');
end

patch([0 0 dy dy], [y_indx(1) y_indx(end)+dx y_indx(end)+dx y_indx(1)],'k','FaceAlpha',0,'LineWidth',2);

yl_indx=[-5:0];
for yy=1:length(yl_indx)
    text(ymin,yl_indx(yy),['10^{' num2str(yl_indx(yy)) '}'],'Fontsize',16);            
end

text(ymin+3,-2.5,Title_Name,'HorizontalAlignment','center','Fontsize',18,'Rotation',270);

axis off;  

NS=length(S);
CC_Risk=interp1(x_risk,C_Risk,risk_measure);

CC_Risk(no_farms,:)=repmat([0.5 0.5 0.5],sum(no_farms),1);

CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;


ax1.Position=[-0.16,-0.15,1.2,1.2];

print(gcf,['Supplementary_Figure_Poultry_Onward.png'],'-dpng','-r300');