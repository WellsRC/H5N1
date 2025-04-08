clear;
clc;

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


[~,~,~,~,Affected_County_Farms,~,~,~,~,~,~,~]= Dairy_Covariates({},{},{});

figure('units','normalized','outerposition',[0.25 0.25 0.4 0.5]);
 ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,0.4,0.6,0.6];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;

subplot('Position',[0.86,0.025,0.01,0.95]);
Title_Name={['Number of outbreaks'],['among dairy farms']};

risk_measure=Affected_County_Farms;
C_Risk=[hex2rgb('#ffffff');
        hex2rgb('#FBF5BB');
        hex2rgb('#FFFF00');
        hex2rgb('#C29545');
        hex2rgb('#FAC898'); % OR
        hex2rgb('#FFA500');
        hex2rgb('#FF7518');
        hex2rgb('#CC5500');
        hex2rgb('#FAA0A0');%R
        hex2rgb('#E30B5C');
        hex2rgb('#D2042D');
        hex2rgb('#880808');
        hex2rgb('#000000')];
 y_indx=unique(risk_measure);
x_risk=y_indx;

 x_indx=y_indx;
 c_indx=linspace(0,1,length(y_indx));
dx=c_indx(2)-c_indx(1);
xlim([0 1]);
ylim([0 1+dx])
ymin=1;
dy=2/(1+sqrt(5));
for ii=1:length(c_indx)
    patch([0 0 dy dy],c_indx(ii)+[dx.*0.8 dx.*0.2 dx.*0.2 dx.*0.8],C_Risk(ii,:),'LineWidth',1.5);
end

for yy=1:length(y_indx)
    text(ymin+1,dx./2+c_indx(yy),[num2str(y_indx(yy),'%2.0f')],'Fontsize',16);            
end

text(ymin+7.5, (1+dx)/2,Title_Name,'HorizontalAlignment','center','Fontsize',18,'Rotation',270);

axis off;  

NS=length(S);
CC_Risk=interp1(x_risk,C_Risk,risk_measure);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;


ax1.Position=[-0.16,-0.15,1.2,1.2];
print(gcf,['Supplementary_Figure_2.png'],'-dpng','-r300');