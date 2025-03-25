clc;
clear;
close all;

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

County_remove=strcmp("AK",US_County.STUSPS) |strcmp("HI",US_County.STUSPS) | strcmp("AS",US_County.STUSPS) | strcmp("GU",US_County.STUSPS) | strcmp("MP",US_County.STUSPS) | strcmp("PR",US_County.STUSPS) | strcmp("VI",US_County.STUSPS);
S=S(~County_remove,:);

NS=length(S);

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);

[~,~,County_Farms,Affected_County_Farms,~,~,~] = Poultry_Covariates({},{},{});
load('Average_Risk_Poultry.mat','w_AIC','outbreak_poultry_farm_County','State_Name');
avg_outbreak_farm_County=outbreak_poultry_farm_County*w_AIC;

avg_outbreak_farm_County(County_Farms==0)=NaN;

Err=(avg_outbreak_farm_County-Affected_County_Farms);

figure('units','normalized','outerposition',[0.25 0.25 0.5 0.5]);
ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;
% ax1.Position=[-0.3,0.4,0.6,0.6];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;


C_Err=[hex2rgb('#b2182b');
hex2rgb('#d6604d');
hex2rgb('#f4a582');
hex2rgb('#fddbc7');
hex2rgb('#f7f7f7');
hex2rgb('#d1e5f0');
hex2rgb('#92c5de');
hex2rgb('#4393c3');
hex2rgb('#2166ac');];

x_Err=[-45 -30 -15 -0.5 0 0.5 15 30 45];


CC_Err=ones(length(S),3);
        
for ii=1:length(S)
    if(~isnan(Err(ii)))
        CC_Err(ii,:)=interp1(x_Err,C_Err,Err(ii));
    else
        CC_Err(ii,:)=[0.7 0.7 0.7];
    end
end

CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Err});
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;