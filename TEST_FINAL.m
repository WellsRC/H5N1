clear;
clc;


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
 load('Total_Spillover_Risk_County_H1N1.mat','avg_localized_transmission_risk_total_County');

 avg_localized_transmission_risk_total_County=(avg_localized_transmission_risk_total_County)./(max(avg_localized_transmission_risk_total_County));
 load('Average_Risk_Swine.mat','avg_susceptible_risk_swine_farm_County');
 avg_susceptible_risk_swine_farm_County=(avg_susceptible_risk_swine_farm_County)./(max(avg_susceptible_risk_swine_farm_County));
load('Underreporting_Dairy_Farms.mat',"Prob_UR",'ur_level');
Underreporting_Rate=0.6;
Underreporting_County=Prob_UR(:,abs(ur_level-Underreporting_Rate)==min(abs(ur_level-Underreporting_Rate)));
Underreporting_County=(Underreporting_County)./(max(Underreporting_County));
load('Surveillance_H5N1_Population.mat','County_Surveillance');
County_Surveillance=(County_Surveillance)./(max(County_Surveillance));

avg_localized_transmission_risk_total_County(isnan(avg_localized_transmission_risk_total_County))=0;
avg_susceptible_risk_swine_farm_County(isnan(avg_susceptible_risk_swine_farm_County))=0;
Underreporting_County(isnan(Underreporting_County))=0;
County_Surveillance(isnan(County_Surveillance))=0;
y=(avg_localized_transmission_risk_total_County+Underreporting_County+County_Surveillance); % +avg_susceptible_risk_swine_farm_County

risk_measure=(y-min(y))./(max(y)-min(y));

cv=prctile(risk_measure,90);
risk_measure(risk_measure>=cv)=1;
risk_measure(risk_measure<cv)=0;

C_Risk=[hex2rgb('#ffffff');
hex2rgb('#000000')];

CC_Risk=0.7.*ones(length(S),3);
CC_Risk(risk_measure==1,:)=repmat(C_Risk(2,:),sum(risk_measure==1),1);
CC_Risk(risk_measure==0,:)=repmat(C_Risk(1,:),sum(risk_measure==0),1);
% CC_Risk=ones(length(S),3);
% for ii=1:length(S)
%     if(~isnan(risk_measure(ii)))
%         CC_Risk(ii,:)=interp1(x_risk,C_Risk,risk_measure(ii));
%     else
%         CC_Risk(ii,:)=[0.7 0.7 0.7];
%     end
% end

figure('units','normalized','outerposition',[0 0.075 1 1]);
 ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;

CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});


            geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
            geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;