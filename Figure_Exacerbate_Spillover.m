function Figure_Exacerbate_Spillover()

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


State_Name=unique({S.STUSPS});
State_Risk=NaN.*zeros(size(State_Name));

load([pwd '/Data/Spillover_Exacerbation.mat'],'US_County_Spillover');
SH=US_County_Spillover.SLAUGHTERHOUSE;
SH=1-exp(-0.16.*SH);
load('Average_Risk_Swine.mat','avg_overall_risk_swine_farm_County','avg_exposure_risk_swine_farm_County','avg_susceptible_risk_swine_farm_County');
avg_overall_risk_farm_County=avg_overall_risk_swine_farm_County;
avg_exposure_risk_farm_County=avg_exposure_risk_swine_farm_County;
avg_susceptible_risk_farm_County=avg_susceptible_risk_swine_farm_County;

for vv=1:4
    switch vv
        case 1
            figure('units','normalized','outerposition',[0 0.075 1 1]);
             ax1=usamap('conus');
        
            framem off; gridm off; mlabel off; plabel off;
            ax1.Position=[-0.3,0.4,0.6,0.6];
            
            states = shaperead('usastatelo', 'UseGeoCoords', true);
            geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;
        case 2
            ax2=usamap('conus');
            
            framem off; gridm off; mlabel off; plabel off;
            ax2.Position=[1.7,0.4,0.6,0.6];
            
            states = shaperead('usastatelo', 'UseGeoCoords', true);
            geoshow(ax2, states,'Facecolor','none','LineWidth',0.5); hold on;
        case 3
            ax3=usamap('conus');
            
            framem off; gridm off; mlabel off; plabel off;
            ax3.Position=[-0.3,-0.1,0.6,0.6];
            
            states = shaperead('usastatelo', 'UseGeoCoords', true);
            geoshow(ax3, states,'Facecolor','none','LineWidth',0.5); hold on;
        case 4
            ax4=usamap('conus');
            
            framem off; gridm off; mlabel off; plabel off;
            ax4.Position=[1.7,-0.1,0.6,0.6];
            
            states = shaperead('usastatelo', 'UseGeoCoords', true);
            geoshow(ax4, states,'Facecolor','none','LineWidth',0.5); hold on;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Risk: colour bar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch vv
        case 1
            subplot('Position',[0.41,0.525,0.01,0.45]);
            Title_Name={'County farms risk of','external','exposure to H5N1'};

            C_Risk=[hex2rgb('#ffffff');
                    hex2rgb('##ffffcc');
                    hex2rgb('#ffeda0');
                    hex2rgb('#fed976');
                    hex2rgb('#feb24c');
                    hex2rgb('#fd8d3c');
                    hex2rgb('#fc4e2a');
                    hex2rgb('#e31a1c');
                    hex2rgb('#bd0026');
                    hex2rgb('#800026');];
            risk_measure=avg_exposure_risk_farm_County;
        case 2
            subplot('Position',[0.885,0.525,0.01,0.45]);
            Title_Name={'Susecptibility of county','farms to H5N1'};

            C_Risk=[hex2rgb('#ffffff');
                    hex2rgb('#fff7fb');
                    hex2rgb('#ece7f2');
                    hex2rgb('#d0d1e6');
                    hex2rgb('#a6bddb');
                    hex2rgb('#74a9cf');
                    hex2rgb('#3690c0');
                    hex2rgb('#0570b0');
                    hex2rgb('#045a8d');
                    hex2rgb('#023858');];

            risk_measure=avg_susceptible_risk_farm_County;
        case 3
            subplot('Position',[0.41,0.025,0.01,0.45]);    
            Title_Name={'County farms overall','risk to H5N1'};

            C_Risk=[hex2rgb('#ffffff');
                    hex2rgb('##fff7f3');
                    hex2rgb('#fde0dd');
                    hex2rgb('#fcc5c0');
                    hex2rgb('#fa9fb5');
                    hex2rgb('#f768a1');
                    hex2rgb('#dd3497');
                    hex2rgb('#ae017e');
                    hex2rgb('#7a0177');
                    hex2rgb('#49006a');];

            risk_measure=avg_overall_risk_farm_County;
        case 4
            subplot('Position',[0.885,0.025,0.01,0.45]);  
            C_Risk=[hex2rgb('#ffffff');
                    hex2rgb('##fff7f3');
                    hex2rgb('#fde0dd');
                    hex2rgb('#fcc5c0');
                    hex2rgb('#fa9fb5');
                    hex2rgb('#f768a1');
                    hex2rgb('#dd3497');
                    hex2rgb('#ae017e');
                    hex2rgb('#7a0177');
                    hex2rgb('#49006a');];

            risk_measure=SH;
    end
    if(vv<4)
        t_upper=risk_measure>prctile(risk_measure,95);
        v_upper=prctile(risk_measure,95);
        t_lower=risk_measure<prctile(risk_measure,5);
        v_lower=prctile(risk_measure,5);
        risk_measure(t_lower)=v_lower;
        risk_measure(t_upper)=v_upper;
        risk_measure=log(risk_measure);
    end
        risk_measure=(risk_measure-min(risk_measure))./(max(risk_measure)-min(risk_measure));
        x_risk=linspace(0,1,size(C_Risk,1));
        c_indx=[1:size(C_Risk,1)];
    
        xlim([0 1]);
        ylim([0 max(c_indx)]);    
        ymin=1;
        dy=2/(1+sqrt(5));
        for ii=1:length(c_indx)
            patch([0 0 dy dy],c_indx(ii)-[1 0 0 1],C_Risk(ii,:),'LineStyle','none');
        end
        
        text(ymin,0,'Lowest','Fontsize',16);
        text(ymin,max(c_indx),'Highest','Fontsize',16);
        text(ymin+2,max(c_indx)./2,Title_Name,'HorizontalAlignment','center','Fontsize',18,'Rotation',270);
        
        axis off;  
    
        text(-40,1,char(64+vv),'FontSize',32,'Units','normalized');
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
    switch vv
        case 1
            geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
            geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;
        case 2
            geoshow(ax2,S,'SymbolSpec',CM,'LineStyle','None'); 
            geoshow(ax2, states,'Facecolor','none','LineWidth',1.5); hold on;
        case 3
            geoshow(ax3,S,'SymbolSpec',CM,'LineStyle','None'); 
            geoshow(ax3, states,'Facecolor','none','LineWidth',1.5); hold on;
        case 4
            geoshow(ax4,S,'SymbolSpec',CM,'LineStyle','None'); 
            geoshow(ax4, states,'Facecolor','none','LineWidth',1.5); hold on;
    end
end


ax1.Position=[-0.075,0.425,0.6,0.6];
ax2.Position=[0.4,0.425,0.6,0.6];
ax3.Position=[-0.075,-0.075,0.6,0.6];
ax4.Position=[0.4,-0.075,0.6,0.6];

print(gcf,['Preliminary_Exacerbate_Spillover_Plot.png'],'-dpng','-r300');
end