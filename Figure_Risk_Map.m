function  Figure_Risk_Map(Var_Plot)
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

if(strcmp(Var_Plot,'Dairy'))
    load('Dairy_Risk_AIC.mat','mle_outbreak_dairy_farm_County','mle_outbreak_risk_dairy_farm_County','mle_spillover_risk_dairy_farm_County','mle_outbreak_dairy_farm_State','mle_spillover_risk_dairy_farm_State');
    mle_outbreak_farm_County=mle_outbreak_dairy_farm_County;
    mle_outbreak_risk_farm_County=mle_outbreak_risk_dairy_farm_County;
    mle_spillover_risk_farm_County=mle_spillover_risk_dairy_farm_County;

    mle_outbreak_farm_State=mle_outbreak_dairy_farm_State;
    mle_spillover_risk_farm_State=mle_spillover_risk_dairy_farm_State;

    farm_type='dairy';
elseif(strcmp(Var_Plot,'Poultry')) 
    load('Poultry_Risk_AIC.mat','mle_outbreak_poultry_farm_County','mle_outbreak_risk_poultry_farm_County','mle_spillover_risk_poultry_farm_County','mle_outbreak_poultry_farm_State','mle_spillover_risk_poultry_farm_State');
    mle_outbreak_farm_County=mle_outbreak_poultry_farm_County;
    mle_outbreak_risk_farm_County=mle_outbreak_risk_poultry_farm_County;
    mle_spillover_risk_farm_County=mle_spillover_risk_poultry_farm_County;

    mle_outbreak_farm_State=mle_outbreak_poultry_farm_State;
    mle_spillover_risk_farm_State=mle_spillover_risk_poultry_farm_State;
    farm_type='poultry';
end

state_nan=~isnan(mle_outbreak_farm_State);
mle_outbreak_farm_State=mle_outbreak_farm_State(state_nan);
mle_spillover_risk_farm_State=mle_spillover_risk_farm_State(state_nan);
State_Name=State_Name(state_nan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Construct plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


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
            ax4=subplot('Position',[0.55,0.20,0.4,0.28]);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Risk: colour bar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch vv
        case 1
            subplot('Position',[0.41,0.525,0.01,0.45]);
            Title_Name={['Expected number of outbreaks'],['among ' farm_type ' farms']};
            
            risk_measure=mle_outbreak_farm_County;
            C_Risk=[hex2rgb('#ffffff');
                    hex2rgb('#f9BA32');
                    hex2rgb('#f16913');
                    hex2rgb('#a50f15');
                    hex2rgb('#000000')];
        case 2
            subplot('Position',[0.885,0.525,0.01,0.45]);
            Title_Name={['Outbreak risk for ' farm_type ' farms']};

            risk_measure=mle_outbreak_risk_farm_County;

            C_Risk=[1.0000    1.0000    1.0000
                1.0000    0.7812    0.4975
                0.9375    0.5859    0.3731
                0.6250    0.3906    0.2487
                0.3125    0.1953    0.1244
                    0         0         0];
        case 3
            subplot('Position',[0.41,0.025,0.01,0.45]);    
            Title_Name={['Spillover risk to humans from ' farm_type ' farms']};


            risk_measure=mle_spillover_risk_farm_County;


            C_Risk=[1.0000    1.0000    1.0000
                    0.9309    0.9309    0.8367
                    0.8563    0.8563    0.6325
                    0.7746    0.6583    0.5164
                    0.6831    0.3651    0.3651
                    0.4082         0         0];
        case 4
            state_filter=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
            State_Name=State_Name(~state_filter);
            risk_measure_outbreak=mle_outbreak_farm_State(~state_filter);
            C_out=[0 0 0];
            risk_measure_spillover=mle_spillover_risk_farm_State(~state_filter);
            C_spill=[0.7746    0.6583    0.5164];
            [~,indx_r]=sort(mean([risk_measure_outbreak risk_measure_spillover],2),'descend');
            risk_measure_outbreak=risk_measure_outbreak(indx_r);
            risk_measure_spillover=risk_measure_spillover(indx_r);
            State_Name=State_Name(indx_r);
    end
    
    if(vv<4)
        text_ub=false;
        
        if(max(risk_measure)-prctile(risk_measure,99)>0.05)
            risk_measure(risk_measure>prctile(risk_measure,99))=prctile(risk_measure,99);
            text_ub=true;
        end

        x_risk=linspace(0,ceil(100.*max(risk_measure))./100,size(C_Risk,1));
        c_indx=linspace(0,ceil(100.*max(risk_measure))./100,251);
        if(length([x_risk(1):0.01:ceil(100.*x_risk(end))./100])<=11)
            y_indx=[x_risk(1):0.01:ceil(100.*x_risk(end))./100];
        elseif(length([x_risk(1):0.05:ceil(20.*x_risk(end))./20])<=11)
            y_indx=[x_risk(1):0.05:ceil(20.*x_risk(end))./20];
        elseif(length([x_risk(1):0.1:ceil(10.*x_risk(end))./10])<=11)
            y_indx=[x_risk(1):0.1:ceil(10.*x_risk(end))./10];
        else
            y_indx=[x_risk(1):0.25:x_risk(end)];
        end


        dx=c_indx(2)-c_indx(1);
        xlim([0 1]);
        ylim([x_risk(1) max(x_risk)+dx]);    
        ymin=1;
        dy=2/(1+sqrt(5));
        for ii=1:length(c_indx)
            patch([0 0 dy dy],c_indx(ii)+[dx 0 0 dx],interp1(x_risk,C_Risk,c_indx(ii)),'LineStyle','none');
        end
        
        for yy=1:length(y_indx)-1
            text(ymin,y_indx(yy),[num2str(y_indx(yy),'%3.2f')],'Fontsize',16);            
        end
        if(text_ub)
            text(ymin,x_risk(end),[num2str(x_risk(end),'%3.2f') '+'],'Fontsize',16);
        else
            text(ymin,x_risk(end),[num2str(x_risk(end),'%3.2f')],'Fontsize',16);
        end
        text(ymin+4.5,[min(c_indx)+max(c_indx)]./2,Title_Name,'HorizontalAlignment','center','Fontsize',18,'Rotation',270);
        
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
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot uptake
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    switch vv
        case 1
            CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
            geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
            geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;
        case 2
            CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
            geoshow(ax2,S,'SymbolSpec',CM,'LineStyle','None'); 
            geoshow(ax2, states,'Facecolor','none','LineWidth',1.5); hold on;
        case 3
            CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
            geoshow(ax3,S,'SymbolSpec',CM,'LineStyle','None'); 
            geoshow(ax3, states,'Facecolor','none','LineWidth',1.5); hold on;
        case 4
            b=bar([1:length(State_Name)],[risk_measure_outbreak],'LineStyle','none');
            b(1).FaceColor = 'flat';
            for ss=1:length(State_Name)
                b(1).CData(ss,:) =C_out;
            end
            box off;
            if(max(risk_measure_outbreak)>300)
                set(gca,'LineWidth',2,'tickdir','out','Fontsize',14,'XTick',[1:length(State_Name)],'XTickLabel',State_Name,'YTick',[0:100:1500]);
            else
                set(gca,'LineWidth',2,'tickdir','out','Fontsize',14,'XTick',[1:length(State_Name)],'XTickLabel',State_Name,'YTick',[0:25:300]);
            end
           
            xlabel('State','FontSize',18);
            yl=ylabel({'State-level outbreaks'},'FontSize',18);
            xlim([0.5 length(State_Name)+0.5]);
            ymax=ceil(max(risk_measure_outbreak(:))./100).*100;
            ylim([0 ymax])

            yyaxis right
            plot([1:length(State_Name)],risk_measure_spillover,'o','color',C_spill,'MarkerFaceColor',C_spill,'LineWidth',2);
            set(gca,'LineWidth',2,'tickdir','out','Fontsize',14,'XTick',[1:length(State_Name)],'YColor',C_spill,'YTick',[0:0.1:1]);
            ylabel({'State-level spillover risk'},'FontSize',18,'Rotation',270);
            text(-0.15,1,'D','FontSize',32,'Units','normalized');
            box off;
            yl.FontSize=18;
    end
end



ax1.Position=[-0.075,0.425,0.6,0.6];
ax2.Position=[0.4,0.425,0.6,0.6];
ax3.Position=[-0.075,-0.075,0.6,0.6];

print(gcf,['Figure_' Var_Plot '_Risk_Plot.png'],'-dpng','-r300');
end

