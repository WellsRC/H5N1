function  Figure_Susceptible_Risk_Map()
% Plot_Variable
close all;
states = shaperead('usastatelo', 'UseGeoCoords', true);
S=shaperead([pwd '/Shapefile/cb_2021_us_county_500k.shp'],'UseGeoCoords',true);
load([pwd '/Data/Data_US_County.mat'],'US_County');
t_f=~isnan(US_County.GINI_2022) & ~strcmp(US_County.STUSPS,'PR');

US_County=readgeotable([pwd '/Shapefile/cb_2021_us_county_500k.shp']);

US_County=US_County(t_f,[2 3 6 7 10 9 12 13]);
Var_Name=US_County.Properties.VariableNames;
Var_Name{end-1}='AREA_LAND';
Var_Name{end}='AREA_WATER';
US_County.Properties.VariableNames=Var_Name;

[~,Indx]=sortrows(US_County,[1 2]);

S=S(Indx);

NS=length(S);


State_Name=unique({S.STUSPS});
State_Risk=NaN.*zeros(size(State_Name));

load([pwd '/Data/Data_US_County.mat'],'US_County');
t_f=~isnan(US_County.GINI_2022) & ~strcmp(US_County.STUSPS,'PR');
US_County=US_County(t_f,:);

load('Average_Risk_Population_COVID.mat','avg_susceptible_risk_population_County_COVID');
load('Average_Risk_Population_H1N1.mat','avg_susceptible_risk_population_County_H1N1');
avg_susceptible_risk_population_County=1-sqrt((1-avg_susceptible_risk_population_County_COVID).*(1-avg_susceptible_risk_population_County_H1N1));

for ss=1:length(State_Name)
    t_state=strcmp(State_Name{ss},US_County.STUSPS);
    c_r=avg_susceptible_risk_population_County(t_state);
    w_c=US_County.POPULATION_SIZE_2022(t_state);
    t_inc=w_c>0 & ~isnan(c_r);
    c_r=c_r(t_inc);
    w_c=w_c(t_inc);
    test=sum((w_c./sum(w_c)).*log(1-c_r));
    if(~isempty(c_r))
        State_Risk(ss)=1-exp(test);
    end
end

state_nan=~isnan(State_Risk);
State_Risk=State_Risk(state_nan);
State_Name=State_Name(state_nan);
[State_Risk,Indx_Risk]=sort(State_Risk,'ascend');
State_Name=State_Name(Indx_Risk);

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
            ax4=subplot('Position',[0.55,0.08,0.44,0.4]);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Risk: colour bar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    switch vv
        case 1
            subplot('Position',[0.41,0.525,0.01,0.45]);
            Title_Name={'Susceptible risk','(COVID-based)'};
            
            risk_measure=avg_susceptible_risk_population_County_COVID;
        case 2
            subplot('Position',[0.885,0.525,0.01,0.45]);
            Title_Name={'Susceptible risk','(H1N1-based)'};
           

            risk_measure=avg_susceptible_risk_population_County_H1N1;
        case 3
            subplot('Position',[0.41,0.025,0.01,0.45]);    
            Title_Name={'Susceptible risk','(COVID- and H1N1 -based)'};


            risk_measure=avg_susceptible_risk_population_County;
        case 4

            risk_measure=State_Risk;
    end
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
    
    if(vv<4)
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
    else
        CC_Risk=ones(length(State_Name),3);
        
        for ii=1:length(State_Name)
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
            b=bar([1:length(State_Name)],risk_measure);
            b.FaceColor = 'flat';
            for ss=1:length(State_Name)
                b.CData(ss,:) =CC_Risk(ss,:);
            end
            box off;
            set(gca,'LineWidth',2,'tickdir','out','YTick',[0 1],'YTickLabel',{'Lowest','Highest'},'Fontsize',16,'XTick',[1:length(State_Name)],'XTickLabel',State_Name);
            xlabel('State','FontSize',18);
            ylabel({'State overall','susceptible risk'},'FontSize',18);
            xlim([0.5 length(State_Name)+0.5]);
            ylim([0 1.01])
            text(-0.15,1,'D','FontSize',32,'Units','normalized');
    end
end


ax1.Position=[-0.075,0.425,0.6,0.6];
ax2.Position=[0.4,0.425,0.6,0.6];
ax3.Position=[-0.075,-0.075,0.6,0.6];

print(gcf,['Preliminary_Susceptible_Risk_Plot.png'],'-dpng','-r300');
end

