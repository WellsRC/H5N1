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


[~,~,~,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,~,~] = Poultry_Covariates({},{});


fprintf([' number of outbreaks :' num2str(sum(Affected_County_Farms),'%3.1f') '\n']);
fprintf([' number of counties affected :' num2str(sum(Affected_County_Farms>0),'%3.1f') '\n']);
fprintf(['Maximum number of outbreaks within a county:' num2str(max(Affected_County_Farms),'%3.1f') '\n']);

p_no=100.*mean(Affected_County_Farms(County_Farms>0)==0);

fprintf(['Percent of counties with no outbreak:' num2str(p_no,'%3.1f') '%% \n']);

State_O=state_weight_matrix*Affected_County_Farms;

fprintf(['Maximum number of outbreaks within a state:' num2str(max(State_O),'%3.1f') '\n']);
fprintf(['Minumum number of outbreaks within a state:' num2str(min(State_O),'%3.1f') '\n']);
fprintf(['Median number of outbreaks within a state:' num2str(median(State_O),'%3.1f') '\n']);


for vv=1:2
     switch vv
        case 1
            figure('units','normalized','outerposition',[0 0.075 1 1]);
             ax1=usamap('conus');
        
            framem off; gridm off; mlabel off; plabel off;
            ax1.Position=[-0.3,0.4,0.6,0.6];
            
            states = shaperead('usastatelo', 'UseGeoCoords', true);
            geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;

            subplot('Position',[0.41,0.525,0.01,0.45]);
            Title_Name={['Number of outbreaks'],['among poultry farms']};
            
            risk_measure=Affected_County_Farms;
            C_Risk=[hex2rgb('#ffffff');
                    hex2rgb('#f9BA32');
                    hex2rgb('#f16913');
                    hex2rgb('#a50f15');
                    hex2rgb('#000000')];
            y_indx=[0:5:max(risk_measure)];

            x_risk=linspace(0,max(risk_measure),size(C_Risk,1));
            c_indx=linspace(0,max(risk_measure),251);
            dx=c_indx(2)-c_indx(1);
            xlim([0 1]);
            ylim([0 max(x_risk)+dx]);    
            ymin=1;
            dy=2/(1+sqrt(5));
            for ii=1:length(c_indx)
                patch([0 0 dy dy],c_indx(ii)+[dx 0 0 dx],interp1(x_risk,C_Risk,c_indx(ii)),'LineStyle','none');
            end
            
            for yy=1:length(y_indx)
                text(ymin,y_indx(yy),[num2str(y_indx(yy),'%2.0f')],'Fontsize',16);            
            end
            
            text(ymin+4.5,[min(c_indx)+max(c_indx)]./2,Title_Name,'HorizontalAlignment','center','Fontsize',18,'Rotation',270);
            
            axis off;  
            patch([0 0 dy dy],[0 max(x_risk)+dx max(x_risk)+dx 0],'k','FaceAlpha',0,'LineWidth',1.5);
        
            text(-40,1,char(64+vv),'FontSize',32,'Units','normalized');
         case 2
            ax2=usamap('conus');
            
            framem off; gridm off; mlabel off; plabel off;
            ax2.Position=[1.7,0.4,0.6,0.6];
            
            states = shaperead('usastatelo', 'UseGeoCoords', true);
            geoshow(ax2, states,'Facecolor','none','LineWidth',0.5); hold on;

            subplot('Position',[0.885,0.525,0.01,0.45]);
            risk_measure=State_Spillover_Events;
            Title_Name={['Human infections linked'],['to poultry farms']};
            y_indx=unique(risk_measure);
            C_Risk=[1.0000    1.0000    1.0000
                    % 0.9309    0.9309    0.8367
                    0.8563    0.8563    0.6325
                    % 0.7746    0.6583    0.5164
                    0.6831    0.3651    0.3651
                    0.4082         0         0];

            x_risk=y_indx;
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
                text(ymin,dx./2+c_indx(yy),[num2str(y_indx(yy),'%2.0f')],'Fontsize',16);            
            end
            
            text(ymin+4.5, (1+dx)/2,Title_Name,'HorizontalAlignment','center','Fontsize',18,'Rotation',270);
            
            axis off;  
        
            text(-40,1,char(64+vv),'FontSize',32,'Units','normalized');
     end
    
     
    switch vv
        case 1
            NS=length(S);
            CC_Risk=interp1(x_risk,C_Risk,risk_measure);
            CC_Risk(County_Farms==0,:)=repmat([0.5 0.5 0.5],sum(County_Farms==0),1);
            CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
            geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
            geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;
        case 2
            
            NS=length(states);
            CC_Risk=ones(NS,3);
            for kk=1:NS
                tf=strcmp(states(kk).Name,State_Name);
                if(sum(tf)>0)
                    CC_Risk(kk,:)=interp1(x_risk,C_Risk,risk_measure(tf));
                end
            end
            CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
            geoshow(ax2,states,'SymbolSpec',CM,'LineWidth',1.5); 
    end
end


[~,~,~,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,~,~,~,~,~]= Dairy_Covariates({},{},{});

fprintf(['number of dairy outbreaks: ' num2str(sum(Affected_State_Farms)) '\n']);
fprintf(['Number of states affected: ' num2str(sum(Affected_State_Farms>0)) '\n']);
fprintf(['Maximum number of dairy outbreaks within a state:' num2str(max(Affected_State_Farms(Affected_State_Farms>0)),'%3.1f') '\n']);
fprintf(['Minumum number of dairy outbreaks within a state:' num2str(min(Affected_State_Farms(Affected_State_Farms>0)),'%3.1f') '\n']);
fprintf(['Median number of  dairy outbreaks within a state:' num2str(median(Affected_State_Farms(Affected_State_Farms>0)),'%3.1f') '\n']);

for vv=1:2
     switch vv
        case 1
            ax3=usamap('conus');
            
            framem off; gridm off; mlabel off; plabel off;
            ax3.Position=[-0.3,-0.1,0.6,0.6];
            
            states = shaperead('usastatelo', 'UseGeoCoords', true);
            geoshow(ax3, states,'Facecolor','none','LineWidth',0.5); hold on;

            subplot('Position',[0.41,0.025,0.01,0.45]);   
            Title_Name={['Number of outbreaks'],['among dairy farms']};
            
            risk_measure=Affected_State_Farms;
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
                text(ymin,dx./2+c_indx(yy),[num2str(y_indx(yy),'%2.0f')],'Fontsize',16);            
            end
            
            text(ymin+4.5, (1+dx)/2,Title_Name,'HorizontalAlignment','center','Fontsize',18,'Rotation',270);
            
            axis off;  
         case 2
            ax4=usamap('conus');
            
            framem off; gridm off; mlabel off; plabel off;
            ax4.Position=[1.7,-0.1,0.6,0.6];
            
            states = shaperead('usastatelo', 'UseGeoCoords', true);
            geoshow(ax4, states,'Facecolor','none','LineWidth',0.5); hold on;

            subplot('Position',[0.885,0.025,0.01,0.45]);
            risk_measure=State_Spillover_Events;
            Title_Name={['Human infections linked'],['to dairy farms']};
            
            C_Risk=[1.0000    1.0000    1.0000
                    0.9309    0.9309    0.8367
                    % 0.8563    0.8563    0.6325
                    0.7746    0.6583    0.5164
                    % 0.6831    0.3651    0.3651
                    0.4082         0         0];
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
                text(ymin,dx./2+c_indx(yy),[num2str(y_indx(yy),'%2.0f')],'Fontsize',16);            
            end
            
            text(ymin+4.5,(1+dx)/2,Title_Name,'HorizontalAlignment','center','Fontsize',18,'Rotation',270);
            
            axis off; 
     end
   

    text(-40,1,char(66+vv),'FontSize',32,'Units','normalized');


    NS=length(states);
    CC_Risk=ones(NS,3);
    for kk=1:NS
        tf=strcmp(states(kk).Name,State_Name);
        if(sum(tf)>0)
            CC_Risk(kk,:)=interp1(x_risk,C_Risk,risk_measure(tf));
        end
    end
    CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});

     switch vv
        case 1
            
            geoshow(ax3,states,'SymbolSpec',CM,'LineWidth',1.5); 
        case 2
            
            geoshow(ax4,states,'SymbolSpec',CM,'LineWidth',1.5); 
    end
end
ax1.Position=[-0.075,0.425,0.6,0.6];
ax2.Position=[0.4,0.425,0.6,0.6];
ax3.Position=[-0.075,-0.075,0.6,0.6];
ax4.Position=[0.4,-0.075,0.6,0.6];

print(gcf,['Supplementary_Figure_1.png'],'-dpng','-r300');