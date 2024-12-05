function  Figure_Surveillance_Map()
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

load('Surveillance_H5N1_Population_H1N1.mat','County_Surveillance','State_Surveillance','State_Name');


state_nan=~isnan(State_Surveillance);
State_Surveillance=State_Surveillance(state_nan);
State_Name=State_Name(state_nan);
[State_Surveillance,Indx_Risk]=sort(State_Surveillance,'descend');
State_Name=State_Name(Indx_Risk);


State_Name(strcmp(State_Name,'District of Columbia'))='D.C.';
figure('units','normalized','outerposition',[0 0.075 1 0.5]);
 ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,0.4,0.6,0.6];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;
       
C_Risk=[hex2rgb('#ffffff');
        hex2rgb('#ffffe5');
hex2rgb('#fff7bc');
hex2rgb('#fee391');
hex2rgb('#fec44f');
hex2rgb('#fe9929');
hex2rgb('#ec7014');
hex2rgb('#cc4c02');
hex2rgb('#993404');
hex2rgb('#662506');];

risk_measure=County_Surveillance;
t_upper=risk_measure>prctile(risk_measure,95);
v_upper=prctile(risk_measure,95);
t_lower=risk_measure<prctile(risk_measure,5);
v_lower=prctile(risk_measure,5);
risk_measure(t_lower)=v_lower;
risk_measure(t_upper)=v_upper;
risk_measure=(risk_measure)./(min(risk_measure));

subplot('Position',[0.38,0.065,0.01,0.915]);

    x_risk=linspace(1,max(max(risk_measure),1.25),size(C_Risk,1));
    c_indx=linspace(1,max(max(risk_measure),1.25),251);
    y_indx=[1.05:0.05:max(max(risk_measure),1.25)];

xlim([0 1]);
ylim([min(c_indx) max(c_indx)]);    
ymin=1;
dy=2/(1+sqrt(5));
for ii=1:length(c_indx)
    patch([0 0 dy dy],c_indx(ii)+[1 0 0 1],interp1(x_risk,C_Risk,c_indx(ii)),'LineStyle','none');
end

    text(ymin+2.5,1,{'Baseline','lowest'},'Fontsize',16,'HorizontalAlignment','center');

for yy=1:length(y_indx)
        text(ymin,y_indx(yy),[num2str(y_indx(yy))],'Fontsize',16);
end
text(ymin+4.5,mean([min(c_indx) max(c_indx)]),{'Fold-increase in likelihood of','undetected H5N1 cases'},'HorizontalAlignment','center','Fontsize',18,'Rotation',270);

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

risk_measure=(State_Surveillance)./(min(State_Surveillance));

if(ceil(max(risk_measure))>=10)
    x_risk=linspace(0,(max(log10(risk_measure))),size(C_Risk,1));
    CC_Risk=ones(length(State_Name),3);

    for ii=1:length(State_Name)
        if(~isnan(risk_measure(ii)))
            CC_Risk(ii,:)=interp1(x_risk,C_Risk,log10(risk_measure(ii)));
        else
            CC_Risk(ii,:)=[0.7 0.7 0.7];
        end
    end
else
    x_risk=linspace(1,(max(risk_measure)),size(C_Risk,1));
    CC_Risk=ones(length(State_Name),3);

    for ii=1:length(State_Name)
        if(~isnan(risk_measure(ii)))
            CC_Risk(ii,:)=interp1(x_risk,C_Risk,risk_measure(ii));
        else
            CC_Risk(ii,:)=[0.7 0.7 0.7];
        end
    end
end


subplot('Position',[0.52,0.45,0.47,0.52]);

b=bar([1:length(State_Name)],risk_measure);
b.FaceColor = 'flat';
for ss=1:length(State_Name)
    b.CData(ss,:) =CC_Risk(ss,:);
end
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',16,'XTick',[1:length(State_Name)],'XTickLabel',State_Name,'YTick',[0.9:0.1:1.5]);
if(max(risk_measure)>=10)
    set(gca,'yscale','log');
    ylim([0.5 10.^ceil(max(log10(risk_measure)))]);
else                
    ylim([0.95 ceil(max((risk_measure)))]);
end
xlabel('State','FontSize',18);
ylabel({'Fold-increase in','likelihood of undetected','H5N1 cases'},'FontSize',18);
xlim([0.5 length(State_Name)+0.5]);
ylim([0.95 1.5]);

text(-0.145,0.97,'B','FontSize',32,'Units','normalized');

ax1.Position=[-0.445,-0.2,1.3,1.3];


print(gcf,['Figure_Surveillance.png'],'-dpng','-r300');
end

