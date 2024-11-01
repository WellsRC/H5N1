function  Figure_Underreporting_Map()

Underreporting_Rate=0.6;
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


load('Underreporting_Dairy_Farms.mat',"Prob_UR",'ur_level');

State_Name=unique({S.STUSPS});
State_Risk=NaN.*zeros(size(State_Name));

load([pwd '/Data/Data_US_County.mat'],'US_County');
Underreporting_County=Prob_UR(:,abs(ur_level-Underreporting_Rate)==min(abs(ur_level-Underreporting_Rate)));

figure('units','normalized','outerposition',[0 0.075 1 0.5]);
 ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,0.4,0.6,0.6];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;
       
C_Risk=[hex2rgb('#ffffff');
        hex2rgb('#fff5eb');
hex2rgb('#fee6ce');
hex2rgb('#fdd0a2');
hex2rgb('#fdae6b');
hex2rgb('#fd8d3c');
hex2rgb('#f16913');
hex2rgb('#d94801');
hex2rgb('#a63603');
hex2rgb('#7f2704');];

risk_measure=Underreporting_County;
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
text(ymin+2,max(c_indx)./2,{'Probability of underreporting','among counties'},'HorizontalAlignment','center','Fontsize',18,'Rotation',270);

axis off;  

text(-38,0.98,'A','FontSize',32,'Units','normalized');
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
PR_UR=zeros(length(State_Name),length(ur_level));
for uu=1:length(ur_level)
    Underreporting_County=Prob_UR(:,uu);
    for ss=1:length(State_Name)
        t_state=strcmp(State_Name{ss},US_County.STUSPS);
        c_r=Underreporting_County(t_state);
        w_c=US_County.TOTAL_DAIRY_OPERATIONS(t_state);
        t_inc=w_c>0 & ~isnan(c_r);
        c_r=c_r(t_inc);
        w_c=w_c(t_inc);
        PR_UR(ss,uu)=sum(w_c(:).*c_r(:))./sum(w_c(:));
    end
end
state_nan=~isnan(mean(PR_UR,2));
PR_UR=PR_UR(state_nan,:);
State_Name=State_Name(state_nan);

State_Risk=zeros(length(State_Name),1);
for ss=1:length(State_Name)
    State_Risk(ss)=PR_UR(ss,:)*ur_level(:);
end

[~,Indx_Risk]=sort(State_Risk,'descend');
State_Name=State_Name(Indx_Risk);
PR_UR=PR_UR(Indx_Risk,:);


subplot('Position',[0.52,0.18,0.47,0.81]);

for uu=1:length(ur_level)
    for ss=1:length(State_Name)
        patch(ss+[-0.5 -0.5 0.5 0.5],100.*(ur_level(uu)+[-0.025 0.025 0.025 -0.025]),interp1(linspace(min(PR_UR(:)),max(PR_UR(:)),size(C_Risk,1)),C_Risk,PR_UR(ss,uu)),'LineStyle','none');
    end
end

xlim([0.5 length(State_Name)+0.5])
ylim([2.5 97.5])
set(gca,'XTick',[1:length(State_Name)],'XTickLabel',State_Name,'YTick',100.*ur_level,'Fontsize',16,'LineWidth',2,'TickDir','out', 'YDir','reverse');
ytickformat('percentage');

ylabel('National underreporting rate','Fontsize',18)
xlabel('State','Fontsize',18)


ax1.Position=[-0.445,-0.2,1.3,1.3];


print(gcf,['Preliminary_Underreporting_Plot.png'],'-dpng','-r300');
end

