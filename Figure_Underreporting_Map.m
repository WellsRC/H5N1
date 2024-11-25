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

State_Name=unique({S.STATE_NAME});
State_Risk=NaN.*zeros(size(State_Name));

load([pwd '/Data/Data_US_County.mat'],'US_County');
Underreporting_County=Prob_UR(:,abs(ur_level-Underreporting_Rate)==min(abs(ur_level-Underreporting_Rate)));

figure('units','normalized','outerposition',[0 0.075 1 0.55]);
 ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,0.4,0.6,0.6];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;
       
C_Risk=[hex2rgb('#ffffff');
        hex2rgb('#f7fcfd');
hex2rgb('#e0ecf4');
hex2rgb('#bfd3e6');
hex2rgb('#9ebcda');
hex2rgb('#8c96c6');
hex2rgb('#8c6bb1');
hex2rgb('#88419d');
hex2rgb('#810f7c');
hex2rgb('#4d004b');];



subplot('Position',[0.38,0.05,0.01,0.9]);

risk_measure=Underreporting_County;
t_upper=risk_measure>prctile(risk_measure,95);
v_upper=prctile(risk_measure,95);
t_lower=risk_measure<prctile(risk_measure,5);
v_lower=prctile(risk_measure,5);
risk_measure(t_lower)=v_lower;
risk_measure(t_upper)=v_upper;
risk_measure=risk_measure./v_lower;

if(ceil(max(risk_measure))>=10)
    risk_measure=log10(risk_measure);
    x_risk=linspace(0,max(risk_measure),size(C_Risk,1));
    c_indx=linspace(0,max(risk_measure),251);
    y_indx=[1:floor(max(risk_measure(risk_measure>0)))];
else
    x_risk=linspace(1,max(max(risk_measure),2),size(C_Risk,1));
    c_indx=linspace(1,max(max(risk_measure),2),251);
    y_indx=[2:ceil(10.*max(risk_measure))./10];
end

xlim([0 1]);
ylim([min(c_indx) max(c_indx)]);    
ymin=1;
dy=2/(1+sqrt(5));
for ii=1:length(c_indx)
    patch([0 0 dy dy],c_indx(ii)+[1 0 0 1],interp1(x_risk,C_Risk,c_indx(ii)),'LineStyle','none');
end

if(min(y_indx)==2)
    text(ymin+2.5,1,{'Baseline','lowest'},'Fontsize',16,'HorizontalAlignment','center');
else
    text(ymin+2.5,0,{'Baseline','lowest'},'Fontsize',16,'HorizontalAlignment','center');
end
for yy=1:length(y_indx)
    if(min(y_indx)==2)
        text(ymin,y_indx(yy),[num2str(y_indx(yy))],'Fontsize',16);
    else
        text(ymin,y_indx(yy),['10^' num2str(y_indx(yy))],'Fontsize',16);
    end
end

text(ymin+4.5,mean([min(c_indx) max(c_indx)]),{'Fold-increase in' ,'probability of underreporting'},'HorizontalAlignment','center','Fontsize',18,'Rotation',270);

axis off;  

text(-37.3,1.01,'A','FontSize',32,'Units','normalized');
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
        t_state=strcmp(State_Name{ss},US_County.STATE_NAME);
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

PR_UR=PR_UR./min(PR_UR(:));

subplot('Position',[0.52,0.445,0.47,0.54]);

for uu=1:length(ur_level)
    for ss=1:length(State_Name)
        patch(ss+[-0.5 -0.5 0.5 0.5],100.*(ur_level(uu)+[-0.05 0.05 0.05 -0.05]),interp1(linspace(0,max(PR_UR(:)),size(C_Risk,1)),C_Risk,PR_UR(ss,uu)),'LineStyle','none');
    end
end

xlim([0.5 length(State_Name)+0.5])
ylim([5 100])
set(gca,'XTick',[1:length(State_Name)],'XTickLabel',State_Name,'YTick',100.*ur_level(1:2:end),'Fontsize',14,'LineWidth',2,'TickDir','out', 'YDir','reverse');
ytickformat('percentage');

ylabel({'National','underreporting rate'},'Fontsize',18)
xl=xlabel('State','Fontsize',18);
xl.Position=xl.Position+[0 2 0];

text(-0.15,0.95,'B','FontSize',32,'Units','normalized');

subplot('Position',[0.52,0.085,0.44,0.04]);
   
ymin=-1;
dy=2/(1+sqrt(5));
for ii=1:length(c_indx)
    patch(max(c_indx)-[c_indx(ii)-[1 0 0 1]],[0 0 dy dy],interp1(x_risk,C_Risk,c_indx(ii)),'LineStyle','none');
end

ylim([0 1]);
xl=xlim;
% xlim([0 max(c_indx)]); 
text(max(xl),ymin,{'Baseline lowest'},'Fontsize',16,'HorizontalAlignment','center');
% for yy=1:length(y_indx)

yy=length(y_indx);
        if(min(y_indx)==2)
            % text(max(c_indx)-y_indx(yy),ymin,[num2str(y_indx(yy))],'Fontsize',16);
            text(0,ymin,[num2str(max(c_indx))],'Fontsize',16);
        else
            % text(max(c_indx)-y_indx(yy),ymin,['10^' num2str(y_indx(yy))],'Fontsize',16);
            text(0,ymin,['10^{' num2str(max(c_indx),2) '}'],'Fontsize',16);
        end
% end
text(mean(xl),ymin,{'Fold-increase in probability of underreporting'},'HorizontalAlignment','center','Fontsize',18);

axis off;  
ax1.Position=[-0.405,-0.15,1.2,1.2];


print(gcf,['Figure_Underreporting.png'],'-dpng','-r300');
end

