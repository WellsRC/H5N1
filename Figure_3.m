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

load('Onward_Transmission.mat','mle_County_dairy_onward','mle_State_dairy_onward','mle_County_onward','mle_State_onward');

no_farms=isnan(mle_County_onward);
avg_onward_transmission_County=mle_County_onward;
avg_onward_transmission_Dairy_County=mle_County_dairy_onward;
avg_onward_transmission_State=mle_State_onward;
avg_onward_transmission_Dairy_State=mle_State_dairy_onward;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Figure 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
figure('units','normalized','outerposition',[0.1 0.15 0.8 0.7]);
 ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,0.4,0.6,0.6];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;

subplot('Position',[0.61,0.035,0.02,0.94]);

Title_Name={'Likelihood of onward transmission'};

risk_measure=log10(avg_onward_transmission_County);
C_Risk=[hex2rgb('#fcfbfd');
hex2rgb('#efedf5');
hex2rgb('#dadaeb');
hex2rgb('#bcbddc');
hex2rgb('#9e9ac8');
hex2rgb('#807dba');
hex2rgb('#6a51a3');
hex2rgb('#54278f');
hex2rgb('#3f007d');];

y_indx=linspace(-5,0,size(C_Risk,1));
x_risk=y_indx;

 x_indx=y_indx;
 c_indx=linspace(y_indx(1),y_indx(end),1001);
dx=c_indx(2)-c_indx(1);
xlim([0 1]);
ylim([y_indx(1) y_indx(end)+dx])
ymin=0.75;
dy=2/(1+sqrt(5));
for ii=1:length(c_indx)
    patch([0 0 dy dy],c_indx(ii)+[dx -dx -dx dx],interp1(x_risk,C_Risk,c_indx(ii)),'LineStyle','none');
end

patch([0 0 dy dy], [y_indx(1) y_indx(end)+dx y_indx(end)+dx y_indx(1)],'k','FaceAlpha',0,'LineWidth',2);

yl_indx=[y_indx(1):y_indx(end)];
for yy=1:length(yl_indx)
    text(ymin,yl_indx(yy),['10^{' num2str(yl_indx(yy)) '}'],'Fontsize',14);            
end

text(ymin+1.65,-2.5,Title_Name,'HorizontalAlignment','center','Fontsize',16,'Rotation',270);

axis off;  

NS=length(S);
CC_Risk=interp1(x_risk,C_Risk,risk_measure);

CC_Risk(no_farms,:)=repmat([0.5 0.5 0.5],sum(no_farms),1);

CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;

subplot('Position',[0.75,0.09,0.24,0.90]);


[~,indxs]=sort(avg_onward_transmission_State,'ascend');
barh(State_Name(indxs),avg_onward_transmission_State(indxs),'FaceColor',hex2rgb('#3f007d'),'LineStyle','none');
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',11,'XTick',[0:0.1:0.5],'XMinorTick','on');

xlabel({'Likelihood of onward transmission'},'Fontsize',16);
text(-0.35,1.02,'B','Fontsize',28,'Units','normalized','VerticalAlignment','top','HorizontalAlignment','center');

ax1.Position=[-0.29,-0.12,1.2,1.2];
text(-3.075,1.02,'A','Fontsize',28,'Units','normalized','VerticalAlignment','top','HorizontalAlignment','center');

print(gcf,['Figure_3.png'],'-dpng','-r300');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % Supplmental Figure
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


figure('units','normalized','outerposition',[0.1 0.15 0.8 0.7]);
 ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,0.4,0.6,0.6];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;

subplot('Position',[0.61,0.035,0.02,0.94]);

Title_Name={'Relative contribution of dairy to onward transmission'};

risk_measure=(avg_onward_transmission_Dairy_County);
C_Risk=[hex2rgb('#ffffff');
hex2rgb('#f0f0f0');
hex2rgb('#d9d9d9');
hex2rgb('#bdbdbd');
hex2rgb('#969696');
hex2rgb('#737373');
hex2rgb('#525252');
hex2rgb('#252525');
hex2rgb('#000000');];

y_indx=linspace(0,1,size(C_Risk,1));
x_risk=y_indx;

 x_indx=y_indx;
 c_indx=linspace(y_indx(1),y_indx(end),1001);
dx=c_indx(2)-c_indx(1);
xlim([0 1]);
ylim([y_indx(1) y_indx(end)+dx])
ymin=0.75;
dy=2/(1+sqrt(5));
for ii=1:length(c_indx)
    patch([0 0 dy dy],c_indx(ii)+[dx -dx -dx dx],interp1(x_risk,C_Risk,c_indx(ii)),'LineStyle','none');
end

patch([0 0 dy dy], [y_indx(1) y_indx(end)+dx y_indx(end)+dx y_indx(1)],'k','FaceAlpha',0,'LineWidth',2);

yl_indx=[0:0.1:1];
for yy=1:length(yl_indx)
    text(ymin,yl_indx(yy),[num2str(yl_indx(yy))],'Fontsize',14);            
end

text(ymin+1.65,0.5,Title_Name,'HorizontalAlignment','center','Fontsize',16,'Rotation',270);

axis off;  

NS=length(S);
CC_Risk=interp1(x_risk,C_Risk,risk_measure);

CC_Risk(no_farms,:)=repmat([0.5 0.5 0.5],sum(no_farms),1);

CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;

subplot('Position',[0.75,0.09,0.24,0.90]);


[~,indxs]=sort(avg_onward_transmission_State,'ascend');
b=barh(State_Name(indxs),[avg_onward_transmission_State(indxs).*avg_onward_transmission_Dairy_State(indxs) avg_onward_transmission_State(indxs).*(1-avg_onward_transmission_Dairy_State(indxs))],'stacked','LineStyle','none');
b(1).FaceColor=[0 0 0];
b(2).FaceColor=[0.7 0.7 0.7];

for ii=1:length(indxs)
    p=avg_onward_transmission_Dairy_State(indxs(ii));
    text(avg_onward_transmission_State(indxs(ii))+0.01,ii,[num2str(100.*p,'%3.1f') '%'],'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',10);
end
box off;
set(gca,'LineWidth',2,'tickdir','out','Fontsize',11,'XTick',[0:0.1:0.5],'XMinorTick','on');
legend({'Dairy','Poultry'},'Location','southeast','Fontsize',14);
legend boxoff;
xlabel({'Likelihood of onward transmission'},'Fontsize',16);
text(-0.35,1.02,'B','Fontsize',28,'Units','normalized','VerticalAlignment','top','HorizontalAlignment','center');

ax1.Position=[-0.29,-0.12,1.2,1.2];
text(-3.075,1.02,'A','Fontsize',28,'Units','normalized','VerticalAlignment','top','HorizontalAlignment','center');

print(gcf,['Supplemental_Figure_Dairy_Proportion.png'],'-dpng','-r300');