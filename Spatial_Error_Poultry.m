% clc;
% clear;
% close all;
% 
% % Plot_Variable
% close all;
% states = shaperead('usastatelo', 'UseGeoCoords', true);
% S=shaperead([pwd '/Shapefile/cb_2021_us_county_500k.shp'],'UseGeoCoords',true);
% 
% US_County=readgeotable([pwd '/Shapefile/cb_2021_us_county_500k.shp']);
% US_County=US_County(:,[2 3 6 7 10 9 12 13]);
% Var_Name=US_County.Properties.VariableNames;
% Var_Name{end-1}='AREA_LAND';
% Var_Name{end}='AREA_WATER';
% US_County.Properties.VariableNames=Var_Name;
% 
% [US_County,Indx]=sortrows(US_County,[1 2]);
% 
% S=S(Indx);
% 
% County_remove=strcmp("AK",US_County.STUSPS) |strcmp("HI",US_County.STUSPS) | strcmp("AS",US_County.STUSPS) | strcmp("GU",US_County.STUSPS) | strcmp("MP",US_County.STUSPS) | strcmp("PR",US_County.STUSPS) | strcmp("VI",US_County.STUSPS);
% S=S(~County_remove,:);
% 
% NS=length(S);
% 
% load([pwd '/Data/Data_US_County.mat'],'US_County');
% 
% State_Name=unique(US_County.STATE_NAME);
% state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
% State_Name=State_Name(~state_remove);
% 
% [~,~,~,County_Farms,Affected_County_Farms,~,~,~] = Poultry_Covariates({},{},{});
% load('Average_Risk_Poultry.mat','w_AIC','outbreak_poultry_farm_County','State_Name');
% avg_outbreak_farm_County=outbreak_poultry_farm_County*w_AIC;
% 
% avg_outbreak_farm_County(County_Farms==0)=NaN;
% 
% Err=(Affected_County_Farms-avg_outbreak_farm_County);
% 
% figure('units','normalized','outerposition',[0.25 0.25 0.45 0.475]);
% ax1=usamap('conus');
% 
% framem off; gridm off; mlabel off; plabel off;
% 
% ax1.Position=[-1.25,-0.2,1.3,1.3];
% 
% states = shaperead('usastatelo', 'UseGeoCoords', true);
% geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;
% 
% 
% C_Err=[hex2rgb('#b2182b');
% hex2rgb('#d6604d');
% hex2rgb('#f4a582');
% hex2rgb('#fddbc7');
% hex2rgb('#f7f7f7');
% hex2rgb('#d1e5f0');
% hex2rgb('#92c5de');
% hex2rgb('#4393c3');
% hex2rgb('#2166ac');];
% 
% x_Err=[-45 -30 -15 -0.5 0 0.5 15 30 45];
% 
% 
% CC_Err=ones(length(S),3);
% 
% for ii=1:length(S)
%     if(~isnan(Err(ii)))
%         CC_Err(ii,:)=interp1(x_Err,C_Err,Err(ii));
%     else
%         CC_Err(ii,:)=[0.7 0.7 0.7];
%     end
% end
% 
% CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Err});
% geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
% geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;
% 
% subplot('Position',[0.8 0.03 0.03 0.94])
% 
% yy=linspace(-45,45,1001);
% for ii=1:1000
%     patch([0 0 1 1],[yy(ii) yy(ii+1) yy(ii+1) yy(ii)],interp1(x_Err,C_Err,yy(ii)),'LineStyle','none'); hold on;
% end
% 
% for er=-45:15:45
%     text(1.5,er,num2str(er),'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','middle');
% end
% 
% text(3,0,'Residual','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',270)
% 
% patch([0 0 1 1],[-45 45 45 -45],'k','LineWidth',2,'FaceAlpha',0);
% axis off;
% 
% ax1.Position=[-0.25,-0.2,1.3,1.3];

% W=zeros(length(S));
% 
% for ii=1:length(S)
%     % p_ii=polyshape(S(ii).Lat,S(ii).Lon,'Simplify',false);
%     for jj=(ii+1):length(S)
%         [~,q]=inpolygon(S(jj).Lat,S(jj).Lon,S(ii).Lat,S(ii).Lon);
%         if(sum(q)>0)
%             W(ii,jj)=1;
%             W(jj,ii)=1;
%         end
%     end
% end

Wt=W(~isnan(Err),:);
Wt=Wt(:,~isnan(Err));
z=Err(~isnan(Err))-mean(Err(~isnan(Err)));
n=length(z);

Wf=Wt;
for ss=1:length(z)
    Wf(ss,ss)=1;
end

I=zeros(length(z),1);
E=zeros(length(z),1);
V=zeros(length(z),1);
for ii=1:length(I)
    m2=sum(z.^2)./n;
    m4=sum(z.^4)./n;
    I(ii)=z(ii).*(Wf(ii,:)*z)./m2;
    E(ii)=-sum(Wf(ii,:))/(n-1);
    zt=z(~ismember(1:n,ii));
    b2=m4./(m2^2);
    wi2=sum(Wt(ii,:).^2);
    wikh=sum(Wt(ii,:))^2;
    V(ii)=wi2*(n-b2)/(n-1)+wikh*(2*b2-n)/((n-1)*(n-2))-E(ii)^2;
end


Z_stat=(I-E)./sqrt(V);

ZZ=real(Z_stat);
pv=normcdf(ZZ);
pv(ZZ>0)=1-normcdf(ZZ(ZZ>0));

p_Err=NaN.*zeros(size(Err));
c=1;
for cc=1:length(Err)
    if(~isnan(Err(cc)))
        p_Err(cc)=pv(c);
        c=c+1;
    end
end

figure('units','normalized','outerposition',[0.25 0.25 0.45 0.475]);
ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;

ax1.Position=[-1.25,-0.2,1.3,1.3];

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;

CC_Err=ones(length(S),3);

for ii=1:length(S)
    if(~isnan(p_Err(ii)) && p_Err(ii)<0.05 && Err(ii)>0)
        CC_Err(ii,:)=[1 0 0];
    elseif(~isnan(p_Err(ii)) && p_Err(ii)<0.05 && Err(ii)<0)
        CC_Err(ii,:)=[0 0 1];
    elseif(~isnan(p_Err(ii)))
        CC_Err(ii,:)=[1 1 1];
    else
        CC_Err(ii,:)=[0.7 0.7 0.7];
    end
end

CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Err});
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;

ax1.Position=[-0.25,-0.2,1.3,1.3];