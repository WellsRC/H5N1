clc;
clear;
close all;

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


[~,~,~,County_Farms,Affected_County_Farms,~,~,~] = Poultry_Covariates({},{});

W=zeros(length(S));

for ii=1:length(S)
    % p_ii=polyshape(S(ii).Lat,S(ii).Lon,'Simplify',false);
    for jj=(ii+1):length(S)
        [~,q]=inpolygon(S(jj).Lat,S(jj).Lon,S(ii).Lat,S(ii).Lon);
        if(sum(q)>0)
            W(ii,jj)=1;
            W(jj,ii)=1;
        end
    end
end

Wt=W(County_Farms>0,:);
Wt=Wt(:,County_Farms>0);

z=Affected_County_Farms(County_Farms>0)-mean(Affected_County_Farms(County_Farms>0));
n=length(z);


% https://www.paulamoraga.com/book-spatial/spatial-autocorrelation.html#local-morans-i
% https://manuals.pqstat.pl/en:przestrzenpl:lokalpl:lmoranpl
 % https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-statistics/h-local-morans-i-additional-math.htm
GI=0;
for ii=1:length(z)
    GI=GI+z(ii).*(Wt(ii,:)*z);    
end
GI=GI./sum(z.^2);

GI=length(z).*GI./sum(Wt(:));

x=mean(Affected_County_Farms(County_Farms>0));

I=zeros(length(z),1);
E=zeros(length(z),1);
V=zeros(length(z),1);
for ii=1:length(I)
    I(ii)=n.*z(ii).*(Wt(ii,:)*z)./sum(z.^2);
    E(ii)=-sum(Wt(ii,:))/(n-1);
    zt=z(~ismember(1:n,ii));
    b2=(n-1).*sum(zt.^4)./(sum(zt.^2)^2);
    wi2=sum(Wt(ii,:).^2);
    wikh=sum(Wt(ii,:))^2;
    V(ii)=wi2*(n-b2)/(n-1)+wikh*(2*b2-n)/((n-1)*(n-2))-E(ii)^2;
end


Z_stat=(I-E)./sqrt(V);

% ZZ=real(Z_stat);
pv=normcdf(Z_stat);
pv(Z_stat>0)=1-normcdf(Z_stat(Z_stat>0));

p_Err=NaN.*zeros(size(County_Farms));
c=1;
CCz=NaN.*zeros(size(County_Farms));
for cc=1:length(County_Farms)
    if(County_Farms(cc)>0)
        p_Err(cc)=pv(c);
        CCz(cc)=I(c);
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


CH=[hex2rgb('#d73027');
hex2rgb('#f46d43');
hex2rgb('#fdae61');
hex2rgb('#fee090');
hex2rgb('#ffffbf');
hex2rgb('#e0f3f8');
hex2rgb('#abd9e9');
hex2rgb('#74add1');
hex2rgb('#4575b4');];


for ii=1:length(S)
    if(~isnan(p_Err(ii)) && p_Err(ii)<0.05)
        CC_Err(ii,:)=interp1(linspace(-1,1,9),CH,CCz(ii));
        
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