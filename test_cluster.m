clear;
clc;
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

load('Average_Risk_Dairy.mat','outbreak_risk_dairy_farm_County','w_AIC','spillover_risk_dairy_farm_County')


clearvars -except S outbreak_risk_dairy_farm_County w_AIC spillover_risk_dairy_farm_County

NS=length(S);

County_Centroid=zeros(NS,2);
for cc=1:NS
    [County_Centroid(cc,1),County_Centroid(cc,2)]=centroid(polyshape(S(cc).Lat,S(cc).Lon));
end


D_County=zeros(NS,NS);
for cc=1:NS
    D_County(cc,:)=(deg2sm(distance(County_Centroid(cc,1),County_Centroid(cc,2),County_Centroid(:,1),County_Centroid(:,2))));
end

MR=zeros(NS,1);
for cc=1:NS
    MR(cc)=median(D_County(cc,:));
end

R_max=ceil(median(MR));

L=zeros(NS,1);
R_c=zeros(NS,1);
HR=zeros(length(w_AIC),NS);
opts=optimoptions('fmincon','FunctionTolerance',10^(-12),'MaxFunctionEvaluations',10^6,'MaxIterations',10^6,'UseParallel',false);

    risk_measure=outbreak_risk_dairy_farm_County*w_AIC;
    % risk_measure=risk_measure;
    C=sum(risk_measure);
    N=NS;
    L0=C.*log(C)+(N-C).*log(N-C)-N.*log(N);
    parfor ii=1:NS
        r0=linspace(-3,log10(R_max),51);
        Lt=zeros(51,1);
        for rr=1:51
            Lt(rr)=Objective_High_Risk_Regions(r0(rr),risk_measure,D_County(ii,:));
        end
        [r_est,L(ii)]=fmincon(@(x)Objective_High_Risk_Regions(x,risk_measure,D_County(ii,:)),min(r0(Lt==min(Lt))),[],[],[],[],-6,log10(R_max),[],opts);
        R_c(ii)=10.^r_est;
    end
    L=-L;
    
    lambda_R=exp(L)./exp(L0);
    
MC_lambda_R=zeros(length(lambda_R),999);
tot_risk_measure=sum(risk_measure);

risk_measure_MC=zeros(length(lambda_R),999);

for mm=1:999
    r=rand(length(lambda_R),1);
    r(risk_measure==0)=0;
    r=tot_risk_measure.*r./sum(r);
    risk_measure_MC(:,mm)=r;
    parfor ii=1:NS
        r0=linspace(-3,log10(R_max),51);
        Lt=zeros(51,1);
        for rr=1:51
            Lt(rr)=Objective_High_Risk_Regions(r0(rr),r,D_County(ii,:));
        end
        [r_est,L(ii)]=fmincon(@(x)Objective_High_Risk_Regions(x,r,D_County(ii,:)),min(r0(Lt==min(Lt))),[],[],[],[],-6,log10(R_max),[],opts);
    end
    L=-L;
    
    MC_lambda_R(:,mm)=exp(L)./exp(L0);
end

p_alpha=zeros(3138,1);
for ii=1:3138
    p_alpha(ii)=1-mean(lambda_R(ii)>MC_lambda_R(ii,:));
end

ff=find(p_alpha<0.01 & lambda_R>median(lambda_R(p_alpha<0.01))); % Find the significant ones and those with the larger lambda_R

figure('units','normalized','outerposition',[0 0.075 1 1]);
ax1=usamap('conus');

framem off; gridm off; mlabel off; plabel off;

states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax1, states,'Facecolor','none','LineWidth',0.5); hold on;


CC_Risk=ones(NS,3);
m=median(risk_measure(risk_measure>0));
ub=prctile(risk_measure(risk_measure>0),99);
xx=linspace(m,ub,1001);
L=zeros(1001,1);
for ii=1:1001
    x_opt=double(risk_measure>xx(ii));
    L(ii)=-Objective_High_Risk_Regions(x_opt,risk_measure);
end
hr=median(xx(L==max(L)));
CC_Risk(risk_measure>hr,:)=repmat([1 0 0],sum(risk_measure>hr),1);

CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Risk});
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); 
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); hold on;