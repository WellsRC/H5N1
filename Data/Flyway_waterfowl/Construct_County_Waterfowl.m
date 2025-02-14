clear;
T=readtable('BandData_2014Feb12.csv');
Spec=unique(T.SPEC);
M=[T.RNODELONG T.RNODELAT];
R=unique(M,'rows');
D=zeros(size(R,1));
for ii=1:size(R,1)
    D(ii,:)=deg2km(distance(R(ii,2),R(ii,1),R(:,2),R(:,1)));    
end

mD=zeros(size(R,1),1);
for ii=1:size(R,1)
    temp_D=D(ii,:);
    temp_D=temp_D(temp_D>0);
    mD(ii)=min(temp_D);
end

r=median(mD);

s=lsqnonlin(@(x)norminv(0.99,0,x)-r,100,0,r);

A=zeros(size(D));

for ii=1:size(R,1)
    A(ii,:)=normpdf(D(ii,:),0,s)./normpdf(0,0,s);
end

Counts=zeros(length(R),1);

temp_pwd=pwd;
temp_pwd=temp_pwd(1:end-length('Data/Flyway_waterfowl'));

S=shaperead([temp_pwd 'Shapefile/cb_2021_us_county_500k.shp'],'UseGeoCoords',true);

CS=zeros(length(S),2);
for ii=1:length(S)
    polyin = polyshape(S(ii).Lon,S(ii).Lat);
    [CS(ii,1),CS(ii,2)] = centroid(polyin);
end

DS=zeros(length(S),size(R,1));

for ii=1:size(R,1)
    DS(:,ii)=deg2km(distance(R(ii,2),R(ii,1),CS(:,2),CS(:,1)));    
end


Water_Fowl.STUSPS={S.STUSPS};
Water_Fowl.NAME={S.NAME};
Water_Fowl.GEOID={S.GEOID};

for ss=1:length(Spec)
    for ii=1:length(R)
        tf = R(ii,1)==M(:,1) & R(ii,2)==M(:,2) & strcmp(T.SPEC,Spec{ss});
        Counts(ii)=sum(tf);
    end
    scale_c=A\Counts;

    WF_County=zeros(length(S),1);
    
    for ii=1:length(S)
        WF_County(ii)=Waterfowl_Surface(DS(ii,:),scale_c,s);
    end

    WF_County(WF_County<0)=0;
    if(strcmp(Spec{ss},'AGWT'))
        Water_Fowl.American_Green_Winged_Teal=WF_County;
    elseif(strcmp(Spec{ss},'CAGO'))
        Water_Fowl.Canada_Goose=WF_County;
    elseif(strcmp(Spec{ss},'MALL'))
        Water_Fowl.Mallard=WF_County;
    elseif(strcmp(Spec{ss},'NOPI'))
        Water_Fowl.Northern_Pintail=WF_County;
    end
end

save('County_Waterfowl_Flyway.mat','Water_Fowl');
