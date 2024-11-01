clear;
clc;

temp_dir=pwd;
temp_dir=temp_dir(1:(length(temp_dir)-length('Data')));
S=shaperead([temp_dir '/Shapefile/cb_2021_us_county_500k.shp'],'UseGeoCoords',true);
Migratory_Bird=readgeotable([temp_dir '/Shapefile/cb_2021_us_county_500k.shp']);
Migratory_Bird=Migratory_Bird(:,[2 3 6 7 10 9 12 13]);
Var_Name=Migratory_Bird.Properties.VariableNames;
Var_Name{end-1}='AREA_LAND';
Var_Name{end}='AREA_WATER';

lat_county=zeros(length(S),1);
lon_county=zeros(length(S),1);

Data=readtable('Migratory_Bird_Flyways.csv');
D_Lon=[Data.BNODELONG;Data.RNODELONG];
D_Lat=[Data.BNODELAT;Data.RNODELAT];

Z_Count=zeros(length(S),1);
parfor cc=1:length(S)
    polyin = polyshape(S(cc).Lon,S(cc).Lat);
    [lon_county(cc),lat_county(cc)] = centroid(polyin);
    [p_in,p_on]=inpolygon(D_Lon,D_Lat,S(cc).Lon,S(cc).Lat);
    Z_Count(cc)=sum(p_in)+sum(p_on);
end

Z_t=log(Z_Count(Z_Count>0));
lon_t=lon_county(Z_Count>0);
lat_t=lat_county(Z_Count>0);


Vq = scatteredInterpolant(lon_t(:),lat_t(:),Z_t(:));

Test_V = exp(Vq(lon_county(:),lat_county(:))');

risk_measure_migratory_bird=log(Test_V);
tf=risk_measure_migratory_bird>prctile(risk_measure_migratory_bird,97.5);
vf=prctile(risk_measure_migratory_bird,97.5);

tg=risk_measure_migratory_bird<prctile(risk_measure_migratory_bird,2.5);
vg=prctile(risk_measure_migratory_bird,2.5);

risk_measure_migratory_bird(tf)=vf;
risk_measure_migratory_bird(tg)=vg;

Migratory_Bird.risk_measure_migratory_bird=risk_measure_migratory_bird(:);
save('Migratory_Bird_Assesment.mat','Migratory_Bird');
