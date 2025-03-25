clear;

temp_pwd=pwd;
temp_pwd=temp_pwd(1:end-length('Data/North_America_Flyway'));

S=shaperead([temp_pwd 'Shapefile/cb_2021_us_county_500k.shp'],'UseGeoCoords',true);

CS=zeros(length(S),2);
for ii=1:length(S)
    polyin = polyshape(S(ii).Lon,S(ii).Lat,'Simplify',false);
    [CS(ii,1),CS(ii,2)] = centroid(polyin);
end

F=shaperead(['WaterfowlFlyways.shp'],'UseGeoCoords',true);

County_Flyway=cell(1,length(S));

for ii=1:length(F)
    p=inpolygon(CS(:,1),CS(:,2),F(ii).Lon,F(ii).Lat);
    County_Flyway(p)={F(ii).NAME};
end

Flyway_Water_Fowl.STUSPS={S.STUSPS};
Flyway_Water_Fowl.NAME={S.NAME};
Flyway_Water_Fowl.GEOID={S.GEOID};
Flyway_Water_Fowl.FLYWAY=County_Flyway;

save('North_American_Flyway_County.mat','Flyway_Water_Fowl');


