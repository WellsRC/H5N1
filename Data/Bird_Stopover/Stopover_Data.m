clear;
clc;

temp_pwd=pwd;
temp_pwd=temp_pwd(1:(length(temp_pwd)-length('Data\Bird_Stopover\')));
S=shaperead([temp_pwd '/Shapefile/cb_2021_us_county_500k.shp'],'UseGeoCoords',true);
US_County_Stopover=readgeotable([temp_pwd '/Shapefile/cb_2021_us_county_500k.shp']);
US_County_Stopover=US_County_Stopover(:,[2 3 6 7 10 9 12 13]);

Spring_Intensity=zeros(height(S),1);
Fall_Intensity=zeros(height(S),1);

[A_spring,R_spring] = readgeoraster('spring_stopover_2500_v9_265_class.tif');
info_spring = georasterinfo('spring_stopover_2500_v9_265_class.tif');
A_spring=double(squeeze(A_spring(:)));

x_spring=linspace(R_spring.XWorldLimits(1),R_spring.XWorldLimits(2),size(A_spring,2));
y_spring=flip(linspace(R_spring.YWorldLimits(1),R_spring.YWorldLimits(2),size(A_spring,1)));
[X_spring,Y_spring] = meshgrid(x_spring,y_spring);
proj_spring = info_spring.CoordinateReferenceSystem;
[lat_spring,lon_spring] = projinv(proj_spring,X_spring(:),Y_spring(:));

lat_spring=lat_spring(~isnan(A_spring));
lon_spring=lon_spring(~isnan(A_spring));
A_spring=A_spring(~isnan(A_spring));

for ss=1:length(S)
    [p_in,p_on]=inpolygon(lon_spring,lat_spring,S(ss).Lon,S(ss).Lat);
    Spring_Intensity(ss)=mean(A_spring(p_on | p_in));
end

[A_fall,R_fall] = readgeoraster('fall_stopover_2500_v9_265_class.tif');
info_fall = georasterinfo('fall_stopover_2500_v9_265_class.tif');
A_fall=double(squeeze(A_fall(:)));

x_fall=linspace(R_fall.XWorldLimits(1),R_fall.XWorldLimits(2),size(A_fall,2));
y_fall=flip(linspace(R_fall.YWorldLimits(1),R_fall.YWorldLimits(2),size(A_fall,1)));
[X_fall,Y_fall] = meshgrid(x_fall,y_fall);
proj_fall = info_fall.CoordinateReferenceSystem;
[lat_fall,lon_fall] = projinv(proj_fall,X_fall(:),Y_fall(:));

lat_fall=lat_fall(~isnan(A_fall));
lon_fall=lon_fall(~isnan(A_fall));
A_fall=A_fall(~isnan(A_fall));

for ss=1:length(S)
    [p_in,p_on]=inpolygon(lon_fall,lat_fall,S(ss).Lon,S(ss).Lat);
    Fall_Intensity(ss)=mean(A_fall(p_on | p_in));
end

US_County_Stopover.Spring_Stopover=Spring_Intensity;
US_County_Stopover.Fall_Stopover=Fall_Intensity;

save('County_Level_Stopover.mat','US_County_Stopover');


