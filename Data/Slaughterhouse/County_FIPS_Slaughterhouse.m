clear;
clc;
temp_pwd = pwd;
temp_pwd=temp_pwd(1:(length(temp_pwd)-length('Data\Slaughterhouse')));

S = shaperead([temp_pwd 'Shapefile/cb_2021_us_county_500k.shp'],'UseGeoCoords',true);

Z = readtable('ZIP-COUNTY-FIPS_2017-06.csv');
SH = readtable('slaughter_productionData.xlsx');

SH_long=NaN.*zeros(height(SH),1);
SH_lat=NaN.*zeros(height(SH),1);
County_FIPS=NaN.*zeros(height(SH),1);

test=zeros(height(SH),1);
for ii=1:height(SH)    
    temp=geoCode([ SH.Address{ii} ', ' SH.City{ii} ', ' SH.State{ii} ',' num2str(SH.Zip(ii)) ],'osm'); 
    if(~isnan(temp(1)))
        SH_lat(ii)=temp(1);
        SH_long(ii)=temp(2);
    else
        t_f=Z.ZIP==SH.Zip(ii);
        if(sum(t_f)==1)
            County_FIPS(ii)=Z.STCOUNTYFP(t_f);
        else
            temp=geoCode([ SH.City{ii} ', ' SH.State{ii} ',' num2str(SH.Zip(ii)) ],'osm'); 
            SH_lat(ii)=temp(1);
            SH_long(ii)=temp(2);
        end
    end
end

County_FIPS_temp=NaN.*zeros(size(County_FIPS));

for ii=1:length(S)
    [p_in,p_on]=inpolygon(SH_long,SH_lat,S(ii).Lon,S(ii).Lat);
    County_FIPS_temp(p_in|p_on)=str2double([S(ii).GEOID]);
end

County_FIPS(~isnan(SH_long))=County_FIPS_temp(~isnan(SH_long));
save('Location_Slaughterhouse.mat','County_FIPS');
