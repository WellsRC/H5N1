clear;
clc;

temp_dir=pwd;
temp_dir=temp_dir(1:(length(temp_dir)-length('Data/Dairy_Network')));
US_Dairy_County=readgeotable([temp_dir '/Shapefile/cb_2021_us_county_500k.shp']);

US_Dairy_County=US_Dairy_County(:,[7 10 9 6]);

US_Dairy_County.GEOID=str2double(US_Dairy_County.GEOID);

Dairy_Transport_i_to_j=zeros(height(US_Dairy_County));

for nn=1:1000
    Net_Num=readtable(['FromToKernelGenhwgaall' num2str(nn) '.txt']);
    for ii=1:height(US_Dairy_County)
        t_out=Net_Num.Var1==US_Dairy_County.GEOID(ii);
        if(sum(t_out)>0)
            Net_Num_trim=Net_Num(t_out,:);
            for jj=1:height(Net_Num_trim)
                t_in=Net_Num_trim.Var2(jj)==US_Dairy_County.GEOID;
                if(sum(t_in)>0)
                    Dairy_Transport_i_to_j(ii,t_in)=Dairy_Transport_i_to_j(ii,t_in)+Net_Num_trim.Var3(jj);
                end
            end
        end
    end
end

Dairy_Transport_i_to_j=Dairy_Transport_i_to_j./1000;
US_Dairy_County.Dairy_Transport_i_to_j=Dairy_Transport_i_to_j;


save('Dairy_County_Network.mat','US_Dairy_County');
