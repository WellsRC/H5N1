clear;
clc;
load('Data_US_County.mat','US_County');

US_County_Spillover=US_County(:,1:12);
clear US_County;
temp_county_fp=str2double(US_County_Spillover.COUNTYFP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Migratory Birds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd '\Slaughterhouse\Location_Slaughterhouse.mat']);

temp_GEOID=str2double(US_County_Spillover.GEOID);
X_temp=zeros(height(US_County_Spillover),1);
for cc=1:length(temp_GEOID)
    t_find=County_FIPS==temp_GEOID(cc);
    X_temp(cc)=sum(t_find);    
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'SLAUGHTERHOUSE'};
US_County_Spillover=[US_County_Spillover T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inferred Cattle Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Swine\Inferred_Hogs_Inventory_County_2022.csv']);

X_temp=zeros(height(US_County_Spillover),1);

for cc=1:height(US_County_Spillover)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County_Spillover.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Value(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'HOG_INVENTORY'};
US_County_Spillover=[US_County_Spillover T_temp];

save('Spillover_Exacerbation.mat','US_County_Spillover');