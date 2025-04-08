clear;
clc;

temp_dir=pwd;
temp_dir=temp_dir(1:(length(temp_dir)-length('Data')));
US_County=readgeotable([temp_dir '/Shapefile/cb_2021_us_county_500k.shp']);
US_County=US_County(:,[2 3 6 7 10 9 12 13]);
Var_Name=US_County.Properties.VariableNames;
Var_Name{end-1}='AREA_LAND';
Var_Name{end}='AREA_WATER';
US_County.Properties.VariableNames=Var_Name;

[US_County,Indx]=sortrows(US_County,[1 2]);

County_remove=strcmp("AK",US_County.STUSPS) | strcmp("HI",US_County.STUSPS) | strcmp("AS",US_County.STUSPS) | strcmp("GU",US_County.STUSPS) | strcmp("MP",US_County.STUSPS) | strcmp("PR",US_County.STUSPS) | strcmp("VI",US_County.STUSPS);
US_County=US_County(~County_remove,:);

US_County_Dairy_to_Human=US_County(:,1:end-2);
US_County_Poultry_to_Human=US_County(:,1:end-2);
temp_county_fp=str2double(US_County.COUNTYFP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Temp=readtable([pwd '\Temperature\County_Average_Temperature_2024.csv']);

temp_ID=cell(height(US_County),1);

for ii=1:length(temp_ID)
    temp_ID{ii}=[US_County.STUSPS{ii} '-' US_County.COUNTYFP{ii}];
end


X_temp=NaN.*zeros(height(US_County),1);
for cc=1:height(US_County)
    t_find= strcmpi(Temp.State,US_County.STATE_NAME{cc}) & strcmpi(Temp.ID,temp_ID{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Temp.Value(t_find);
    else
         t_find= strcmpi(Temp.State,US_County.STATE_NAME{cc});
         if(sum(t_find)>0)
            X_temp(cc)=mean(Temp.Value(cc));
         end
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'TEMP'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Light_Intensity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd '\Bird_Stopover\County_Level_Stopover.mat'],'US_County_Stopover');
X_temp=zeros(height(US_County),1);
for cc=1:height(US_County_Stopover)
    t_find=strcmpi(US_County_Stopover.STUSPS{cc},US_County.STUSPS) & strcmpi(US_County_Stopover.NAME{cc},US_County.NAME) & strcmpi(US_County_Stopover.GEOID{cc},US_County.GEOID);
    if(sum(t_find)>0)
        X_temp(t_find)=(US_County_Stopover.Spring_Stopover(cc)+US_County_Stopover.Fall_Stopover(cc))./2;
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'LIGHT_INT'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Waterfowl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd '\Flyway_waterfowl\County_Waterfowl_Flyway.mat'],'Water_Fowl');
X_temp=zeros(height(US_County),4);

for cc=1:length(Water_Fowl.STUSPS)
    t_find=strcmpi(Water_Fowl.STUSPS{cc},US_County.STUSPS) & strcmpi(Water_Fowl.NAME{cc},US_County.NAME) & strcmpi(Water_Fowl.GEOID{cc},US_County.GEOID);
    if(sum(t_find)>0)
        X_temp(t_find,1)=Water_Fowl.Mallard(cc);
        X_temp(t_find,2)=Water_Fowl.Canada_Goose(cc);
        X_temp(t_find,3)=Water_Fowl.American_Green_Winged_Teal(cc);
        X_temp(t_find,4)=Water_Fowl.Northern_Pintail(cc);
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'MALLARD','CANADA_GOOSE','AGW_TEAL','NORTH_PINTAIL'};
US_County=[US_County T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Waterfowl Flyway
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd '\North_America_Flyway\North_American_Flyway_County.mat'],'Flyway_Water_Fowl');

FL={'Atlantic Flyway','Mississippi Flyway','Pacific Flyway','Central Flyway'};
X_temp=zeros(height(US_County),4);
for cc=1:length(Flyway_Water_Fowl.STUSPS)
    t_find=strcmpi(Flyway_Water_Fowl.STUSPS{cc},US_County.STUSPS) & strcmpi(Flyway_Water_Fowl.NAME{cc},US_County.NAME) & strcmpi(Flyway_Water_Fowl.GEOID{cc},US_County.GEOID);
    if(sum(t_find)>0)
        X_temp(t_find,strcmp(Flyway_Water_Fowl.FLYWAY{cc},FL))=1;
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'ATLANTIC_FLYWAY','MISSISSIPPI_FLYWAY','PACIFIC_FLYWAY','CENTRAL_FLYWAY'};
US_County=[US_County T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Milk operations with Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Dairy\County_Operations_with_Inventory_Cattle_Milk_2022.csv']);
Test=unique(Data.DomainCategory);
Test=Test([1 2 4 6 3 5 7]);

X_temp=zeros(height(US_County),length(Test));

for tt=1:length(Test)
    for cc=1:height(US_County)
        t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc}) & strcmp(Data.DomainCategory,Test{tt});
        if(sum(t_find)>0)
            X_temp(cc,tt)=Data.Value(t_find);
        end
    end
end

X_temp=[X_temp sum(X_temp,2)];
T_temp=array2table(X_temp);
for jj=1:length(Test)
    T_temp.Properties.VariableNames{jj}=['CATTLE ' Test{jj}];
end
T_temp.Properties.VariableNames{end}='TOTAL_DAIRY_OPERATIONS';
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inferred Cattle Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Dairy\Inferred_County_Milk_Cow_Inventory_2022.csv']);

X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Value(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'CATTLE_INVENTORY'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connectivity: Dairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd '\Dairy_Network\Dairy_County_Network.mat'],'US_Dairy_County');
US_Dairy_County=US_Dairy_County(Indx,:); % Need to reorder the rows relative to the indx shuffle 
US_Dairy_County=US_Dairy_County(~County_remove,:); % Need to reorder the rows relative to the indx shuffle 

T_temp=US_Dairy_County(:,end);
test=table2array(T_temp);
test=test(:,Indx); % Need to reorder the columns relative to the indx shuffle 
test=test(:,~County_remove);

T_temp=table(test);
T_temp.Properties.VariableNames={'CONNECT_DAIRY'};
US_County=[US_County T_temp];
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Bovine Testing
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Data=readtable([pwd '\Cattle\Cattle_Inventory_State_2022.csv']);
% 
% X_temp=zeros(height(US_County),1);
% 
% for cc=1:height(US_County)
%     t_find=strcmpi(Data.State,US_County.STATE_NAME{cc});
%     if(sum(t_find)>0)
%         X_temp(cc)=Data.Value(t_find);
%     end
% end
% 
% Data=readtable([pwd '\Cattle\Bovine_TB_Test_State_2022.csv']);
% X_temp2=zeros(height(US_County),1);
% 
% for cc=1:height(US_County)
%     t_find=strcmpi(Data.State,US_County.STATE_NAME{cc});
%     if(sum(t_find)>0)
%         X_temp2(cc)=Data.Value(t_find);
%     end
% end
% 
% X_temp=X_temp2./X_temp;
% 
% T_temp=array2table(X_temp);
% T_temp.Properties.VariableNames={'SURVIELLANCE_PROXY'};
% US_County=[US_County T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State H5N1 among dairy farms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\H5N1_Outbreaks\HPAI_Dairy_Farms.xlsx'],'Sheet','Dairy');

X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(Data.State,US_County.STATE_NAME{cc}) & strcmpi(Data.County,US_County.NAME{cc}) & Data.CountySuppressed==0;
    X_temp(cc)=sum(t_find);
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'DAIRY_HPAI_OUTBREAK_KNOWN'};
US_County=[US_County T_temp];


U_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(Data.State,US_County.STATE_NAME{cc}) & Data.CountySuppressed==1;
    t_opr=strcmpi(US_County.STATE_NAME,US_County.STATE_NAME{cc});
    temp_TOTAL=sum(US_County.TOTAL_DAIRY_OPERATIONS(t_opr));
    U_temp(cc,1)=sum(t_find).*US_County.TOTAL_DAIRY_OPERATIONS(cc)./temp_TOTAL;
end
T_temp=array2table(U_temp);
T_temp.Properties.VariableNames={'DAIRY_HPAI_OUTBREAK_UNKNOWN'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inferred Broiler Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Poultry\Inferred_Broiler_Inventory_County_2022.csv']);

X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Value(t_find);
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'BROILER_INVENTORY'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inferred Layer Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Poultry\Inferred_Layers_Inventory_County_2022.csv']);
X_temp=zeros(height(US_County),1);
for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=X_temp(cc)+Data.Value(t_find);
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'LAYER_INVENTORY'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pullet Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Poultry\Inferred_Pullet_Inventory_County_2022.csv']);
X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=X_temp(cc)+Data.Value(t_find);
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'PULLET_INVENTORY'};
US_County=[US_County T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turkey Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data=readtable([pwd '\Poultry\Inferred_Turkey_Inventory_County_2022.csv']);
X_temp=zeros(height(US_County),1);
for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=X_temp(cc)+Data.Value(t_find);
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'TURKEY_INVENTORY'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Poulty operations with Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Poultry\County_Operations_with_Inventory_Poultry_2022.csv']);

X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Value(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'POULTRY_OPR_w_INVENTORY'};
US_County=[US_County T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Broiler operations with Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Poultry\County_Level_Broiler_Operations_with_Inventory_2022.csv']);

X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Value(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'BROILER_OPR_w_INVENTORY'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer operations with Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Poultry\County_Level_Layer_Operations_with_Inventory_2022.csv']);

X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Value(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'LAYER_OPR_w_INVENTORY'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pullet operations with Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Poultry\County_Level_Pullet_Operations_with_Inventory_2022.csv']);

X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Value(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'PULLET_OPR_w_INVENTORY'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turkey operations with Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Poultry\County_Turkey_Farms_with_Inventory_2022.csv']);

X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Value(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'TURKEY_OPR_w_INVENTORY'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% County H5N1 among poultry farms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\H5N1_Outbreaks\HPAI_Poulty_Farms.xlsx'],'Sheet','Poultry');

X_temp=zeros(height(US_County),1);
% Total_Stratified_HPAI_Outbreak=max(Data.Total_Indicator,1);
for cc=1:height(US_County)
    t_find=strcmpi(Data.State,US_County.STATE_NAME{cc}) & strcmpi(Data.CountyName,US_County.NAME{cc});
    % X_temp(cc,2)=sum(Data.Pullet_Indicator(t_find)./Total_Stratified_HPAI_Outbreak(t_find));
    % X_temp(cc,3)=sum(Data.Broiler_Indicator(t_find)./Total_Stratified_HPAI_Outbreak(t_find));
    % X_temp(cc,4)=sum(Data.Layer_Indicator(t_find)./Total_Stratified_HPAI_Outbreak(t_find));
    % X_temp(cc,5)=sum(Data.Turkey_Indicator(t_find)./Total_Stratified_HPAI_Outbreak(t_find));
    X_temp(cc,1)=sum(t_find); %-sum(X_temp(cc,2:5),2);
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'POULTRY_HPAI_OUTBREAK'};%,'PULLET_HPAI_OUTBREAK','BROILER_HPAI_OUTBREAK','LAYER_HPAI_OUTBREAK','TURKEY_HPAI_OUTBREAK'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H5N1 Cases Among Humans: Dairy Connected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\H5N1_Outbreaks\H5N1_State_Human.csv']);
Data=Data(Data.Dairy>0,:);
Spillover=zeros(height(US_County),1);

for ii=1:height(Data)
    t_state=strcmpi(Data.StateName{ii},US_County.STATE_NAME) & (US_County.DAIRY_HPAI_OUTBREAK_UNKNOWN+US_County.DAIRY_HPAI_OUTBREAK_KNOWN>0);
    w_county=US_County.DAIRY_HPAI_OUTBREAK_UNKNOWN+US_County.DAIRY_HPAI_OUTBREAK_KNOWN;
    
    Spillover(t_state)=Data.Dairy(ii).*w_county(t_state)./sum(w_county(t_state));
end

T_temp=[array2table(Spillover)];
T_temp.Properties.VariableNames={'SPILLOVER_DAIRY'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H5N1 Cases Among Humans: Poultry Connected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\H5N1_Outbreaks\H5N1_State_Human.csv']);
Data=Data(Data.Poultry>0,:);
Spillover=zeros(height(US_County),1);


OB=US_County.POULTRY_HPAI_OUTBREAK; %+US_County.PULLET_HPAI_OUTBREAK+US_County.BROILER_HPAI_OUTBREAK+US_County.LAYER_HPAI_OUTBREAK+US_County.TURKEY_HPAI_OUTBREAK;
for ii=1:height(Data)
    t_state=strcmpi(Data.StateName{ii},US_County.STATE_NAME) & OB>0;
    w_county=OB;    
    Spillover(t_state)=Data.Poultry(ii).*w_county(t_state)./sum(w_county(t_state));
end

T_temp=[array2table(Spillover)];
T_temp.Properties.VariableNames={'SPILLOVER_POULTRY'};
US_County=[US_County T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('Data_US_County.mat','US_County');