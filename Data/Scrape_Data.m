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

County_remove=strcmp("HI",US_County.STUSPS) | strcmp("AS",US_County.STUSPS) | strcmp("GU",US_County.STUSPS) | strcmp("MP",US_County.STUSPS) | strcmp("PR",US_County.STUSPS) | strcmp("VI",US_County.STUSPS);
US_County=US_County(~County_remove,:);

US_County_Dairy_to_Human=US_County(:,1:end-2);
US_County_Poultry_to_Human=US_County(:,1:end-2);
temp_county_fp=str2double(US_County.COUNTYFP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Migratory Birds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd '\Flyway\Migratory_Bird_Assesment.mat'],'Migratory_Bird');
X_temp=zeros(height(US_County),1);
for cc=1:height(Migratory_Bird)
    t_find=strcmpi(Migratory_Bird.STUSPS{cc},US_County.STUSPS) & strcmpi(Migratory_Bird.NAME{cc},US_County.NAME) & strcmpi(Migratory_Bird.GEOID{cc},US_County.GEOID);
    if(sum(t_find)>0)
        X_temp(t_find)=Migratory_Bird.risk_measure_migratory_bird(cc);
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'MIGRATORY_BIRD'};
US_County=[US_County T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% County H5N1 among migratory birds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\H5N1_Outbreaks\HPAI Detections in Wild Birds.csv']);

X_temp=zeros(height(US_County),3);
for cc=1:height(US_County)
    t_find=strcmpi(Data.State,US_County.STATE_NAME{cc}) & strcmpi(Data.County,US_County.NAME{cc}) & Data.CollectionDate<datetime('January 1, 2023') & Data.CollectionDate>=datetime('January 1, 2022');
    X_temp(cc,1)=sum(t_find);
    t_find=strcmpi(Data.State,US_County.STATE_NAME{cc}) & strcmpi(Data.County,US_County.NAME{cc}) & Data.CollectionDate<datetime('January 1, 2024') & Data.CollectionDate>=datetime('January 1, 2023');
    X_temp(cc,2)=sum(t_find);
    t_find=strcmpi(Data.State,US_County.STATE_NAME{cc}) & strcmpi(Data.County,US_County.NAME{cc}) & Data.CollectionDate<datetime('January 1, 2025') & Data.CollectionDate>=datetime('January 1, 2024');
    X_temp(cc,3)=sum(t_find);
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'HPAI_2022_MIGRATORY_BIRDS','HPAI_2023_MIGRATORY_BIRDS','HPAI_2024_MIGRATORY_BIRDS'};
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State H5N1 among dairy farms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\H5N1_Outbreaks\HPAI_Dairy_Farms.xlsx'],'Sheet','Dairy');

X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(Data.State,US_County.STATE_NAME{cc}) & strcmpi(Data.County,US_County.NAME{cc}) & Data.CountySuppressed==0;
    t_supp=strcmpi(Data.State,US_County.STATE_NAME{cc}) & Data.CountySuppressed==1 & US_County.TOTAL_DAIRY_OPERATIONS(cc)>0;
    if(sum(t_find)>0)
        X_temp(cc)=sum(t_find);
    elseif(sum(t_supp)>0)
        X_temp(cc)=NaN;
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'DAIRY_HPAI_OUTBREAK_KNOWN'};
US_County=[US_County T_temp];


U_temp=NaN.*zeros(height(US_County),2);

for cc=1:height(US_County)
    if(isnan(X_temp(cc)))
        t_find=strcmpi(Data.State,US_County.STATE_NAME{cc}) & Data.CountySuppressed==1;
        t_opr=strcmpi(US_County.STATE_NAME,US_County.STATE_NAME{cc}) & isnan(X_temp);
        temp_TOTAL=sum(US_County.TOTAL_DAIRY_OPERATIONS(t_opr));
        U_temp(cc,1)=sum(t_find);
        U_temp(cc,2)=US_County.TOTAL_DAIRY_OPERATIONS(cc)./temp_TOTAL;
    end
end
T_temp=array2table(U_temp);
T_temp.Properties.VariableNames={'DAIRY_HPAI_REMAIN_STATE','STATE_REM_WEIGHT'};
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
% Rooster Inventory - County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data=readtable([pwd '\Poultry\Inferred_Rooster_Inventory_County_2022.csv']);
X_temp=zeros(height(US_County),1);
for cc=1:height(US_County)
    t_find=temp_county_fp(cc)==Data.CountyANSI & strcmpi(Data.State,US_County.STATE_NAME{cc});
    if(sum(t_find)>0)
        X_temp(cc)=X_temp(cc)+Data.Value(t_find);
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'ROOSTER_INVENTORY'};
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

X_temp=zeros(height(US_County),5);
Total_Stratified_HPAI_Outbreak=max(Data.Total_Indicator,1);
for cc=1:height(US_County)
    t_find=strcmpi(Data.State,US_County.STATE_NAME{cc}) & strcmpi(Data.CountyName,US_County.NAME{cc});
    X_temp(cc,2)=sum(Data.Pullet_Indicator(t_find)./Total_Stratified_HPAI_Outbreak(t_find));
    X_temp(cc,3)=sum(Data.Broiler_Indicator(t_find)./Total_Stratified_HPAI_Outbreak(t_find));
    X_temp(cc,4)=sum(Data.Layer_Indicator(t_find)./Total_Stratified_HPAI_Outbreak(t_find));
    X_temp(cc,5)=sum(Data.Turkey_Indicator(t_find)./Total_Stratified_HPAI_Outbreak(t_find));
    X_temp(cc,1)=sum(t_find)-sum(X_temp(cc,2:5),2);
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'POULTRY_HPAI_OUTBREAK','PULLET_HPAI_OUTBREAK','BROILER_HPAI_OUTBREAK','LAYER_HPAI_OUTBREAK','TURKEY_HPAI_OUTBREAK'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GINI Index: 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Gini_Index\ACSDT5Y2010.B19083-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

X_temp=NaN.*zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Estimate__GiniIndex(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'GINI_2010'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GINI Index: 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Gini_Index\ACSDT5Y2020.B19083-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

X_temp=NaN.*zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Estimate__GiniIndex(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'GINI_2020'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GINI Index: 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Gini_Index\ACSDT5Y2022.B19083-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

X_temp=NaN.*zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Estimate__GiniIndex(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'GINI_2022'};
US_County=[US_County T_temp];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill in GINI Index gaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_fit= ~isnan(US_County.GINI_2022) & ~isnan(US_County.GINI_2020) & ~isnan(US_County.GINI_2010);
f_fill_2010= ~isnan(US_County.GINI_2022) & ~isnan(US_County.GINI_2020) & isnan(US_County.GINI_2010);
f_fill_2020= ~isnan(US_County.GINI_2022) & isnan(US_County.GINI_2020) & isnan(US_County.GINI_2010);
f_fill_2022= isnan(US_County.GINI_2022) & ~isnan(US_County.GINI_2020) & ~isnan(US_County.GINI_2010);
% 2010
X=[US_County.GINI_2020 US_County.GINI_2022];
Y=US_County.GINI_2010;
mdl_f=fitlm(X(f_fit,:),Y(f_fit));
US_County.GINI_2010(f_fill_2010)=predict(mdl_f,X(f_fill_2010,:));

% 2020
X=[US_County.GINI_2010 US_County.GINI_2022];
Y=US_County.GINI_2020;
mdl_f=fitlm(X(f_fit,:),Y(f_fit));
US_County.GINI_2020(f_fill_2020)=predict(mdl_f,X(f_fill_2020,:));

% 2022
X=[US_County.GINI_2010 US_County.GINI_2020];
Y=US_County.GINI_2022;
mdl_f=fitlm(X(f_fit,:),Y(f_fit));
US_County.GINI_2022(f_fill_2022)=predict(mdl_f,X(f_fill_2022,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Household income: 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Household_Income\ACSST5Y2010.S1901-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

X_temp=NaN.*zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Households__Estimate__MedianIncome_dollars_(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'HOUSEHOLD_INCOME_2010'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Household income: 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Household_Income\ACSST5Y2020.S1901-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

X_temp=NaN.*zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Estimate__Households__MedianIncome_dollars_(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'HOUSEHOLD_INCOME_2020'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Household income: 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Household_Income\ACSST5Y2022.S1901-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

X_temp=NaN.*zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.Estimate__Households__MedianIncome_dollars_(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'HOUSEHOLD_INCOME_2022'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill in HOUSEHOLD_INCOME gaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_fit= ~isnan(US_County.HOUSEHOLD_INCOME_2022) & ~isnan(US_County.HOUSEHOLD_INCOME_2020) & ~isnan(US_County.HOUSEHOLD_INCOME_2010);
f_fill_2010= ~isnan(US_County.HOUSEHOLD_INCOME_2022) & ~isnan(US_County.HOUSEHOLD_INCOME_2020) & isnan(US_County.HOUSEHOLD_INCOME_2010);
f_fill_2020= ~isnan(US_County.HOUSEHOLD_INCOME_2022) & isnan(US_County.HOUSEHOLD_INCOME_2020) & isnan(US_County.HOUSEHOLD_INCOME_2010);
f_fill_2022= isnan(US_County.HOUSEHOLD_INCOME_2022) & ~isnan(US_County.HOUSEHOLD_INCOME_2020) & ~isnan(US_County.HOUSEHOLD_INCOME_2010);
% 2010
X=[US_County.HOUSEHOLD_INCOME_2020 US_County.HOUSEHOLD_INCOME_2022];
Y=US_County.HOUSEHOLD_INCOME_2010;
mdl_f=fitlm(X(f_fit,:),Y(f_fit));
US_County.HOUSEHOLD_INCOME_2010(f_fill_2010)=predict(mdl_f,X(f_fill_2010,:));

% 2020
X=[US_County.HOUSEHOLD_INCOME_2010 US_County.HOUSEHOLD_INCOME_2022];
Y=US_County.HOUSEHOLD_INCOME_2020;
mdl_f=fitlm(X(f_fit,:),Y(f_fit));
US_County.HOUSEHOLD_INCOME_2020(f_fill_2020)=predict(mdl_f,X(f_fill_2020,:));

% 2022
X=[US_County.HOUSEHOLD_INCOME_2010 US_County.HOUSEHOLD_INCOME_2020];
Y=US_County.HOUSEHOLD_INCOME_2022;
mdl_f=fitlm(X(f_fit,:),Y(f_fit));
US_County.HOUSEHOLD_INCOME_2022(f_fill_2022)=predict(mdl_f,X(f_fill_2022,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Education: 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Education\ACSST5Y2010.S1501-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

X_temp=NaN.*zeros(height(US_County),1);
Est_Ed=(Data.Total__Estimate__Population18To24Years.*Data.Total__Estimate__LessThanHighSchoolGraduate+Data.Total__Estimate__Population25YearsAndOver.*(Data.Total__Estimate__9thTo12thGrade_NoDiploma+Data.Total__Estimate__LessThan9thGrade))./(Data.Total__Estimate__Population18To24Years+Data.Total__Estimate__Population25YearsAndOver);
for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Est_Ed(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'EDUCATION_2010'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Education: 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Education\ACSST5Y2020.S1501-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

Est_Ed=100.*(Data.Estimate__Total__AGEBYEDUCATIONALATTAINMENT__Population18To24_1+Data.Estimate__Total__AGEBYEDUCATIONALATTAINMENT__Population25Year_1+Data.Estimate__Total__AGEBYEDUCATIONALATTAINMENT__Population25Year_2)./(Data.Estimate__Total__AGEBYEDUCATIONALATTAINMENT__Population18To24Ye+Data.Estimate__Total__AGEBYEDUCATIONALATTAINMENT__Population25YearsA);
X_temp=NaN.*zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Est_Ed(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'EDUCATION_2020'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Education: 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Education\ACSST5Y2022.S1501-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end
Est_Ed=100.*(Data.Estimate__Total__AGEBYEDUCATIONALATTAINMENT__Population18To24_1+Data.Estimate__Total__AGEBYEDUCATIONALATTAINMENT__Population25Year_1+Data.Estimate__Total__AGEBYEDUCATIONALATTAINMENT__Population25Year_2)./(Data.Estimate__Total__AGEBYEDUCATIONALATTAINMENT__Population18To24Ye+Data.Estimate__Total__AGEBYEDUCATIONALATTAINMENT__Population25YearsA);
X_temp=NaN.*zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Est_Ed(t_find);
    end
end
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'EDUCATION_2022'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill in EDUCATION gaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_fit= ~isnan(US_County.EDUCATION_2022) & ~isnan(US_County.EDUCATION_2020) & ~isnan(US_County.EDUCATION_2010);
f_fill_2010= ~isnan(US_County.EDUCATION_2022) & ~isnan(US_County.EDUCATION_2020) & isnan(US_County.EDUCATION_2010);
f_fill_2020= ~isnan(US_County.EDUCATION_2022) & isnan(US_County.EDUCATION_2020) & isnan(US_County.EDUCATION_2010);
f_fill_2022= isnan(US_County.EDUCATION_2022) & ~isnan(US_County.EDUCATION_2020) & ~isnan(US_County.EDUCATION_2010);
% 2010
X=[US_County.EDUCATION_2020 US_County.EDUCATION_2022];
Y=US_County.EDUCATION_2010;
mdl_f=fitlm(X(f_fit,:),Y(f_fit));
US_County.EDUCATION_2010(f_fill_2010)=predict(mdl_f,X(f_fill_2010,:));

% 2020
X=[US_County.EDUCATION_2010 US_County.EDUCATION_2022];
Y=US_County.EDUCATION_2020;
mdl_f=fitlm(X(f_fit,:),Y(f_fit));
US_County.EDUCATION_2020(f_fill_2020)=predict(mdl_f,X(f_fill_2020,:));

% 2022
X=[US_County.EDUCATION_2010 US_County.EDUCATION_2020];
Y=US_County.EDUCATION_2022;
mdl_f=fitlm(X(f_fit,:),Y(f_fit));
US_County.EDUCATION_2022(f_fill_2022)=predict(mdl_f,X(f_fill_2022,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population size: 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Age_Sex\ACSST5Y2010.S0101-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

Under_15=Data.Total__Estimate__AGE__Under5Years+Data.Total__Estimate__AGE__5To9Years+Data.Total__Estimate__AGE__10To14Years;
Age_15_64=Data.Total__Estimate__SELECTEDAGECATEGORIES__15To44Years+Data.Total__Estimate__AGE__45To49Years+Data.Total__Estimate__AGE__50To54Years+Data.Total__Estimate__AGE__55To59Years+Data.Total__Estimate__AGE__60To64Years;
Over_64=Data.Total__Estimate__SELECTEDAGECATEGORIES__65YearsAndOver;
X_temp=NaN.*zeros(height(US_County),4);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc,1)=Data.Total__Estimate__TotalPopulation(t_find);
        X_temp(cc,2)=Under_15(t_find);
        X_temp(cc,3)=Age_15_64(t_find);
        X_temp(cc,4)=Over_64(t_find);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill in gaps for 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X=[US_County.HOUSEHOLD_INCOME_2010 US_County.GINI_2010 US_County.EDUCATION_2010];
for ii=1:4
    Y=X_temp(:,ii);

    t_nan=~isnan(sum(X,2)) & ~isnan(Y);
    t_need=~isnan(sum(X,2)) & isnan(Y);
    mdl_fit=fitlm(X(t_nan,:),Y(t_nan));
    X_temp(t_need,ii)=predict(mdl_fit,X(t_need,:));
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'POPULATION_SIZE_2010','AGE_UNDER_15_2010','AGE_15_64_2010','AGE_65_OLDER_2010'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population size: 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Age_Sex\ACSST5Y2020.S0101-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

Under_15=100.*(Data.Estimate__Total__TotalPopulation__AGE__Under5Years+Data.Estimate__Total__TotalPopulation__AGE__5To9Years+Data.Estimate__Total__TotalPopulation__AGE__10To14Years)./Data.Estimate__Total__TotalPopulation;
Age_15_64=Data.Estimate__Percent__TotalPopulation__SELECTEDAGECATEGORIES__15_1+100.*(Data.Estimate__Total__TotalPopulation__AGE__45To49Years+Data.Estimate__Total__TotalPopulation__AGE__50To54Years+Data.Estimate__Total__TotalPopulation__AGE__55To59Years+Data.Estimate__Total__TotalPopulation__AGE__60To64Years)./Data.Estimate__Total__TotalPopulation;
Over_64=Data.Estimate__Percent__TotalPopulation__SELECTEDAGECATEGORIES__65YearsAndOver;

X_temp=NaN.*zeros(height(US_County),4);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc,1)=Data.Estimate__Total__TotalPopulation(t_find);
        X_temp(cc,2)=Under_15(t_find);
        X_temp(cc,3)=Age_15_64(t_find);
        X_temp(cc,4)=Over_64(t_find);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill in gaps for 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X=[US_County.POPULATION_SIZE_2010 US_County.AGE_UNDER_15_2010 US_County.AGE_15_64_2010 US_County.AGE_65_OLDER_2010];
for ii=1:4
    Y=X_temp(:,ii);

    t_nan=~isnan(sum(X,2)) & ~isnan(Y);
    t_need=~isnan(sum(X,2)) & isnan(Y);
    mdl_fit=fitlm(X(t_nan,:),Y(t_nan));
    X_temp(t_need,ii)=predict(mdl_fit,X(t_need,:));
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'POPULATION_SIZE_2020','AGE_UNDER_15_2020','AGE_15_64_2020','AGE_65_OLDER_2020'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population size: 2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Age_Sex\ACSST5Y2022.S0101-Data.csv']);
temp_ID=Data.Geography;
County_ID=cell(size(temp_ID));
for ii=1:length(temp_ID)
    temp=temp_ID{ii};
    County_ID{ii}=temp(end-4:end);
end

Under_15=100.*(Data.Estimate__Total__TotalPopulation__AGE__Under5Years+Data.Estimate__Total__TotalPopulation__AGE__5To9Years+Data.Estimate__Total__TotalPopulation__AGE__10To14Years)./Data.Estimate__Total__TotalPopulation;
Age_15_64=Data.Estimate__Percent__TotalPopulation__SELECTEDAGECATEGORIES__15_1+100.*(Data.Estimate__Total__TotalPopulation__AGE__45To49Years+Data.Estimate__Total__TotalPopulation__AGE__50To54Years+Data.Estimate__Total__TotalPopulation__AGE__55To59Years+Data.Estimate__Total__TotalPopulation__AGE__60To64Years)./Data.Estimate__Total__TotalPopulation;
Over_64=Data.Estimate__Percent__TotalPopulation__SELECTEDAGECATEGORIES__65YearsAndOver;

X_temp=NaN.*zeros(height(US_County),4);

for cc=1:height(US_County)
    t_find=strcmpi(County_ID,US_County.GEOID{cc});
    if(sum(t_find)>0)
        X_temp(cc,1)=Data.Estimate__Total__TotalPopulation(t_find);
        X_temp(cc,2)=Under_15(t_find);
        X_temp(cc,3)=Age_15_64(t_find);
        X_temp(cc,4)=Over_64(t_find);
    end
end

X=[US_County.POPULATION_SIZE_2010 US_County.POPULATION_SIZE_2020];
Y=X_temp(:,1);
t_nan=~isnan(sum(X,2)) & ~isnan(Y);
t_need=~isnan(sum(X,2)) & isnan(Y);
mdl_fit=fitlm(X(t_nan,:),Y(t_nan));
X_temp(t_need,1)=predict(mdl_fit,X(t_need,:));

X=[US_County.AGE_UNDER_15_2010 US_County.AGE_UNDER_15_2020 US_County.AGE_15_64_2010 US_County.AGE_15_64_2020 US_County.AGE_65_OLDER_2010 US_County.AGE_65_OLDER_2020];
Y=X_temp(:,2);
t_nan=~isnan(sum(X,2)) & ~isnan(Y);
t_need=~isnan(sum(X,2)) & isnan(Y);
mdl_fit=fitlm(X(t_nan,:),Y(t_nan));
X_temp(t_need,2)=predict(mdl_fit,X(t_need,:));

Y=X_temp(:,3);
t_nan=~isnan(sum(X,2)) & ~isnan(Y);
t_need=~isnan(sum(X,2)) & isnan(Y);
mdl_fit=fitlm(X(t_nan,:),Y(t_nan));
X_temp(t_need,3)=predict(mdl_fit,X(t_need,:));

Y=X_temp(:,4);
t_nan=~isnan(sum(X,2)) & ~isnan(Y);
t_need=~isnan(sum(X,2)) & isnan(Y);
mdl_fit=fitlm(X(t_nan,:),Y(t_nan));
X_temp(t_need,4)=predict(mdl_fit,X(t_need,:));

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'POPULATION_SIZE_2022','AGE_UNDER_15_2022','AGE_15_64_2022','AGE_65_OLDER_2022'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rural Population: 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Rural_Population\Rural_Population_County_2020.xlsx']);

X_temp=NaN.*zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(Data.STATE,US_County.STATEFP{cc}) & strcmpi(Data.COUNTY,US_County.COUNTYFP{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.POP_URB(t_find)./Data.POP_COU(t_find);
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'PERCENT_URBAN_2020'};
US_County=[US_County T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rural Population: 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\Rural_Population\Rural_Population_County_2010.xls']);


X_temp=NaN.*zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=strcmpi(Data.STATE,US_County.STATEFP{cc}) & strcmpi(Data.COUNTY,US_County.COUNTYFP{cc});
    if(sum(t_find)>0)
        X_temp(cc)=Data.POP_URBAN(t_find)./Data.POP_COU(t_find);
    end
end
Y=X_temp;
X=[US_County.PERCENT_URBAN_2020 log(US_County.POPULATION_SIZE_2020) log(US_County.POPULATION_SIZE_2010)];

t_nan=~isnan(sum(X,2)) & ~isnan(Y);
t_need=~isnan(sum(X,2)) & isnan(Y);
mdl_fit=fitlm(X(t_nan,:),Y(t_nan));
X_temp(t_need)=predict(mdl_fit,X(t_need,:));

X_temp(X_temp>1)=1;
X_temp(X_temp<0)=0;
T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'PERCENT_URBAN_2010'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connection Population: 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=shaperead([temp_dir '/Shapefile/cb_2021_us_county_500k.shp'],'UseGeoCoords',true);
S=S(Indx);
S=S(~County_remove);

lat_county=zeros(height(US_County),1);
lon_county=zeros(height(US_County),1);
for cc=1:height(US_County)
    t_indx=strcmp({S.GEOID},US_County.GEOID(cc));
    polyin = polyshape(S(t_indx).Lon,S(t_indx).Lat);
    [lon_county(cc),lat_county(cc)] = centroid(polyin);
end

GM=zeros(length(S));
P=US_County.POPULATION_SIZE_2010;
P(isnan(P))=0;
for ii=1:length(S)
    if(P(ii)>0)
        temp_g=P(ii).*P./(deg2sm(distance(lat_county(ii),lon_county(ii),lat_county,lon_county)));
        GM(ii,:)=temp_g;
        GM(ii,ii)=0;
    end
end
X_temp=log(mean(GM,2));
X_temp(isinf(X_temp))=0;

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'CONNECT_POP_2010'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connection Population: 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GM=zeros(length(S));
P=US_County.POPULATION_SIZE_2020;
P(isnan(P))=0;
for ii=1:length(S)
    if(P(ii)>0)
        temp_g=P(ii).*P./(deg2sm(distance(lat_county(ii),lon_county(ii),lat_county,lon_county)));
        GM(ii,:)=temp_g;
        GM(ii,ii)=0;
    end
end
X_temp=log(mean(GM,2));
X_temp(isinf(X_temp))=0;

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'CONNECT_POP_2020'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connection Population: 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GM=zeros(length(S));
P=US_County.POPULATION_SIZE_2022;
P(isnan(P))=0;
for ii=1:length(S)
    if(P(ii)>0)
        temp_g=P(ii).*P./(deg2sm(distance(lat_county(ii),lon_county(ii),lat_county,lon_county)));
        GM(ii,:)=temp_g;
        GM(ii,ii)=0;
    end
end
X_temp=log(mean(GM,2));
X_temp(isinf(X_temp))=0;

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'CONNECT_POP_2022'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COVID-19 Deaths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\COVID-19\AH_County_of_Residence_COVID-19_Deaths_Counts__2020_Provisional.csv']);

X_temp=zeros(height(US_County),1);

for cc=1:height(US_County)
    t_find=Data.FipsCode==str2double(US_County.GEOID{cc});
    
    if(sum(t_find)>0)
        
        if(isnan(Data.COVID_19Deaths(t_find)))
            temp_v=9; % the upper bound of what is supressed;
        else
            temp_v=Data.COVID_19Deaths(t_find);
        end
        X_temp(cc)=temp_v;
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'COVID_DEATHS'};
US_County=[US_County T_temp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COVID STRINGENCY INDEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\COVID-19\OxCGRTUS_timeseries_all.xlsx']);
Avg_STR=mean(table2array(Data(:,6:371)),2);
ST=cell(size(Avg_STR));
for ii=1:length(ST)
    temp=Data.region_code{ii};
    ST{ii}=temp(end-1:end);
end
X_temp=zeros(height(US_County),1);
for cc=1:height(Data)
    t_find=strcmpi(ST{cc},US_County.STUSPS);
    if(sum(t_find)>0)
        X_temp(t_find)=Avg_STR(cc);
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'COVID_STRI'};
US_County=[US_County T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H1N1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\H1N1\h1n1-070609.xlsx'],'Sheet','US_H1N1');
X_temp=zeros(height(US_County),1);

for cc=1:height(Data)
    t_find=strcmpi(Data.state{cc},US_County.STUSPS) & strcmpi(Data.county{cc},US_County.NAME);
    if(sum(t_find)>0 && ~isnan(Data.confirmed(cc)))
        X_temp(t_find)=Data.confirmed(cc)./sum(t_find);
    end
end

T_temp=array2table(X_temp);
T_temp.Properties.VariableNames={'H1N1_CASES'};
US_County=[US_County T_temp];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H5N1 Cases Among Humans: Dairy Connected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\H5N1_Outbreaks\H5N1_State_Human.csv']);
Data=Data(Data.Dairy>0,:);
N_Samp=10^3;
Spillover=zeros(height(US_County),N_Samp);
Spillover(~isnan(US_County.DAIRY_HPAI_REMAIN_STATE),:)=NaN;
Spillover_Remaining_State=zeros(height(US_County),N_Samp);
Spillover_Remaining_State(isnan(US_County.DAIRY_HPAI_REMAIN_STATE),:)=NaN;

for ii=1:height(Data)
    t_state=strcmpi(Data.StateName{ii},US_County.STATE_NAME) & ~isnan(US_County.DAIRY_HPAI_OUTBREAK_KNOWN);
    w_county=US_County.DAIRY_HPAI_OUTBREAK_KNOWN;
    t_state_unknown=strcmpi(Data.StateName{ii},US_County.STATE_NAME) & ~isnan(US_County.DAIRY_HPAI_REMAIN_STATE);
    w_county(~t_state)=0;
    w_county=cumsum(w_county)./(sum(w_county)+sum(unique(US_County.DAIRY_HPAI_REMAIN_STATE(t_state_unknown))));
    for nn=1:Data.Dairy(ii)
        for ss=1:N_Samp
            indx_f=find(rand(1)<=w_county,1);
            if(~isempty(indx_f))
                Spillover(indx_f,ss)=Spillover(indx_f,ss)+1;
            else
                Spillover_Remaining_State(t_state_unknown,ss)=Spillover_Remaining_State(t_state_unknown,ss)+1;
            end
        end
    end
end

Samples=[array2table(Spillover) array2table(Spillover_Remaining_State)];
US_County_Dairy_to_Human=[US_County_Dairy_to_Human Samples];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H5N1 Cases Among Humans: Poultry Connected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data=readtable([pwd '\H5N1_Outbreaks\H5N1_State_Human.csv']);
Data=Data(Data.Poultry>0,:);
N_Samp=10^3;
Spillover=zeros(height(US_County),N_Samp);

for ii=1:height(Data)
    t_state=strcmpi(Data.StateName{ii},US_County.STATE_NAME);
    w_county=US_County.POULTRY_HPAI_OUTBREAK;
    w_county(~t_state)=0;
    w_county=cumsum(w_county)./sum(w_county);
    for nn=1:Data.Poultry(ii)
        for ss=1:N_Samp
            indx_f=find(rand(1)<=w_county,1);
            Spillover(indx_f,ss)=Spillover(indx_f,ss)+1;
        end
    end
end

Samples=array2table(Spillover);
US_County_Poultry_to_Human=[US_County_Poultry_to_Human Samples];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('Data_US_County.mat','US_County','US_County_Dairy_to_Human','US_County_Poultry_to_Human');