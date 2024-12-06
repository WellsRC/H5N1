function [X_County,Y_County,County_Farms,Affected_County_Farms,State_Spillover_Matrix,State_Spillover_Events,Remainaing_Affected_State_Farms,Remainaing_Total_State_Farms,state_weight_hpai_matrix,Dairy_Network,logic_par] = Dairy_Covariates(H5N1_Variable,Farm_Variables,Stratified_Operations_Variables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logic_par=false(1,25);
logic_par(25)=true; % Hyper-paramter
logic_par(1)=true; % Constant
logic_par(20)=true; % Constant

load([pwd '/Data/Data_US_County.mat'],'US_County');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% HPAI Outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
County_Farms=US_County.TOTAL_DAIRY_OPERATIONS;
Affected_County_Farms = US_County.DAIRY_HPAI_OUTBREAK_KNOWN;
Dairy_HPAI_Remain_State_Level=unique(US_County.STATE_NAME(US_County.DAIRY_HPAI_OUTBREAK_UNKNOWN>0));
Remainaing_Affected_State_Farms = zeros(size(Dairy_HPAI_Remain_State_Level));
Remainaing_Total_State_Farms = zeros(size(Dairy_HPAI_Remain_State_Level));
state_weight_hpai_matrix=zeros(length(Dairy_HPAI_Remain_State_Level),height(US_County));

for ss=1:length(Dairy_HPAI_Remain_State_Level)
    t_find = strcmp(US_County.STATE_NAME,Dairy_HPAI_Remain_State_Level{ss});
    Remainaing_Affected_State_Farms(ss)= sum(US_County.DAIRY_HPAI_OUTBREAK_UNKNOWN(t_find));
    Remainaing_Total_State_Farms(ss)= sum(US_County.TOTAL_DAIRY_OPERATIONS(t_find));
    state_weight_hpai_matrix(ss,t_find)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Spillover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dairy_Spillover_State=unique(US_County.STATE_NAME(US_County.SPILLOVER_DAIRY>0));
State_Spillover_Events=zeros(size(Dairy_Spillover_State));
State_Spillover_Matrix=zeros(length(Dairy_Spillover_State),height(US_County));

for ss=1:length(Dairy_Spillover_State)
    t_state=strcmp(Dairy_Spillover_State{ss},US_County.STATE_NAME);
    State_Spillover_Events(ss)=sum(US_County.SPILLOVER_DAIRY(t_state)); 
    State_Spillover_Matrix(ss,t_state)=1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Exposure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


Y_County=zeros(length(H5N1_Variable),height(US_County));
for yy=1:length(H5N1_Variable)
    if(strcmp(H5N1_Variable{yy},'Migratory_Birds_H5N1'))
        Y_County(yy,:)=US_County.HPAI_2022_MIGRATORY_BIRDS+US_County.HPAI_2023_MIGRATORY_BIRDS+US_County.HPAI_2024_MIGRATORY_BIRDS;
        logic_par(21)=true;
    elseif(strcmp(H5N1_Variable{yy},'Migratory_Bird_Density'))
        Y_County(yy,:)=US_County.MIGRATORY_BIRD;
        logic_par(22)=true;
    elseif(strcmp(H5N1_Variable{yy},'Detection_H5N1_Birds'))
        Y_County(yy,:)=min(US_County.HPAI_2022_MIGRATORY_BIRDS+US_County.HPAI_2023_MIGRATORY_BIRDS+US_County.HPAI_2024_MIGRATORY_BIRDS,1);
        logic_par(23)=true;
     elseif(strcmp(H5N1_Variable{yy},'Detection_H5N1_Bird_Density'))
        Y_County(yy,:)=min(US_County.HPAI_2022_MIGRATORY_BIRDS+US_County.HPAI_2023_MIGRATORY_BIRDS+US_County.HPAI_2024_MIGRATORY_BIRDS,1).*US_County.MIGRATORY_BIRD;
        logic_par(24)=true;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Suseptbility  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if(strcmp(Stratified_Operations_Variables,'All'))
    X_County=zeros(length(Farm_Variables)+7,height(US_County));
elseif(~isempty(Stratified_Operations_Variables))
    X_County=zeros(length(Farm_Variables)+2,height(US_County));
else
    X_County=zeros(length(Farm_Variables),height(US_County));
end

Dairy_Network={};
for ff=1:length(Farm_Variables)
    if(strcmp(Farm_Variables{ff},'Inventory'))
        X_County(ff,:)=log(US_County.CATTLE_INVENTORY);
        logic_par(2)=true;
    elseif(strcmp(Farm_Variables{ff},'Poultry_Operations'))
        X_County(ff,:)=US_County.POULTRY_OPR_w_INVENTORY;
        logic_par(3)=true;
    elseif(strcmp(Farm_Variables{ff},'Connectivity'))
       X_County(ff,:)=Inf.*ones(size(US_County.CATTLE_INVENTORY));
       Dairy_Network=US_County.CONNECT_DAIRY;
       Dairy_Network=Dairy_Network-diag(diag(Dairy_Network));
       logic_par(4)=true;
    end

end

if strcmp(Stratified_Operations_Variables,'All')
    X_County(length(Farm_Variables)+1,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (1 TO 9 HEAD)");
    X_County(length(Farm_Variables)+2,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (10 TO 19 HEAD)");
    X_County(length(Farm_Variables)+3,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (20 TO 49 HEAD)");
    X_County(length(Farm_Variables)+4,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (50 TO 99 HEAD)");
    X_County(length(Farm_Variables)+5,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (100 TO 199 HEAD)");
    X_County(length(Farm_Variables)+6,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (200 TO 499 HEAD)");
    X_County(length(Farm_Variables)+7,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (500 OR MORE HEAD)");
    logic_par(13:19)=true;
elseif strcmp(Stratified_Operations_Variables,'Cattle_Inventory_50')
    X_County(length(Farm_Variables)+1,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (1 TO 9 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (10 TO 19 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (20 TO 49 HEAD)");
    X_County(length(Farm_Variables)+2,:)=County_Farms'-X_County(length(Farm_Variables)+1,:);
    logic_par(5:6)=true;
elseif strcmp(Stratified_Operations_Variables,'Cattle_Inventory_100')
    X_County(length(Farm_Variables)+1,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (1 TO 9 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (10 TO 19 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (20 TO 49 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (50 TO 99 HEAD)");
    X_County(length(Farm_Variables)+2,:)=County_Farms'-X_County(length(Farm_Variables)+1,:);
    logic_par(7:8)=true;
elseif strcmp(Stratified_Operations_Variables,'Cattle_Inventory_200')
    X_County(length(Farm_Variables)+1,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (1 TO 9 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (10 TO 19 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (20 TO 49 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (50 TO 99 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (100 TO 199 HEAD)");
    X_County(length(Farm_Variables)+2,:)=County_Farms'-X_County(length(Farm_Variables)+1,:);
    logic_par(9:10)=true;
elseif strcmp(Stratified_Operations_Variables,'Cattle_Inventory_500')
    X_County(length(Farm_Variables)+2,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (500 OR MORE HEAD)");
    X_County(length(Farm_Variables)+1,:)=County_Farms'-X_County(length(Farm_Variables)+2,:);
    logic_par(11:12)=true;
end

end

