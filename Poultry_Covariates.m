function [F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_par] = Poultry_Covariates(H5N1_Variable,Farm_Variables,Stratified_Inventory_Variables)
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([pwd '/Data/Data_US_County.mat'],'US_County');


State_Names=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Names,"Alaska") | strcmp(State_Names,"District of Columbia");
State_Names=State_Names(~state_remove);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% HPAI Outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Affected_County_Farms = US_County.POULTRY_HPAI_OUTBREAK;
County_Farms=US_County.POULTRY_OPR_w_INVENTORY;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Spillover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Poultry_Spillover_State=unique(US_County.STATE_NAME);
state_remove=strcmp(Poultry_Spillover_State,"Alaska") | strcmp(Poultry_Spillover_State,"District of Columbia");
Poultry_Spillover_State=Poultry_Spillover_State(~state_remove);

State_Spillover_Events=zeros(size(Poultry_Spillover_State));

for ss=1:length(Poultry_Spillover_State)
    t_state=strcmp(Poultry_Spillover_State{ss},US_County.STATE_NAME);
    State_Spillover_Events(ss)=round(sum(US_County.SPILLOVER_POULTRY(t_state))); 
end

state_weight_matrix=zeros(length(State_Names),height(US_County));


for ss=1:length(State_Names)
    t_find = strcmp(US_County.STATE_NAME,State_Names{ss});
    state_weight_matrix(ss,t_find)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Flyway
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
F_County=[US_County.ATLANTIC_FLYWAY US_County.MISSISSIPPI_FLYWAY US_County.PACIFIC_FLYWAY US_County.CENTRAL_FLYWAY]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Suseptbility  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

if(strcmp(Stratified_Inventory_Variables,'All'))
    X_County=zeros(length(Farm_Variables)+3,height(US_County));
elseif(~isempty(Stratified_Inventory_Variables))
    X_County=zeros(length(Farm_Variables)+1,height(US_County));
else
    X_County=zeros(length(Farm_Variables),height(US_County));
end

logic_par=false(1,8);

for ff=1:length(Farm_Variables)
    if(strcmp(Farm_Variables{ff},'Turkey_Operations'))
       X_County(ff,:)=US_County.TURKEY_OPR_w_INVENTORY;
       logic_par(1)=true;
    elseif(strcmp(Farm_Variables{ff},'Broiler_Operations'))
       X_County(ff,:)=US_County.BROILER_OPR_w_INVENTORY;
       logic_par(2)=true;
    elseif(strcmp(Farm_Variables{ff},'Layer_Operations'))
       X_County(ff,:)=US_County.LAYER_OPR_w_INVENTORY;
       logic_par(3)=true;
    elseif(strcmp(Farm_Variables{ff},'Pullet_Operations'))
       X_County(ff,:)=US_County.PULLET_OPR_w_INVENTORY;
       logic_par(4)=true;
    end
end

if strcmp(Stratified_Inventory_Variables,'All')
    temp_p=US_County.PULLET_INVENTORY;
    X_County(length(Farm_Variables)+1,:)=log10(1+temp_p);

    temp_p=US_County.LAYER_INVENTORY;
    X_County(length(Farm_Variables)+2,:)=log10(1+temp_p);

    temp_p=US_County.BROILER_INVENTORY;
    X_County(length(Farm_Variables)+3,:)=log10(1+temp_p);

    logic_par(6:8)=true;
elseif strcmp(Stratified_Inventory_Variables,'Total_Inventory')
    temp_p=US_County.BROILER_INVENTORY+US_County.ROOSTER_INVENTORY+US_County.PULLET_INVENTORY+US_County.LAYER_INVENTORY;
        X_County(length(Farm_Variables)+1,:)=log10(1+temp_p);
    logic_par(5)=true;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Exposure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% X_County=zeros(length(H5N1_Variable).*size(X_County,1),height(US_County));
P_County=zeros(length(H5N1_Variable),height(US_County));

logic_exp=false(5,1);
for yy=1:length(H5N1_Variable)
    if(strcmp(H5N1_Variable{yy},'Light_Intensity'))
        Y_County=US_County.LIGHT_INT';
        logic_exp(1)=true;
    elseif(strcmp(H5N1_Variable{yy},'Waterfowl_Mallard'))
        Y_County=log10(1+US_County.MALLARD)';
        logic_exp(2)=true;
    elseif(strcmp(H5N1_Variable{yy},'Waterfowl_Canada_Goose'))
        Y_County=log10(1+US_County.CANADA_GOOSE)';
        logic_exp(3)=true;
    elseif(strcmp(H5N1_Variable{yy},'Waterfowl_AGW_Teal'))
        Y_County=log10(1+US_County.AGW_TEAL)';
        logic_exp(4)=true;
     elseif(strcmp(H5N1_Variable{yy},'Waterfowl_N_Pintail'))
        Y_County=log10(1+US_County.NORTH_PINTAIL)';
        logic_exp(5)=true;
    end
    P_County(yy,:)=Y_County;
    % for xx=1:size(X_County_temp,1)
    %     X_County(yy+(xx-1).*length(H5N1_Variable),:)=X_County_temp(xx,:).*Y_County;
    % end
end

logic_par=[logic_exp(:); logic_par(:)];
end

