function [F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_par] = Dairy_Covariates(H5N1_Variable,Farm_Variables,Stratified_Operations_Variables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logic_par=false(5,18);
logic_connect=false(5,18);

load([pwd '/Data/Data_US_County.mat'],'US_County');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% HPAI Outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
County_Farms=US_County.TOTAL_DAIRY_OPERATIONS;
Affected_County_Farms = US_County.DAIRY_HPAI_OUTBREAK_KNOWN;

State_Names=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Names,"Alaska") | strcmp(State_Names,"District of Columbia");
State_Names=State_Names(~state_remove);

Affected_State_Farms = zeros(size(State_Names));
state_weight_matrix=zeros(length(State_Names),height(US_County));


for ss=1:length(State_Names)
    t_find = strcmp(US_County.STATE_NAME,State_Names{ss});
    Affected_State_Farms(ss)= round(sum(US_County.DAIRY_HPAI_OUTBREAK_UNKNOWN(t_find))+sum(US_County.DAIRY_HPAI_OUTBREAK_KNOWN(t_find))); % Need the state level cases otherwise we will be underestimating state-level outbreaks in our likelihood function
    state_weight_matrix(ss,t_find)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Spillover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dairy_Spillover_State=unique(US_County.STATE_NAME);
state_remove=strcmp(Dairy_Spillover_State,"Alaska") | strcmp(Dairy_Spillover_State,"District of Columbia");
Dairy_Spillover_State=Dairy_Spillover_State(~state_remove);
State_Spillover_Events=zeros(size(Dairy_Spillover_State));

for ss=1:length(Dairy_Spillover_State)
    t_state=strcmp(Dairy_Spillover_State{ss},US_County.STATE_NAME);
    State_Spillover_Events(ss)=round(sum(US_County.SPILLOVER_DAIRY(t_state))); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Flyway
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
F_County=[US_County.ATLANTIC_FLYWAY US_County.MISSISSIPPI_FLYWAY US_County.PACIFIC_FLYWAY US_County.CENTRAL_FLYWAY]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Suseptbility  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if(strcmp(Stratified_Operations_Variables,'All'))
    X_County_temp=zeros(length(Farm_Variables)+7,height(US_County));
elseif(~isempty(Stratified_Operations_Variables))
    X_County_temp=zeros(length(Farm_Variables)+2,height(US_County));
else
    X_County_temp=zeros(length(Farm_Variables),height(US_County));
end

logic_par_t=false(1,18);
Dairy_Network={};
logic_conn_included=false;
for ff=1:length(Farm_Variables)
    if(strcmp(Farm_Variables{ff},'Inventory'))
        X_County_temp(ff,:)=log10(US_County.CATTLE_INVENTORY);
        logic_par_t(1)=true;
    elseif(strcmp(Farm_Variables{ff},'Poultry_Operations'))
        X_County_temp(ff,:)=log10(1+US_County.POULTRY_OPR_w_INVENTORY);
        logic_par_t(2)=true;
    elseif(strcmp(Farm_Variables{ff},'Connectivity'))
       X_County_temp(ff,:)=ones(size(US_County.CATTLE_INVENTORY));
       Dairy_Network=US_County.CONNECT_DAIRY;
       Dairy_Network=Dairy_Network-diag(diag(Dairy_Network));
       logic_par_t(3)=true;
       logic_connect(:,3)=true(size(logic_connect(:,3)));
       logic_conn_included=true;
    end

end

if strcmp(Stratified_Operations_Variables,'Cattle_Inventory_50')
    X_County_temp(length(Farm_Variables)+1,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (1 TO 9 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (10 TO 19 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (20 TO 49 HEAD)");
    X_County_temp(length(Farm_Variables)+2,:)=County_Farms'-X_County_temp(length(Farm_Variables)+1,:);
    logic_par_t(4:5)=true;
elseif strcmp(Stratified_Operations_Variables,'Cattle_Inventory_100')
    X_County_temp(length(Farm_Variables)+1,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (1 TO 9 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (10 TO 19 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (20 TO 49 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (50 TO 99 HEAD)");
    X_County_temp(length(Farm_Variables)+2,:)=County_Farms'-X_County_temp(length(Farm_Variables)+1,:);
    logic_par_t(6:7)=true;
elseif strcmp(Stratified_Operations_Variables,'Cattle_Inventory_200')
    X_County_temp(length(Farm_Variables)+1,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (1 TO 9 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (10 TO 19 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (20 TO 49 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (50 TO 99 HEAD)")+US_County.("CATTLE INVENTORY OF MILK COWS: (100 TO 199 HEAD)");
    X_County_temp(length(Farm_Variables)+2,:)=County_Farms'-X_County_temp(length(Farm_Variables)+1,:);
    logic_par_t(8:9)=true;
elseif strcmp(Stratified_Operations_Variables,'Cattle_Inventory_500')
    X_County_temp(length(Farm_Variables)+2,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (500 OR MORE HEAD)");
    X_County_temp(length(Farm_Variables)+1,:)=County_Farms'-X_County_temp(length(Farm_Variables)+2,:);
    logic_par_t(10:11)=true;
elseif strcmp(Stratified_Operations_Variables,'All')
    X_County_temp(length(Farm_Variables)+1,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (1 TO 9 HEAD)");
    X_County_temp(length(Farm_Variables)+2,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (10 TO 19 HEAD)");
    X_County_temp(length(Farm_Variables)+3,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (20 TO 49 HEAD)");
    X_County_temp(length(Farm_Variables)+4,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (50 TO 99 HEAD)");
    X_County_temp(length(Farm_Variables)+5,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (100 TO 199 HEAD)");
    X_County_temp(length(Farm_Variables)+6,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (200 TO 499 HEAD)");
    X_County_temp(length(Farm_Variables)+7,:)=US_County.("CATTLE INVENTORY OF MILK COWS: (500 OR MORE HEAD)");
    logic_par_t(12:18)=true;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Exposure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
X_County=zeros(length(H5N1_Variable).*size(X_County_temp,1),height(US_County));

logic_exp=false(6,1);
if(~logic_conn_included)
    P_County=zeros(length(H5N1_Variable),height(US_County));
    logic_connect_p=false(length(H5N1_Variable),1);
else
    P_County=zeros(length(H5N1_Variable)+1,height(US_County));
    logic_connect_p=false(length(H5N1_Variable)+1,1);
    logic_connect_p(end)=logic_conn_included;
    logic_exp(6)=true;
end

for yy=1:length(H5N1_Variable)
    if(strcmp(H5N1_Variable{yy},'Light_Intensity'))
        Y_County=US_County.LIGHT_INT';
        logic_par(1,:)=logic_par_t;
        logic_exp(1)=true;
    elseif(strcmp(H5N1_Variable{yy},'Waterfowl_Mallard'))
        Y_County=log10(1+US_County.MALLARD)';
        logic_par(2,:)=logic_par_t;
        logic_exp(2)=true;
    elseif(strcmp(H5N1_Variable{yy},'Waterfowl_Canada_Goose'))
        Y_County=log10(1+US_County.CANADA_GOOSE)';
        logic_par(3,:)=logic_par_t;
        logic_exp(3)=true;
    elseif(strcmp(H5N1_Variable{yy},'Waterfowl_AGW_Teal'))
        Y_County=log10(1+US_County.AGW_TEAL)';
        logic_par(4,:)=logic_par_t;
        logic_exp(4)=true;
     elseif(strcmp(H5N1_Variable{yy},'Waterfowl_N_Pintail'))
        Y_County=log10(1+US_County.NORTH_PINTAIL)';
        logic_par(5,:)=logic_par_t;
        logic_exp(5)=true;
    end
    P_County(yy,:)=Y_County;
    for xx=1:size(X_County_temp,1)
        X_County(yy+(xx-1).*length(H5N1_Variable),:)=X_County_temp(xx,:).*Y_County;
    end
end


logic_connect=logic_connect(logic_par);
logic_par=[logic_exp(:); logic_par(:)];


end

