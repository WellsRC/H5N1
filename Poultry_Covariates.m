function [X_County,Y_County,County_Farms,Affected_County_Farms,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,Spillover,indx,logic_par] = Poultry_Covariates(H5N1_Variable,Farm_Variables,Stratified_Inventory_Variables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logic_par=false(1,15);
logic_par(15)=true; % Hyper-paramter
logic_par(1)=true; % Constant
logic_par(10)=true; % Constant

load([pwd '/Data/Data_US_County.mat'],'US_County','US_County_Poultry_to_Human');

Spillover=table2array(US_County_Poultry_to_Human(:,7:end));

Affected_County_Farms = US_County.POULTRY_HPAI_OUTBREAK;
County_Farms=US_County.POULTRY_OPR_w_INVENTORY;

Pullet_Farms=US_County.PULLET_OPR_w_INVENTORY;
Layer_Farms=US_County.LAYER_OPR_w_INVENTORY;
Turkey_Farms=US_County.TURKEY_OPR_w_INVENTORY;
Broiler_Farms=US_County.BROILER_OPR_w_INVENTORY;

HPAI_Pullet_Farms=US_County.PULLET_HPAI_OUTBREAK;
HPAI_Layer_Farms=US_County.LAYER_HPAI_OUTBREAK;
HPAI_Turkey_Farms=US_County.TURKEY_HPAI_OUTBREAK;
HPAI_Broiler_Farms=US_County.BROILER_HPAI_OUTBREAK;

indx.Pullet=[];
indx.Layer=[];
indx.Broiler=[];
indx.Turkey=[];

Y_County=zeros(length(H5N1_Variable),height(US_County));
for yy=1:length(H5N1_Variable)
    if(strcmp(H5N1_Variable{yy},'Migratory_Birds_H5N1'))
        Y_County(yy,:)=US_County.HPAI_2022_MIGRATORY_BIRDS+US_County.HPAI_2023_MIGRATORY_BIRDS+US_County.HPAI_2024_MIGRATORY_BIRDS;
        logic_par(11)=true;
    elseif(strcmp(H5N1_Variable{yy},'Migratory_Bird_Density'))
        Y_County(yy,:)=US_County.MIGRATORY_BIRD;
        logic_par(12)=true;
    elseif(strcmp(H5N1_Variable{yy},'Detection_H5N1_Birds'))
        Y_County(yy,:)=min(US_County.HPAI_2022_MIGRATORY_BIRDS+US_County.HPAI_2023_MIGRATORY_BIRDS+US_County.HPAI_2024_MIGRATORY_BIRDS,1);
        logic_par(13)=true;
     elseif(strcmp(H5N1_Variable{yy},'Detection_H5N1_Bird_Density'))
        Y_County(yy,:)=min(US_County.HPAI_2022_MIGRATORY_BIRDS+US_County.HPAI_2023_MIGRATORY_BIRDS+US_County.HPAI_2024_MIGRATORY_BIRDS,1).*US_County.MIGRATORY_BIRD;
        logic_par(14)=true;
    end
end

if(strcmp(Stratified_Inventory_Variables,'All'))
    X_County=zeros(length(Farm_Variables)+3,height(US_County));
elseif(~isempty(Stratified_Inventory_Variables))
    X_County=zeros(length(Farm_Variables)+1,height(US_County));
else
    X_County=zeros(length(Farm_Variables),height(US_County));
end

for ff=1:length(Farm_Variables)
    if(strcmp(Farm_Variables{ff},'Turkey_Operations'))
       X_County(ff,:)=US_County.TURKEY_OPR_w_INVENTORY;
       logic_par(2)=true;
       indx.Turkey=[indx.Turkey ff];
    elseif(strcmp(Farm_Variables{ff},'Broiler_Operations'))
       X_County(ff,:)=US_County.BROILER_OPR_w_INVENTORY;
       logic_par(3)=true;
       indx.Broiler=[indx.Broiler ff];
    elseif(strcmp(Farm_Variables{ff},'Layer_Operations'))
       X_County(ff,:)=US_County.LAYER_OPR_w_INVENTORY;
       logic_par(4)=true;
       indx.Layer=[indx.Layer ff];
    elseif(strcmp(Farm_Variables{ff},'Pullet_Operations'))
       X_County(ff,:)=US_County.PULLET_OPR_w_INVENTORY;
       logic_par(5)=true;
       indx.Pullet=[indx.Pullet ff];
    end
end

if strcmp(Stratified_Inventory_Variables,'All')
    temp_p=US_County.PULLET_INVENTORY;
    temp_p(temp_p==0)=exp(-1);
    X_County(length(Farm_Variables)+1,:)=log(temp_p);
    indx.Pullet=[indx.Pullet length(Farm_Variables)+1];

    temp_p=US_County.LAYER_INVENTORY;
    temp_p(temp_p==0)=exp(-1);
    X_County(length(Farm_Variables)+2,:)=log(temp_p);
    indx.Layer=[indx.Layer length(Farm_Variables)+2];

    temp_p=US_County.BROILER_INVENTORY;
    temp_p(temp_p==0)=exp(-1);
    X_County(length(Farm_Variables)+3,:)=log(temp_p);
    indx.Broiler=[indx.Broiler length(Farm_Variables)+3];

    logic_par(7:9)=true;
elseif strcmp(Stratified_Inventory_Variables,'Total_Inventory')
    temp_p=US_County.BROILER_INVENTORY+US_County.ROOSTER_INVENTORY+US_County.PULLET_INVENTORY+US_County.LAYER_INVENTORY;
    temp_p(temp_p==0)=exp(-1);
    X_County(length(Farm_Variables)+1,:)=log(temp_p);
    indx.Layer=[indx.Layer length(Farm_Variables)+1];
    indx.Broiler=[indx.Broiler length(Farm_Variables)+1];
    indx.Pullet=[indx.Pullet length(Farm_Variables)+1];
    logic_par(6)=true;
end

end

