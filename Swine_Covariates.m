function [X_County,Y_County,logic_par] = Swine_Covariates(H5N1_Variable)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logic_par=false(1,25);
logic_par(25)=false; % Hyper-paramter NOT NEEDED
logic_par(1)=false; % Constant
logic_par(20)=true; % Constant

load([pwd '/Data/Spillover_Exacerbation.mat'],'US_County_Spillover');

Y_County=zeros(length(H5N1_Variable),height(US_County_Spillover));
for yy=1:length(H5N1_Variable)
    if(strcmp(H5N1_Variable{yy},'Migratory_Birds_H5N1'))
        Y_County(yy,:)=US_County_Spillover.HPAI_2022_MIGRATORY_BIRDS+US_County_Spillover.HPAI_2023_MIGRATORY_BIRDS+US_County_Spillover.HPAI_2024_MIGRATORY_BIRDS;
        logic_par(21)=true;
    elseif(strcmp(H5N1_Variable{yy},'Migratory_Bird_Density'))
        Y_County(yy,:)=US_County_Spillover.MIGRATORY_BIRD;
        logic_par(22)=true;
    elseif(strcmp(H5N1_Variable{yy},'Detection_H5N1_Birds'))
        Y_County(yy,:)=min(US_County_Spillover.HPAI_2022_MIGRATORY_BIRDS+US_County_Spillover.HPAI_2023_MIGRATORY_BIRDS+US_County_Spillover.HPAI_2024_MIGRATORY_BIRDS,1);
        logic_par(23)=true;
     elseif(strcmp(H5N1_Variable{yy},'Detection_H5N1_Bird_Density'))
        Y_County(yy,:)=min(US_County_Spillover.HPAI_2022_MIGRATORY_BIRDS+US_County_Spillover.HPAI_2023_MIGRATORY_BIRDS+US_County_Spillover.HPAI_2024_MIGRATORY_BIRDS,1).*US_County_Spillover.MIGRATORY_BIRD;
        logic_par(24)=true;
    end
end

X_County=log(US_County_Spillover.HOG_INVENTORY);
end

