function [X_County_H5N1,Y_County_H5N1,Population_H5N1,logic_par_H5N1] = Population_Covariates_H5N1(Exposure_Variables,Suceptibility_Variables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


logic_par_H5N1=false(1,12);
logic_par_H5N1(12)=true; % Hyper-paramter
logic_par_H5N1(1)=true; % Constant
logic_par_H5N1(8)=true; % Constant

load([pwd '/Data/Data_US_County.mat'],'US_County');
US_County=US_County(~isnan(US_County.GINI_2022) & ~strcmp(US_County.STUSPS,'PR'),:);

Population_H5N1=US_County.POPULATION_SIZE_2022;


Y_County_H5N1=zeros(length(Exposure_Variables),height(US_County));
for yy=1:length(Exposure_Variables)
    if(strcmp(Exposure_Variables{yy},'Connectivity'))
        Y_County_H5N1(yy,:)=US_County.CONNECT_POP_2022;
        logic_par_H5N1(9)=true;
    elseif(strcmp(Exposure_Variables{yy},'Population'))
        Y_County_H5N1(yy,:)=log(US_County.POPULATION_SIZE_2022);
        logic_par_H5N1(10)=true;
    elseif(strcmp(Exposure_Variables{yy},'Urban'))
        temp_p=US_County.PERCENT_URBAN_2020;
        temp_p(temp_p<=0)=10.^(-16);
        temp_p(temp_p>=1)=1-10.^(-16);
        Y_County_H5N1(yy,:)=log(temp_p./(1-temp_p));
        logic_par_H5N1(11)=true;
    end
end

X_County_H5N1=zeros(length(Suceptibility_Variables),height(US_County));
for yy=1:length(Suceptibility_Variables)
    if(strcmp(Suceptibility_Variables{yy},'Education'))
        temp_p=US_County.EDUCATION_2022./100;
        temp_p(temp_p==0)=10.^(-16);
        temp_p(temp_p==1)=1-10.^(-16);
        X_County_H5N1(yy,:)=log(temp_p./(1-temp_p));
        logic_par_H5N1(2)=true;
    elseif(strcmp(Suceptibility_Variables{yy},'Gini'))
        X_County_H5N1(yy,:)=log(US_County.GINI_2022./(1-US_County.GINI_2022));
        logic_par_H5N1(3)=true;
    elseif(strcmp(Suceptibility_Variables{yy},'Income'))
        X_County_H5N1(yy,:)=log(US_County.HOUSEHOLD_INCOME_2022);
        logic_par_H5N1(4)=true;
    elseif(strcmp(Suceptibility_Variables{yy},'Density'))
        X_County_H5N1(yy,:)=log(US_County.POPULATION_SIZE_2022./double(US_County.AREA_LAND));
        logic_par_H5N1(5)=true;
    elseif(strcmp(Suceptibility_Variables{yy},'Eldery'))
        temp_p=US_County.AGE_65_OLDER_2022;
        temp_p(temp_p==0)=10.^(-16);
        temp_p(temp_p==1)=1-10.^(-16);
        X_County_H5N1(yy,:)=log(temp_p./(1-temp_p));
        logic_par_H5N1(6)=true;
    elseif(strcmp(Suceptibility_Variables{yy},'Under_15'))
        temp_p=US_County.AGE_UNDER_15_2022;
        temp_p(temp_p==0)=10.^(-16);
        temp_p(temp_p==1)=1-10.^(-16);
        X_County_H5N1(yy,:)=log(temp_p./(1-temp_p));
        logic_par_H5N1(7)=true;
    end
end

end

