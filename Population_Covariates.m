function [X_County_H1N1,Y_County_H1N1,Cases_County_H1N1,Population_H1N1,logic_par_H1N1,X_County_COVID,Y_County_COVID,Deaths_County_COVID,Population_COVID,logic_par_COVID] = Population_Covariates(Exposure_Variables,Suceptibility_Variables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


logic_par_H1N1=false(1,12);
logic_par_H1N1(12)=true; % Hyper-paramter
logic_par_H1N1(1)=true; % Constant
logic_par_H1N1(8)=true; % Constant


logic_par_COVID=false(1,13);
logic_par_COVID(13)=true; % Hyper-paramter
logic_par_COVID(1)=true; % Constant and COVID STRI
logic_par_COVID(8:9)=true; % Constant

load([pwd '/Data/Data_US_County.mat'],'US_County');
US_County=US_County(~isnan(US_County.GINI_2022) & ~strcmp(US_County.STUSPS,'PR'),:);

Population_H1N1=US_County.POPULATION_SIZE_2010;
Population_COVID=US_County.POPULATION_SIZE_2020;
Cases_County_H1N1=round(US_County.H1N1_CASES);
Deaths_County_COVID=round(US_County.COVID_DEATHS);


Y_County_H1N1=zeros(length(Exposure_Variables),height(US_County));
Y_County_COVID=zeros(length(Exposure_Variables)+1,height(US_County));
for yy=1:length(Exposure_Variables)
    if(strcmp(Exposure_Variables{yy},'Connectivity'))
        Y_County_H1N1(yy,:)=US_County.CONNECT_POP_2010;
        logic_par_H1N1(9)=true;
        Y_County_COVID(yy+1,:)=US_County.CONNECT_POP_2020;
        logic_par_COVID(10)=true;
    elseif(strcmp(Exposure_Variables{yy},'Population'))
        Y_County_H1N1(yy,:)=log(US_County.POPULATION_SIZE_2010);
        logic_par_H1N1(10)=true;
        Y_County_COVID(yy+1,:)=log(US_County.POPULATION_SIZE_2020);
        logic_par_COVID(11)=true;
    elseif(strcmp(Exposure_Variables{yy},'Urban'))
        temp_p=US_County.PERCENT_URBAN_2010;
        temp_p(temp_p<=0)=10.^(-16);
        temp_p(temp_p>=1)=1-10.^(-16);
        Y_County_H1N1(yy,:)=log(temp_p./(1-temp_p));
        logic_par_H1N1(11)=true;

        temp_p=US_County.PERCENT_URBAN_2020;
        temp_p(temp_p<=0)=10.^(-16);
        temp_p(temp_p>=1)=1-10.^(-16);
        Y_County_COVID(yy+1,:)=log(temp_p./(1-temp_p));
        logic_par_COVID(12)=true;
    end
end
Y_County_COVID(1,:)=US_County.COVID_STRI;

X_County_H1N1=zeros(length(Suceptibility_Variables),height(US_County));
X_County_COVID=zeros(length(Suceptibility_Variables),height(US_County));
for yy=1:length(Suceptibility_Variables)
    if(strcmp(Suceptibility_Variables{yy},'Education'))
        temp_p=US_County.EDUCATION_2010./100;
        temp_p(temp_p==0)=10.^(-16);
        temp_p(temp_p==1)=1-10.^(-16);
        X_County_H1N1(yy,:)=log(temp_p./(1-temp_p));
        logic_par_H1N1(2)=true;

        temp_p=US_County.EDUCATION_2020./100;
        temp_p(temp_p==0)=10.^(-16);
        temp_p(temp_p==1)=1-10.^(-16);
        X_County_COVID(yy,:)=log(temp_p./(1-temp_p));
        logic_par_COVID(2)=true;
    elseif(strcmp(Suceptibility_Variables{yy},'Gini'))
        X_County_H1N1(yy,:)=log(US_County.GINI_2010./(1-US_County.GINI_2010));
        logic_par_H1N1(3)=true;

        X_County_COVID(yy,:)=log(US_County.GINI_2020./(1-US_County.GINI_2020));
        logic_par_COVID(3)=true;
    elseif(strcmp(Suceptibility_Variables{yy},'Income'))
        X_County_H1N1(yy,:)=log(US_County.HOUSEHOLD_INCOME_2010);
        logic_par_H1N1(4)=true;

        tt=US_County.HOUSEHOLD_INCOME_2020;
        tt1=US_County.HOUSEHOLD_INCOME_2010;
        tt2=US_County.HOUSEHOLD_INCOME_2022;
        tt(isnan(tt))=pchip([2010 2022],[tt1(isnan(tt)) tt2(isnan(tt))],2020);
        X_County_COVID(yy,:)=log(tt);
        logic_par_COVID(4)=true;
    elseif(strcmp(Suceptibility_Variables{yy},'Density'))
        X_County_H1N1(yy,:)=log(US_County.POPULATION_SIZE_2010./double(US_County.AREA_LAND));
        logic_par_H1N1(5)=true;

        X_County_COVID(yy,:)=log(US_County.POPULATION_SIZE_2020./double(US_County.AREA_LAND));
        logic_par_COVID(5)=true;
    elseif(strcmp(Suceptibility_Variables{yy},'Eldery'))
        temp_p=US_County.AGE_65_OLDER_2010;
        temp_p(temp_p==0)=10.^(-16);
        temp_p(temp_p==1)=1-10.^(-16);
        X_County_H1N1(yy,:)=log(temp_p./(1-temp_p));
        logic_par_H1N1(6)=true;

        temp_p=US_County.AGE_65_OLDER_2020;
        temp_p(temp_p==0)=10.^(-16);
        temp_p(temp_p==1)=1-10.^(-16);
        X_County_COVID(yy,:)=log(temp_p./(1-temp_p));
        logic_par_COVID(6)=true;
    elseif(strcmp(Suceptibility_Variables{yy},'Under_15'))
        temp_p=US_County.AGE_UNDER_15_2010;
        temp_p(temp_p==0)=10.^(-16);
        temp_p(temp_p==1)=1-10.^(-16);
        X_County_H1N1(yy,:)=log(temp_p./(1-temp_p));
        logic_par_H1N1(7)=true;

        temp_p=US_County.AGE_UNDER_15_2020;
        temp_p(temp_p==0)=10.^(-16);
        temp_p(temp_p==1)=1-10.^(-16);
        X_County_COVID(yy,:)=log(temp_p./(1-temp_p));
        logic_par_COVID(7)=true;
    end
end

end

