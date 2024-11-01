clear;
% parpool(6)
load([pwd '/Data/Data_US_County.mat'],'US_County');
US_County=US_County(~isnan(US_County.GINI_2022) & ~strcmp(US_County.STUSPS,'PR'),:);
load('Population_Models_Fit_COVID.mat',"par_est_COVID","w_AIC_COVID");

% Define the variables to loop through
Exposure_Variables_v={'Connectivity','Population','Urban'};
Suceptibility_Variables_v={'Education','Gini','Income','Density','Eldery','Under_15'};
bin_suc=dec2bin([0:2^9-1]',9)=='1';
[~,srt_indx]=sort(sum(bin_suc,2),'ascend');
bin_suc=bin_suc(srt_indx,:);

susceptible_risk_population_County_COVID=zeros(height(US_County),length(par_est_COVID));

for mm=1:size(bin_suc,1)
    if(sum(bin_suc(mm,1:3))>0)
        Exposure_Variables=Exposure_Variables_v(bin_suc(mm,1:3));
    else
        Exposure_Variables={};
    end
    if(sum(bin_suc(mm,4:9))>0)
        Suceptibility_Variables=Suceptibility_Variables_v(bin_suc(mm,4:9));
    else
        Suceptibility_Variables={};
    end
    % Obtain the covariates and data needed for the model
    [~,~,~,~,~,X_County_COVID,Y_County_COVID,Deaths_County_COVID,Population_COVID,~] = Population_Covariates(Exposure_Variables,Suceptibility_Variables);

    x=par_est_COVID{mm};
    beta_x=x(1:(1+size(X_County_COVID,1)));

    beta_y=x((1+length(beta_x)):(end-1));

    sigma_ln=10.^x(end);
    susceptible_risk_population_County_COVID(:,mm) = Risk_Assesment_Farms(beta_x,X_County_COVID);                
end

avg_susceptible_risk_population_County_COVID=susceptible_risk_population_County_COVID*w_AIC_COVID;

save('Average_Risk_Population_COVID.mat','avg_susceptible_risk_population_County_COVID','susceptible_risk_population_County_COVID','w_AIC_COVID');

load('Population_Models_Fit_H1N1.mat',"par_est_H1N1","w_AIC_H1N1");

% Define the variables to loop through
Exposure_Variables_v={'Connectivity','Population','Urban'};
Suceptibility_Variables_v={'Education','Gini','Income','Density','Eldery','Under_15'};
bin_suc=dec2bin([0:2^9-1]',9)=='1';
[~,srt_indx]=sort(sum(bin_suc,2),'ascend');
bin_suc=bin_suc(srt_indx,:);

susceptible_risk_population_County_H1N1=zeros(height(US_County),length(par_est_H1N1));

for mm=1:size(bin_suc,1)
    if(sum(bin_suc(mm,1:3))>0)
        Exposure_Variables=Exposure_Variables_v(bin_suc(mm,1:3));
    else
        Exposure_Variables={};
    end
    if(sum(bin_suc(mm,4:9))>0)
        Suceptibility_Variables=Suceptibility_Variables_v(bin_suc(mm,4:9));
    else
        Suceptibility_Variables={};
    end
    % Obtain the covariates and data needed for the model
    [X_County_H1N1,Y_County_H1N1,Deaths_County_H1N1,Population_H1N1,~,~,~,~,~,~] = Population_Covariates(Exposure_Variables,Suceptibility_Variables);
    
    x=par_est_H1N1{mm};
    beta_x=x(1:(1+size(X_County_H1N1,1)));

    beta_y=x((1+length(beta_x)):(end-1));
    
    sigma_ln=10.^x(end);
    susceptible_risk_population_County_H1N1(:,mm) = Risk_Assesment_Farms(beta_x,X_County_H1N1);                
end

avg_susceptible_risk_population_County_H1N1=susceptible_risk_population_County_H1N1*w_AIC_H1N1;

save('Average_Risk_Population_H1N1.mat','avg_susceptible_risk_population_County_H1N1','susceptible_risk_population_County_H1N1','w_AIC_H1N1');