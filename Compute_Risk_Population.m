clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COVID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd '/Data/Data_US_County.mat'],'US_County');
load('Population_Models_Fit_COVID.mat',"par_est_COVID","w_AIC_COVID");

State_Name=unique(US_County.STATE_NAME);

% Define the variables to loop through
Exposure_Variables_v={'Connectivity','Population','Urban'};
Suceptibility_Variables_v={'Education','Gini','Income','Density','Eldery','Under_15'};
bin_suc=dec2bin([0:2^9-1]',9)=='1';
[~,srt_indx]=sort(sum(bin_suc,2),'ascend');
bin_suc=bin_suc(srt_indx,:);

susceptible_risk_population_County_COVID=zeros(height(US_County),length(par_est_COVID));
susceptible_risk_population_State_COVID=zeros(length(State_Name),length(par_est_COVID));

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

    for ss=1:length(State_Name)
        t_state=strcmp(State_Name{ss},US_County.STATE_NAME);

        c_r=susceptible_risk_population_County_COVID(t_state,mm);
        w_c=US_County.POPULATION_SIZE_2022(t_state);
        t_inc=w_c>0 & ~isnan(c_r);
        c_r=c_r(t_inc);
        w_c=w_c(t_inc);
        susceptible_risk_population_State_COVID(ss,mm)=1-nthroot(prod((1-c_r).^w_c),sum(w_c));
    end
end

avg_susceptible_risk_population_County_COVID=susceptible_risk_population_County_COVID*w_AIC_COVID;
avg_susceptible_risk_population_State_COVID=susceptible_risk_population_State_COVID*w_AIC_COVID;

save('Average_Risk_Population_COVID.mat','avg_susceptible_risk_population_County_COVID','susceptible_risk_population_County_COVID','w_AIC_COVID','avg_susceptible_risk_population_State_COVID','susceptible_risk_population_State_COVID','State_Name');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H1N1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Population_Models_Fit_H1N1.mat',"par_est_H1N1","w_AIC_H1N1");

% Define the variables to loop through
Exposure_Variables_v={'Connectivity','Population','Urban'};
Suceptibility_Variables_v={'Education','Gini','Income','Density','Eldery','Under_15'};
bin_suc=dec2bin([0:2^9-1]',9)=='1';
[~,srt_indx]=sort(sum(bin_suc,2),'ascend');
bin_suc=bin_suc(srt_indx,:);

susceptible_risk_population_County_H1N1=zeros(height(US_County),length(par_est_H1N1));
susceptible_risk_population_State_H1N1=zeros(length(State_Name),length(par_est_H1N1));

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
    [~,~,~,~,~,X_County_H1N1,Y_County_H1N1,Deaths_County_H1N1,Population_H1N1,~] = Population_Covariates(Exposure_Variables,Suceptibility_Variables);

    x=par_est_H1N1{mm};
    beta_x=x(1:(1+size(X_County_H1N1,1)));

    beta_y=x((1+length(beta_x)):(end-1));

    sigma_ln=10.^x(end);
    susceptible_risk_population_County_H1N1(:,mm) = Risk_Assesment_Farms(beta_x,X_County_H1N1);                

    for ss=1:length(State_Name)
        t_state=strcmp(State_Name{ss},US_County.STATE_NAME);

        c_r=susceptible_risk_population_County_H1N1(t_state,mm);
        w_c=US_County.POPULATION_SIZE_2022(t_state);
        t_inc=w_c>0 & ~isnan(c_r);
        c_r=c_r(t_inc);
        w_c=w_c(t_inc);
        susceptible_risk_population_State_H1N1(ss,mm)=1-nthroot(prod((1-c_r).^w_c),sum(w_c));
    end
end

avg_susceptible_risk_population_County_H1N1=susceptible_risk_population_County_H1N1*w_AIC_H1N1;
avg_susceptible_risk_population_State_H1N1=susceptible_risk_population_State_H1N1*w_AIC_H1N1;

save('Average_Risk_Population_H1N1.mat','avg_susceptible_risk_population_County_H1N1','susceptible_risk_population_County_H1N1','w_AIC_H1N1','avg_susceptible_risk_population_State_H1N1','susceptible_risk_population_State_H1N1','State_Name');
