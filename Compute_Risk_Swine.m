clear;

load([pwd '/Data/Spillover_Exacerbation.mat'],'US_County_Spillover');
load('Dairy_Models_Fit.mat',"par_est","w_AIC","Dairy_Model");

overall_risk_swine_farm_County=zeros(height(US_County_Spillover),length(par_est));
exposure_risk_swine_farm_County=zeros(height(US_County_Spillover),length(par_est));
susceptible_risk_swine_farm_County=zeros(height(US_County_Spillover),length(par_est));
spillover_risk_swine_farm_County=zeros(height(US_County_Spillover),length(par_est));

for mm=1:height(Dairy_Model)
    [X_County,Y_County,logic_par] = Swine_Covariates(Dairy_Model.Model_H5N1{mm});
    x=par_est{mm};
    x=x(1:end-1);
    beta_y=x((length(x)-(sum(logic_par)-1)):length(x));
    if(length(beta_y)>1)
        beta_y(2:end)=10.^beta_y(2:end);
    end
    exposure_risk_swine_farm_County(:,mm) = Risk_Assesment_Farms(beta_y,Y_County);
    exposure_risk_swine_farm_County(isinf(X_County),mm)=0;

    temp_s= X_County;                
    temp_s(isinf(X_County))=NaN;

    susceptible_risk_swine_farm_County(:,mm)=(temp_s-min(temp_s))./(max(temp_s)-min(temp_s));
    susceptible_risk_swine_farm_County(isnan(temp_s),mm)=0;

    overall_risk_swine_farm_County(:,mm) = Risk_Assesment_Farms(beta_y,Y_County).*(susceptible_risk_swine_farm_County(:,mm)');
    overall_risk_swine_farm_County(isinf(X_County),mm)=0;
end

avg_overall_risk_swine_farm_County=overall_risk_swine_farm_County*w_AIC;
avg_exposure_risk_swine_farm_County=exposure_risk_swine_farm_County*w_AIC;
avg_susceptible_risk_swine_farm_County=susceptible_risk_swine_farm_County*w_AIC;
avg_spillover_risk_swine_farm_County=spillover_risk_swine_farm_County*w_AIC;

load('Poultry_Models_Fit.mat',"par_est","w_AIC","Poultry_Model");

overall_risk_swine_farm_County=zeros(height(US_County_Spillover),length(par_est));
exposure_risk_swine_farm_County=zeros(height(US_County_Spillover),length(par_est));
susceptible_risk_swine_farm_County=zeros(height(US_County_Spillover),length(par_est));
spillover_risk_swine_farm_County=zeros(height(US_County_Spillover),length(par_est));

for mm=1:height(Dairy_Model)
    [X_County,Y_County,logic_par] = Swine_Covariates(Dairy_Model.Model_H5N1{mm});
    x=par_est{mm};
    x=x(1:end-1);
    beta_y=x((length(x)-(sum(logic_par)-1)):length(x));
    if(length(beta_y)>1)
        beta_y(2:end)=10.^beta_y(2:end);
    end
    exposure_risk_swine_farm_County(:,mm) = Risk_Assesment_Farms(beta_y,Y_County);
    exposure_risk_swine_farm_County(isinf(X_County),mm)=0;

    temp_s= X_County;                
    temp_s(isinf(X_County))=NaN;

    susceptible_risk_swine_farm_County(:,mm)=(temp_s-min(temp_s))./(max(temp_s)-min(temp_s));
    susceptible_risk_swine_farm_County(isnan(temp_s),mm)=0;

    overall_risk_swine_farm_County(:,mm) = Risk_Assesment_Farms(beta_y,Y_County).*(susceptible_risk_swine_farm_County(:,mm)');
    overall_risk_swine_farm_County(isinf(X_County),mm)=0;
end

avg_overall_risk_swine_farm_County=(avg_overall_risk_swine_farm_County+overall_risk_swine_farm_County*w_AIC)./2;
avg_exposure_risk_swine_farm_County=(avg_exposure_risk_swine_farm_County+exposure_risk_swine_farm_County*w_AIC)./2;
avg_susceptible_risk_swine_farm_County=(avg_susceptible_risk_swine_farm_County+susceptible_risk_swine_farm_County*w_AIC)./2;
avg_spillover_risk_swine_farm_County=(avg_spillover_risk_swine_farm_County+spillover_risk_swine_farm_County*w_AIC)./2;

avg_overall_risk_swine_farm_County(isinf(X_County))=NaN;
avg_exposure_risk_swine_farm_County(isinf(X_County))=NaN;
avg_susceptible_risk_swine_farm_County(isinf(X_County))=NaN;
avg_spillover_risk_swine_farm_County(isinf(X_County))=NaN;

save('Average_Risk_Swine.mat','avg_overall_risk_swine_farm_County','avg_exposure_risk_swine_farm_County','avg_susceptible_risk_swine_farm_County');