clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');
load('Poultry_Models_Fit.mat',"par_est","w_AIC",'Poultry_Model');

overall_risk_poultry_farm_County=zeros(height(US_County),length(par_est));
exposure_risk_poultry_farm_County=zeros(height(US_County),length(par_est));
susceptible_risk_poultry_farm_County=zeros(height(US_County),length(par_est));
spillover_risk_poultry_farm_County=zeros(height(US_County),length(par_est));

for mm=1:height(Poultry_Model)
    
    [X_County,Y_County,County_Farms,Affected_County_Farms] = Poultry_Covariates(Poultry_Model.Model_H5N1{mm},Poultry_Model.Model_Farm{mm},Poultry_Model.Model_Stratified_Chicken_Inventory{mm});
    x=par_est{mm};
    beta_x=x(1:(1+size(X_County,1)));
    if(length(beta_x)>1)
        beta_x(2:end)=10.^beta_x(2:end);
    end
    beta_y=x((1+length(beta_x)):(end-1));
    if(length(beta_y)>1)
        beta_y(2:end)=10.^beta_y(2:end);
    end
    
    r_nbin=10.^x(end);
    
    overall_risk_poultry_farm_County(:,mm) = Risk_Assesment_Farms(beta_y,Y_County).*Risk_Assesment_Farms(beta_x,X_County);                
    overall_risk_poultry_farm_County(County_Farms==0,mm)=0;

    exposure_risk_poultry_farm_County(:,mm) = Risk_Assesment_Farms(beta_y,Y_County);
    exposure_risk_poultry_farm_County(County_Farms==0,mm)=0;

    susceptible_risk_poultry_farm_County(:,mm) = Risk_Assesment_Farms(beta_x,X_County);                
    susceptible_risk_poultry_farm_County(County_Farms==0,mm)=0;

    spillover_risk_poultry_farm_County(:,mm)=1-nbincdf(0,r_nbin,1-overall_risk_poultry_farm_County(:,mm));
    spillover_risk_poultry_farm_County(County_Farms==0,mm)=0;
end

avg_overall_risk_poultry_farm_County=overall_risk_poultry_farm_County*w_AIC;
avg_exposure_risk_poultry_farm_County=exposure_risk_poultry_farm_County*w_AIC;
avg_susceptible_risk_poultry_farm_County=susceptible_risk_poultry_farm_County*w_AIC;
avg_spillover_risk_poultry_farm_County=spillover_risk_poultry_farm_County*w_AIC;

avg_overall_risk_poultry_farm_County(County_Farms==0)=NaN;
avg_exposure_risk_poultry_farm_County(County_Farms==0)=NaN;
avg_susceptible_risk_poultry_farm_County(County_Farms==0)=NaN;
avg_spillover_risk_poultry_farm_County(County_Farms==0)=NaN;

save('Average_Risk_Poultry.mat','avg_overall_risk_poultry_farm_County','avg_exposure_risk_poultry_farm_County','avg_susceptible_risk_poultry_farm_County','overall_risk_poultry_farm_County','exposure_risk_poultry_farm_County','susceptible_risk_poultry_farm_County','w_AIC','avg_spillover_risk_poultry_farm_County','spillover_risk_poultry_farm_County');