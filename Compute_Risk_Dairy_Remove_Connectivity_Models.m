clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');
load('Dairy_Models_Fit.mat',"par_est","w_AIC","Dairy_Model");

overall_risk_dairy_farm_County=zeros(height(US_County),length(par_est));
exposure_risk_dairy_farm_County=zeros(height(US_County),length(par_est));
susceptible_risk_dairy_farm_County=zeros(height(US_County),length(par_est));
spillover_risk_dairy_farm_County=zeros(height(US_County),length(par_est));

Model_Keep=true(size(w_AIC));

for mm=1:height(Dairy_Model)
    [X_County,Y_County,County_Farms,Affected_County_Farms,~,~,~,~,~,Dairy_Network,logic_par] = Dairy_Covariates(Dairy_Model.Model_H5N1{mm},Dairy_Model.Model_Farm{mm},Dairy_Model.Model_Stratified_Operations{mm});
    if(logic_par(4))
        Model_Keep(mm)=false(1);
    else
        x=par_est{mm};
        beta_x=x(1:(1+size(X_County,1)));
        if(length(beta_x)>1)
            beta_x(2:end)=10.^beta_x(2:end);        
        end
        beta_y=x((1+length(beta_x)):(end-1));
        if(length(beta_y)>1)
            beta_y(2:end)=10.^beta_y(2:end);
        end
        
        t_connect=sum(isinf(X_County),2)==size(X_County,2);
        if(~isempty(Dairy_Network))
            beta_x_temp=beta_x([1 1+find(~t_connect)']);
            r_farm_temp = Risk_Assesment_Farms(beta_y,Y_County).*Risk_Assesment_Farms(beta_x_temp,X_County(~t_connect,:));
            r_farm_temp(County_Farms==0)=0;
            X_County(t_connect,:)=(r_farm_temp(:)'*Dairy_Network);
        end
        
        r_nbin=10.^x(end);
        
        overall_risk_dairy_farm_County(:,mm) = Risk_Assesment_Farms(beta_y,Y_County).*Risk_Assesment_Farms(beta_x,X_County);                
        overall_risk_dairy_farm_County(County_Farms==0,mm)=0;
    
        exposure_risk_dairy_farm_County(:,mm) = Risk_Assesment_Farms(beta_y,Y_County);
        exposure_risk_dairy_farm_County(County_Farms==0,mm)=0;
    
        susceptible_risk_dairy_farm_County(:,mm) = Risk_Assesment_Farms(beta_x,X_County);                
        susceptible_risk_dairy_farm_County(County_Farms==0,mm)=0;
    
        spillover_risk_dairy_farm_County(:,mm)=1-nbincdf(0,r_nbin,1-overall_risk_dairy_farm_County(:,mm));
        spillover_risk_dairy_farm_County(County_Farms==0,mm)=0;
    end
end

overall_risk_dairy_farm_County=overall_risk_dairy_farm_County(:,Model_Keep);
spillover_risk_dairy_farm_County=spillover_risk_dairy_farm_County(:,Model_Keep);
susceptible_risk_dairy_farm_County=susceptible_risk_dairy_farm_County(:,Model_Keep);
exposure_risk_dairy_farm_County=exposure_risk_dairy_farm_County(:,Model_Keep);

w_AIC=w_AIC(Model_Keep);
w_AIC=w_AIC./sum(w_AIC);

avg_overall_risk_dairy_farm_County=overall_risk_dairy_farm_County*w_AIC;
avg_exposure_risk_dairy_farm_County=exposure_risk_dairy_farm_County*w_AIC;
avg_susceptible_risk_dairy_farm_County=susceptible_risk_dairy_farm_County*w_AIC;
avg_spillover_risk_dairy_farm_County=spillover_risk_dairy_farm_County*w_AIC;

avg_overall_risk_dairy_farm_County(County_Farms==0)=NaN;
avg_exposure_risk_dairy_farm_County(County_Farms==0)=NaN;
avg_susceptible_risk_dairy_farm_County(County_Farms==0)=NaN;
avg_spillover_risk_dairy_farm_County(County_Farms==0)=NaN;

save('Average_Risk_Dairy_Remove_Connectivity_Models.mat','avg_overall_risk_dairy_farm_County','avg_exposure_risk_dairy_farm_County','avg_susceptible_risk_dairy_farm_County','overall_risk_dairy_farm_County','exposure_risk_dairy_farm_County','susceptible_risk_dairy_farm_County','w_AIC','avg_spillover_risk_dairy_farm_County','spillover_risk_dairy_farm_County');