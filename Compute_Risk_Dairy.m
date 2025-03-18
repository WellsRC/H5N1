clear;
load([pwd '/Data/Data_US_County.mat'],'US_County');
load('Dairy_Models_Fit.mat',"par_est","w_AIC","Dairy_Model");

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);

exposure_dairy_farm_County=zeros(height(US_County),length(par_est));

outbreak_dairy_farm_County=zeros(height(US_County),length(par_est));
outbreak_risk_dairy_farm_County=zeros(height(US_County),length(par_est));

spillover_dairy_farm_County=zeros(height(US_County),length(par_est));
spillover_risk_dairy_farm_County=zeros(height(US_County),length(par_est));

outbreak_dairy_farm_State=zeros(length(State_Name),length(par_est));
outbreak_risk_dairy_farm_State=zeros(length(State_Name),length(par_est));
spillover_dairy_farm_State=zeros(length(State_Name),length(par_est));
spillover_risk_dairy_farm_State=zeros(length(State_Name),length(par_est));

post_outbreak_dairy_farm_State=zeros(length(State_Name),2501,length(par_est));
post_spillover_dairy_farm_State=zeros(length(State_Name),101,length(par_est));



for mm=1:height(Dairy_Model)
    [X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_par]= Dairy_Covariates(Dairy_Model.Model_H5N1{mm},Dairy_Model.Model_Farm{mm},Dairy_Model.Model_Stratified_Operations{mm});
    x=par_est{mm};
    no_farms=County_Farms==0;
    
    beta_x=x([1 (1+2+size(P_County,1)):(2+size(P_County,1)+size(X_County,1))]);
    if(length(beta_x)>1)
        beta_x(2:end)=10.^beta_x(2:end);
    end
    
    beta_p=x([2:(2+size(P_County,1))]);
    if(length(beta_x)>1)
        beta_p(2:end)=-10.^beta_p(2:end);
    end
    
    kappa_spillover=10.^x(end);
    
    
    if(~isempty(Dairy_Network))
        beta_x_temp=beta_x([1 1+find(~logic_connect)']);
        mu_farm_temp = Risk_Assesment_Farms(beta_x_temp,X_County(~logic_connect,:));
        mu_farm_temp(County_Farms==0)=0;
    
        beta_p_temp=beta_p([1 1+find(~logic_connect_p)']);
        p_inf_County_temp=Zero_Inflation(beta_p_temp,P_County(~logic_connect_p,:));
        p_inf_County_temp(County_Farms==0)=1;
    
        temp_r=(1-p_inf_County_temp(:)).*mu_farm_temp(:);
        X_County(logic_connect,:)=repmat((temp_r(:)'*Dairy_Network),sum(logic_connect),1);
        P_County(logic_connect_p,:)=repmat((temp_r(:)'*Dairy_Network),sum(logic_connect_p),1);
    end
    
    p_inf_County=Zero_Inflation(beta_p,P_County);
    p_inf_County(County_Farms==0)=1;

    mu_farm_County = Risk_Assesment_Farms(beta_x,X_County);
    mu_farm_County(County_Farms==0)=0;
    
    mu_farm_State=zeros(size(State_Spillover_Events));
    for ss=1:length(mu_farm_State)
        temp_county=(1-p_inf_County(:)).*mu_farm_County(:); 
        mu_farm_State(ss)=state_weight_matrix(ss,:)*temp_county; 
    end
     
    exposure_dairy_farm_County(:,mm)=(1-p_inf_County(:));

    outbreak_dairy_farm_County(:,mm)=(1-p_inf_County(:)).*mu_farm_County(:);
    outbreak_dairy_farm_State(:,mm)=mu_farm_State;

    spillover_dairy_farm_County(:,mm)=kappa_spillover.*outbreak_dairy_farm_County(:,mm);
    spillover_dairy_farm_State(:,mm)=kappa_spillover.*outbreak_dairy_farm_State(:,mm);

    post_outbreak_dairy_farm_State(:,1,mm)=poisspdf(0,mu_farm_State(:));
    for ii=1:2500
        post_outbreak_dairy_farm_State(:,1+ii,mm)=poisspdf(ii,mu_farm_State(:));
    end
    
    post_spillover_dairy_farm_State(:,1,mm)=poisspdf(0,spillover_dairy_farm_State(:,mm));
    for ii=1:100
        post_spillover_dairy_farm_State(:,1+ii,mm)=poisspdf(ii,spillover_dairy_farm_State(:,mm));
    end

    outbreak_risk_dairy_farm_County(:,mm)=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:)));
    t_nan=isnan(outbreak_risk_dairy_farm_County(:,mm));
    outbreak_risk_dairy_farm_County(t_nan,mm)=0;

    outbreak_risk_dairy_farm_State(:,mm)=1-(poisspdf(0,mu_farm_State(:)));
    t_nan=isnan(outbreak_risk_dairy_farm_State(:,mm));
    outbreak_risk_dairy_farm_State(t_nan,mm)=0;


    spillover_risk_dairy_farm_County(:,mm)=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County(:)));
    spillover_risk_dairy_farm_State(:,mm)=1-(p_inf_State(:)+(1-p_inf_State(:)).*poisspdf(0,kappa_spillover.*mu_farm_State(:)));

    outbreak_dairy_farm_County(no_farms,mm)=NaN;
    
    exposure_dairy_farm_County(no_farms,mm)=NaN;

    spillover_dairy_farm_County(no_farms,mm)=NaN;
    
    outbreak_risk_dairy_farm_County(no_farms,mm)=NaN;
    
    spillover_risk_dairy_farm_County(no_farms,mm)=NaN;
end

save('Average_Risk_Dairy.mat','post_spillover_dairy_farm_State','post_outbreak_dairy_farm_State','w_AIC','State_Name','exposure_dairy_farm_County','outbreak_dairy_farm_County','outbreak_risk_dairy_farm_County','spillover_dairy_farm_County','spillover_risk_dairy_farm_County','outbreak_dairy_farm_State','outbreak_risk_dairy_farm_State','spillover_dairy_farm_State','spillover_risk_dairy_farm_State');