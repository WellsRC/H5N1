clear;
load([pwd '/Data/Data_US_County.mat'],'US_County');
load('Dairy_Models_Fit.mat',"par_est","w_AIC","Dairy_Model");

State_Name=unique(US_County.STATE_NAME);

phi_farm_County=zeros(height(US_County),length(par_est));
outbreak_risk_dairy_farm_County=zeros(height(US_County),length(par_est));
spillover_risk_dairy_farm_County=zeros(height(US_County),length(par_est));

phi_farm_State=zeros(length(State_Name),length(par_est));
outbreak_risk_dairy_farm_State=zeros(length(State_Name),length(par_est));
spillover_risk_dairy_farm_State=zeros(length(State_Name),length(par_est));

for mm=1:height(Dairy_Model)
    [X_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect,logic_par]= Dairy_Covariates(Dairy_Model.Model_H5N1{mm},Dairy_Model.Model_Farm{mm},Dairy_Model.Model_Stratified_Operations{mm});
    x=par_est{mm};
    
    beta_x=x(1:(1+size(X_County,1)));
    if(length(beta_x)>1)
        beta_x(2:end)=10.^beta_x(2:end);
    end
    k_outbreak=10.^x(end-1);
    k_spillover=10.^x(end);
    
    
    if(~isempty(Dairy_Network))
        beta_x_temp=beta_x([1 1+find(~logic_connect)']);
        r_farm_temp = Risk_Assesment_Farms(beta_x_temp,X_County(~logic_connect,:));
        r_farm_temp(County_Farms==0)=0;
        X_County(logic_connect,:)=X_County(logic_connect,:).*repmat((r_farm_temp(:)'*Dairy_Network),sum(logic_connect),1);
    end
        
    phi_farm_County(:,mm) = Risk_Assesment_Farms(beta_x,X_County);  
    outbreak_risk_dairy_farm_County(:,mm) = 1-nbinpdf(0,k_outbreak,1-phi_farm_County(:,mm));             
    spillover_risk_dairy_farm_County(:,mm)=1-nbinpdf(0,k_spillover,1-outbreak_risk_dairy_farm_County(:,mm));
    
    phi_farm_County(County_Farms==0,mm)=NaN;
    outbreak_risk_dairy_farm_County(County_Farms==0,mm)=NaN;
    spillover_risk_dairy_farm_County(County_Farms==0,mm)=NaN;


    
    for ss=1:length(State_Name)
        t_state=strcmp(State_Name{ss},US_County.STATE_NAME);

        c_r=phi_farm_County(t_state,mm);
        w_c=US_County.TOTAL_DAIRY_OPERATIONS(t_state);
        t_inc=w_c>0 & ~isnan(c_r);
        c_r=c_r(t_inc);
        phi_farm_State(ss,mm)=1-prod((1-c_r));
    end
    
    outbreak_risk_dairy_farm_State(:,mm) = 1-nbinpdf(0,k_outbreak,1-phi_farm_State(:,mm));             
    spillover_risk_dairy_farm_State(:,mm)=1-nbinpdf(0,k_spillover,1-outbreak_risk_dairy_farm_State(:,mm));
    
    phi_farm_County(County_Farms==0,mm)=0;
    outbreak_risk_dairy_farm_County(County_Farms==0,mm)=0;
    spillover_risk_dairy_farm_County(County_Farms==0,mm)=0;
end

avg_phi_farm_County=phi_farm_County*w_AIC;
avg_outbreak_risk_dairy_farm_County=outbreak_risk_dairy_farm_County*w_AIC;
avg_spillover_risk_dairy_farm_County=spillover_risk_dairy_farm_County*w_AIC;

avg_phi_farm_State=phi_farm_State*w_AIC;
avg_outbreak_risk_dairy_farm_State=outbreak_risk_dairy_farm_State*w_AIC;
avg_spillover_risk_dairy_farm_State=spillover_risk_dairy_farm_State*w_AIC;

avg_phi_farm_County(County_Farms==0)=NaN;
avg_outbreak_risk_dairy_farm_County(County_Farms==0)=NaN;
avg_spillover_risk_dairy_farm_County(County_Farms==0)=NaN;

save('Average_Risk_Dairy.mat','State_Name','w_AIC','phi_farm_County','phi_farm_State','outbreak_risk_dairy_farm_County','outbreak_risk_dairy_farm_State','spillover_risk_dairy_farm_County','spillover_risk_dairy_farm_State','avg_phi_farm_County','avg_phi_farm_State','avg_outbreak_risk_dairy_farm_County','avg_outbreak_risk_dairy_farm_State','avg_spillover_risk_dairy_farm_County','avg_spillover_risk_dairy_farm_State');