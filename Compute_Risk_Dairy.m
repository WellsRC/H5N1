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

onward_transmission_dairy_farm_County=zeros(height(US_County),length(par_est));
onward_transmission_dairy_farm_State=zeros(length(State_Name),length(par_est));

post_outbreak_dairy_farm_State=zeros(length(State_Name),2501,length(par_est));
post_spillover_dairy_farm_State=zeros(length(State_Name),101,length(par_est));

k_onward_transmission=2.69;
R0=0.05;
p_no_onward_transmission=nbinpdf(0,k_onward_transmission,k_onward_transmission./(k_onward_transmission+R0));

for mm=1:height(Dairy_Model)
    [F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_par,logic_temperature_indx]= Dairy_Covariates(Dairy_Model.Model_H5N1{mm},Dairy_Model.Model_Farm{mm},Dairy_Model.Model_Stratified_Operations{mm});
    x=par_est{mm};
    no_farms=County_Farms==0;
    
    indx_pinf=[5:(8+size(P_County,1))];

    beta_x=x(~ismember(1:length(x),[indx_pinf length(x)]));
    temp_beta=x(~ismember(1:length(x),[indx_pinf length(x)]));
    if(length(beta_x)>4)
            beta_x(5:end)=10.^beta_x(5:end);
        if(~isempty(logic_temperature_indx))
            beta_x(logic_temperature_indx+4)=-temp_beta(logic_temperature_indx+4);
        end
    end
    
    beta_p=x(indx_pinf);
    if(length(beta_p)>4)
        beta_p(5:end)=-10.^beta_p(5:end);
    end
    
    kappa_spillover=10.^x(end);
    
    
    if(~isempty(Dairy_Network))
        beta_x_temp=beta_x([1:4 4+find(~logic_connect)]);
        mu_farm_temp = Risk_Assesment_Farms(beta_x_temp,[F_County; X_County(~logic_connect,:)]);
        mu_farm_temp(County_Farms==0)=0;
    
        beta_p_temp=beta_p([1:4 4+find(~logic_connect_p)]);
        p_inf_County_temp=Zero_Inflation(beta_p_temp,[F_County; P_County(~logic_connect_p,:)]);
        p_inf_County_temp(County_Farms==0)=1;
    
        temp_r=(1-p_inf_County_temp(:)).*mu_farm_temp(:);
        X_County(logic_connect,:)=repmat((temp_r(:)'*Dairy_Network),sum(logic_connect),1);
        P_County(logic_connect_p,:)=repmat((temp_r(:)'*Dairy_Network),sum(logic_connect_p),1);
    end
        
    p_inf_County=Zero_Inflation(beta_p,[F_County; P_County]);
    p_inf_County(County_Farms==0)=1;
    
    mu_farm_County = Risk_Assesment_Farms(beta_x,[F_County; X_County]);
    mu_farm_County(County_Farms==0)=0;
    
    mu_farm_State=zeros(size(State_Spillover_Events));
    k_State=zeros(size(State_Spillover_Events));
    k_spill_State=zeros(size(State_Spillover_Events));
    
    p_temp=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:));
    p_temp_spill=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County(:));
    
    r=10.^linspace(-3,4,501);
        
    w_state=zeros(size(State_Spillover_Events));
    County_w_Farms=County_Farms>0;
    
    for ss=1:length(mu_farm_State)
    
        temp_county=(1-p_inf_County(:)).*mu_farm_County(:); 
        mu_farm_State(ss)=state_weight_matrix(ss,:)*temp_county;   
        w_state(ss)=state_weight_matrix(ss,:)*County_w_Farms(:);
        
        p_zero_county=prod(p_temp(state_weight_matrix(ss,:)==1));
        if(p_zero_county==0)
            p_zero_county=10^(-64);
        end
        p_temp_state=r.*log(r./(r+mu_farm_State(ss)));
        [p_temp_state,ia]=unique(p_temp_state);
        rt=r(ia);
        rt=rt(~isinf(p_temp_state) & ~isnan(p_temp_state));
        p_temp_state=p_temp_state(~isinf(p_temp_state) & ~isnan(p_temp_state));
        k_State(ss)=interp1(p_temp_state,rt,log(p_zero_county),"pchip");
    
        p_zero_county=prod(p_temp_spill(state_weight_matrix(ss,:)==1));
        if(p_zero_county==0)
            p_zero_county=10^(-64);
        end
        p_temp_state=r.*log(r./(r+mu_farm_State(ss)));
        [p_temp_state,ia]=unique(p_temp_state);
        rt=r(ia);
        rt=rt(~isinf(p_temp_state) & ~isnan(p_temp_state));
        p_temp_state=p_temp_state(~isinf(p_temp_state) & ~isnan(p_temp_state));
        k_spill_State(ss)=interp1(p_temp_state,rt,log(p_zero_county),"pchip");
        
    end
     
    exposure_dairy_farm_County(:,mm)=(1-p_inf_County(:));

    outbreak_dairy_farm_County(:,mm)=(1-p_inf_County(:)).*mu_farm_County(:);
    outbreak_dairy_farm_State(:,mm)=mu_farm_State;

    spillover_dairy_farm_County(:,mm)=kappa_spillover.*outbreak_dairy_farm_County(:,mm);
    spillover_dairy_farm_State(:,mm)=kappa_spillover.*outbreak_dairy_farm_State(:,mm);

    for ii=0:2500
        post_outbreak_dairy_farm_State(:,1+ii,mm)=nbinpdf(ii,k_State(:),k_State(:)./(k_State(:)+mu_farm_State(:)));
    end
    
    for ii=0:100
        if(ii>0)
            onward_transmission_dairy_farm_County(:,mm)=onward_transmission_dairy_farm_County(:,mm)+(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:)).*(1-p_no_onward_transmission.^ii);
            onward_transmission_dairy_farm_State(:,mm)=onward_transmission_dairy_farm_State(:,mm)+nbinpdf(ii,k_spill_State(:),k_spill_State(:)./(k_spill_State(:)+spillover_dairy_farm_State(:,mm))).*(1-p_no_onward_transmission.^ii);
        end
        post_spillover_dairy_farm_State(:,1+ii,mm)=nbinpdf(ii,k_spill_State(:),k_spill_State(:)./(k_spill_State(:)+spillover_dairy_farm_State(:,mm))); 
    end

    outbreak_risk_dairy_farm_County(:,mm)=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:)));
    t_nan=isnan(outbreak_risk_dairy_farm_County(:,mm));
    outbreak_risk_dairy_farm_County(t_nan,mm)=0;

    outbreak_risk_dairy_farm_State(:,mm)=1-nbinpdf(0,k_State(:),k_State(:)./(k_State(:)+mu_farm_State(:)));
    t_nan=isnan(outbreak_risk_dairy_farm_State(:,mm));
    outbreak_risk_dairy_farm_State(t_nan,mm)=0;


    spillover_risk_dairy_farm_County(:,mm)=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County(:)));
    spillover_risk_dairy_farm_State(:,mm)=1-nbinpdf(0,k_spill_State(:),k_spill_State(:)./(k_spill_State(:)+spillover_dairy_farm_State(:,mm))); 

    outbreak_dairy_farm_County(no_farms,mm)=NaN;
    
    exposure_dairy_farm_County(no_farms,mm)=NaN;

    spillover_dairy_farm_County(no_farms,mm)=NaN;
    
    outbreak_risk_dairy_farm_County(no_farms,mm)=NaN;
    
    spillover_risk_dairy_farm_County(no_farms,mm)=NaN;
end

save('Average_Risk_Dairy.mat','onward_transmission_dairy_farm_County','onward_transmission_dairy_farm_State','post_spillover_dairy_farm_State','post_outbreak_dairy_farm_State','w_AIC','State_Name','exposure_dairy_farm_County','outbreak_dairy_farm_County','outbreak_risk_dairy_farm_County','spillover_dairy_farm_County','spillover_risk_dairy_farm_County','outbreak_dairy_farm_State','outbreak_risk_dairy_farm_State','spillover_dairy_farm_State','spillover_risk_dairy_farm_State');