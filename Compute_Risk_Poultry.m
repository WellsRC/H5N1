clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');
load('Poultry_Models_Refined_Fit.mat',"par_est","w_AIC",'Poultry_Model');

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);

potential_outbreak_poultry_farm_County=zeros(height(US_County),length(par_est));

outbreak_poultry_farm_County=zeros(height(US_County),length(par_est));
outbreak_risk_poultry_farm_County=zeros(height(US_County),length(par_est));

spillover_poultry_farm_County=zeros(height(US_County),length(par_est));
spillover_risk_poultry_farm_County=zeros(height(US_County),length(par_est));

outbreak_poultry_farm_State=zeros(length(State_Name),length(par_est));
outbreak_risk_poultry_farm_State=zeros(length(State_Name),length(par_est));
spillover_poultry_farm_State=zeros(length(State_Name),length(par_est));
spillover_risk_poultry_farm_State=zeros(length(State_Name),length(par_est));

onward_transmission_poultry_farm_County=zeros(height(US_County),length(par_est));
onward_transmission_poultry_farm_State=zeros(length(State_Name),length(par_est));

post_outbreak_poultry_farm_State=zeros(length(State_Name),1001,length(par_est));
post_spillover_poultry_farm_State=zeros(length(State_Name),76,length(par_est));

post_spillover_poultry_farm_County=zeros(height(US_County),76,length(par_est));
post_outbreak_poultry_farm_County_CI=zeros(height(US_County),5,length(par_est));
post_spillover_poultry_farm_County_CI=zeros(height(US_County),5,length(par_est));

par_spillover=zeros(1,length(par_est));

k_onward_transmission=2.69;
R0=0.05;
p_no_onward_transmission=nbinpdf(0,k_onward_transmission,k_onward_transmission./(k_onward_transmission+R0));

for mm=1:length(par_est)
    
    [F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_par,logic_temperature] = Poultry_Covariates(Poultry_Model.Model_H5N1{mm},Poultry_Model.Model_Farm{mm});

    x=par_est{mm};  
    no_farms=County_Farms==0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    indx_pinf=[5:(8+size(P_County,1))];

    beta_x=x(~ismember(1:length(x),[indx_pinf length(x)]));
    if(length(beta_x)>4)
        beta_x(5:end)=10.^beta_x(5:end);    
    end
    
    beta_p=x(indx_pinf);
    if(length(beta_p)>4)
        if(logic_temperature)
            beta_p(5:end-1)=-10.^beta_p(5:end-1);    
        else
            beta_p(5:end)=-10.^beta_p(5:end);
        end
    end
    
    kappa_spillover=10.^x(end);
    par_spillover(mm)=kappa_spillover;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Probability of extra zeros and avg for poisson
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

    p_inf_County=Zero_Inflation(beta_p,[F_County; P_County]);
    p_inf_County(County_Farms==0)=1;
    
    mu_farm_County = Risk_Assesment_Farms(beta_x,[F_County; X_County]);
    
    mu_farm_County(County_Farms==0)=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Probability of extra zeros and avg State-Level
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    mu_farm_State=zeros(size(State_Spillover_Events));
    
    k_state=zeros(size(State_Spillover_Events));
    p_nb_state=zeros(size(State_Spillover_Events));
    z_state=zeros(size(State_Spillover_Events));
    
    k_state_spill=zeros(size(State_Spillover_Events));
    z_state_spill=zeros(size(State_Spillover_Events));
    p_nb_state_spill=zeros(size(State_Spillover_Events));
    
    Affected_State_Farms=zeros(size(State_Spillover_Events));
    
    p_temp=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:));
    p_temp_spill=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County(:));
    
    w_state=zeros(size(State_Spillover_Events));
    County_w_Farms=County_Farms>0;
    
    L_Nan=false;
    for ss=1:length(w_state) 
        mu_farm_State(ss)=state_weight_matrix(ss,:)*((1-p_inf_County(:)).*mu_farm_County(:)); 
        w_state(ss)=state_weight_matrix(ss,:)*County_w_Farms(:);
    
        mu_County_NB=state_weight_matrix(ss,:)*((1-p_inf_County(:)).*mu_farm_County(:)); 
        var_County_NB=state_weight_matrix(ss,:)*((1-p_inf_County(:)).*(mu_farm_County(:)+p_inf_County(:).*mu_farm_County(:).^2)); 
        p_zero_county=prod(p_temp(state_weight_matrix(ss,:)==1));
        if(p_zero_county==0)
            p_zero_county=10^(-64);
        end
        
        [z_state(ss),k_state(ss),p_nb_state(ss),L_Nan1]=Estimate_ZIF_NB(mu_County_NB,var_County_NB,p_zero_county);
    
    
        mu_County_NB=state_weight_matrix(ss,:)*((1-p_inf_County(:)).*kappa_spillover.*mu_farm_County(:)); 
        var_County_NB=state_weight_matrix(ss,:)*((1-p_inf_County(:)).*(kappa_spillover.*mu_farm_County(:)+p_inf_County(:).*(kappa_spillover.*mu_farm_County(:)).^2)); 
        p_zero_county=prod(p_temp_spill(state_weight_matrix(ss,:)==1));
        if(p_zero_county==0)
            p_zero_county=10^(-64);
        end
        [z_state_spill(ss),k_state_spill(ss),p_nb_state_spill(ss),L_Nan2]=Estimate_ZIF_NB(mu_County_NB,var_County_NB,p_zero_county);
    end

    potential_outbreak_poultry_farm_County(:,mm)=mu_farm_County(:);
    outbreak_poultry_farm_County(:,mm)=(1-p_inf_County(:)).*mu_farm_County(:);

    outbreak_poultry_farm_State(:,mm)=mu_farm_State(:);

    spillover_poultry_farm_County(:,mm)=kappa_spillover.*outbreak_poultry_farm_County(:,mm);
    spillover_poultry_farm_State(:,mm)=kappa_spillover.*outbreak_poultry_farm_State(:,mm);

    post_outbreak_poultry_farm_County_CI(:,1,mm)=poissinv((0.025-p_inf_County(:))./(1-p_inf_County(:)),mu_farm_County(:));
    post_outbreak_poultry_farm_County_CI(:,2,mm)=poissinv((0.25-p_inf_County(:))./(1-p_inf_County(:)),mu_farm_County(:));
    post_outbreak_poultry_farm_County_CI(:,3,mm)=poissinv((0.5-p_inf_County(:))./(1-p_inf_County(:)),mu_farm_County(:));
    post_outbreak_poultry_farm_County_CI(:,4,mm)=poissinv((0.75-p_inf_County(:))./(1-p_inf_County(:)),mu_farm_County(:));
    post_outbreak_poultry_farm_County_CI(:,5,mm)=poissinv((0.975-p_inf_County(:))./(1-p_inf_County(:)),mu_farm_County(:));

    lb_z=p_inf_County>=0.025;
    ub_z=p_inf_County>=0.975;
    post_outbreak_poultry_farm_County_CI(lb_z,1,mm)=0;
    post_outbreak_poultry_farm_County_CI(ub_z,2,mm)=0;

    post_spillover_poultry_farm_County_CI(:,1,mm)=poissinv((0.025-p_inf_County(:))./(1-p_inf_County(:)),kappa_spillover.*mu_farm_County(:));
    post_spillover_poultry_farm_County_CI(:,2,mm)=poissinv((0.25-p_inf_County(:))./(1-p_inf_County(:)),kappa_spillover.*mu_farm_County(:));
    post_spillover_poultry_farm_County_CI(:,3,mm)=poissinv((0.5-p_inf_County(:))./(1-p_inf_County(:)),kappa_spillover.*mu_farm_County(:));
    post_spillover_poultry_farm_County_CI(:,4,mm)=poissinv((0.75-p_inf_County(:))./(1-p_inf_County(:)),kappa_spillover.*mu_farm_County(:));
    post_spillover_poultry_farm_County_CI(:,5,mm)=poissinv((0.975-p_inf_County(:))./(1-p_inf_County(:)),kappa_spillover.*mu_farm_County(:));

    lb_z=p_inf_County>=0.025;
    ub_z=p_inf_County>=0.975;
    post_spillover_poultry_farm_County_CI(lb_z,1,mm)=0;
    post_spillover_poultry_farm_County_CI(ub_z,2,mm)=0;

    for ii=0:1000
        if(ii==0)
            post_outbreak_poultry_farm_State(:,1+ii,mm)=z_state(:)+(1-z_state(:)).*nbinpdf(ii,k_state(:),p_nb_state(:));
        else
            post_outbreak_poultry_farm_State(:,1+ii,mm)=(1-z_state(:)).*nbinpdf(ii,k_state(:),p_nb_state(:));
        end
    end
    
    for ii=0:75
        if(ii>0)
            onward_transmission_poultry_farm_County(:,mm)=onward_transmission_poultry_farm_County(:,mm)+(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:)).*(1-p_no_onward_transmission.^ii);
            onward_transmission_poultry_farm_State(:,mm)=onward_transmission_poultry_farm_State(:,mm)+(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:)).*(1-p_no_onward_transmission.^ii);
        end
        if(ii==0)
            post_spillover_poultry_farm_State(:,1+ii,mm)=z_state_spill(:)+(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:));
            post_spillover_poultry_farm_County(:,1+ii,mm)=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
        else
            post_spillover_poultry_farm_State(:,1+ii,mm)=(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:));
            post_spillover_poultry_farm_County(:,1+ii,mm)=(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
        end
    end

    outbreak_risk_poultry_farm_County(:,mm)=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:)));
    t_nan=isnan(outbreak_risk_poultry_farm_County(:,mm));
    outbreak_risk_poultry_farm_County(t_nan,mm)=0;

    outbreak_risk_poultry_farm_State(:,mm)=1-(z_state(:)+(1-z_state(:)).*nbinpdf(0,k_state(:),p_nb_state(:)));
    t_nan=isnan(outbreak_risk_poultry_farm_State(:,mm));
    outbreak_risk_poultry_farm_State(t_nan,mm)=0;

    spillover_risk_poultry_farm_County(:,mm)=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County(:)));
    spillover_risk_poultry_farm_State(:,mm)=1-(z_state_spill(:)+(1-z_state_spill(:)).*nbinpdf(0,k_state_spill(:),p_nb_state_spill(:)));
end

save('Poultry_Risk_Models.mat','par_spillover','no_farms','post_spillover_poultry_farm_County','potential_outbreak_poultry_farm_County','post_outbreak_poultry_farm_County_CI','post_spillover_poultry_farm_County_CI','onward_transmission_poultry_farm_State','onward_transmission_poultry_farm_County','post_spillover_poultry_farm_State','post_outbreak_poultry_farm_State','w_AIC','State_Name','outbreak_poultry_farm_County','outbreak_risk_poultry_farm_County','spillover_poultry_farm_County','spillover_risk_poultry_farm_County','outbreak_poultry_farm_State','outbreak_risk_poultry_farm_State','spillover_poultry_farm_State','spillover_risk_poultry_farm_State');



