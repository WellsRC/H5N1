clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);

load('Uncertainty_AIC_Poultry_Model.mat','L_samp','par_samp','Poultry_Model','par_mle','L_mle');

% Only model that is being run
[F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_par,logic_temperature] = Poultry_Covariates(Poultry_Model.Model_H5N1{1},Poultry_Model.Model_Farm{1});
no_farms=County_Farms==0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% MLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
N_Samp=1;

mle_onward_transmission_poultry_farm_County=zeros(height(US_County),N_Samp);
mle_onward_transmission_poultry_farm_State=zeros(length(State_Name),N_Samp);

mle_post_outbreak_poultry_farm_State=zeros(length(State_Name),1001);
mle_post_spillover_poultry_farm_State=zeros(length(State_Name),76);

mle_post_spillover_poultry_farm_County=zeros(height(US_County),76);
mle_post_outbreak_poultry_farm_County_CI=zeros(height(US_County),5);
mle_post_spillover_poultry_farm_County_CI=zeros(height(US_County),5);

k_onward_transmission=2.69;
R0=0.05;
p_no_onward_transmission=nbinpdf(0,k_onward_transmission,k_onward_transmission./(k_onward_transmission+R0));

x=par_mle;  
    
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
mle_par_spillover=kappa_spillover;

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

mle_potential_outbreak_poultry_farm_County=mu_farm_County(:);
mle_outbreak_poultry_farm_County=(1-p_inf_County(:)).*mu_farm_County(:);

mle_outbreak_poultry_farm_State=mu_farm_State(:);

mle_spillover_poultry_farm_County=kappa_spillover.*mle_outbreak_poultry_farm_County;
mle_spillover_poultry_farm_State=kappa_spillover.*mle_outbreak_poultry_farm_State;

mle_post_outbreak_poultry_farm_County_CI(:,1)=poissinv((0.025-p_inf_County(:))./(1-p_inf_County(:)),mu_farm_County(:));
mle_post_outbreak_poultry_farm_County_CI(:,2)=poissinv((0.25-p_inf_County(:))./(1-p_inf_County(:)),mu_farm_County(:));
mle_post_outbreak_poultry_farm_County_CI(:,3)=poissinv((0.5-p_inf_County(:))./(1-p_inf_County(:)),mu_farm_County(:));
mle_post_outbreak_poultry_farm_County_CI(:,4)=poissinv((0.75-p_inf_County(:))./(1-p_inf_County(:)),mu_farm_County(:));
mle_post_outbreak_poultry_farm_County_CI(:,5)=poissinv((0.975-p_inf_County(:))./(1-p_inf_County(:)),mu_farm_County(:));

mle_post_outbreak_poultry_farm_County_CI(p_inf_County>=0.025,1)=0;
mle_post_outbreak_poultry_farm_County_CI(p_inf_County>=0.25,2)=0;
mle_post_outbreak_poultry_farm_County_CI(p_inf_County>=0.5,3)=0;
mle_post_outbreak_poultry_farm_County_CI(p_inf_County>=0.75,4)=0;
mle_post_outbreak_poultry_farm_County_CI(p_inf_County>=0.975,5)=0;

mle_post_spillover_poultry_farm_County_CI(:,1)=poissinv((0.025-p_inf_County(:))./(1-p_inf_County(:)),kappa_spillover.*mu_farm_County(:));
mle_post_spillover_poultry_farm_County_CI(:,2)=poissinv((0.25-p_inf_County(:))./(1-p_inf_County(:)),kappa_spillover.*mu_farm_County(:));
mle_post_spillover_poultry_farm_County_CI(:,3)=poissinv((0.5-p_inf_County(:))./(1-p_inf_County(:)),kappa_spillover.*mu_farm_County(:));
mle_post_spillover_poultry_farm_County_CI(:,4)=poissinv((0.75-p_inf_County(:))./(1-p_inf_County(:)),kappa_spillover.*mu_farm_County(:));
mle_post_spillover_poultry_farm_County_CI(:,5)=poissinv((0.975-p_inf_County(:))./(1-p_inf_County(:)),kappa_spillover.*mu_farm_County(:));

mle_post_spillover_poultry_farm_County_CI(p_inf_County>=0.025,1)=0;
mle_post_spillover_poultry_farm_County_CI(p_inf_County>=0.25,2)=0;
mle_post_spillover_poultry_farm_County_CI(p_inf_County>=0.5,3)=0;
mle_post_spillover_poultry_farm_County_CI(p_inf_County>=0.75,4)=0;
mle_post_spillover_poultry_farm_County_CI(p_inf_County>=0.975,5)=0;

for ii=0:1000
    if(ii==0)
        mle_post_outbreak_poultry_farm_State(:,1+ii)=z_state(:)+(1-z_state(:)).*nbinpdf(ii,k_state(:),p_nb_state(:));
    else
        mle_post_outbreak_poultry_farm_State(:,1+ii)=(1-z_state(:)).*nbinpdf(ii,k_state(:),p_nb_state(:));
    end
end

for ii=0:75
    if(ii>0)
        mle_onward_transmission_poultry_farm_County=mle_onward_transmission_poultry_farm_County+(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:)).*(1-p_no_onward_transmission.^ii);
        mle_onward_transmission_poultry_farm_State=mle_onward_transmission_poultry_farm_State+(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:)).*(1-p_no_onward_transmission.^ii);
    end
    if(ii==0)
        mle_post_spillover_poultry_farm_State(:,1+ii)=z_state_spill(:)+(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:));
        mle_post_spillover_poultry_farm_County(:,1+ii)=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
    else
        mle_post_spillover_poultry_farm_State(:,1+ii)=(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:));
        mle_post_spillover_poultry_farm_County(:,1+ii)=(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
    end
end

mle_outbreak_risk_poultry_farm_County=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:)));
t_nan=isnan(mle_outbreak_risk_poultry_farm_County);
mle_outbreak_risk_poultry_farm_County(t_nan)=0;

mle_outbreak_risk_poultry_farm_State=1-(z_state(:)+(1-z_state(:)).*nbinpdf(0,k_state(:),p_nb_state(:)));
t_nan=isnan(mle_outbreak_risk_poultry_farm_State);
mle_outbreak_risk_poultry_farm_State(t_nan)=0;

mle_spillover_risk_poultry_farm_County=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County(:)));
mle_spillover_risk_poultry_farm_State=1-(z_state_spill(:)+(1-z_state_spill(:)).*nbinpdf(0,k_state_spill(:),p_nb_state_spill(:)));

mle_outbreak_poultry_farm_County(no_farms)=NaN;

mle_potntial_outbreak_poultry_farm_County(no_farms)=NaN;

mle_spillover_poultry_farm_County(no_farms)=NaN;

mle_outbreak_risk_poultry_farm_County(no_farms)=NaN;

mle_spillover_risk_poultry_farm_County(no_farms)=NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Uncertainty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
N_Samp=length(L_samp);

potential_outbreak_poultry_farm_County=zeros(height(US_County),N_Samp);

outbreak_poultry_farm_County=zeros(height(US_County),N_Samp);
outbreak_risk_poultry_farm_County=zeros(height(US_County),N_Samp);

spillover_poultry_farm_County=zeros(height(US_County),N_Samp);
spillover_risk_poultry_farm_County=zeros(height(US_County),N_Samp);

outbreak_poultry_farm_State=zeros(length(State_Name),N_Samp);
outbreak_risk_poultry_farm_State=zeros(length(State_Name),N_Samp);
spillover_poultry_farm_State=zeros(length(State_Name),N_Samp);
spillover_risk_poultry_farm_State=zeros(length(State_Name),N_Samp);

onward_transmission_poultry_farm_County=zeros(height(US_County),N_Samp);
onward_transmission_poultry_farm_State=zeros(length(State_Name),N_Samp);

par_spillover=zeros(1,N_Samp);

for mm=1:N_Samp       

    x=par_samp(mm,:);  
    
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
    
    for ii=0:75
        if(ii>0)
            onward_transmission_poultry_farm_County(:,mm)=onward_transmission_poultry_farm_County(:,mm)+(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:)).*(1-p_no_onward_transmission.^ii);
            onward_transmission_poultry_farm_State(:,mm)=onward_transmission_poultry_farm_State(:,mm)+(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:)).*(1-p_no_onward_transmission.^ii);
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

    outbreak_poultry_farm_County(no_farms,mm)=NaN;
    
    potential_outbreak_poultry_farm_County(no_farms,mm)=NaN;

    spillover_poultry_farm_County(no_farms,mm)=NaN;
    
    outbreak_risk_poultry_farm_County(no_farms,mm)=NaN;
    
    spillover_risk_poultry_farm_County(no_farms,mm)=NaN;
end

min_outbreaks=min(outbreak_poultry_farm_County(~no_farms,:),[],1);
min_outbreaks_95=prctile(min_outbreaks,[2.5 97.5]);
max_outbreaks=max(outbreak_poultry_farm_County(~no_farms,:),[],1);
max_outbreaks_95=prctile(max_outbreaks,[2.5 97.5]);
median_outbreaks=median(outbreak_poultry_farm_County(~no_farms,:),1);
median_outbreaks_95=prctile(median_outbreaks,[2.5 97.5]);

County_outbreak_risk_over_50=sum(outbreak_risk_poultry_farm_County>0.5,1);
County_outbreak_risk_over_50_95=prctile(County_outbreak_risk_over_50,[2.5 97.5]);
potential_outbreak_95_prctl=prctile(potential_outbreak_poultry_farm_County(~no_farms,:),95,1);
potential_outbreak_95_prctl_95=prctile(potential_outbreak_95_prctl,[2.5 97.5]);

temp_C=repmat(potential_outbreak_95_prctl,size(potential_outbreak_poultry_farm_County,1),1);
temp_C_over_95p=potential_outbreak_poultry_farm_County>temp_C;

per_risk_under_10=zeros(1,N_Samp);
for kk=1:N_Samp
    temp_r=outbreak_risk_poultry_farm_County(temp_C_over_95p(:,kk),kk)<0.1;
    per_risk_under_10(kk)=mean(temp_r);
end

per_risk_under_10_95=prctile(per_risk_under_10,[2.5 97.5]);

potential_outbreak_poultry_farm_County_95=prctile(potential_outbreak_poultry_farm_County,[2.5 97.5],2);

outbreak_poultry_farm_County_95=prctile(outbreak_poultry_farm_County,[2.5 97.5],2);
outbreak_risk_poultry_farm_County_95=prctile(outbreak_risk_poultry_farm_County,[2.5 97.5],2);

spillover_poultry_farm_County_95=prctile(spillover_poultry_farm_County,[2.5 97.5],2);
spillover_risk_poultry_farm_County_95=prctile(spillover_risk_poultry_farm_County,[2.5 97.5],2);

outbreak_poultry_farm_State_95=prctile(outbreak_poultry_farm_State,[2.5 97.5],2);
outbreak_risk_poultry_farm_State_95=prctile(outbreak_risk_poultry_farm_State,[2.5 97.5],2);
spillover_poultry_farm_State_95=prctile(spillover_poultry_farm_State,[2.5 97.5],2);
spillover_risk_poultry_farm_State_95=prctile(spillover_risk_poultry_farm_State,[2.5 97.5],2);

onward_transmission_poultry_farm_County_95=prctile(onward_transmission_poultry_farm_County,[2.5 97.5],2);
onward_transmission_poultry_farm_State_95=prctile(onward_transmission_poultry_farm_State,[2.5 97.5],2);

par_spillover_95=prctile(par_spillover,[2.5 97.5]);


save('Poultry_Risk_AIC.mat','mle_potential_outbreak_poultry_farm_County','min_outbreaks_95','max_outbreaks_95','median_outbreaks_95','per_risk_under_10_95','potential_outbreak_95_prctl_95','County_outbreak_risk_over_50_95','potential_outbreak_poultry_farm_County_95','outbreak_poultry_farm_County_95','outbreak_risk_poultry_farm_County_95','spillover_poultry_farm_County_95','spillover_risk_poultry_farm_County_95','outbreak_poultry_farm_State_95','outbreak_risk_poultry_farm_State_95','spillover_poultry_farm_State_95','spillover_risk_poultry_farm_State_95','onward_transmission_poultry_farm_County_95','onward_transmission_poultry_farm_State_95','par_spillover_95','no_farms','mle_outbreak_poultry_farm_County','mle_outbreak_risk_poultry_farm_County','mle_spillover_poultry_farm_County','mle_spillover_risk_poultry_farm_County','mle_outbreak_poultry_farm_State','mle_outbreak_risk_poultry_farm_State','mle_spillover_poultry_farm_State','mle_spillover_risk_poultry_farm_State','mle_onward_transmission_poultry_farm_County','mle_onward_transmission_poultry_farm_State','mle_post_outbreak_poultry_farm_State','mle_post_spillover_poultry_farm_State','mle_post_spillover_poultry_farm_County','mle_post_outbreak_poultry_farm_County_CI','mle_post_spillover_poultry_farm_County_CI','mle_par_spillover');



