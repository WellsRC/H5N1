clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');
load('Poultry_Models_Fit.mat',"par_est","w_AIC",'Poultry_Model');

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);

outbreak_poultry_farm_County=zeros(height(US_County),length(par_est));
outbreak_risk_poultry_farm_County=zeros(height(US_County),length(par_est));

outbreak_layer_farm_County=zeros(height(US_County),length(par_est));
outbreak_pullet_farm_County=zeros(height(US_County),length(par_est));
outbreak_broiler_farm_County=zeros(height(US_County),length(par_est));
outbreak_turkey_farm_County=zeros(height(US_County),length(par_est));

outbreak_risk_layer_farm_County=zeros(height(US_County),length(par_est));
outbreak_risk_pullet_farm_County=zeros(height(US_County),length(par_est));
outbreak_risk_broiler_farm_County=zeros(height(US_County),length(par_est));
outbreak_risk_turkey_farm_County=zeros(height(US_County),length(par_est));

spillover_poultry_farm_County=zeros(height(US_County),length(par_est));
spillover_risk_poultry_farm_County=zeros(height(US_County),length(par_est));

outbreak_poultry_farm_State=zeros(length(State_Name),length(par_est));
outbreak_risk_poultry_farm_State=zeros(length(State_Name),length(par_est));
spillover_poultry_farm_State=zeros(length(State_Name),length(par_est));
spillover_risk_poultry_farm_State=zeros(length(State_Name),length(par_est));

outbreak_layer_farm_State=zeros(length(State_Name),length(par_est));
outbreak_pullet_farm_State=zeros(length(State_Name),length(par_est));
outbreak_broiler_farm_State=zeros(length(State_Name),length(par_est));
outbreak_turkey_farm_State=zeros(length(State_Name),length(par_est));

outbreak_risk_layer_farm_State=zeros(length(State_Name),length(par_est));
outbreak_risk_pullet_farm_State=zeros(length(State_Name),length(par_est));
outbreak_risk_broiler_farm_State=zeros(length(State_Name),length(par_est));
outbreak_risk_turkey_farm_State=zeros(length(State_Name),length(par_est));


post_outbreak_poultry_farm_State=zeros(length(State_Name),2501,length(par_est));
post_spillover_poultry_farm_State=zeros(length(State_Name),101,length(par_est));

for mm=1:length(par_est)
    
    [X_County,P_County,County_Farms,Affected_County_Farms_Unknown,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,state_weight_matrix,State_Spillover_Events,logic_par] = Poultry_Covariates(Poultry_Model.Model_H5N1{mm},Poultry_Model.Model_Farm{mm},Poultry_Model.Model_Stratified_Chicken_Inventory{mm});
    x=par_est{mm};
    
    no_farms=County_Farms==0;
    no_layer_farms=Layer_Farms==0;
    no_pullet_farms=Pullet_Farms==0;
    no_broiler_farms=Broiler_Farms==0;
    no_turkey_farms=Turkey_Farms==0;

    
    beta_x=x([1 (1+2+size(P_County,1)):(2+size(P_County,1)+size(X_County,1))]);
    if(length(beta_x)>1)
        beta_x(2:end)=10.^beta_x(2:end);
    end
    
    beta_p=x([2:(2+size(P_County,1))]);
    if(length(beta_x)>1)
        beta_p(2:end)=-10.^beta_p(2:end);
    end
    
    f_outbreak_turkey=10.^x(end-3);
    f_outbreak_pullet=10.^x(end-2);
    f_outbreak_broiler=10.^x(end-1);
    kappa_spillover=10.^x(end);
    
    p_inf_County=Zero_Inflation(beta_p,P_County);


    p_inf_County_Baseline=p_inf_County;
    p_inf_County_Layer=p_inf_County;
    p_inf_County_Turkey=p_inf_County;
    p_inf_County_Pullet=p_inf_County;
    p_inf_County_Broiler=p_inf_County;
    
    p_inf_County_Layer(Layer_Farms==0)=1;
    p_inf_County_Turkey(Turkey_Farms==0)=1;
    p_inf_County_Pullet(Pullet_Farms==0)=1;
    p_inf_County_Broiler(Broiler_Farms==0)=1;

    mu_farm_County_Basline = Risk_Assesment_Farms(beta_x,X_County);
    mu_farm_County_Layer = Risk_Assesment_Farms(beta_x,X_County);
    
    mu_farm_County_Turkey = f_outbreak_turkey.*mu_farm_County_Basline;
    mu_farm_County_Pullet = f_outbreak_pullet.*mu_farm_County_Basline;
    mu_farm_County_Broiler = f_outbreak_broiler.*mu_farm_County_Basline;
    
    mu_farm_County_Basline(County_Farms==0)=0;
    
    dX=County_Farms-Layer_Farms-Turkey_Farms-Pullet_Farms-Broiler_Farms;
    mu_farm_County_Basline(dX<=0)=0;
    p_inf_County_Baseline(dX<=0)=1;
    
    mu_farm_County_Layer(Layer_Farms==0)=0;
    mu_farm_County_Turkey(Turkey_Farms==0)=0;
    mu_farm_County_Pullet(Pullet_Farms==0)=0;
    mu_farm_County_Broiler(Broiler_Farms==0)=0;
    
    mu_farm_State_Baseline=zeros(size(State_Spillover_Events));
    mu_farm_State_Layer=zeros(size(State_Spillover_Events));
    mu_farm_State_Turkey=zeros(size(State_Spillover_Events));
    mu_farm_State_Pullet=zeros(size(State_Spillover_Events));
    mu_farm_State_Broiler=zeros(size(State_Spillover_Events));
    
    for ss=1:length(mu_farm_State_Layer)
        mu_farm_State_Baseline(ss)=state_weight_matrix(ss,:)*((1-p_inf_County_Baseline(:)).*mu_farm_County_Basline(:));
        mu_farm_State_Layer(ss)=state_weight_matrix(ss,:)*((1-p_inf_County_Layer(:)).*mu_farm_County_Layer(:));
        mu_farm_State_Turkey(ss)=state_weight_matrix(ss,:)*((1-p_inf_County_Turkey(:)).*mu_farm_County_Turkey(:));
        mu_farm_State_Pullet(ss)=state_weight_matrix(ss,:)*((1-p_inf_County_Pullet(:)).*mu_farm_County_Pullet(:));
        mu_farm_State_Broiler(ss)=state_weight_matrix(ss,:)*((1-p_inf_County_Broiler(:)).*mu_farm_County_Broiler(:));
    end
    
    mu_farm_County_Total=mu_farm_County_Layer+mu_farm_County_Turkey+mu_farm_County_Pullet+mu_farm_County_Broiler+mu_farm_County_Basline;
    mu_farm_State_Total=mu_farm_State_Layer+mu_farm_State_Turkey+mu_farm_State_Pullet+mu_farm_State_Broiler+mu_farm_State_Baseline;

    outbreak_poultry_farm_County(:,mm)=(1-p_inf_County(:)).*mu_farm_County_Total(:);
    outbreak_layer_farm_County(:,mm)=(1-p_inf_County_Layer(:)).*mu_farm_County_Layer(:);
    outbreak_pullet_farm_County(:,mm)=(1-p_inf_County_Pullet(:)).*mu_farm_County_Pullet(:);
    outbreak_broiler_farm_County(:,mm)=(1-p_inf_County_Broiler(:)).*mu_farm_County_Broiler(:);
    outbreak_turkey_farm_County(:,mm)=(1-p_inf_County_Turkey(:)).*mu_farm_County_Turkey(:);

    outbreak_poultry_farm_State(:,mm)=mu_farm_State_Total(:);
    outbreak_layer_farm_State(:,mm)=mu_farm_State_Layer(:);
    outbreak_pullet_farm_State(:,mm)=mu_farm_State_Pullet(:);
    outbreak_broiler_farm_State(:,mm)=mu_farm_State_Broiler(:);
    outbreak_turkey_farm_State(:,mm)=mu_farm_State_Turkey(:);

    spillover_poultry_farm_County(:,mm)=kappa_spillover.*outbreak_poultry_farm_County(:,mm);
    spillover_poultry_farm_State(:,mm)=kappa_spillover.*outbreak_poultry_farm_State(:,mm);

    post_outbreak_poultry_farm_State(:,1,mm)=poisspdf(0,mu_farm_State_Total(:));
    for ii=1:2500
        post_outbreak_poultry_farm_State(:,1+ii,mm)=poisspdf(ii,mu_farm_State_Total(:));
    end
    
    post_spillover_poultry_farm_State(:,1,mm)=poisspdf(0,spillover_poultry_farm_State(:,mm));
    for ii=1:100
        post_outbreak_poultry_farm_State(:,1+ii,mm)=poisspdf(ii,spillover_poultry_farm_State(:,mm));
    end

    outbreak_risk_poultry_farm_County(:,mm)=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County_Total(:)));
    t_nan=isnan(outbreak_risk_poultry_farm_County(:,mm));
    outbreak_risk_poultry_farm_County(t_nan,mm)=0;

    outbreak_risk_layer_farm_County(:,mm)=1-(p_inf_County_Layer(:)+(1-p_inf_County_Layer(:)).*poisspdf(0,mu_farm_County_Layer(:)));
    t_nan=isnan(outbreak_risk_layer_farm_County(:,mm));
    outbreak_risk_layer_farm_County(t_nan,mm)=0;

    outbreak_risk_pullet_farm_County(:,mm)=1-(p_inf_County_Pullet(:)+(1-p_inf_County_Pullet(:)).*poisspdf(0,mu_farm_County_Pullet(:)));
    t_nan=isnan(outbreak_risk_pullet_farm_County(:,mm));
    outbreak_risk_pullet_farm_County(t_nan,mm)=0;

    outbreak_risk_broiler_farm_County(:,mm)=1-(p_inf_County_Broiler(:)+(1-p_inf_County_Broiler(:)).*poisspdf(0,mu_farm_County_Broiler(:)));
    t_nan=isnan(outbreak_risk_broiler_farm_County(:,mm));
    outbreak_risk_broiler_farm_County(t_nan,mm)=0;

    outbreak_risk_turkey_farm_County(:,mm)=1-(p_inf_County_Turkey(:)+(1-p_inf_County_Turkey(:)).*poisspdf(0,mu_farm_County_Turkey(:)));
    t_nan=isnan(outbreak_risk_turkey_farm_County(:,mm));
    outbreak_risk_turkey_farm_County(t_nan,mm)=0;

    outbreak_risk_poultry_farm_State(:,mm)=1-(poisspdf(0,mu_farm_State_Total(:)));
    t_nan=isnan(outbreak_risk_poultry_farm_State(:,mm));
    outbreak_risk_poultry_farm_State(t_nan,mm)=0;

    outbreak_risk_layer_farm_State(:,mm)=1-(poisspdf(0,mu_farm_State_Layer(:)));
    t_nan=isnan(outbreak_risk_layer_farm_State(:,mm));
    outbreak_risk_layer_farm_State(t_nan,mm)=0;

    outbreak_risk_pullet_farm_State(:,mm)=1-(poisspdf(0,mu_farm_State_Pullet(:)));
    t_nan=isnan(outbreak_risk_pullet_farm_State(:,mm));
    outbreak_risk_pullet_farm_State(t_nan,mm)=0;

    outbreak_risk_broiler_farm_State(:,mm)=1-(poisspdf(0,mu_farm_State_Broiler(:)));
    t_nan=isnan(outbreak_risk_broiler_farm_State(:,mm));
    outbreak_risk_broiler_farm_State(t_nan,mm)=0;

    outbreak_risk_turkey_farm_State(:,mm)=1-(poisspdf(0,mu_farm_State_Turkey(:)));
    t_nan=isnan(outbreak_risk_turkey_farm_State(:,mm));
    outbreak_risk_turkey_farm_State(t_nan,mm)=0;

    spillover_risk_poultry_farm_County(:,mm)=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County_Total(:)));
    spillover_risk_poultry_farm_State(:,mm)=1-(poisspdf(0,kappa_spillover.*mu_farm_State_Total(:)));

    outbreak_poultry_farm_County(no_farms,mm)=NaN;
    outbreak_layer_farm_County(no_layer_farms,mm)=NaN;
    outbreak_pullet_farm_County(no_pullet_farms,mm)=NaN;
    outbreak_broiler_farm_County(no_broiler_farms,mm)=NaN;
    outbreak_turkey_farm_County(no_turkey_farms,mm)=NaN;
    
    spillover_poultry_farm_County(no_farms,mm)=NaN;
    
    outbreak_risk_poultry_farm_County(no_farms,mm)=NaN;
    outbreak_risk_layer_farm_County(no_layer_farms,mm)=NaN;
    outbreak_risk_pullet_farm_County(no_pullet_farms,mm)=NaN;
    outbreak_risk_broiler_farm_County(no_broiler_farms,mm)=NaN;
    outbreak_risk_turkey_farm_County(no_turkey_farms,mm)=NaN;
    
    spillover_risk_poultry_farm_County(no_farms,mm)=NaN;
end

save('Average_Risk_Poultry.mat','post_spillover_poultry_farm_State','post_outbreak_poultry_farm_State','w_AIC','State_Name','outbreak_poultry_farm_County','outbreak_risk_poultry_farm_County','outbreak_layer_farm_County','outbreak_pullet_farm_County','outbreak_broiler_farm_County','outbreak_turkey_farm_County','outbreak_risk_layer_farm_County','outbreak_risk_pullet_farm_County','outbreak_risk_broiler_farm_County','outbreak_risk_turkey_farm_County','spillover_poultry_farm_County','spillover_risk_poultry_farm_County','outbreak_poultry_farm_State','outbreak_risk_poultry_farm_State','spillover_poultry_farm_State','spillover_risk_poultry_farm_State','outbreak_layer_farm_State','outbreak_pullet_farm_State','outbreak_broiler_farm_State','outbreak_turkey_farm_State','outbreak_risk_layer_farm_State','outbreak_risk_pullet_farm_State','outbreak_risk_broiler_farm_State','outbreak_risk_turkey_farm_State');



