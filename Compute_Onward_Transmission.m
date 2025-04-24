clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');

State_Name=unique(US_County.STATE_NAME);
state_remove=strcmp(State_Name,"Alaska") | strcmp(State_Name,"District of Columbia");
State_Name=State_Name(~state_remove);

load('Uncertainty_AIC_Dairy_Model.mat','par_samp','Dairy_Model','par_mle');
dairy_par_samp=par_samp;

% Only model that is being run
[F_County_Dairy,X_County_Dairy,P_County_Dairy,County_Farms_Dairy,~,State_Spillover_Events_dairy,Affected_State_Farms_Dairy,state_weight_matrix_Dairy,Dairy_Network,logic_connect_Dairy,logic_connect_p_Dairy,~,logic_temperature_Dairy]= Dairy_Covariates(Dairy_Model.Model_H5N1{1},Dairy_Model.Model_Farm{1},Dairy_Model.Model_Stratified_Operations{1});
no_diary_farms=County_Farms_Dairy==0;

k_onward_transmission=2.69;
R0=0.05;
p_no_onward_transmission=nbinpdf(0,k_onward_transmission,k_onward_transmission./(k_onward_transmission+R0));

post_spillover_dairy_farm_State=zeros(length(State_Name),101);
post_spillover_dairy_farm_County=zeros(height(US_County),101);

x=par_mle;  

indx_pinf=[5:(8+size(P_County_Dairy,1))];

beta_x=x(~ismember(1:length(x),[indx_pinf length(x)]));
if(length(beta_x)>4)
    beta_x(5:end)=10.^beta_x(5:end);
end
    
beta_p=x(indx_pinf);
if(length(beta_p)>4)
    if(logic_temperature_Dairy & isempty(Dairy_Network))
        beta_p(5:end-1)=-10.^beta_p(5:end-1);    
    elseif(logic_temperature_Dairy & ~isempty(Dairy_Network))
        beta_p(5:end-2)=-10.^beta_p(5:end-2);
        beta_p(end)=-10.^beta_p(end);
    else
        beta_p(5:end)=-10.^beta_p(5:end);
    end
end


kappa_spillover=10.^x(end);

if(~isempty(Dairy_Network))
    beta_x_temp=beta_x([1:4 4+find(~logic_connect_Dairy)]);
    
    mu_farm_temp = Risk_Assesment_Farms(beta_x_temp,[F_County_Dairy; X_County_Dairy(~logic_connect_Dairy,:)]);
    mu_farm_temp(County_Farms_Dairy==0)=0;

    beta_p_temp=beta_p([1:4 4+find(~logic_connect_p_Dairy)]);
    p_inf_County_temp=Zero_Inflation(beta_p_temp,[F_County_Dairy; P_County_Dairy(~logic_connect_p_Dairy,:)]);
    p_inf_County_temp(County_Farms_Dairy==0)=1;

    temp_r=(1-p_inf_County_temp(:)).*mu_farm_temp(:);

    temp_r=(temp_r(:)'*Dairy_Network);
    tt=min(temp_r(temp_r>0));
    temp_r(temp_r(:)==0 & County_Farms_Dairy>0)=tt./10;

    X_County_Dairy(logic_connect_Dairy,:)=log10(temp_r);
    P_County_Dairy(logic_connect_p_Dairy,:)=log10(temp_r);
end
    
p_inf_County=Zero_Inflation(beta_p,[F_County_Dairy; P_County_Dairy]);
p_inf_County(County_Farms_Dairy==0)=1;

mu_farm_County = Risk_Assesment_Farms(beta_x,[F_County_Dairy; X_County_Dairy]);
mu_farm_County(County_Farms_Dairy==0)=0;

mu_farm_State=zeros(size(State_Spillover_Events_dairy));

k_state=zeros(size(State_Spillover_Events_dairy));
p_nb_state=zeros(size(State_Spillover_Events_dairy));
z_state=zeros(size(State_Spillover_Events_dairy));

k_spill_State=zeros(size(State_Spillover_Events_dairy));
z_state_spill=zeros(size(State_Spillover_Events_dairy));
p_nb_state_spill=zeros(size(State_Spillover_Events_dairy));

p_temp=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:));
p_temp_spill=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County(:));

r=10.^linspace(-3,4,501);
    
w_state=zeros(size(State_Spillover_Events_dairy));
County_w_Farms=County_Farms_Dairy>0;

for ss=1:length(w_state)
    w_state(ss)=state_weight_matrix_Dairy(ss,:)*County_w_Farms(:);
    mu_farm_State(ss)=state_weight_matrix_Dairy(ss,:)*((1-p_inf_County(:)).*mu_farm_County(:));

    mu_County_NB=state_weight_matrix_Dairy(ss,:)*((1-p_inf_County(:)).*mu_farm_County(:)); 
    var_County_NB=state_weight_matrix_Dairy(ss,:)*((1-p_inf_County(:)).*(mu_farm_County(:)+p_inf_County(:).*mu_farm_County(:).^2)); 
    p_zero_county=prod(p_temp(state_weight_matrix_Dairy(ss,:)==1));
    if(p_zero_county==0)
        p_zero_county=10^(-64);
    end
    
    [z_state(ss),k_state(ss),p_nb_state(ss),L_Nan1]=Estimate_ZIF_NB(mu_County_NB,var_County_NB,p_zero_county);


    mu_County_NB=state_weight_matrix_Dairy(ss,:)*((1-p_inf_County(:)).*kappa_spillover.*mu_farm_County(:)); 
    var_County_NB=state_weight_matrix_Dairy(ss,:)*((1-p_inf_County(:)).*(kappa_spillover.*mu_farm_County(:)+p_inf_County(:).*(kappa_spillover.*mu_farm_County(:)).^2)); 
    p_zero_county=prod(p_temp_spill(state_weight_matrix_Dairy(ss,:)==1));
    if(p_zero_county==0)
        p_zero_county=10^(-64);
    end
    [z_state_spill(ss),k_spill_State(ss),p_nb_state_spill(ss),L_Nan2]=Estimate_ZIF_NB(mu_County_NB,var_County_NB,p_zero_county);

    
end
 

for ii=0:100
    if(ii==0)
        post_spillover_dairy_farm_County(:,1+ii)=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
        post_spillover_dairy_farm_State(:,1+ii)=z_state_spill(:)+(1-z_state_spill(:)).*nbinpdf(ii,k_spill_State(:),p_nb_state_spill(:)); 
    else
        post_spillover_dairy_farm_State(:,1+ii)=(1-z_state_spill(:)).*nbinpdf(ii,k_spill_State(:),p_nb_state_spill(:)); 
        post_spillover_dairy_farm_County(:,1+ii)=(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Poultry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('Uncertainty_AIC_Poultry_Model.mat','par_samp','Poultry_Model','par_mle');
poultry_par_samp=par_samp;
% Only model that is being run
[F_County,X_County,P_County,County_Farms,~,state_weight_matrix,State_Spillover_Events,~,logic_temperature] = Poultry_Covariates(Poultry_Model.Model_H5N1{1},Poultry_Model.Model_Farm{1});
no_farms=County_Farms==0;

post_spillover_poultry_farm_State=zeros(length(State_Name),76);
post_spillover_poultry_farm_County=zeros(height(US_County),76);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
x=par_mle; 
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

for ii=0:75
    if(ii==0)
        post_spillover_poultry_farm_State(:,1+ii)=z_state_spill(:)+(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:));
        post_spillover_poultry_farm_County(:,1+ii)=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
    else
        post_spillover_poultry_farm_State(:,1+ii)=(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:));
        post_spillover_poultry_farm_County(:,1+ii)=(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Onward transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mle_County_onward=zeros(length(mu_farm_County),1);
mle_State_onward=zeros(length(mu_farm_State),1);
mle_State_dairy_onward=zeros(length(mu_farm_State),1);
mle_County_dairy_onward=zeros(length(mu_farm_County),1);

for pp=0:size(post_spillover_poultry_farm_County,2)-1
    for dd=0:size(post_spillover_dairy_farm_County,2)-1
        mle_County_onward=mle_County_onward+post_spillover_poultry_farm_County(:,pp+1).*post_spillover_dairy_farm_County(:,dd+1).*(1-p_no_onward_transmission.^(pp+dd));
        mle_County_dairy_onward=mle_County_dairy_onward+post_spillover_poultry_farm_County(:,pp+1).*post_spillover_dairy_farm_County(:,dd+1).*(1-p_no_onward_transmission.^(dd));
        mle_State_onward=mle_State_onward+post_spillover_poultry_farm_State(:,pp+1).*post_spillover_dairy_farm_State(:,dd+1).*(1-p_no_onward_transmission.^(pp+dd));
        mle_State_dairy_onward=mle_State_dairy_onward+post_spillover_poultry_farm_State(:,pp+1).*post_spillover_dairy_farm_State(:,dd+1).*(1-p_no_onward_transmission.^(dd));
    end
end

mle_State_dairy_onward=mle_State_dairy_onward./mle_State_onward;
mle_County_dairy_onward=mle_County_dairy_onward./mle_County_onward;

County_onward=zeros(length(mle_County_onward),size(dairy_par_samp,1));
County_dairy_onward=zeros(length(mle_County_onward),size(dairy_par_samp,1));
State_onward=zeros(length(mle_State_onward),size(dairy_par_samp,1));
State_dairy_onward=zeros(length(mle_State_onward),size(dairy_par_samp,1));

for mm=1:size(dairy_par_samp,1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Dairy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    x=dairy_par_samp(mm,:);  
    
    indx_pinf=[5:(8+size(P_County_Dairy,1))];
    
    beta_x=x(~ismember(1:length(x),[indx_pinf length(x)]));
    if(length(beta_x)>4)
        beta_x(5:end)=10.^beta_x(5:end);
    end
        
    beta_p=x(indx_pinf);
    if(length(beta_p)>4)
        if(logic_temperature_Dairy & isempty(Dairy_Network))
            beta_p(5:end-1)=-10.^beta_p(5:end-1);    
        elseif(logic_temperature_Dairy & ~isempty(Dairy_Network))
            beta_p(5:end-2)=-10.^beta_p(5:end-2);
            beta_p(end)=-10.^beta_p(end);
        else
            beta_p(5:end)=-10.^beta_p(5:end);
        end
    end
    
    
    kappa_spillover=10.^x(end);
    
    if(~isempty(Dairy_Network))
        beta_x_temp=beta_x([1:4 4+find(~logic_connect_Dairy)]);
        
        mu_farm_temp = Risk_Assesment_Farms(beta_x_temp,[F_County_Dairy; X_County_Dairy(~logic_connect_Dairy,:)]);
        mu_farm_temp(County_Farms_Dairy==0)=0;
    
        beta_p_temp=beta_p([1:4 4+find(~logic_connect_p_Dairy)]);
        p_inf_County_temp=Zero_Inflation(beta_p_temp,[F_County_Dairy; P_County_Dairy(~logic_connect_p_Dairy,:)]);
        p_inf_County_temp(County_Farms_Dairy==0)=1;
    
        temp_r=(1-p_inf_County_temp(:)).*mu_farm_temp(:);
    
        temp_r=(temp_r(:)'*Dairy_Network);
        tt=min(temp_r(temp_r>0));
        temp_r(temp_r(:)==0 & County_Farms_Dairy>0)=tt./10;
    
        X_County_Dairy(logic_connect_Dairy,:)=log10(temp_r);
        P_County_Dairy(logic_connect_p_Dairy,:)=log10(temp_r);
    end
        
    p_inf_County=Zero_Inflation(beta_p,[F_County_Dairy; P_County_Dairy]);
    p_inf_County(County_Farms_Dairy==0)=1;
    
    mu_farm_County = Risk_Assesment_Farms(beta_x,[F_County_Dairy; X_County_Dairy]);
    mu_farm_County(County_Farms_Dairy==0)=0;
    
    mu_farm_State=zeros(size(State_Spillover_Events_dairy));
    
    k_state=zeros(size(State_Spillover_Events_dairy));
    p_nb_state=zeros(size(State_Spillover_Events_dairy));
    z_state=zeros(size(State_Spillover_Events_dairy));
    
    k_spill_State=zeros(size(State_Spillover_Events_dairy));
    z_state_spill=zeros(size(State_Spillover_Events_dairy));
    p_nb_state_spill=zeros(size(State_Spillover_Events_dairy));
    
    p_temp=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:));
    p_temp_spill=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County(:));
    
    r=10.^linspace(-3,4,501);
        
    w_state=zeros(size(State_Spillover_Events_dairy));
    County_w_Farms=County_Farms_Dairy>0;
    
    for ss=1:length(w_state)
        w_state(ss)=state_weight_matrix_Dairy(ss,:)*County_w_Farms(:);
        mu_farm_State(ss)=state_weight_matrix_Dairy(ss,:)*((1-p_inf_County(:)).*mu_farm_County(:));
    
        mu_County_NB=state_weight_matrix_Dairy(ss,:)*((1-p_inf_County(:)).*mu_farm_County(:)); 
        var_County_NB=state_weight_matrix_Dairy(ss,:)*((1-p_inf_County(:)).*(mu_farm_County(:)+p_inf_County(:).*mu_farm_County(:).^2)); 
        p_zero_county=prod(p_temp(state_weight_matrix_Dairy(ss,:)==1));
        if(p_zero_county==0)
            p_zero_county=10^(-64);
        end
        
        [z_state(ss),k_state(ss),p_nb_state(ss),L_Nan1]=Estimate_ZIF_NB(mu_County_NB,var_County_NB,p_zero_county);
    
    
        mu_County_NB=state_weight_matrix_Dairy(ss,:)*((1-p_inf_County(:)).*kappa_spillover.*mu_farm_County(:)); 
        var_County_NB=state_weight_matrix_Dairy(ss,:)*((1-p_inf_County(:)).*(kappa_spillover.*mu_farm_County(:)+p_inf_County(:).*(kappa_spillover.*mu_farm_County(:)).^2)); 
        p_zero_county=prod(p_temp_spill(state_weight_matrix_Dairy(ss,:)==1));
        if(p_zero_county==0)
            p_zero_county=10^(-64);
        end
        [z_state_spill(ss),k_spill_State(ss),p_nb_state_spill(ss),L_Nan2]=Estimate_ZIF_NB(mu_County_NB,var_County_NB,p_zero_county);
    
        
    end
     
    
    for ii=0:100
        if(ii==0)
            post_spillover_dairy_farm_County(:,1+ii)=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
            post_spillover_dairy_farm_State(:,1+ii)=z_state_spill(:)+(1-z_state_spill(:)).*nbinpdf(ii,k_spill_State(:),p_nb_state_spill(:)); 
        else
            post_spillover_dairy_farm_State(:,1+ii)=(1-z_state_spill(:)).*nbinpdf(ii,k_spill_State(:),p_nb_state_spill(:)); 
            post_spillover_dairy_farm_County(:,1+ii)=(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Poultry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    x=poultry_par_samp(mm,:); 
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
    
    for ii=0:75
        if(ii==0)
            post_spillover_poultry_farm_State(:,1+ii)=z_state_spill(:)+(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:));
            post_spillover_poultry_farm_County(:,1+ii)=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
        else
            post_spillover_poultry_farm_State(:,1+ii)=(1-z_state_spill(:)).*nbinpdf(ii,k_state_spill(:),p_nb_state_spill(:));
            post_spillover_poultry_farm_County(:,1+ii)=(1-p_inf_County(:)).*poisspdf(ii,kappa_spillover.*mu_farm_County(:));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Onward transmission
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    for pp=0:size(post_spillover_poultry_farm_County,2)-1
        for dd=0:size(post_spillover_dairy_farm_County,2)-1
             County_onward(:,mm)= County_onward(:,mm)+post_spillover_poultry_farm_County(:,pp+1).*post_spillover_dairy_farm_County(:,dd+1).*(1-p_no_onward_transmission.^(pp+dd));
             County_dairy_onward(:,mm)= County_dairy_onward(:,mm)+post_spillover_poultry_farm_County(:,pp+1).*post_spillover_dairy_farm_County(:,dd+1).*(1-p_no_onward_transmission.^(dd));
            State_onward(:,mm)=State_onward(:,mm)+post_spillover_poultry_farm_State(:,pp+1).*post_spillover_dairy_farm_State(:,dd+1).*(1-p_no_onward_transmission.^(pp+dd));
            State_dairy_onward(:,mm)=State_dairy_onward(:,mm)+post_spillover_poultry_farm_State(:,pp+1).*post_spillover_dairy_farm_State(:,dd+1).*(1-p_no_onward_transmission.^(dd));
        end
    end
    County_dairy_onward(:,mm)=County_dairy_onward(:,mm)./County_onward(:,mm);
    State_dairy_onward(:,mm)=State_dairy_onward(:,mm)./State_onward(:,mm);

end

mle_County_onward(no_farms & no_diary_farms)=NaN;
County_onward(no_farms & no_diary_farms,:)=NaN;

min_County_onward=min(County_onward(~(no_farms & no_diary_farms),:),[],1);
max_County_onward=max(County_onward(~(no_farms & no_diary_farms),:),[],1);
median_County_onward=median(County_onward(~(no_farms & no_diary_farms),:),1);

min_County_onward_95=prctile(min_County_onward,[2.5 97.5]);
max_County_onward_95=prctile(max_County_onward,[2.5 97.5]);
median_County_onward_95=prctile(median_County_onward,[2.5 97.5]);

min_County_dairy_onward=min(County_dairy_onward(~(no_farms & no_diary_farms),:),[],1);
max_County_dairy_onward=max(County_dairy_onward(~(no_farms & no_diary_farms),:),[],1);
median_County_dairy_onward=max(County_dairy_onward(~(no_farms & no_diary_farms),:),[],1);

min_County_dairy_onward_95=prctile(min_County_dairy_onward,[2.5 97.5]);
max_County_dairy_onward_95=prctile(max_County_dairy_onward,[2.5 97.5]);
median_County_dairy_onward_95=prctile(median_County_dairy_onward,[2.5 97.5]);

County_onward_95=prctile(County_onward,[2.5 97.5],2);
State_onward_95=prctile(State_onward,[2.5 97.5],2);

County_dairy_onward_95=prctile(County_dairy_onward,[2.5 97.5],2);
State_dairy_onward_95=prctile(State_dairy_onward,[2.5 97.5],2);

save('Onward_Transmission.mat','mle_County_dairy_onward','County_dairy_onward_95','State_dairy_onward_95','mle_State_dairy_onward','min_County_dairy_onward_95','median_County_dairy_onward_95','max_County_dairy_onward_95','mle_County_onward','County_onward_95','State_onward_95','mle_State_onward','min_County_onward_95','median_County_onward_95','max_County_onward_95');