function F = Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature)

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


if(~isempty(Dairy_Network))
    beta_x_temp=beta_x([1:4 4+find(~logic_connect)]);
    
    mu_farm_temp = Risk_Assesment_Farms(beta_x_temp,[F_County; X_County(~logic_connect,:)]);
    mu_farm_temp(County_Farms==0)=0;

    beta_p_temp=beta_p([1:4 4+find(~logic_connect_p)]);
    p_inf_County_temp=Zero_Inflation(beta_p_temp,[F_County; P_County(~logic_connect_p,:)]);
    p_inf_County_temp(County_Farms==0)=1;

    temp_r=(1-p_inf_County_temp(:)).*mu_farm_temp(:);
    X_County(logic_connect,:)=(temp_r(:)'*Dairy_Network);
    P_County(logic_connect_p,:)=(temp_r(:)'*Dairy_Network);
end

delta_Affected_State=Affected_State_Farms-state_weight_matrix*Affected_County_Farms;
upper_County_Unknown=(delta_Affected_State')*state_weight_matrix;
upper_County_Unknown=upper_County_Unknown(:);

p_inf_County=Zero_Inflation(beta_p,[F_County; P_County]);
p_inf_County(County_Farms==0)=1;

mu_farm_County = Risk_Assesment_Farms(beta_x,[F_County; X_County]);
mu_farm_County(County_Farms==0)=0;

k_state=zeros(size(State_Spillover_Events));
p_nb_state=zeros(size(State_Spillover_Events));
z_state=zeros(size(State_Spillover_Events));

k_state_spill=zeros(size(State_Spillover_Events));
z_state_spill=zeros(size(State_Spillover_Events));
p_nb_state_spill=zeros(size(State_Spillover_Events));

p_temp=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:));
p_temp_spill=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County(:)); % Use the "TRUE" county risk as aspects are not observed if lack of surviellance 

w_state=zeros(size(State_Spillover_Events));
County_w_Farms=County_Farms>0;
L_Nan=false;
for ss=1:length(w_state)

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

    L_Nan=L_Nan |L_Nan1 | L_Nan2;
    Affected_State_Farms(ss)=state_weight_matrix(ss,:)*(Affected_County_Farms(:));
    if(L_Nan)
        break;
    end
end
if (~L_Nan)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % County
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k_outbreak=mu_farm_County(:);
    L_County=log((1-p_inf_County(:)).*poisspdf(Affected_County_Farms(:),k_outbreak(:)));
    L_County(Affected_County_Farms==0)=log(p_inf_County(Affected_County_Farms==0)+(1-p_inf_County(Affected_County_Farms==0)).*poisspdf(Affected_County_Farms(Affected_County_Farms==0),k_outbreak(Affected_County_Farms==0)));
    
    cdf=(1-p_inf_County(:)).*poisscdf(Affected_County_Farms(:)+upper_County_Unknown(:),k_outbreak(:))-(1-p_inf_County(:)).*poisscdf(Affected_County_Farms(:),k_outbreak(:)); % do not need the addition of p_inf_County alone since it cancels out in subtration because computing cdf in each 
    L_County(upper_County_Unknown>0)=log(cdf(upper_County_Unknown>0));
    
    L_County=L_County(County_Farms>0 & ~isnan(Affected_County_Farms));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % State outbreaks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    L_State=log((1-z_state(:)).*nbinpdf(Affected_State_Farms(:),k_state(:),p_nb_state(:)));
    L_State(Affected_State_Farms==0)=log(z_state(Affected_State_Farms==0)+(1-z_state(Affected_State_Farms==0)).*nbinpdf(Affected_State_Farms(Affected_State_Farms==0),k_state(Affected_State_Farms==0),p_nb_state(Affected_State_Farms==0)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spillover
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    L_Spillover_State=log((1-z_state_spill(:)).*nbinpdf(State_Spillover_Events(:),k_state_spill(:),p_nb_state_spill(:)));
    L_Spillover_State(State_Spillover_Events(:)==0)=log(z_state_spill(State_Spillover_Events(:)==0)+(1-z_state_spill(State_Spillover_Events(:)==0)).*nbinpdf(State_Spillover_Events(State_Spillover_Events(:)==0),k_state_spill(State_Spillover_Events(:)==0),p_nb_state_spill(State_Spillover_Events(:)==0)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Objective function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F=-sum(L_County) -sum(w_state(:).*L_Spillover_State(:)) -sum(w_state(:).*L_State(:));
else
    F=NaN;
end
end

