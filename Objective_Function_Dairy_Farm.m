function F = Objective_Function_Dairy_Farm(x,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p)

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

delta_Affected_State=Affected_State_Farms-state_weight_matrix*Affected_County_Farms;
upper_County_Unknown=(delta_Affected_State')*state_weight_matrix;

p_inf_County=Zero_Inflation(beta_p,P_County);
p_inf_County(County_Farms==0)=1;

mu_farm_County = Risk_Assesment_Farms(beta_x,X_County);
mu_farm_County(County_Farms==0)=0;

mu_farm_state=zeros(size(State_Spillover_Events));
for ss=1:length(mu_farm_state)
    temp_county=(1-p_inf_County(:)).*mu_farm_County(:); 
    mu_farm_state(ss)=state_weight_matrix(ss,:)*temp_county;
end

k_outbreak_county=mu_farm_County(:);
k_outbreak_state=mu_farm_state(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_County=log((1-p_inf_County(:)).*poisspdf(Affected_County_Farms(:),k_outbreak_county(:)));
L_County(Affected_County_Farms==0)=log(p_inf_County(Affected_County_Farms==0)+(1-p_inf_County(Affected_County_Farms==0)).*poisspdf(Affected_County_Farms(Affected_County_Farms==0),k_outbreak_county(Affected_County_Farms==0)));

cdf=(1-p_inf_County(:)).*poisscdf(Affected_County_Farms(:)+upper_County_Unknown(:),k_outbreak_county(:))-(1-p_inf_County(:)).*poisscdf(Affected_County_Farms(:),k_outbreak_county(:));
L_County(upper_County_Unknown>0)=log(cdf(upper_County_Unknown>0));

L_County=L_County(County_Farms>0 & ~isnan(Affected_County_Farms));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_State=log(poisspdf(Affected_State_Farms(County_Suppressed_State),k_outbreak_state(County_Suppressed_State)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spillover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_Spillover_State=log(poisspdf(State_Spillover_Events(:),mu_farm_state(:).*kappa_spillover));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F=-sum(L_County)-sum(L_State)-sum(L_Spillover_State);
end

