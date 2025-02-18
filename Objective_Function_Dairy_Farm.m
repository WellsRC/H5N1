function F = Objective_Function_Dairy_Farm(x,X_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect)

beta_x=x(1:(1+size(X_County,1)));
if(length(beta_x)>1)
    beta_x(2:end)=10.^beta_x(2:end);
end
p_outbreak=10.^x(end-1);
kappa_spillover=10.^x(end);


if(~isempty(Dairy_Network))
    beta_x_temp=beta_x([1 1+find(~logic_connect)']);
    mu_farm_temp = Risk_Assesment_Farms(beta_x_temp,X_County(~logic_connect,:));
    mu_farm_temp(County_Farms==0)=0;
    X_County(logic_connect,:)=X_County(logic_connect,:).*repmat((mu_farm_temp(:)'*Dairy_Network),sum(logic_connect),1);
end

mu_farm_County = Risk_Assesment_Farms(beta_x,X_County);
mu_farm_County(County_Farms==0)=0;

mu_farm_state=zeros(size(State_Spillover_Events));
for ss=1:length(mu_farm_state)
    temp_county=state_weight_matrix(ss,:).*mu_farm_County; 
    mu_farm_state(ss)=sum(temp_county);
end

k_outbreak_county=mu_farm_County.*p_outbreak./(1-p_outbreak);
k_outbreak_state=mu_farm_state.*p_outbreak./(1-p_outbreak);

L_County=log(nbinpdf(Affected_County_Farms(:),k_outbreak_county(:),p_outbreak));
L_County=L_County(County_Farms>0 & ~isnan(Affected_County_Farms) & County_Nonsuppressed);
L_State=log(nbinpdf(Affected_State_Farms(County_Suppressed_State),k_outbreak_state(County_Suppressed_State),p_outbreak));

L_Spillover_State=log(poisspdf(State_Spillover_Events(:),mu_farm_state(:).*kappa_spillover));

F=-sum(L_County)-sum(L_State)-sum(L_Spillover_State);
end

