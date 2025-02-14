function F = Objective_Function_Dairy_Farm(x,X_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect)

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

phi_farm_County = Risk_Assesment_Farms(beta_x,X_County);
phi_farm_County(County_Farms==0)=0;

phi_farm_state=zeros(size(State_Spillover_Events));
for ss=1:length(phi_farm_state)
    temp_county=state_weight_matrix(ss,:).*log((1-phi_farm_County)); 
    temp_state=exp(sum(temp_county));
    phi_farm_state(ss)=1-temp_state;
end


L_County=log(nbinpdf(Affected_County_Farms(:),k_outbreak,1-phi_farm_County(:)));
L_County=L_County(County_Farms>0 & ~isnan(Affected_County_Farms) & County_Nonsuppressed);
L_State=log(nbinpdf(Affected_State_Farms(County_Suppressed_State),k_outbreak,1-phi_farm_state(County_Suppressed_State)));

Risk_Outbreak_State=1-nbinpdf(0,k_outbreak,1-phi_farm_state);
L_Spillover_State=log(nbinpdf(State_Spillover_Events(:),k_spillover,1-Risk_Outbreak_State(:)));

F=-sum(L_County)-sum(L_State)-sum(L_Spillover_State);
end

