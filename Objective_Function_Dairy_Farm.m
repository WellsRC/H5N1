function F = Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,State_Spillover_Matrix,State_Spillover_Events,Remainaing_Affected_State_Farms,Remainaing_Total_State_Farms,state_weight_hpai_matrix,Dairy_Network)

beta_x=x(1:(1+size(X_County,1)));
if(length(beta_x)>1)
    beta_x(2:end)=10.^beta_x(2:end);
end
beta_y=x((1+length(beta_x)):(end-1));
if(length(beta_y)>1)
    beta_y(2:end)=10.^beta_y(2:end);
end

r_nbin=10.^x(end);


t_connect=sum(isinf(X_County),2)==size(X_County,2);
if(~isempty(Dairy_Network))
    beta_x_temp=beta_x([1 1+find(~t_connect)']);
    r_farm_temp = Risk_Assesment_Farms(beta_y,Y_County).*Risk_Assesment_Farms(beta_x_temp,X_County(~t_connect,:));
    r_farm_temp(County_Farms==0)=0;
    X_County(t_connect,:)=(r_farm_temp(:)'*Dairy_Network);
end

r_farm_County = Risk_Assesment_Farms(beta_y,Y_County).*Risk_Assesment_Farms(beta_x,X_County);
r_farm_County(County_Farms==0)=0;

r_farm_state=zeros(size(Remainaing_Affected_State_Farms));
for ss=1:length(r_farm_state)
    temp_county=state_weight_hpai_matrix(ss,:).*log((1-r_farm_County)); 
    temp_state=exp(sum(temp_county));
    r_farm_state(ss)=1-temp_state;
end

L_Spillover_State=zeros(size(State_Spillover_Events));
for ss=1:length(State_Spillover_Events)
    temp_county=State_Spillover_Matrix(ss,:).*log((1-r_farm_County)); 
    temp_state=exp(sum(temp_county));
    r_spillover_state=1-temp_state;
    L_Spillover_State(ss)=nbinpdf(round(State_Spillover_Events(ss)),r_nbin,1-r_spillover_state);
end

L_County=Affected_County_Farms(:).*log(r_farm_County(:))+(County_Farms(:)-Affected_County_Farms(:)).*log(1-r_farm_County(:));
L_County=L_County(County_Farms>0 & ~isnan(Affected_County_Farms));
L_State=Remainaing_Affected_State_Farms(:).*log(r_farm_state(:))+(Remainaing_Total_State_Farms(:)-Remainaing_Affected_State_Farms(:)).*log(1-r_farm_state(:));


F=-sum(L_County)-sum(L_State)-sum(L_Spillover_State);
end

