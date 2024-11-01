function F = Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network)

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
    temp_county=nthroot((1-r_farm_County).^(Remainaing_State_Farms(ss).*state_weight_matrix(ss,:)),Remainaing_State_Farms(ss));
    r_farm_state(ss)=1-prod(temp_county);
end

L_Spillover=nbinpdf(County_Spillover,r_nbin,repmat(1-r_farm_County(:),1,size(County_Spillover,2)));
L_Spillover=mean(L_Spillover,2);

L_County=Affected_County_Farms(:).*log(r_farm_County(:))+(County_Farms(:)-Affected_County_Farms(:)).*log(1-r_farm_County(:));

L_County=L_County(County_Farms>0 & ~isnan(Affected_County_Farms));
L_Spillover=L_Spillover(County_Farms>0 & ~isnan(Affected_County_Farms));


L_State=Remainaing_Affected_State_Farms(:).*log(r_farm_state(:))+(Remainaing_State_Farms(:)-Remainaing_Affected_State_Farms(:)).*log(1-r_farm_state(:));
L_Spillover_State=nbinpdf(Remaining_State_Spillover,r_nbin,repmat(1-r_farm_state(:),1,size(Remaining_State_Spillover,2)));
L_Spillover_State=mean(L_Spillover_State,2);

F=-sum(L_County)-sum(L_Spillover)-sum(L_State)-sum(L_Spillover_State);
end

