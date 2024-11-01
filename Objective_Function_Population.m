function F = Objective_Function_Population(x,X_County,Y_County,Cases_County,Population)

beta_x=x(1:(1+size(X_County,1)));

beta_y=x((1+length(beta_x)):(end-1));

sigma_ln=10.^x(end);

r_County = Risk_Assesment_Farms(beta_y,Y_County).*Risk_Assesment_Farms(beta_x,X_County);
r_County=r_County(:);
z_temp=Cases_County(:)./Population(:);
z_temp=log(z_temp./(1-z_temp));
z_censor=1./Population(:);
z_censor=log(z_censor./(1-z_censor));

m_temp=log(r_County(:)./(1-r_County(:)));

L_County=log(normpdf(z_temp(:),m_temp(:),sigma_ln));
L_County(Cases_County==0)=log(normcdf(real(z_censor(Cases_County==0)),real(m_temp(Cases_County==0)),sigma_ln));

F=-sum(L_County);
if(isnan(F))
    F=10^(20);
end
end

