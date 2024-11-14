function F = Objective_Function_Poultry_Farm(x,X_County,Y_County,N_County,A_County,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,Spillover,indx)

beta_x=x(1:(1+size(X_County,1)));
if(length(beta_x)>1)
    beta_x(2:end)=10.^beta_x(2:end);
end
beta_y=x((1+length(beta_x)):(end-1));
if(length(beta_y)>1)
    beta_y(2:end)=10.^beta_y(2:end);
end

N_County_temp=N_County-HPAI_Pullet_Farms-HPAI_Layer_Farms-HPAI_Turkey_Farms-HPAI_Broiler_Farms; % Need to remove the farms that are affected because this is the set of farms that we are uncertain about

r_nbin=10.^x(end);

r_farm_County = Risk_Assesment_Farms(beta_y,Y_County).*Risk_Assesment_Farms(beta_x,X_County);
r_farm_temp=r_farm_County;
r_farm_temp(N_County_temp==0)=0;
r_farm_County(N_County==0)=0;

L_County=A_County(:).*log(r_farm_temp(:))+(N_County_temp(:)-A_County(:)).*log(1-r_farm_temp(:));
L_County=L_County(N_County_temp>0);

r_Pullet_County = Risk_Assesment_Farms(beta_x([1 1+indx.Pullet]),X_County(indx.Pullet,:));
r_Pullet_County(Pullet_Farms==0)=0;
L_HPAI_Pullet=HPAI_Pullet_Farms(:).*log(r_Pullet_County(:))+(Pullet_Farms(:)-HPAI_Pullet_Farms(:)).*log(1-r_Pullet_County(:));
L_HPAI_Pullet=L_HPAI_Pullet(Pullet_Farms>0);

r_Layer_County = Risk_Assesment_Farms(beta_x([1 1+indx.Layer]),X_County(indx.Layer,:));
r_Layer_County(Layer_Farms==0)=0;
L_HPAI_Layer=HPAI_Layer_Farms(:).*log(r_Layer_County(:))+(Layer_Farms(:)-HPAI_Layer_Farms(:)).*log(1-r_Layer_County(:))
L_HPAI_Layer=L_HPAI_Layer(Layer_Farms>0);

r_Broiler_County = Risk_Assesment_Farms(beta_x([1 1+indx.Broiler]),X_County(indx.Broiler,:));
r_Broiler_County(Broiler_Farms==0)=0;
L_HPAI_Broiler=HPAI_Broiler_Farms(:).*log(r_Broiler_County(:))+(Broiler_Farms(:)-HPAI_Broiler_Farms(:)).*log(1-r_Broiler_County(:));
L_HPAI_Broiler=L_HPAI_Broiler(Broiler_Farms>0);

r_Turkey_County = Risk_Assesment_Farms(beta_x([1 1+indx.Turkey]),X_County(indx.Turkey,:));
r_Turkey_County(Turkey_Farms==0)=0;
L_HPAI_Turkey=HPAI_Turkey_Farms(:).*log(r_Turkey_County(:))+(Turkey_Farms(:)-HPAI_Turkey_Farms(:)).*log(1-r_Turkey_County(:));
L_HPAI_Turkey=L_HPAI_Turkey(Turkey_Farms>0);

L_Spillover=nbinpdf(Spillover,r_nbin,repmat(1-r_farm_County(:),1,size(Spillover,2)));
L_Spillover=mean(L_Spillover,2);


L_Spillover=L_Spillover(N_County>0);

F=-sum(L_County)-sum(L_Spillover) - sum(L_HPAI_Pullet(:)) - sum(L_HPAI_Layer(:)) - sum(L_HPAI_Broiler(:)) - sum(L_HPAI_Turkey(:));
end

