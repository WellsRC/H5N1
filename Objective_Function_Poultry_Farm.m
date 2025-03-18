function F = Objective_Function_Poultry_Farm(x,X_County,P_County,County_Farms,Affected_County_Farms_Unknown,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,state_weight_matrix,State_Spillover_Events)


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total outbreaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_HPAI_Pullet_Farms=HPAI_Pullet_Farms;
temp_HPAI_Pullet_Farms(isnan(temp_HPAI_Pullet_Farms))=0;

temp_HPAI_Layer_Farms=HPAI_Layer_Farms;
temp_HPAI_Layer_Farms(isnan(temp_HPAI_Layer_Farms))=0;

temp_HPAI_Broiler_Farms=HPAI_Broiler_Farms;
temp_HPAI_Broiler_Farms(isnan(temp_HPAI_Broiler_Farms))=0;

temp_HPAI_Turkey_Farms=HPAI_Turkey_Farms;
temp_HPAI_Turkey_Farms(isnan(temp_HPAI_Turkey_Farms))=0;

Tot_Affected_Farms=Affected_County_Farms_Unknown(:);
Tot_Affected_Farms(isnan(Tot_Affected_Farms))=0;

Tot_Affected_Farms=Tot_Affected_Farms+temp_HPAI_Pullet_Farms+temp_HPAI_Layer_Farms+temp_HPAI_Broiler_Farms+temp_HPAI_Turkey_Farms;

k_outbreak=mu_farm_County_Total(:);
L_County=log((1-p_inf_County(:)).*poisspdf(Tot_Affected_Farms(:),k_outbreak(:)));
L_County(Tot_Affected_Farms==0)=log(p_inf_County(Tot_Affected_Farms==0)+(1-p_inf_County(Tot_Affected_Farms==0)).*poisspdf(Tot_Affected_Farms(Tot_Affected_Farms==0),k_outbreak(Tot_Affected_Farms==0)));

L_County=L_County(County_Farms>0 & ~isnan(k_outbreak));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pullet farms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp_unknown=Affected_County_Farms_Unknown(:);
temp_unknown(isnan(temp_unknown) | Pullet_Farms==0)=0;
k_outbreak=mu_farm_County_Pullet(:);

L_HPAI_Pullet=log((1-p_inf_County_Pullet(:)).*poisspdf(HPAI_Pullet_Farms(:),k_outbreak(:)));
t_frac=HPAI_Pullet_Farms(:) ~= round(HPAI_Pullet_Farms(:));
L_HPAI_Pullet(t_frac)=log((1-p_inf_County_Pullet(t_frac)).*poisscdf(ceil(HPAI_Pullet_Farms(t_frac))+temp_unknown(t_frac),k_outbreak(t_frac))-(1-p_inf_County_Pullet(t_frac)).*poisscdf(floor(HPAI_Pullet_Farms(t_frac)),k_outbreak(t_frac)));
L_HPAI_Pullet(HPAI_Pullet_Farms==0)=log(p_inf_County_Pullet(HPAI_Pullet_Farms==0)+(1-p_inf_County_Pullet(HPAI_Pullet_Farms==0)).*poisspdf(HPAI_Pullet_Farms(HPAI_Pullet_Farms==0),k_outbreak(HPAI_Pullet_Farms==0)));

cdf=(1-p_inf_County_Pullet(:)).*poisscdf(HPAI_Pullet_Farms(:)+temp_unknown(:),k_outbreak(:))-((1-p_inf_County_Pullet(:)).*poisscdf(HPAI_Pullet_Farms(:),k_outbreak(:)));
L_HPAI_Pullet(temp_unknown>0 & ~t_frac)=log(cdf(temp_unknown>0 & ~t_frac)); 

L_HPAI_Pullet=L_HPAI_Pullet(Pullet_Farms>0 & ~isnan(HPAI_Pullet_Farms));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer farms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_unknown=Affected_County_Farms_Unknown(:);
temp_unknown(isnan(temp_unknown) | Layer_Farms==0)=0;
k_outbreak=mu_farm_County_Layer(:);

L_HPAI_Layer=log((1-p_inf_County_Layer(:)).*poisspdf(HPAI_Layer_Farms(:),k_outbreak(:)));
L_HPAI_Layer(HPAI_Layer_Farms==0)=log(p_inf_County_Layer(HPAI_Layer_Farms==0)+(1-p_inf_County_Layer(HPAI_Layer_Farms==0)).*poisspdf(HPAI_Layer_Farms(HPAI_Layer_Farms==0),k_outbreak(HPAI_Layer_Farms==0)));

cdf=(1-p_inf_County_Layer(:)).*poisscdf(HPAI_Layer_Farms(:)+temp_unknown(:),k_outbreak(:))-((1-p_inf_County_Layer(:)).*poisscdf(HPAI_Layer_Farms(:),k_outbreak(:)));
L_HPAI_Layer(temp_unknown>0)=log(cdf(temp_unknown>0)); 

L_HPAI_Layer=L_HPAI_Layer(Layer_Farms>0 & ~isnan(HPAI_Layer_Farms));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Broiler farms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp_unknown=Affected_County_Farms_Unknown(:);
temp_unknown(isnan(temp_unknown) | Broiler_Farms==0)=0;
k_outbreak=mu_farm_County_Broiler(:);

L_HPAI_Broiler=log((1-p_inf_County_Broiler(:)).*poisspdf(HPAI_Broiler_Farms(:),k_outbreak(:)));
t_frac=HPAI_Broiler_Farms(:) ~= round(HPAI_Broiler_Farms(:));
L_HPAI_Broiler(t_frac)=log((1-p_inf_County_Broiler(t_frac)).*poisscdf(ceil(HPAI_Broiler_Farms(t_frac))+temp_unknown(t_frac),k_outbreak(t_frac))-(1-p_inf_County_Broiler(t_frac)).*poisscdf(floor(HPAI_Broiler_Farms(t_frac)),k_outbreak(t_frac)));
L_HPAI_Broiler(HPAI_Broiler_Farms==0)=log(p_inf_County_Broiler(HPAI_Broiler_Farms==0)+(1-p_inf_County_Broiler(HPAI_Broiler_Farms==0)).*poisspdf(HPAI_Broiler_Farms(HPAI_Broiler_Farms==0),k_outbreak(HPAI_Broiler_Farms==0)));

cdf=(1-p_inf_County_Broiler(:)).*poisscdf(HPAI_Broiler_Farms(:)+temp_unknown(:),k_outbreak(:))-((1-p_inf_County_Broiler(:)).*poisscdf(HPAI_Broiler_Farms(:),k_outbreak(:)));
L_HPAI_Broiler(temp_unknown>0 & ~t_frac)=log(cdf(temp_unknown>0 & ~t_frac)); 
L_HPAI_Broiler=L_HPAI_Broiler(Broiler_Farms>0 & ~isnan(HPAI_Broiler_Farms));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turkey farms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_unknown=Affected_County_Farms_Unknown(:);
temp_unknown(isnan(temp_unknown) | Turkey_Farms==0)=0;
k_outbreak=mu_farm_County_Turkey(:);
L_HPAI_Turkey=log((1-p_inf_County_Turkey(:)).*poisspdf(HPAI_Turkey_Farms(:),k_outbreak(:)));
L_HPAI_Turkey(HPAI_Turkey_Farms==0)=log(p_inf_County_Turkey(HPAI_Turkey_Farms==0)+(1-p_inf_County_Turkey(HPAI_Turkey_Farms==0)).*poisspdf(HPAI_Turkey_Farms(HPAI_Turkey_Farms==0),k_outbreak(HPAI_Turkey_Farms==0)));

cdf=(1-p_inf_County_Turkey(:)).*poisscdf(HPAI_Turkey_Farms(:)+temp_unknown(:),k_outbreak(:))-((1-p_inf_County_Turkey(:)).*poisscdf(HPAI_Turkey_Farms(:),k_outbreak(:)));
L_HPAI_Layer(temp_unknown>0)=log(cdf(temp_unknown>0)); 

L_HPAI_Turkey=L_HPAI_Turkey(Turkey_Farms>0 & ~isnan(HPAI_Turkey_Farms));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spillover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_Spillover_State=log(poisspdf(State_Spillover_Events(:),mu_farm_State_Total(:).*kappa_spillover));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F=-sum(L_County) -sum(L_Spillover_State) -sum(L_HPAI_Pullet(:)) -sum(L_HPAI_Layer(:)) -sum(L_HPAI_Broiler(:)) -sum(L_HPAI_Turkey(:));
end

