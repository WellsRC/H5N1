function F = Objective_Function_Poultry_Farm(x,X_County,County_Farms,Affected_County_Farms_Unknown,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,state_weight_matrix,State_Spillover_Events)

beta_x=x(1:(1+size(X_County,1)));
if(length(beta_x)>1)
    beta_x(2:end)=10.^beta_x(2:end);
end

f_outbreak_turkey=10.^x(end-4);
f_outbreak_pullet=10.^x(end-3);
f_outbreak_broiler=10.^x(end-2);
p_outbreak=10.^x(end-1);
kappa_spillover=10.^x(end);

mu_farm_County_Basline = Risk_Assesment_Farms(beta_x,X_County);
mu_farm_County_Layer = Risk_Assesment_Farms(beta_x,X_County);

mu_farm_County_Turkey = f_outbreak_turkey.*mu_farm_County_Basline;
mu_farm_County_Pullet = f_outbreak_pullet.*mu_farm_County_Basline;
mu_farm_County_Broiler = f_outbreak_broiler.*mu_farm_County_Basline;

mu_farm_County_Basline(County_Farms==0)=0;
dX=County_Farms-Layer_Farms-Turkey_Farms-Pullet_Farms-Broiler_Farms;
mu_farm_County_Basline(dX<=0)=0;
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
    mu_farm_State_Baseline(ss)=sum(state_weight_matrix(ss,:).*mu_farm_County_Basline);
    mu_farm_State_Layer(ss)=sum(state_weight_matrix(ss,:).*mu_farm_County_Layer);
    mu_farm_State_Turkey(ss)=sum(state_weight_matrix(ss,:).*mu_farm_County_Turkey);
    mu_farm_State_Pullet(ss)=sum(state_weight_matrix(ss,:).*mu_farm_County_Pullet);
    mu_farm_State_Broiler(ss)=sum(state_weight_matrix(ss,:).*mu_farm_County_Broiler);
end

mu_farm_County_Total=mu_farm_County_Layer+mu_farm_County_Turkey+mu_farm_County_Pullet+mu_farm_County_Broiler+mu_farm_County_Basline;
mu_farm_State_Total=mu_farm_State_Layer+mu_farm_State_Turkey+mu_farm_State_Pullet+mu_farm_State_Broiler+mu_farm_State_Baseline;


k_outbreak=mu_farm_County_Total(:).*p_outbreak./(1-p_outbreak);
L_County=log(nbinpdf(Affected_County_Farms_Unknown(:),k_outbreak(:),p_outbreak));
L_County=L_County(County_Farms>0 & ~isnan(Affected_County_Farms_Unknown) & ~isnan(k_outbreak));

k_outbreak=mu_farm_County_Pullet(:).*p_outbreak./(1-p_outbreak);
L_HPAI_Pullet=log(nbinpdf(HPAI_Pullet_Farms(:),k_outbreak(:),p_outbreak));
t_frac=HPAI_Pullet_Farms(:) ~= round(HPAI_Pullet_Farms(:));
L_HPAI_Pullet(t_frac)=log(nbincdf(ceil(HPAI_Pullet_Farms(t_frac)),k_outbreak(t_frac),p_outbreak)-nbincdf(floor(HPAI_Pullet_Farms(t_frac)),k_outbreak(t_frac),p_outbreak));
L_HPAI_Pullet=L_HPAI_Pullet(Pullet_Farms>0 & ~isnan(HPAI_Pullet_Farms));

k_outbreak=mu_farm_County_Layer(:).*p_outbreak./(1-p_outbreak);
L_HPAI_Layer=log(nbinpdf(HPAI_Layer_Farms(:),k_outbreak(:),p_outbreak));
L_HPAI_Layer=L_HPAI_Layer(Layer_Farms>0 & ~isnan(HPAI_Layer_Farms));

k_outbreak=mu_farm_County_Broiler(:).*p_outbreak./(1-p_outbreak);
L_HPAI_Broiler=log(nbinpdf(HPAI_Broiler_Farms(:),k_outbreak(:),p_outbreak));

t_frac=HPAI_Broiler_Farms(:) ~= round(HPAI_Broiler_Farms(:));
L_HPAI_Broiler(t_frac)=log(nbincdf(ceil(HPAI_Broiler_Farms(t_frac)),k_outbreak(t_frac),p_outbreak)-nbincdf(floor(HPAI_Broiler_Farms(t_frac)),k_outbreak(t_frac),p_outbreak));
L_HPAI_Broiler=L_HPAI_Broiler(Broiler_Farms>0 & ~isnan(HPAI_Broiler_Farms));

k_outbreak=mu_farm_County_Turkey(:).*p_outbreak./(1-p_outbreak);
L_HPAI_Turkey=log(nbinpdf(HPAI_Turkey_Farms(:),k_outbreak(:),p_outbreak));
L_HPAI_Turkey=L_HPAI_Turkey(Turkey_Farms>0 & ~isnan(HPAI_Turkey_Farms));

L_Spillover_State=log(poisspdf(State_Spillover_Events(:),mu_farm_State_Total(:).*kappa_spillover));

F=-sum(L_County) -sum(L_Spillover_State) -sum(L_HPAI_Pullet(:)) -sum(L_HPAI_Layer(:)) -sum(L_HPAI_Broiler(:)) -sum(L_HPAI_Turkey(:));
end

