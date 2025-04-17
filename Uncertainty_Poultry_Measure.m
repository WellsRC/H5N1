clear;
clc;

load('Poultry_Models_Refined_Fit.mat',"par_est","w_AIC",'Poultry_Model');

Poultry_Model=Poultry_Model(w_AIC==max(w_AIC),:);
par_est=par_est{w_AIC==max(w_AIC)};

[F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_par,logic_temperature] = Poultry_Covariates(Poultry_Model.Model_H5N1{1},Poultry_Model.Model_Farm{1});

if(logic_temperature)
    lb=[-7.*ones(1,size(F_County,1)) -7.*ones(1,size(F_County,1)) -16.*ones(1,size(P_County,1)-1) -0.5 -16.*ones(1,size(X_County,1))  -3];
    ub=[7.*ones(1,size(F_County,1)) 7.*ones(1,size(F_County,1))  log10(5).*ones(1,size(P_County,1)-1) 0.5 log10(5).*ones(1,size(X_County,1)) 0];
else
    lb=[-7.*ones(1,size(F_County,1)) -7.*ones(1,size(F_County,1)) -16.*ones(1,size(P_County,1)) -16.*ones(1,size(X_County,1)) -3];
    ub=[7.*ones(1,size(F_County,1)) 7.*ones(1,size(F_County,1))  log10(5).*ones(1,size(P_County,1)) log10(5).*ones(1,size(X_County,1)) 0];
end

L_samp=zeros(1+10^6,1);
par_samp=zeros(1+10^6,length(par_est));

L_samp(1)=Poultry_Model.L;
par_samp(1,:)=par_est;

dx=abs(par_est).*0.25;

x_lb=par_est-dx;
x_lb(x_lb<lb)=lb(x_lb<lb);
x_ub=par_est+dx;
x_ub(x_ub>ub)=ub(x_ub>ub);

par_samp(2:end,:)=repmat(x_lb,10^6,1)+repmat((x_ub-x_lb),10^6,1).*rand(10^6,length(par_est));

parfor ii=2:length(L_samp)
    L_samp(ii)=-Objective_Function_Poultry_Farm(par_samp(ii,:),F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature);
end

w=exp(L_samp-L_samp(1))./sum(exp(L_samp-L_samp(1)));

samp_indx=zeros(10^3,1);
r=rand(10^3,1);

wc=cumsum(w);

for jj=1:10^3
    samp_indx(jj)=find(r(jj)<=wc,1,"first");
end

par_mle=par_samp(L==max(L),:);
par_samp=par_samp(samp_indx,:);

save('Uncertainty_AIC_Poultry_Model.mat','par_mle','par_samp','Poultry_Model')
