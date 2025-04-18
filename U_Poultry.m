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


x_lb=[-1 -0.9 -0.04 -0.44 1.114 1.2298 1.095 1.3719 0.1273 -1.8334 -0.4569 -2];
x_ub=[-0.7 -0.75 -0.015 -0.39 1.117 1.2302 1.0985 1.3729 0.1274 -1.8327 -0.4567 -1.5];
L_samp=zeros(1,12*200);
par_samp=0.*repmat(par_est,12.*200,1);
opt=optimoptions('patternsearch','UseParallel',true,'FunctionTolerance',10^(-8),'PlotFcn',[],'MeshTolerance',10^(-12));
for indx=1:length(x_lb)    
    x0=flip(linspace(x_lb(indx),par_est(indx),101));
    x0=x0(2:end);
    for nn=1:100
        if(nn==1)
            p0=par_est(~ismember([1:length(par_est)],indx));
        else
            p0=par_temp;
        end
        [par_temp,L_samp(nn+200.*(indx-1))]=patternsearch(@(y)Objective_Function_Poultry_Farm_LP(indx,x0(nn),y,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),p0,[],[],[],[],lb(~ismember([1:length(par_est)],indx)),ub(~ismember([1:length(par_est)],indx)),[],opt);
        par_samp(nn+200.*(indx-1),indx)=x0(nn);
        par_samp(nn+200.*(indx-1),~ismember([1:length(par_est)],indx))=par_temp;
    end

    x0=linspace(par_est(indx),x_ub(indx),101);
    x0=x0(2:end);
    for nn=101:200
        if(nn==1)
            p0=par_est(~ismember([1:length(par_est)],indx));
        else
            p0=par_temp;
        end
        [par_temp,L_samp(nn+200.*(indx-1))]=patternsearch(@(y)Objective_Function_Poultry_Farm_LP(indx,x0(nn-100),y,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),p0,[],[],[],[],lb(~ismember([1:length(par_est)],indx)),ub(~ismember([1:length(par_est)],indx)),[],opt);
        par_samp(nn+200.*(indx-1),indx)=x0(nn-100);
        par_samp(nn+200.*(indx-1),~ismember([1:length(par_est)],indx))=par_temp;
    end
end

w=exp(L_samp-L_samp(1))./sum(exp(L_samp-L_samp(1)));

samp_indx=zeros(10^3,1);
r=rand(10^3,1);

wc=cumsum(w);

for jj=1:10^3
    samp_indx(jj)=find(r(jj)<=wc,1,"first");
end

par_mle=par_samp(L_samp==max(L_samp),:);
par_samp=par_samp(samp_indx,:);

save('Uncertainty_AIC_Poultry_Model.mat','par_mle','par_samp','Poultry_Model')
