clear;
clc;

load('Poultry_Models_Fit.mat',"par_est","w_AIC",'Poultry_Model','L','AIC');

mm=find(w_AIC==max(w_AIC));
mle_L=L(mm);
mle_AIC=AIC(mm);
[F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_par] = Poultry_Covariates(Poultry_Model.Model_H5N1{mm},Poultry_Model.Model_Farm{mm},Poultry_Model.Model_Stratified_Chicken_Inventory{mm});
x=par_est{mm}; 

x(x<-16)=-Inf;
L_test=-Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events);
AIC_test=aicbic(L_test,sum(~isinf(x)));