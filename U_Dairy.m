% parpool(48)
load('Dairy_Models_Refined_Fit.mat',"par_est","w_AIC",'Dairy_Model');

Dairy_Model=Dairy_Model(w_AIC==max(w_AIC),:);

L_mle=Dairy_Model.L;
par_mle=par_est{w_AIC==max(w_AIC)};

[F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_par,logic_temperature]= Dairy_Covariates(Dairy_Model.Model_H5N1{1},Dairy_Model.Model_Farm{1},Dairy_Model.Model_Stratified_Operations{1});


if(logic_temperature && isempty(Dairy_Network))
    lb=[-7.*ones(1,size(F_County,1)) -7.*ones(1,size(F_County,1)) -16.*ones(1,size(P_County,1)-1) -0.5 -16.*ones(1,size(X_County,1))  -3];
    ub=[7.*ones(1,size(F_County,1)) 7.*ones(1,size(F_County,1))  log10(5).*ones(1,size(P_County,1)-1) 0.5 log10(5).*ones(1,size(X_County,1)) 0];
elseif(logic_temperature && ~isempty(Dairy_Network))
    lb=[-7.*ones(1,size(F_County,1)) -7.*ones(1,size(F_County,1)) -16.*ones(1,size(P_County,1)-2) -0.5 -16 -16.*ones(1,size(X_County,1))  -3];
    ub=[7.*ones(1,size(F_County,1)) 7.*ones(1,size(F_County,1))  log10(5).*ones(1,size(P_County,1)-2) 0.5 1 log10(5).*ones(1,size(X_County,1)) 0];
else
    lb=[-7.*ones(1,size(F_County,1)) -7.*ones(1,size(F_County,1)) -16.*ones(1,size(P_County,1)) -16.*ones(1,size(X_County,1)) -3];
    ub=[7.*ones(1,size(F_County,1)) 7.*ones(1,size(F_County,1))  log10(5).*ones(1,size(P_County,1)) log10(5).*ones(1,size(X_County,1)) 0];
end


x_lb=[-2.3 -1.1 0.05 -0.035 0.5 2.4 0.53 9.806*10^(-3) -0.2676436 -5.66 -6.45 -1.76897 -16 -5.18 -16 -16 -1.2010021 -1.4827886 -1.95];
x_ub=[-2 -0.8 0.18 -0.03 0.97 2.5 0.62 9.8067*10^(-3) -0.26764335 -5.57 -5.5 -1.76894 -5 -5.08 -5.5 -6.5 -1.2010014 -1.4827874 -1.57];


dx_lb=par_mle-x_lb;
dx_ub=x_ub-par_mle;

dx=2.*min(dx_ub,dx_lb);
N_Samp=5000;
L_samp=zeros(N_Samp,1);
par_samp=zeros(N_Samp,length(x_lb));
cc=0;
while(cc<N_Samp)
    p_old=par_mle;
    burn_in=50;
    for kk=1:1000
        rv=randperm(length(x_lb));
        for zz=1:length(x_lb)
            r_int=rv(zz);
            dp=dx(r_int).*(rand(1)-0.5);
            p_new=p_old;
            p_new(r_int)=p_new(r_int)+dp;
            if(p_new(r_int)>=lb(r_int) && p_new(r_int)<=ub(r_int))              
                Lt=-Objective_Function_Dairy_Farm(p_new,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature);
                if(rand(1)<=exp(Lt-L_mle)/(exp(Lt-L_mle)+1)) && burn_in<=0
                    cc=cc+1;
                    L_samp(cc)=Lt;
                    par_samp(cc,:)=p_new;
                    p_old=p_new;
                    if(cc==N_Samp)
                        break;
                    end
                elseif(rand(1)<=exp(Lt-L_mle)/(exp(Lt-L_mle)+1))
                    p_old=p_new;
                    burn_in=burn_in-1;
                end
            end
        end
        if(cc==N_Samp)
            break;
        end
    end

end
save('Uncertainty_AIC_Dairy_Model.mat','L_samp','par_samp','Dairy_Model','par_mle','L_mle')



