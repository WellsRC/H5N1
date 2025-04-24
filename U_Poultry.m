% parpool(48)
load('Poultry_Models_Refined_Fit.mat',"par_est","w_AIC",'Poultry_Model');

Poultry_Model=Poultry_Model(w_AIC==max(w_AIC),:);

L_mle=Poultry_Model.L;
par_mle=par_est{w_AIC==max(w_AIC)};

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
                Lt=-Objective_Function_Poultry_Farm(p_new,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature);
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
save('Uncertainty_AIC_Poultry_Model.mat','L_samp','par_samp','Poultry_Model','par_mle','L_mle')
