function [par_est,L,AIC]=Optimize_Dairy_Farm_Risk(X_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect,x0)
lb=[-6 -18.*ones(1,size(X_County,1))  -4 -6];
ub=[ 0 -2.*ones(1,size(X_County,1)) 0 1];
options=optimoptions("surrogateopt","UseParallel",false,'MaxFunctionEvaluations',250.*length(lb),'PlotFcn',[]); 
[par_est]=surrogateopt(@(x)Objective_Function_Dairy_Farm(x,X_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect),lb,ub,options);

options=optimoptions('patternsearch','MaxFunctionEvaluations',250.*length(lb),'PlotFcn',[],'UseParallel',false,'FunctionTolerance',10^(-8),'MaxIterations',125.*length(lb),'MeshTolerance',10^(-9),'StepTolerance',10^(-9));
lb=[-7 -12.*ones(1,size(X_County,1))  -6 -8];
ub=[ 1.5   ones(1,size(X_County,1))   0  2];
[par_temp,f_temp]=patternsearch(@(x)Objective_Function_Dairy_Farm(x,X_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect),par_est,[],[],[],[],lb,ub,[],options);

lb=[-50 -32.*ones(1,size(X_County,1)) -32 -32];
ub=[ 3  4.*ones(1,size(X_County,1))   0  3];
options=optimoptions('fmincon','UseParallel',false,'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',10^4,'MaxIterations',10^4,'OptimalityTolerance',10^(-16),'StepTolerance',10^(-16));
[par_est,fval_final]=fmincon(@(x)Objective_Function_Dairy_Farm(x,X_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect),par_temp,[],[],[],[],lb,ub,[],options);


if(f_temp<fval_final)
    fval_final=f_temp;
    par_est=par_temp;
end

if(~isempty(x0))
    lb=[-50 -64.*ones(1,size(X_County,1)) -64 -64];
    ub=[ 3  5.*ones(1,size(X_County,1))   0    3];

    options=optimoptions("surrogateopt","UseParallel",false,'MaxFunctionEvaluations',100.*length(lb),'PlotFcn',[],'InitialPoints',x0);
    [par_temp]=surrogateopt(@(x)Objective_Function_Dairy_Farm(x,X_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect),lb,ub,options);

    options=optimoptions('patternsearch','MaxFunctionEvaluations',100.*length(lb),'PlotFcn',[],'UseParallel',false,'FunctionTolerance',10^(-8),'MaxIterations',125.*length(lb),'MeshTolerance',10^(-9),'StepTolerance',10^(-9));
    [par_temp,r_temp]=patternsearch(@(x)Objective_Function_Dairy_Farm(x,X_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect),par_temp,[],[],[],[],lb,ub,[],options);

    options=optimoptions('fmincon','UseParallel',false,'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',10^4,'MaxIterations',10^4,'OptimalityTolerance',10^(-16),'StepTolerance',10^(-16));
    [par_final,r_final]=fmincon(@(x)Objective_Function_Dairy_Farm(x,X_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect),par_temp,[],[],[],[],lb,ub,[],options);
    if(r_temp<r_final)
        r_final=r_temp;
        par_final=par_temp;    
    end
    if(r_final<fval_final)
        fval_final=r_final;
        par_est=par_final;
    end
end

L=-fval_final;
AIC=aicbic(L,length(lb));
end