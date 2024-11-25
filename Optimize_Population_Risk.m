function [par_est,L,AIC]=Optimize_Population_Risk(X_County,Y_County,Cases_County,Population,x0)

lb=[-10 -10.*ones(1,size(X_County,1))  -5 -10.*ones(1,size(Y_County,1)) -2];
ub=[ 0  10.*ones(1,size(X_County,1))  -2  10.*ones(1,size(Y_County,1))  1];
options=optimoptions("surrogateopt","UseParallel",false,'MaxFunctionEvaluations',250.*length(lb),'PlotFcn',[]); 
[par_est_old]=surrogateopt(@(x)Objective_Function_Population(x,X_County,Y_County,Cases_County,Population),lb,ub,options);

lb=[-15 -50.*ones(1,size(X_County,1))  -7 -50.*ones(1,size(Y_County,1)) -5];
ub=[ 10  50.*ones(1,size(X_County,1))   3  50.*ones(1,size(Y_County,1))  2];
options=optimoptions('patternsearch','MaxFunctionEvaluations',250.*length(lb),'PlotFcn',[],'UseParallel',false,'FunctionTolerance',10^(-8),'MaxIterations',125.*length(lb),'MeshTolerance',10^(-9),'StepTolerance',10^(-9));
[par_temp,f_temp]=patternsearch(@(x)Objective_Function_Population(x,X_County,Y_County,Cases_County,Population),par_est_old,[],[],[],[],lb,ub,[],options);

lb=[-25 -100.*ones(1,size(X_County,1))  -16 -100.*ones(1,size(Y_County,1)) -32];
ub=[ 25  100.*ones(1,size(X_County,1))   10  100.*ones(1,size(Y_County,1))  3];
options=optimoptions('fmincon','UseParallel',false,'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',10^4,'MaxIterations',10^4,'OptimalityTolerance',10^(-16),'StepTolerance',10^(-16));
[par_est,fval_final]=fmincon(@(x)Objective_Function_Population(x,X_County,Y_County,Cases_County,Population),par_temp,[],[],[],[],lb,ub,[],options);

if(f_temp<fval_final)
    fval_final=f_temp;
    par_est=par_temp;
end

if(~isempty(x0))
    lb=[-50 -125.*ones(1,size(X_County,1))  -50 -125.*ones(1,size(Y_County,1)) -64];
    ub=[ 50  125.*ones(1,size(X_County,1))   50  125.*ones(1,size(Y_County,1)) 4];

    options=optimoptions("surrogateopt","UseParallel",false,'MaxFunctionEvaluations',100.*length(lb),'PlotFcn',[],'InitialPoints',x0);
    [par_temp]=surrogateopt(@(x)Objective_Function_Population(x,X_County,Y_County,Cases_County,Population),lb,ub,options);

    options=optimoptions('patternsearch','MaxFunctionEvaluations',100.*length(lb),'PlotFcn',[],'UseParallel',false,'FunctionTolerance',10^(-8),'MaxIterations',125.*length(lb),'MeshTolerance',10^(-9),'StepTolerance',10^(-9));
    [par_temp,r_temp]=patternsearch(@(x)Objective_Function_Population(x,X_County,Y_County,Cases_County,Population),par_temp,[],[],[],[],lb,ub,[],options);

    options=optimoptions('fmincon','UseParallel',false,'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',10^4,'MaxIterations',10^4,'OptimalityTolerance',10^(-16),'StepTolerance',10^(-16));
    [par_final,r_final]=fmincon(@(x)Objective_Function_Population(x,X_County,Y_County,Cases_County,Population),par_temp,[],[],[],[],lb,ub,[],options);
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