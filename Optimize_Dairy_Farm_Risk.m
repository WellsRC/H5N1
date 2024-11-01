function [par_est,L,AIC]=Optimize_Dairy_Farm_Risk(X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network,x0)

if(isempty(x0))
    lb=[-6 -4.*ones(1,size(X_County,1))  -5 -3.*ones(1,size(Y_County,1)) -2];
    ub=[ 0  zeros(1,size(X_County,1))    -2  zeros(1,size(Y_County,1)) 2];
    options=optimoptions("surrogateopt","UseParallel",false,'MaxFunctionEvaluations',250.*length(lb),'PlotFcn',[]); 
    [par_est]=surrogateopt(@(x)Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network),lb,ub,options);
    
    options=optimoptions('patternsearch','MaxFunctionEvaluations',250.*length(lb),'PlotFcn',[],'UseParallel',false,'FunctionTolerance',10^(-8),'MaxIterations',125.*length(lb),'MeshTolerance',10^(-9),'StepTolerance',10^(-9));
    lb=[-7 -9.*ones(1,size(X_County,1))  -9 -6.*ones(1,size(Y_County,1)) -4];
    ub=[ 3  2.*ones(1,size(X_County,1))  -1  ones(1,size(Y_County,1)) 3];
    [par_temp,f_temp]=patternsearch(@(x)Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network),par_est,[],[],[],[],lb,ub,[],options);

    lb=[-50 -32.*ones(1,size(X_County,1)) -50 -32.*ones(1,size(Y_County,1)) -32];
    ub=[ 50  4.*ones(1,size(X_County,1))   50   2.*ones(1,size(Y_County,1)) 4];
    options=optimoptions('fmincon','UseParallel',false,'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',10^4,'MaxIterations',10^4,'OptimalityTolerance',10^(-16),'StepTolerance',10^(-16));
    [par_est,fval_final]=fmincon(@(x)Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network),par_temp,[],[],[],[],lb,ub,[],options);

else
    lb=[-6 -4.*ones(1,size(X_County,1))  -5 -3.*ones(1,size(Y_County,1)) -2];
    ub=[ 0  zeros(1,size(X_County,1))    -2  zeros(1,size(Y_County,1)) 2];
    options=optimoptions("surrogateopt","UseParallel",false,'MaxFunctionEvaluations',100.*length(lb),'PlotFcn',[]); 
    [x_temp]=surrogateopt(@(x)Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network),lb,ub,options);

    lb=[-7 -9.*ones(1,size(X_County,1))  -9 -6.*ones(1,size(Y_County,1)) -4];
    ub=[ 3  2.*ones(1,size(X_County,1))  -1  ones(1,size(Y_County,1)) 3];

    x_samp=repmat(lb,250,1)+repmat(ub-lb,250,1).*lhsdesign(250,length(lb));

    options=optimoptions("surrogateopt","UseParallel",false,'MaxFunctionEvaluations',100.*length(lb),'PlotFcn',[],'InitialPoints',repmat(x_temp,100,1).*(1+0.02.*(rand(100,length(lb))-0.5)));
    [x_temp]=surrogateopt(@(x)Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network),lb,ub,options);

    lb=[-50 -64.*ones(1,size(X_County,1)) -50 -64.*ones(1,size(Y_County,1)) -32];
    ub=[ 50  4.*ones(1,size(X_County,1))   50   2.*ones(1,size(Y_County,1)) 4];
    options=optimoptions('fmincon','UseParallel',false,'FunctionTolerance',10^(-6),'MaxFunctionEvaluations',10^3,'MaxIterations',10^3,'OptimalityTolerance',10^(-6),'StepTolerance',10^(-6));
    x_new=x0;
    for xx=1:size(x0,1)
        [x_new(xx,:)]=fmincon(@(x)Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network),x0(xx,:),[],[],[],[],lb,ub,[],options);
    end
    x0=[x0;x_new;x_samp;x_temp];

    options=optimoptions("surrogateopt","UseParallel",false,'MaxFunctionEvaluations',75.*length(lb)+size(x0,1),'PlotFcn',[],'InitialPoints',x0); 
    [par_est]=surrogateopt(@(x)Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network),lb,ub,options);

    options=optimoptions('patternsearch','MaxFunctionEvaluations',125.*length(lb),'PlotFcn',[],'UseParallel',false,'FunctionTolerance',10^(-8),'MaxIterations',125.*length(lb),'MeshTolerance',10^(-9),'StepTolerance',10^(-9));
    [par_temp,f_temp]=patternsearch(@(x)Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network),par_est,[],[],[],[],lb,ub,[],options);

    options=optimoptions('fmincon','UseParallel',false,'FunctionTolerance',10^(-16),'MaxFunctionEvaluations',10^4,'MaxIterations',10^4,'OptimalityTolerance',10^(-16),'StepTolerance',10^(-16));
    [par_est,fval_final]=fmincon(@(x)Objective_Function_Dairy_Farm(x,X_County,Y_County,County_Farms,Affected_County_Farms,County_Spillover,Remaining_State_Spillover,Remainaing_Affected_State_Farms,Remainaing_State_Farms,state_weight_matrix,Dairy_Network),par_temp,[],[],[],[],lb,ub,[],options);
end

if(f_temp<fval_final)
    fval_final=f_temp;
    par_est=par_temp;
end
L=-fval_final;
AIC=aicbic(L,length(lb));
end