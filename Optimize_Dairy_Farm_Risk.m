function [par_est,L,AIC]=Optimize_Dairy_Farm_Risk(F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,x0)


lb=[-10.*ones(1,size(F_County,1)) -10.*ones(1,size(F_County,1)) -32.*ones(1,size(P_County,1)) -32.*ones(1,size(X_County,1)) -4];
ub=[15.*ones(1,size(F_County,1)) 10.*ones(1,size(F_County,1))  ones(1,size(P_County,1)) ones(1,size(X_County,1)) 0];

opts=optimoptions("ga","PlotFcn",[],"MaxGenerations",300,"FunctionTolerance",10^(-6),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250);
[par_est,f0]=ga(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p),length(lb),[],[],[],[],lb,ub,[],[],opts);

if(isnan(f0) || isinf(f0))
    par_est=[-4.70452471258611	1.50823595682272	2.63439274388045	2.41195276928562	0.271267253730793	4.84860581083441	2.48088057729165	3.50344040761827 -32.*ones(1,size(P_County,1)) -32.*ones(1,size(X_County,1)) -1.32804699856517];
    opts=optimoptions("ga","PlotFcn",[],"MaxGenerations",300,"FunctionTolerance",10^(-6),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',par_est);
    [par_est,~]=ga(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p),length(lb),[],[],[],[],lb,ub,[],[],opts);

end
PS_opts=optimoptions('patternsearch','UseParallel',false,'FunctionTolerance',10^(-9),'MaxIterations',750,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
[par_est,fval_final]=patternsearch(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p),par_est,[],[],[],[],lb,ub,[],PS_opts);


if(~isempty(x0))
    lb=[-10.*ones(1,size(F_County,1)) -10.*ones(1,size(F_County,1)) -32.*ones(1,size(P_County,1)) -32.*ones(1,size(X_County,1)) -4];
    ub=[15.*ones(1,size(F_County,1)) 10.*ones(1,size(F_County,1))  ones(1,size(P_County,1)) ones(1,size(X_County,1)) 0];

    PS_opts=optimoptions('patternsearch','UseParallel',false,'FunctionTolerance',10^(-9),'MaxIterations',375,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
    opts=optimoptions("ga","PlotFcn",[],"MaxGenerations",150,"FunctionTolerance",10^(-6),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',x0);
    [par_final,~]=ga(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p),length(lb),[],[],[],[],lb,ub,[],[],opts);

    [par_final,r_final]=patternsearch(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p),par_final,[],[],[],[],lb,ub,[],PS_opts);

    if(r_final<fval_final)
        fval_final=r_final;
        par_est=par_final;
    end
end

L=-fval_final;
AIC=aicbic(L,length(lb));
end