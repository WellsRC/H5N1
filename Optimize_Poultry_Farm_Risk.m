function [par_est,L,AIC]=Optimize_Poultry_Farm_Risk(F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature,x0)

if(logic_temperature)
    lb=[-10.*ones(1,size(F_County,1)) -15.*ones(1,size(F_County,1)) -32.*ones(1,size(P_County,1)-1) -2 -32.*ones(1,size(X_County,1))  -5];
    ub=[10.*ones(1,size(F_County,1)) 15.*ones(1,size(F_County,1))  ones(1,size(P_County,1) -1) 2 ones(1,size(X_County,1)-1) 0];
else
    lb=[-10.*ones(1,size(F_County,1)) -15.*ones(1,size(F_County,1)) -32.*ones(1,size(P_County,1)) -32.*ones(1,size(X_County,1)) -5];
    ub=[10.*ones(1,size(F_County,1)) 15.*ones(1,size(F_County,1))  ones(1,size(P_County,1)) ones(1,size(X_County,1)) 0];
end

opts=optimoptions("ga","PlotFcn",[],"MaxGenerations",300,"FunctionTolerance",10^(-6),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250);
[par_est,f_0]=ga(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),length(lb),[],[],[],[],lb,ub,[],[],opts);

if(isnan(f_0) || isinf(f_0))
    if(logic_temperature)
        par_est=[1.57166766095570	1.79837598258649	1.66761874608149	1.73989227385420	3.18179460249133	2.84161860181550	0.536912520430342	3.54270138326663 -32.*ones(1,size(P_County,1)-1) 0 -32.*ones(1,size(X_County,1)) -1.81684280638710];
    else
        par_est=[1.57166766095570	1.79837598258649	1.66761874608149	1.73989227385420	3.18179460249133	2.84161860181550	0.536912520430342	3.54270138326663 -32.*ones(1,size(P_County,1)) -32.*ones(1,size(X_County,1)) -1.81684280638710];
    end

    opts=optimoptions("ga","PlotFcn",[],"MaxGenerations",200,"FunctionTolerance",10^(-6),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',par_est);
    [par_est,~]=ga(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),length(lb),[],[],[],[],lb,ub,[],[],opts);
end
PS_opts=optimoptions('patternsearch',"PlotFcn",[],'UseParallel',false,'FunctionTolerance',10^(-9),'MaxIterations',750,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
[par_est,fval_final]=patternsearch(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),par_est,[],[],[],[],lb,ub,[],PS_opts);


if(~isempty(x0))

    PS_opts=optimoptions('patternsearch','UseParallel',false,'FunctionTolerance',10^(-9),'MaxIterations',375,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
    opts=optimoptions("ga","PlotFcn",[],"MaxGenerations",150,"FunctionTolerance",10^(-6),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',x0);

    [par_final,~]=ga(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),length(lb),[],[],[],[],lb,ub,[],[],opts);
    [par_final,r_final]=patternsearch(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),par_final,[],[],[],[],lb,ub,[],PS_opts);

    if(r_final<fval_final)
        fval_final=r_final;
        par_est=par_final;
    end
end

L=-fval_final;
AIC=aicbic(L,length(lb));

end