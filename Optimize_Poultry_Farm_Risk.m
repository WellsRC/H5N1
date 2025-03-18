function [par_est,L,AIC]=Optimize_Poultry_Farm_Risk(X_County,P_County,County_Farms,Affected_County_Farms_Unknown,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,state_weight_matrix,State_Spillover_Events,x0)

lb=[-10 -10 -16.*ones(1,size(P_County,1)) -16.*ones(1,size(X_County,1)) -4 -4 -4  -5];
ub=[15 15  ones(1,size(P_County,1)) ones(1,size(X_County,1)) 2 2 2 0];

opts=optimoptions("ga","PlotFcn",[],"MaxGenerations",300,"FunctionTolerance",10^(-6),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250);
[par_est,~]=ga(@(x)Objective_Function_Poultry_Farm(x,X_County,P_County,County_Farms,Affected_County_Farms_Unknown,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,state_weight_matrix,State_Spillover_Events),length(lb),[],[],[],[],lb,ub,[],[],opts);

if(isnan(par_est))
    par_est=[-0.743945316736287	1.14165153184578 -15.95.*ones(1,size(P_County,1)) -15.95.*ones(1,size(X_County,1)) 0.627163824461746	-0.0742059036982640	-0.0139474941565318	-1.92300407928467];
end
PS_opts=optimoptions('patternsearch','UseParallel',false,'FunctionTolerance',10^(-9),'MaxIterations',750,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
[par_est,fval_final]=patternsearch(@(x)Objective_Function_Poultry_Farm(x,X_County,P_County,County_Farms,Affected_County_Farms_Unknown,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,state_weight_matrix,State_Spillover_Events),par_est,[],[],[],[],lb,ub,[],PS_opts);


if(~isempty(x0))

    lb=[-10 -10 -16.*ones(1,size(P_County,1)) -16.*ones(1,size(X_County,1)) -4 -4 -4  -5];
ub=[15 15  ones(1,size(P_County,1)) ones(1,size(X_County,1)) 2 2 2 0];

    PS_opts=optimoptions('patternsearch','UseParallel',false,'FunctionTolerance',10^(-9),'MaxIterations',375,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
    opts=optimoptions("ga","PlotFcn",[],"MaxGenerations",150,"FunctionTolerance",10^(-6),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',x0);

    [par_final,~]=ga(@(x)Objective_Function_Poultry_Farm(x,X_County,P_County,County_Farms,Affected_County_Farms_Unknown,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,state_weight_matrix,State_Spillover_Events),length(lb),[],[],[],[],lb,ub,[],[],opts);
    [par_final,r_final]=patternsearch(@(x)Objective_Function_Poultry_Farm(x,X_County,P_County,County_Farms,Affected_County_Farms_Unknown,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,state_weight_matrix,State_Spillover_Events),par_final,[],[],[],[],lb,ub,[],PS_opts);

    if(r_final<fval_final)
        fval_final=r_final;
        par_est=par_final;
    end
end

L=-fval_final;
AIC=aicbic(L,length(lb));

end