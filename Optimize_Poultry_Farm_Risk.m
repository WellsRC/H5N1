function [par_est,L,AIC]=Optimize_Poultry_Farm_Risk(F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature,x0)

if(logic_temperature)
    lb=[-7.*ones(1,size(F_County,1)) -7.*ones(1,size(F_County,1)) -16.*ones(1,size(P_County,1)-1) -0.5 -16.*ones(1,size(X_County,1))  -3];
    ub=[7.*ones(1,size(F_County,1)) 7.*ones(1,size(F_County,1))  log10(5).*ones(1,size(P_County,1)-1) 0.5 log10(5).*ones(1,size(X_County,1)) 0];
else
    lb=[-7.*ones(1,size(F_County,1)) -7.*ones(1,size(F_County,1)) -16.*ones(1,size(P_County,1)) -16.*ones(1,size(X_County,1)) -3];
    ub=[7.*ones(1,size(F_County,1)) 7.*ones(1,size(F_County,1))  log10(5).*ones(1,size(P_County,1)) log10(5).*ones(1,size(X_County,1)) 0];
end

if(isempty(x0))
    par_0=zeros(10,length(lb));
    
    opts_psw=optimoptions("particleswarm",'UseParallel',false,"PlotFcn",[],"FunctionTolerance",10^(-6),"MaxIterations",100,'SwarmSize',10.*length(lb));
    for ii=1:10
        [par_0(ii,:),~]=particleswarm(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),length(lb),lb,ub,opts_psw);
    end
    
    opts=optimoptions("ga",'UseParallel',false,"PlotFcn",[],"MaxGenerations",300,"FunctionTolerance",10^(-9),'CrossoverFcn',{@crossoverintermediate},'MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'CrossoverFraction',0.7,'InitialPopulationMatrix',par_0);
    [par_est,f_0]=ga(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),length(lb),[],[],[],[],lb,ub,[],[],opts);
    
    if(isnan(f_0) || isinf(f_0))
        if(logic_temperature)
            par_est=[0.831770408546465	0.855557420402483	1.15922161174091	0.893964923241808	2.08439682499956	1.64089980304958	0.548850994019000	1.78088795061238 -16.*ones(1,size(P_County,1)-1) 0 -16.*ones(1,size(X_County,1)) -1.72268101284425];
        else
            par_est=[0.831770408546465	0.855557420402483	1.15922161174091	0.893964923241808	2.08439682499956	1.64089980304958	0.548850994019000	1.78088795061238 -16.*ones(1,size(P_County,1)) -16.*ones(1,size(X_County,1)) -1.72268101284425];
        end
    
        opts=optimoptions("ga","PlotFcn",[],"MaxGenerations",300,"FunctionTolerance",10^(-9),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',par_est);
        [par_est,~]=ga(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),length(lb),[],[],[],[],lb,ub,[],[],opts);
    end
    
    PS_opts=optimoptions('patternsearch',"PlotFcn",[],'UseParallel',false,'FunctionTolerance',10^(-12),'MaxIterations',750,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
    [par_est,fval_final]=patternsearch(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),par_est,[],[],[],[],lb,ub,[],PS_opts);

else

    sws=max(min(100,10*length(lb)),size(x0,1));
    opts_psw=optimoptions("particleswarm",'UseParallel',false,"PlotFcn",[],"FunctionTolerance",10^(-6),"MaxIterations",100,'InitialPoints',x0,'SwarmSize',sws);

    par_0=zeros(10,length(lb));
    for ii=1:10
        [par_0(ii,:),~]=particleswarm(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),length(lb),lb,ub,opts_psw);
    end

    opts=optimoptions("ga","PlotFcn",[],"MaxGenerations",300,"FunctionTolerance",10^(-9),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',par_0);
    [par_final,~]=ga(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),length(lb),[],[],[],[],lb,ub,[],[],opts);

    PS_opts=optimoptions('patternsearch','UseParallel',false,'FunctionTolerance',10^(-12),'MaxIterations',750,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
    [par_est,fval_final]=patternsearch(@(x)Objective_Function_Poultry_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_temperature),par_final,[],[],[],[],lb,ub,[],PS_opts);

end

L=-fval_final;
AIC=aicbic(L,length(lb));

end