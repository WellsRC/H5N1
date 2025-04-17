function [par_est,L,AIC]=Optimize_Dairy_Farm_Risk(F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature,x0)

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
 
if(isempty(x0))
    par_0=zeros(10,length(lb));
    fv=zeros(10,1);
    opts_psw=optimoptions("particleswarm",'UseParallel',false,"PlotFcn",[],"FunctionTolerance",10^(-6),"MaxIterations",100);
    for ii=1:10
        [par_0(ii,:),fv(ii)]=particleswarm(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature),length(lb),lb,ub,opts_psw);
    end
    
    opts=optimoptions("ga","PlotFcn",[],'UseParallel',false,"MaxGenerations",300,"FunctionTolerance",10^(-9),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',par_0);
    [par_est,f0]=ga(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature),length(lb),[],[],[],[],lb,ub,[],[],opts);
    
    if(isnan(f0)||isinf(f0))
        if(logic_temperature && isempty(Dairy_Network))
            par_est=[-2.22236971565865	0.154804918547251	0.617560829468310	2.54703013158430	0.0973854659412203	2.61965273214736	0.357968834057754	2.16792947795670 -16.*ones(1,size(P_County,1)-1) 0 -16.*ones(1,size(X_County,1)) -1.41773995561889];
        elseif(logic_temperature && ~isempty(Dairy_Network))
            par_est=[-2.22236971565865	0.154804918547251	0.617560829468310	2.54703013158430	0.0973854659412203	2.61965273214736	0.357968834057754	2.16792947795670 -16.*ones(1,size(P_County,1)-2) 0 -16 -16.*ones(1,size(X_County,1)) -1.41773995561889];
        else
            par_est=[-2.22236971565865	0.154804918547251	0.617560829468310	2.54703013158430	0.0973854659412203	2.61965273214736	0.357968834057754	2.16792947795670 -16.*ones(1,size(P_County,1)) -16.*ones(1,size(X_County,1)) -1.41773995561889];
        end
        opts=optimoptions("ga",'UseParallel',false,"PlotFcn",[],"MaxGenerations",300,"FunctionTolerance",10^(-9),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',par_est);
        [par_est,~]=ga(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature),length(lb),[],[],[],[],lb,ub,[],[],opts);                                                                 
    end
    
    PS_opts=optimoptions('patternsearch','UseParallel',false,"PlotFcn",[],'FunctionTolerance',10^(-12),'MaxIterations',750,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
    [par_est,fval_final]=patternsearch(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature),par_est,[],[],[],[],lb,ub,[],PS_opts);
   
else
    sws=max(min(100,10*length(lb)),size(x0,1));
    opts_psw=optimoptions("particleswarm",'UseParallel',false,"PlotFcn",[],"FunctionTolerance",10^(-6),"MaxIterations",100,'InitialPoints',x0,'SwarmSize',sws);
    par_0=zeros(10,length(lb));
    fv=zeros(10,1);
    for ii=1:10
        [par_0(ii,:),fv(ii)]=particleswarm(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature),length(lb),lb,ub,opts_psw);
    end

    opts=optimoptions("ga","PlotFcn",[],'UseParallel',false,"MaxGenerations",300,"FunctionTolerance",10^(-9),'CrossoverFcn','crossoverheuristic','MigrationInterval',25,'SelectionFcn',{@selectiontournament,8},'PopulationSize',250,'InitialPopulationMatrix',par_0);
    [par_final,~]=ga(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature),length(lb),[],[],[],[],lb,ub,[],[],opts);

    PS_opts=optimoptions('patternsearch',"PlotFcn",[],'UseParallel',false,'FunctionTolerance',10^(-12),'MaxIterations',750,'StepTolerance',10^(-9),'MaxFunctionEvaluations',10^4,'Cache','on');
    [par_est,fval_final]=patternsearch(@(x)Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature),par_final,[],[],[],[],lb,ub,[],PS_opts);

end

L=-fval_final;
AIC=aicbic(L,length(lb));
end