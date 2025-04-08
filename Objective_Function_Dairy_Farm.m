function F = Objective_Function_Dairy_Farm(x,F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature)

indx_pinf=[5:(8+size(P_County,1))];

beta_x=x(~ismember(1:length(x),[indx_pinf length(x)]));
if(length(beta_x)>4)
    beta_x(5:end)=10.^beta_x(5:end);
end

beta_p=x(indx_pinf);
if(length(beta_p)>4)
    if(logic_temperature)
        beta_p(5:end-1)=-10.^beta_p(5:end-1);    
    else
        beta_p(5:end)=-10.^beta_p(5:end);
    end
end

kappa_spillover=10.^x(end);


if(~isempty(Dairy_Network))
    beta_x_temp=beta_x([1:4 4+find(~logic_connect)]);
    
    mu_farm_temp = Risk_Assesment_Farms(beta_x_temp,[F_County; X_County(~logic_connect,:)]);
    mu_farm_temp(County_Farms==0)=0;

    beta_p_temp=beta_p([1:4 4+find(~logic_connect_p)]);
    p_inf_County_temp=Zero_Inflation(beta_p_temp,[F_County; P_County(~logic_connect_p,:)]);
    p_inf_County_temp(County_Farms==0)=1;

    temp_r=(1-p_inf_County_temp(:)).*mu_farm_temp(:);
    X_County(logic_connect,:)=(temp_r(:)'*Dairy_Network);
    P_County(logic_connect_p,:)=(temp_r(:)'*Dairy_Network);
end

delta_Affected_State=Affected_State_Farms-state_weight_matrix*Affected_County_Farms;
upper_County_Unknown=(delta_Affected_State')*state_weight_matrix;
upper_County_Unknown=upper_County_Unknown(:);

p_inf_County=Zero_Inflation(beta_p,[F_County; P_County]);
p_inf_County(County_Farms==0)=1;

mu_farm_County = Risk_Assesment_Farms(beta_x,[F_County; X_County]);
mu_farm_County(County_Farms==0)=0;

mu_farm_State=zeros(size(State_Spillover_Events));
k_State=zeros(size(State_Spillover_Events));
k_State_spillover=zeros(size(State_Spillover_Events));

p_temp=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:));
p_temp_spill=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,kappa_spillover.*mu_farm_County(:)); % Use the "TRUE" county risk as aspects are not observed if lack of surviellance 

r=10.^linspace(-3,4,501);
L_Nan=false;

w_state=zeros(size(State_Spillover_Events));
County_w_Farms=County_Farms>0;

for ss=1:length(mu_farm_State)

    temp_county=(1-p_inf_County(:)).*mu_farm_County(:); 
    mu_farm_State(ss)=state_weight_matrix(ss,:)*temp_county;   
    w_state(ss)=state_weight_matrix(ss,:)*County_w_Farms(:);
    
    p_zero_county=prod(p_temp(state_weight_matrix(ss,:)==1));
    if(p_zero_county==0)
        p_zero_county=10^(-64);
    end
    p_temp_state=r.*log(r./(r+mu_farm_State(ss)));
    [p_temp_state,ia]=unique(p_temp_state);
    rt=r(ia);
    rt=rt(~isinf(p_temp_state) & ~isnan(p_temp_state));
    p_temp_state=p_temp_state(~isinf(p_temp_state) & ~isnan(p_temp_state));
    if(length(rt)>=2)
        k_State(ss)=interp1(p_temp_state,rt,log(p_zero_county),"pchip");
    else
        k_State(ss)=1; % Arbitrary since we will be returning a NaN value for the lieklihood
        L_Nan=true;
    end

    p_zero_county=prod(p_temp_spill(state_weight_matrix(ss,:)==1));
    if(p_zero_county==0)
        p_zero_county=10^(-64);
    end
    p_temp_state=r.*log(r./(r+kappa_spillover.*mu_farm_State(ss)));
    [p_temp_state,ia]=unique(p_temp_state);
    rt=r(ia);
    rt=rt(~isinf(p_temp_state) & ~isnan(p_temp_state));
    p_temp_state=p_temp_state(~isinf(p_temp_state) & ~isnan(p_temp_state));
    if(length(rt)>=2)
        k_State_spillover(ss)=interp1(p_temp_state,rt,log(p_zero_county),"pchip");
    else
        k_State_spillover(ss)=1; % Arbitrary since we will be returning a NaN value for the lieklihood
        L_Nan=true;
    end
end

k_outbreak_county=mu_farm_County(:);
k_outbreak_state=mu_farm_State(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_County=log((1-p_inf_County(:)).*poisspdf(Affected_County_Farms(:),k_outbreak_county(:)));
L_County(Affected_County_Farms==0)=log(p_inf_County(Affected_County_Farms==0)+(1-p_inf_County(Affected_County_Farms==0)).*poisspdf(Affected_County_Farms(Affected_County_Farms==0),k_outbreak_county(Affected_County_Farms==0)));

cdf=(1-p_inf_County(:)).*poisscdf(Affected_County_Farms(:)+upper_County_Unknown(:),k_outbreak_county(:))-(1-p_inf_County(:)).*poisscdf(Affected_County_Farms(:),k_outbreak_county(:)); % do not need the addition of p_inf_County alone since it cancels out in subtration because computing cdf in each 
L_County(upper_County_Unknown>0)=log(cdf(upper_County_Unknown>0));

L_County=L_County(County_Farms>0 & ~isnan(Affected_County_Farms));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_State=log(nbinpdf(Affected_State_Farms(:),k_State(:),k_State(:)./(k_State(:)+k_outbreak_state(:))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spillover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_Spillover_State=log(nbinpdf(State_Spillover_Events(:),k_State_spillover(:),k_State_spillover(:)./(k_State_spillover(:)+kappa_spillover.*mu_farm_State(:))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~L_Nan)
    F=-sum(L_County) -sum(w_state(:).*L_Spillover_State(:)) -sum(w_state(:).*L_State(:));
else
    F=NaN;
end
end

