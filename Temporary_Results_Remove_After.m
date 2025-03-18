function [avg_outbreak_farm_County,avg_outbreak_risk_farm_County,avg_spillover_risk_farm_County,avg_outbreak_risk_farm_State,avg_spillover_risk_farm_State,Spillover_Post_State,Outbreak_Post_State]=Temporary_Results_Remove_After()

H5N1_Variable_v={'Light_Intensity','Waterfowl_Mallard','Waterfowl_Canada_Goose','Waterfowl_AGW_Teal','Waterfowl_N_Pintail'};
Farm_Variables_v={'Inventory','Poultry_Operations','Connectivity'};
Stratified_Operations_Variables_v={'Cattle_Inventory_50','Cattle_Inventory_100','Cattle_Inventory_200','Cattle_Inventory_500','All'};

bin_farm=dec2bin([0:2^13-1]',13)=='1';
bin_farm=bin_farm(sum(bin_farm(:,9:13),2)<=1,:);
bin_farm=bin_farm(sum(bin_farm(:,1:5),2)>=1,:);
bin_farm=bin_farm(sum(bin_farm(:,6:13),2)>=1,:);
[~,srt_indx]=sort(sum(bin_farm,2),'ascend');
bin_farm=bin_farm(srt_indx,:);

par_est=cell(size(bin_farm,1),1);
L=zeros(size(bin_farm,1),1);
AIC=zeros(size(bin_farm,1),1);
Model_H5N1=cell(size(bin_farm,1),1);
Model_Farm=cell(size(bin_farm,1),1);
Model_Stratified_Operations=cell(size(bin_farm,1),1);

Sub_Model_Fit=NaN.*zeros(size(bin_farm,1),90);
logic_temp=cell(size(bin_farm,1),1);

ss=2;
    m_indx=find(sum(bin_farm,2)==ss);
    m_start=min(m_indx);
    m_end=max(m_indx);
    
    if(ss>2)
        x0_pot=Sub_Model_Fit(sum(bin_farm,2)<ss,:);
    else
        x0_pot=[];
    end
    
     mm=5;
        if(sum(bin_farm(mm,1:5))>0)
            H5N1_Variable=H5N1_Variable_v(bin_farm(mm,1:5));
        else
            H5N1_Variable={};
        end
        if(sum(bin_farm(mm,6:8))>0)
            Farm_Variables=Farm_Variables_v(bin_farm(mm,6:8));
        else
            Farm_Variables={};
        end
        if(sum(bin_farm(mm,9:13))>0)
            Stratified_Operations_Variables=Stratified_Operations_Variables_v(bin_farm(mm,9:13));
        else
            Stratified_Operations_Variables={};
        end
    
        [X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect,logic_par]= Dairy_Covariates(H5N1_Variable,Farm_Variables,Stratified_Operations_Variables);
        [par_est,~,~]=Optimize_Dairy_Farm_Risk(X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,County_Suppressed_State,County_Nonsuppressed,state_weight_matrix,Dairy_Network,logic_connect,[]);

        x=par_est;

beta_x=x([1 (1+2+size(P_County,1)):(2+size(P_County,1)+size(X_County,1))]);
if(length(beta_x)>1)
    beta_x(2:end)=10.^beta_x(2:end);
end

beta_p=x([2:(2+size(P_County,1))]);
if(length(beta_x)>1)
    beta_p(2:end)=-10.^beta_p(2:end);
end

kappa_spillover=10.^x(end);

if(~isempty(Dairy_Network))
    beta_x_temp=beta_x([1 1+find(~logic_connect)']);
    mu_farm_temp = Risk_Assesment_Farms(beta_x_temp,X_County(~logic_connect,:));
    mu_farm_temp(County_Farms==0)=0;
    X_County(logic_connect,:)=X_County(logic_connect,:).*repmat((mu_farm_temp(:)'*Dairy_Network),sum(logic_connect),1);
end



p_inf_County=Zero_Inflation(beta_p,P_County);

mu_farm_County = Risk_Assesment_Farms(beta_x,X_County);
mu_farm_County(County_Farms==0)=0;

mu_farm_state=zeros(size(State_Spillover_Events));
p_inf_State=zeros(size(State_Spillover_Events));
for ss=1:length(mu_farm_state)
    temp_county=state_weight_matrix(ss,:).*mu_farm_County; 
    mu_farm_state(ss)=sum(temp_county);
    p_inf_State(ss)=prod(p_inf_County(state_weight_matrix(ss,:)==1));
end

Spillover_Post_State=zeros(length(mu_farm_state),26);
Outbreak_Post_State=zeros(length(mu_farm_state),10001);
Spillover_Post_State(:,1)=(p_inf_State(:)+(1-p_inf_State(:)).*poisspdf(0,mu_farm_state(:).*kappa_spillover));
Outbreak_Post_State(:,1)=(p_inf_State(:)+(1-p_inf_State(:)).*poisspdf(0,mu_farm_state(:)));
for ss=1:length(mu_farm_state)
    Spillover_Post_State(ss,2:end)=((1-p_inf_State(ss)).*poisspdf([1:25],mu_farm_state(ss).*kappa_spillover));
    Outbreak_Post_State(ss,2:end)=((1-p_inf_State(ss)).*poisspdf([1:10000],mu_farm_state(ss)));
end

avg_outbreak_farm_County=(1-p_inf_County(:)).*mu_farm_County(:);

avg_outbreak_risk_farm_County=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:)));
avg_outbreak_risk_farm_State=1-(p_inf_State(:)+(1-p_inf_State(:)).*poisspdf(0,mu_farm_state(:)));

avg_spillover_risk_farm_State=1-(p_inf_State(:)+(1-p_inf_State(:)).*poisspdf(0,mu_farm_state(:).*kappa_spillover));
avg_spillover_risk_farm_County=1-(p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_farm_County(:).*kappa_spillover));
end