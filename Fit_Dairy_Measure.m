clear;
clc;
parpool(24);
% % Define the variables to loop through
H5N1_Variable_v={'Light_Intensity','Waterfowl_Mallard','Waterfowl_Canada_Goose','Waterfowl_AGW_Teal','Waterfowl_N_Pintail'};
Farm_Variables_v={'Inventory','Connectivity'};
Stratified_Operations_Variables_v={'Cattle_Inventory_50','Cattle_Inventory_100','Cattle_Inventory_200','Cattle_Inventory_500','All'};

bin_farm=dec2bin([0:2^12-1]',12)=='1';
bin_farm=bin_farm(sum(bin_farm(:,8:12),2)<=1,:);
bin_farm=bin_farm(sum(bin_farm(:,1:5),2)>=1,:);
bin_farm=bin_farm(sum(bin_farm(:,6:12),2)>=1,:);
[~,srt_indx]=sort(sum(bin_farm,2),'ascend');
bin_farm=bin_farm(srt_indx,:);

par_est=cell(size(bin_farm,1),1);
L=zeros(size(bin_farm,1),1);
AIC=zeros(size(bin_farm,1),1);
Model_H5N1=cell(size(bin_farm,1),1);
Model_Farm=cell(size(bin_farm,1),1);
Model_Stratified_Operations=cell(size(bin_farm,1),1);

Sub_Model_Fit=NaN.*zeros(size(bin_farm,1),32);
logic_temp=cell(size(bin_farm,1),1);

for ss=2:8
    m_indx=find(sum(bin_farm,2)==ss);
    m_start=min(m_indx);
    m_end=max(m_indx);
    
    if(ss>2)
        x0_pot=Sub_Model_Fit(sum(bin_farm,2)<ss,:);
    else
        x0_pot=[];
    end
    
    parfor mm=m_start:m_end
        if(sum(bin_farm(mm,1:5))>0)
            H5N1_Variable=H5N1_Variable_v(bin_farm(mm,1:5));
        else
            H5N1_Variable={};
        end
        if(sum(bin_farm(mm,6:7))>0)
            Farm_Variables=Farm_Variables_v(bin_farm(mm,6:7));
        else
            Farm_Variables={};
        end
        if(sum(bin_farm(mm,8:12))>0)
            Stratified_Operations_Variables=Stratified_Operations_Variables_v(bin_farm(mm,8:12));
        else
            Stratified_Operations_Variables={};
        end
    
        [F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_par]= Dairy_Covariates(H5N1_Variable,Farm_Variables,Stratified_Operations_Variables);
        logic_temp{mm}=logic_par;
        if(~isempty(x0_pot))
            lt=[true(8,1); logic_par; true(1)];
            xt=max(abs(x0_pot(:,~lt)),[],2);
            x0=x0_pot(isnan(xt),:);
            x0(isnan(x0))=-32;
            x0=x0(:,lt);
        else
            x0=[];
        end
        [par_est{mm},L(mm),AIC(mm)]=Optimize_Dairy_Farm_Risk(F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,x0);
        Model_H5N1{mm}=H5N1_Variable;
        Model_Farm{mm}=Farm_Variables;
        Model_Stratified_Operations{mm}=Stratified_Operations_Variables;

        save('Dairy_Models_Fit_temp.mat','Sub_Model_Fit',"par_est","L","AIC",'ss');
    end
    for mm= m_start:m_end
        lt=[true(8,1); logic_temp{mm}; true(1)];
        Sub_Model_Fit(mm,lt)=par_est{mm};
    end
    save('Dairy_Models_Fit_temp.mat','Sub_Model_Fit',"par_est","L","AIC",'ss');
end
dAIC=AIC-min(AIC);
w_AIC=exp(-dAIC./2)./sum(exp(-dAIC./2));
Dairy_Model=table(Model_H5N1,Model_Farm,Model_Stratified_Operations,L,AIC,w_AIC);

save('Dairy_Models_Fit.mat',"Dairy_Model","par_est","L","AIC","w_AIC");