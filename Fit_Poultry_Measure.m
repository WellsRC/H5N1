clear;
clc;
parpool(24);
% Define the variables to loop through
H5N1_Variable_v={'Migratory_Birds_H5N1','Migratory_Bird_Density','Detection_H5N1_Birds','Detection_H5N1_Bird_Density'};
Farm_Variables_v={'Turkey_Operations','Broiler_Operations','Layer_Operations','Pullet_Operations'};
Stratified_Chicken_Inventory_Variables_v={'Total_Inventory','All'};
bin_farm=dec2bin([0:2^10-1]',10)=='1';
 
bin_farm=bin_farm(sum(bin_farm(:,9:10),2)<=1,:);
[~,srt_indx]=sort(sum(bin_farm,2),'ascend');
bin_farm=bin_farm(srt_indx,:);

par_est=cell(size(bin_farm,1),1);
L=zeros(size(bin_farm,1),1);
AIC=zeros(size(bin_farm,1),1);
Model_H5N1=cell(size(bin_farm,1),1);
Model_Farm=cell(size(bin_farm,1),1);
Model_Stratified_Chicken_Inventory=cell(size(bin_farm,1),1);

Sub_Model_Fit=NaN.*zeros(size(bin_farm,1),15);
logic_temp=cell(size(bin_farm,1),1);

for ss=0:9
    m_indx=find(sum(bin_farm,2)==ss);
    m_start=min(m_indx);
    m_end=max(m_indx);

    if(ss>0)
        x0_pot=Sub_Model_Fit(sum(bin_farm,2)<ss,:);
    else
        x0_pot=[];
    end

    parfor mm=m_start:m_end
        if(sum(bin_farm(mm,1:4))>0)
            H5N1_Variable=H5N1_Variable_v(bin_farm(mm,1:4));
        else
            H5N1_Variable={};
        end
        if(sum(bin_farm(mm,5:8))>0)
            Farm_Variables=Farm_Variables_v(bin_farm(mm,5:8));
        else
            Farm_Variables={};
        end
        if(sum(bin_farm(mm,9:10))>0)
            Stratified_Chicken_Inventory_Variables=Stratified_Chicken_Inventory_Variables_v(bin_farm(mm,9:10));
        else
            Stratified_Chicken_Inventory_Variables={};
        end

        [X_County,Y_County,County_Farms,Affected_County_Farms,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,State_Spillover_Matrix,State_Spillover_Events,indx,logic_par] = Poultry_Covariates(H5N1_Variable,Farm_Variables,Stratified_Chicken_Inventory_Variables);
        logic_temp{mm}=logic_par;
        if(~isempty(x0_pot))
            xt=max(abs(x0_pot(:,~logic_par)),[],2);
            x0=x0_pot(isnan(xt),:);
            x0(isnan(x0))=-64;
            x0=x0(:,logic_par);
        else
            x0=[];
        end     
        [par_est{mm},L(mm),AIC(mm)]=Optimize_Poultry_Farm_Risk(X_County,Y_County,County_Farms,Affected_County_Farms,Pullet_Farms,Layer_Farms,Turkey_Farms,Broiler_Farms,HPAI_Pullet_Farms,HPAI_Layer_Farms,HPAI_Turkey_Farms,HPAI_Broiler_Farms,State_Spillover_Matrix,State_Spillover_Events,indx,x0);
        Model_H5N1{mm}=H5N1_Variable;
        Model_Farm{mm}=Farm_Variables;
        Model_Stratified_Chicken_Inventory{mm}=Stratified_Chicken_Inventory_Variables;
    end
    for mm= m_start:m_end
        Sub_Model_Fit(mm,logic_temp{mm})=par_est{mm};
    end
end
dAIC=AIC-min(AIC);
w_AIC=exp(-dAIC./2)./sum(exp(-dAIC./2));
Poultry_Model=table(Model_H5N1,Model_Farm,Model_Stratified_Chicken_Inventory,L,AIC,w_AIC);

save('Poultry_Models_Fit_Alternative_Likelihood.mat',"Poultry_Model","par_est","L","AIC","w_AIC");