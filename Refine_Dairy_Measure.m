clear;
clc;
parpool(48);
% % Define the variables to loop through
H5N1_Variable_v={'Light_Intensity','Waterfowl_Mallard','Waterfowl_Canada_Goose','Waterfowl_AGW_Teal','Waterfowl_N_Pintail','Temperature'};
Farm_Variables_v={'Inventory','Connectivity'};
Stratified_Operations_Variables_v={'Total','All'};

bin_farm=dec2bin([0:2^10-1]',10)=='1';
bin_farm=bin_farm(sum(bin_farm(:,9:10),2)<=1,:);
bin_farm=bin_farm(sum(bin_farm(:,1:6),2)>=1,:);
bin_farm=bin_farm(sum(bin_farm(:,7:10),2)>=1,:);
[~,srt_indx]=sort(sum(bin_farm,2),'ascend');
bin_farm=flip(bin_farm(srt_indx,:));


load('Dairy_Models_Fit.mat','par_est','L','AIC');
L_old=flip(L);
AIC_old=flip(AIC);
Sub_Model_Fit=flip(par_est);

par_est=cell(size(bin_farm,1),1);
L=zeros(size(bin_farm,1),1);
AIC=zeros(size(bin_farm,1),1);
Model_H5N1=cell(size(bin_farm,1),1);
Model_Farm=cell(size(bin_farm,1),1);
Model_Stratified_Operations=cell(size(bin_farm,1),1);

logic_temp=cell(size(bin_farm,1),1);
for ss=9:-1:2
    m_indx=find(sum(bin_farm,2)==ss);
    m_start=min(m_indx);
    m_end=max(m_indx);
    

    logic_temp_2=logic_temp;
    
    parfor mm=m_start:m_end
        if(sum(bin_farm(mm,1:6))>0)
            H5N1_Variable=H5N1_Variable_v(bin_farm(mm,1:6));
        else
            H5N1_Variable={};
        end
        if(sum(bin_farm(mm,7:8))>0)
            Farm_Variables=Farm_Variables_v(bin_farm(mm,7:8));
        else
            Farm_Variables={};
        end
        if(sum(bin_farm(mm,9:10))>0)
            Stratified_Operations_Variables=Stratified_Operations_Variables_v(bin_farm(mm,9:10));
        else
            Stratified_Operations_Variables={};
        end
    
        [F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_par,logic_temperature]= Dairy_Covariates(H5N1_Variable,Farm_Variables,Stratified_Operations_Variables);
        logic_temp{mm}=logic_par;
        x0=Sub_Model_Fit{mm,:};
        lt=[true(8,1); logic_par; true(1)];
        if(sum(sum(bin_farm,2)>ss)>0)
            d_Model=min(double(bin_farm(sum(bin_farm,2)>ss,:)) -repmat(double(bin_farm(mm,:)),sum(sum(bin_farm,2)>ss),1),[],2);
            f_mdl_indx=find(sum(bin_farm,2)>ss);
            f_mdl_indx=f_mdl_indx(d_Model>=0);
            for zz=1:length(f_mdl_indx)
                x_temp=zeros(1,length(lt));
                l_temp=[true(8,1); logic_temp_2{f_mdl_indx(zz)}; true(1)];
                x_temp(l_temp)=Sub_Model_Fit{f_mdl_indx(zz),:};
                x0=[x0; x_temp(lt)];
            end
        end
        [par_est{mm},L(mm),AIC(mm)]=Optimize_Dairy_Farm_Risk(F_County,X_County,P_County,County_Farms,Affected_County_Farms,State_Spillover_Events,Affected_State_Farms,state_weight_matrix,Dairy_Network,logic_connect,logic_connect_p,logic_temperature,x0);
        if(L(mm)<L_old(mm))
            par_est{mm}=Sub_Model_Fit{mm,:};
            L(mm)=L_old(mm);
            AIC(mm)=AIC_old(mm);
        end
        Model_H5N1{mm}=H5N1_Variable;
        Model_Farm{mm}=Farm_Variables;
        Model_Stratified_Operations{mm}=Stratified_Operations_Variables;
    end
    save('Dairy_Models_Fit_temp.mat','Sub_Model_Fit',"par_est","L","AIC",'ss');
end
dAIC=AIC-min(AIC);
w_AIC=exp(-dAIC./2)./sum(exp(-dAIC./2));
Dairy_Model=table(Model_H5N1,Model_Farm,Model_Stratified_Operations,L,AIC,w_AIC);

save('Dairy_Models_Refined_Fit.mat',"Dairy_Model","par_est","L","AIC","w_AIC");