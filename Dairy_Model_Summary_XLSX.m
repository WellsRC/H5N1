clear;
clc;

H5N1_Variable_v={'Light_Intensity','Waterfowl_Mallard','Waterfowl_Canada_Goose','Waterfowl_AGW_Teal','Waterfowl_N_Pintail','Connectivity'};
Other_Variables_v={'Inventory','Connectivity','Cattle_Inventory_50_or_less','Cattle_Inventory_over_50','Cattle_Inventory_100_or_less','Cattle_Inventory_over_100','Cattle_Inventory_200_or_less','Cattle_Inventory_over_200','Cattle_Inventory_500_or_less','Cattle_Inventory_over_500','Cattle_Inventory_1_to_9','Cattle_Inventory_10_to_19','Cattle_Inventory_20_to_49','Cattle_Inventory_50_to_99','Cattle_Inventory_100_to_199','Cattle_Inventory_200_to_499','Cattle_Inventory_500_or_more'};

load('Dairy_Models_Fit.mat')

Delta_AIC=AIC-min(AIC);

Zero_Inflated=NaN.*zeros(length(L),7);
Outbreak_Model=NaN.*zeros(length(L),91);

Spillover_per_Outbreak=zeros(length(L),1);

Variable_Names=cell(1,103);
Variable_Names(1)={'Zero_Inflated_Intercept'};
Variable_Names(8)={'Outbreak_Intercept'};
for ii=1:length(H5N1_Variable_v)
    Variable_Names(ii+1)={['Zero_Inflated_' H5N1_Variable_v{ii}]};
    if(ii<length(H5N1_Variable_v))
        for ff=1:length(Other_Variables_v)
            Variable_Names(8+ff+length(Other_Variables_v).*(ii-1))={['Outbreak_' H5N1_Variable_v{ii} '_' Other_Variables_v{ff}]};
        end
    end
end

Variable_Names(99)={'Spillover_per_Outbreak'};
Variable_Names(100)={'Log_Likelihood'};
Variable_Names(101)={'AIC'};
Variable_Names(102)={'Delta_AIC'};
Variable_Names(103)={'AIC_Weight'};

for mm=1:length(L)
    [X_County,P_County,~,~,~,~,~,~,logic_connect,logic_connect_p,logic_par] = Dairy_Covariates(Dairy_Model.Model_H5N1{mm},Dairy_Model.Model_Farm{mm},Dairy_Model.Model_Stratified_Operations{mm});

    x_mle=par_est{mm};
    
    beta_x=x_mle([1 (1+2+size(P_County,1)):(2+size(P_County,1)+size(X_County,1))]);
    if(length(beta_x)>1)
        beta_x(2:end)=10.^beta_x(2:end);
    end
    l_out=[true; logic_par(7:end)];
    Outbreak_Model(mm,l_out)=beta_x;
    
    beta_p=x_mle([2:(2+size(P_County,1))]);
    if(length(beta_x)>1)
        beta_p(2:end)=-10.^beta_p(2:end);
    end
    
    l_p=[true; logic_par(1:6)];
    Zero_Inflated(mm,l_p)=beta_p;

    Spillover_per_Outbreak(mm)=10.^x_mle(end);

end


T=[array2table(Zero_Inflated) array2table(Outbreak_Model) table(Spillover_per_Outbreak,L,AIC,Delta_AIC,w_AIC)];
T.Properties.VariableNames=Variable_Names;

writetable(T,'H5N1_Dairy_Risk_Model.xlsx','Sheet','Estimated_Coefficients');