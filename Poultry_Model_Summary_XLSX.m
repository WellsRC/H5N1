clear;
clc;

H5N1_Variable_v={'Light_Intensity','Waterfowl_Mallard','Waterfowl_Canada_Goose','Waterfowl_AGW_Teal','Waterfowl_N_Pintail'};
Other_Variables_v={'Turkey_Operations','Broiler_Operations','Layer_Operations','Pullet_Operations','Total_Inventory','Pullet_Inventory','Layer_Inventory','Broiler_Inventory'};

load('Poultry_Models_Fit.mat')

Delta_AIC=AIC-min(AIC);

Zero_Inflated=NaN.*zeros(length(L),6);
Outbreak_Model=NaN.*zeros(length(L),41);

Fold_Turkey=zeros(length(L),1);
Fold_Pullet=zeros(length(L),1);
Fold_Broiler=zeros(length(L),1);

Spillover_per_Outbreak=zeros(length(L),1);

Variable_Names=cell(1,55);
Variable_Names(1)={'Zero_Inflated_Intercept'};
Variable_Names(7)={'Outbreak_Intercept'};
for ii=1:length(H5N1_Variable_v)
    Variable_Names(ii+1)={['Zero_Inflated_' H5N1_Variable_v{ii}]};
    for ff=1:length(Other_Variables_v)
        Variable_Names(7+ff+length(Other_Variables_v).*(ii-1))={['Outbreak_' H5N1_Variable_v{ii} '_' Other_Variables_v{ff}]};
    end
end

Variable_Names(48)={'Fold_Turkey'};
Variable_Names(49)={'Fold_Pullet'};
Variable_Names(50)={'Fold_Broiler'};
Variable_Names(51)={'Spillover_per_Outbreak'};
Variable_Names(52)={'Log_Likelihood'};
Variable_Names(53)={'AIC'};
Variable_Names(54)={'Delta_AIC'};
Variable_Names(55)={'AIC_Weight'};

for mm=1:length(L)
    [X_County,P_County,~,~,~,~,~,~,~,~,~,~,~,~,logic_par] = Poultry_Covariates(Poultry_Model.Model_H5N1{mm},Poultry_Model.Model_Farm{mm},Poultry_Model.Model_Stratified_Chicken_Inventory{mm});

    x_mle=par_est{mm};
    
    beta_x=x_mle([1 (1+2+size(P_County,1)):(2+size(P_County,1)+size(X_County,1))]);
    if(length(beta_x)>1)
        beta_x(2:end)=10.^beta_x(2:end);
    end
    l_out=[true; logic_par(6:end)];
    Outbreak_Model(mm,l_out)=beta_x;
    
    beta_p=x_mle([2:(2+size(P_County,1))]);
    if(length(beta_x)>1)
        beta_p(2:end)=-10.^beta_p(2:end);
    end
    
    l_p=[true; logic_par(1:5)];
    Zero_Inflated(mm,l_p)=beta_p;

    Fold_Turkey(mm)=10.^x_mle(end-3);
    Fold_Pullet(mm)=10.^x_mle(end-2);
    Fold_Broiler(mm)=10.^x_mle(end-1);
    Spillover_per_Outbreak(mm)=10.^x_mle(end);

end


T=[array2table(Zero_Inflated) array2table(Outbreak_Model) table(Fold_Turkey,Fold_Pullet,Fold_Broiler,Spillover_per_Outbreak,L,AIC,Delta_AIC,w_AIC)];
T.Properties.VariableNames=Variable_Names;

writetable(T,'H5N1_Poultry_Risk_Model.xlsx','Sheet','Estimated_Coefficients');