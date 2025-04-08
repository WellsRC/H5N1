clear;
clc;

H5N1_Variable_v={'Light_Intensity','Waterfowl_Mallard','Waterfowl_Canada_Goose','Waterfowl_AGW_Teal','Waterfowl_N_Pintail'};
Other_Variables_v={'Turkey_Operations','Turkey_Inventory','Broiler_Operations','Broiler_Inventory','Layer_Operations','Layer_Inventory','Pullet_Operations','Pullet_Inventory','Temperature'};

load('Poultry_Models_Fit.mat')

Delta_AIC=AIC-min(AIC);

Zero_Inflated=NaN.*zeros(length(L),9);
Outbreak_Model=NaN.*zeros(length(L),13);

Fold_Turkey=zeros(length(L),1);
Fold_Pullet=zeros(length(L),1);
Fold_Broiler=zeros(length(L),1);

Spillover_per_Outbreak=zeros(length(L),1);

Variable_Names=cell(1,27);

Variable_Names(1)={'Zero_Inflated_Intercept_Atlantic_Flyway'};
Variable_Names(2)={'Zero_Inflated_Intercept_Missippi_Flyway'};
Variable_Names(3)={'Zero_Inflated_Intercept_Pacific_Flyway'};
Variable_Names(4)={'Zero_Inflated_Intercept_Central_Flyway'};

Variable_Names(10)={'Outbreak_Intercept_Atlantic_Flyway'};
Variable_Names(11)={'Outbreak_Intercept_Missippi_Flyway'};
Variable_Names(12)={'Outbreak_Intercept_Pacific_Flyway'};
Variable_Names(13)={'Outbreak_Intercept_Central_Flyway'};
for ii=1:length(H5N1_Variable_v)
    Variable_Names(ii+4)={['Zero_Inflated_' H5N1_Variable_v{ii}]};
end

for ff=1:length(Other_Variables_v)
    Variable_Names(13+ff)={['Outbreak_' Other_Variables_v{ff}]};
end

Variable_Names(23)={'Spillover_per_Outbreak'};
Variable_Names(24)={'Log_Likelihood'};
Variable_Names(25)={'AIC'};
Variable_Names(26)={'Delta_AIC'};
Variable_Names(27)={'AIC_Weight'};

for mm=1:length(L)
    [F_County,X_County,P_County,~,~,~,~,logic_par,logic_temperature] = Poultry_Covariates(Poultry_Model.Model_H5N1{mm},Poultry_Model.Model_Farm{mm});

    x=par_est{mm};
    
    indx_pinf=[5:(8+size(P_County,1))];

    beta_x=x(~ismember(1:length(x),[indx_pinf length(x)]));
    if(length(beta_x)>4)
        if(logic_temperature)
            indx_betax=[1:(size(X_County,1)-1)];
        else
            indx_betax=[1:(size(X_County,1))];
        end
        beta_x(4+indx_betax)=10.^beta_x(4+indx_betax);    
    end
    
    beta_p=x(indx_pinf);
    if(length(beta_p)>4)
        beta_p(5:end)=-10.^beta_p(5:end);
    end
        
    l_out=[true(4,1); logic_par(6:end)];
    Outbreak_Model(mm,l_out)=beta_x;
    
    
    l_p=[true(4,1); logic_par(1:5)];
    Zero_Inflated(mm,l_p)=beta_p;

    Spillover_per_Outbreak(mm)=10.^x(end);

end


T=[array2table(Zero_Inflated) array2table(Outbreak_Model) table(Spillover_per_Outbreak,L,AIC,Delta_AIC,w_AIC)];
T.Properties.VariableNames=Variable_Names;

writetable(T,'H5N1_Poultry_Risk_Model.xlsx','Sheet','Estimated_Coefficients');