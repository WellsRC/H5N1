clear;
clc;

H5N1_Variable_v={'Light_Intensity','Waterfowl_Mallard','Waterfowl_Canada_Goose','Waterfowl_AGW_Teal','Waterfowl_N_Pintail','Temperature'};
Other_Variables_v={'Turkey_Operations','Turkey_Inventory','Broiler_Operations','Broiler_Inventory','Layer_Operations','Layer_Inventory','Pullet_Operations','Pullet_Inventory'};

load('Poultry_Models_Refined_Fit.mat')
Poultry_Model=flip(Poultry_Model);
AIC=flip(AIC);
L=flip(L);
par_est=flip(par_est);
w_AIC=flip(w_AIC);

Delta_AIC=AIC-min(AIC);

Variable_Included=zeros(length(L),10);

Zero_Inflated=NaN.*zeros(length(L),10);
Outbreak_Model=NaN.*zeros(length(L),12);

Spillover_per_Outbreak=zeros(length(L),1);

Variable_Names=cell(1,38);

Variable_Names(1)={'Model'};

Variable_Names(2)={'Light_Intensity'};
Variable_Names(3)={'Waterfowl_Mallard'};
Variable_Names(4)={'Waterfowl_Canada_Goose'};
Variable_Names(5)={'Waterfowl_AGW_Teal'};
Variable_Names(6)={'Waterfowl_N_Pintail'};
Variable_Names(7)={'Temperature'};
Variable_Names(8)={'Turkey'};
Variable_Names(9)={'Broiler'};
Variable_Names(10)={'Layer'};
Variable_Names(11)={'Pullet'};


Variable_Names(12)={'Zero_Inflated_Intercept_Atlantic_Flyway'};
Variable_Names(13)={'Zero_Inflated_Intercept_Missippi_Flyway'};
Variable_Names(14)={'Zero_Inflated_Intercept_Pacific_Flyway'};
Variable_Names(15)={'Zero_Inflated_Intercept_Central_Flyway'};

Variable_Names(22)={'Outbreak_Intercept_Atlantic_Flyway'};
Variable_Names(23)={'Outbreak_Intercept_Missippi_Flyway'};
Variable_Names(24)={'Outbreak_Intercept_Pacific_Flyway'};
Variable_Names(25)={'Outbreak_Intercept_Central_Flyway'};
for ii=1:length(H5N1_Variable_v)
    Variable_Names(ii+15)={['Zero_Inflated_' H5N1_Variable_v{ii}]};
end

for ff=1:length(Other_Variables_v)
    Variable_Names(25+ff)={['Outbreak_' Other_Variables_v{ff}]};
end

Variable_Names(34)={'Spillover_per_Outbreak'};
Variable_Names(35)={'Log_Likelihood'};
Variable_Names(36)={'AIC'};
Variable_Names(37)={'Delta_AIC'};
Variable_Names(38)={'AIC_Weight'};

Model=cell(length(L),1);

for mm=1:length(L)
    Model{mm}=['Model ' num2str(mm)];
    [F_County,X_County,P_County,~,~,~,~,logic_par,logic_temperature] = Poultry_Covariates(Poultry_Model.Model_H5N1{mm},Poultry_Model.Model_Farm{mm});

    Variable_Included(mm,1:6)=double(ismember(H5N1_Variable_v,Poultry_Model.Model_H5N1{mm}));
    Variable_Included(mm,7)=double(contains('Turkey_Operations',Poultry_Model.Model_Farm{mm}));
    Variable_Included(mm,8)=double(contains('Broiler_Operations',Poultry_Model.Model_Farm{mm}));
    Variable_Included(mm,9)=double(contains('Layer_Operations',Poultry_Model.Model_Farm{mm}));
    Variable_Included(mm,10)=double(contains('Pullet_Operations',Poultry_Model.Model_Farm{mm}));
    x=par_est{mm};
    
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
        
    l_out=[true(4,1); logic_par(7:end)];
    Outbreak_Model(mm,l_out)=beta_x;
    
    
    l_p=[true(4,1); logic_par(1:6)];
    Zero_Inflated(mm,l_p)=beta_p;

    Spillover_per_Outbreak(mm)=10.^x(end);

end


T=[table(Model) array2table(Variable_Included) array2table(Zero_Inflated) array2table(Outbreak_Model) table(Spillover_per_Outbreak,L,AIC,Delta_AIC,w_AIC)];
T.Properties.VariableNames=Variable_Names;

Avg_V=(Variable_Included'*w_AIC)';
Variable_Included=Avg_V;
Model=cell(1,1);
Model{1}='Weighted_Average';

X_temp=table2array(T(:,12:34));
X_temp(isnan(X_temp))=0;

Zero_Inflated=w_AIC'*X_temp(:,1:10);
Outbreak_Model=w_AIC'*X_temp(:,11:22);
Spillover_per_Outbreak=w_AIC'*X_temp(:,23);

L=NaN.*zeros(1,1);
AIC=NaN.*zeros(1,1);
Delta_AIC=NaN.*zeros(1,1);
w_AIC=NaN.*zeros(1,1);

T2=[table(Model) array2table(Variable_Included) array2table(Zero_Inflated) array2table(Outbreak_Model) table(Spillover_per_Outbreak,L,AIC,Delta_AIC,w_AIC)];
T2.Properties.VariableNames=Variable_Names;

T=[T2; T];
writetable(T,'H5N1_Poultry_Risk_Model.xlsx','Sheet','Estimated_Coefficients');