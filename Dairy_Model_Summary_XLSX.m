clear;
clc;

H5N1_Variable_v={'Light_Intensity','Waterfowl_Mallard','Waterfowl_Canada_Goose','Waterfowl_AGW_Teal','Waterfowl_N_Pintail','Temperature','Connectivity'};
Other_Variables_v={'Inventory','Connectivity','Total_Operations','Cattle_Inventory_1_to_9','Cattle_Inventory_10_to_19','Cattle_Inventory_20_to_49','Cattle_Inventory_50_to_99','Cattle_Inventory_100_to_199','Cattle_Inventory_200_to_499','Cattle_Inventory_500_or_more'};

load('Dairy_Models_Fit.mat')

Delta_AIC=AIC-min(AIC);

Variable_Included=zeros(length(L),10);

Zero_Inflated=NaN.*zeros(length(L),11);
Outbreak_Model=NaN.*zeros(length(L),14);

Spillover_per_Outbreak=zeros(length(L),1);

Variable_Names=cell(1,41);

Variable_Names(1)={'Model'};

Variable_Names(2)={'Light_Intensity'};
Variable_Names(3)={'Waterfowl_Mallard'};
Variable_Names(4)={'Waterfowl_Canada_Goose'};
Variable_Names(5)={'Waterfowl_AGW_Teal'};
Variable_Names(6)={'Waterfowl_N_Pintail'};
Variable_Names(7)={'Temperature'};
Variable_Names(8)={'Inventory'};
Variable_Names(9)={'Connectivity'};
Variable_Names(10)={'Total_Operations'};
Variable_Names(11)={'Stratified_Operations'};

Variable_Names(12)={'Zero_Inflated_Intercept_Atlantic_Flyway'};
Variable_Names(13)={'Zero_Inflated_Intercept_Missippi_Flyway'};
Variable_Names(14)={'Zero_Inflated_Intercept_Pacific_Flyway'};
Variable_Names(15)={'Zero_Inflated_Intercept_Central_Flyway'};

Variable_Names(23)={'Outbreak_Intercept_Atlantic_Flyway'};
Variable_Names(24)={'Outbreak_Intercept_Missippi_Flyway'};
Variable_Names(25)={'Outbreak_Intercept_Pacific_Flyway'};
Variable_Names(26)={'Outbreak_Intercept_Central_Flyway'};
for ii=1:length(H5N1_Variable_v)
    Variable_Names(ii+15)={['Zero_Inflated_' H5N1_Variable_v{ii}]};
end

for ff=1:length(Other_Variables_v)
    Variable_Names(26+ff)={['Outbreak_' Other_Variables_v{ff}]};
end

Variable_Names(37)={'Spillover_per_Outbreak'};
Variable_Names(38)={'Log_Likelihood'};
Variable_Names(39)={'AIC'};
Variable_Names(40)={'Delta_AIC'};
Variable_Names(41)={'AIC_Weight'};

Model=cell(length(L),1);

for mm=1:length(L)
    Model{mm}=['Model ' num2str(mm)];
    [F_County,X_County,P_County,~,~,~,~,~,Dairy_Network,logic_connect,logic_connect_p,logic_par,logic_temperature]= Dairy_Covariates(Dairy_Model.Model_H5N1{mm},Dairy_Model.Model_Farm{mm},Dairy_Model.Model_Stratified_Operations{mm});

    Variable_Included(mm,1:6)=double(ismember(H5N1_Variable_v(1:end-1),Dairy_Model.Model_H5N1{mm}));
    Variable_Included(mm,7)=double(contains('Inventory',Dairy_Model.Model_Farm{mm}));
    Variable_Included(mm,8)=double(contains('Connectivity',Dairy_Model.Model_Farm{mm}));
    Variable_Included(mm,9)=double(contains('Total',Dairy_Model.Model_Stratified_Operations{mm}));
    Variable_Included(mm,10)=double(contains('All',Dairy_Model.Model_Stratified_Operations{mm}));

    x=par_est{mm};

    indx_pinf=[5:(8+size(P_County,1))];

   beta_x=x(~ismember(1:length(x),[indx_pinf length(x)]));
    if(length(beta_x)>4)
        beta_x(5:end)=10.^beta_x(5:end);
    end

    beta_p=x(indx_pinf);
    if(length(beta_p)>4)
        if(logic_temperature & isempty(Dairy_Network))
            beta_p(5:end-1)=-10.^beta_p(5:end-1);    
        elseif(logic_temperature & ~isempty(Dairy_Network))
            beta_p(5:end-2)=-10.^beta_p(5:end-2);
            beta_p(end)=-10.^beta_p(end);
        else
            beta_p(5:end)=-10.^beta_p(5:end);
        end
    end

    l_out=[true(4,1); logic_par(8:end)];
    Outbreak_Model(mm,l_out)=beta_x;

    l_p=[true(4,1);  logic_par(1:7)];
    Zero_Inflated(mm,l_p)=beta_p;

    Spillover_per_Outbreak(mm)=10.^x(end);

end


T=[table(Model) array2table(Variable_Included) array2table(Zero_Inflated) array2table(Outbreak_Model) table(Spillover_per_Outbreak,L,AIC,Delta_AIC,w_AIC)];
T.Properties.VariableNames=Variable_Names;


Avg_V=(Variable_Included'*w_AIC)';
Variable_Included=[Avg_V; NaN.*zeros(3,length(Avg_V))];
Model=cell(4,1);
Model{1}='Average';
Model{2}='Median';
Model{3}='2.5 Percentile';
Model{4}='97.5 Percentile';

X_temp=table2array(T(:,12:37));
X_temp(isnan(X_temp))=0;
Zero_Inflated=NaN.*zeros(3,11);
Outbreak_Model=NaN.*zeros(3,14);
Spillover_per_Outbreak=zeros(3,1);

Zero_Inflated(1,:)=w_AIC'*X_temp(:,1:11);
Outbreak_Model(1,:)=w_AIC'*X_temp(:,12:25);
Spillover_per_Outbreak(1)=w_AIC'*X_temp(:,26);

X_Samp=zeros(10^3,size(X_temp,2));
r=rand(10^3,1);
wc=cumsum(w_AIC);

for ii=1:length(r)
    X_Samp(ii,:)=X_temp(find(r(ii)<=wc,1,"first"),:);
end

Zero_Inflated(2,:)=prctile(X_Samp(:,1:11),50);
Zero_Inflated(3,:)=prctile(X_Samp(:,1:11),2.5);
Zero_Inflated(4,:)=prctile(X_Samp(:,1:11),97.5);

Outbreak_Model(2,:)=prctile(X_Samp(:,12:25),50);
Outbreak_Model(3,:)=prctile(X_Samp(:,12:25),2.5);
Outbreak_Model(4,:)=prctile(X_Samp(:,12:25),97.5);

Spillover_per_Outbreak(2:4)=prctile(X_Samp(:,26),[50 2.5 97.5]);

L=NaN.*zeros(4,1);
AIC=NaN.*zeros(4,1);
Delta_AIC=NaN.*zeros(4,1);
w_AIC=NaN.*zeros(4,1);

T2=[table(Model) array2table(Variable_Included) array2table(Zero_Inflated) array2table(Outbreak_Model) table(Spillover_per_Outbreak,L,AIC,Delta_AIC,w_AIC)];
T2.Properties.VariableNames=Variable_Names;

T=[T2; T];

writetable(T,'H5N1_Dairy_Risk_Model.xlsx','Sheet','Estimated_Coefficients');