function T=Table_Risk_Model(Var_Table)

if(strcmp(Var_Table,'Dairy'))
    load('Dairy_Models_Fit.mat',"par_est","w_AIC","Dairy_Model");

    Table_Summary=Dairy_Model;

    % Exposure_H5N1=unique(Table_Summary.Model_H5N1);
    % Susceptible_H5N1=unique(Table_Summary.Model_Farm);
    % Susceptible_Stratification_H5N1=unique(Table_Summary.Model_Stratified_Operations);

    % Inclusion_Weight=cell(length(Exposure_H5N1)+length(Susceptible_H5N1)+length(Susceptible_Stratification_H5N1),1);
    Parameter_Full=zeros(24,length(w_AIC));
    Indicator_Full=zeros(24,length(w_AIC));
    % Weighted_Uncertainty=cell(length(Exposure_H5N1)+length(Susceptible_H5N1)+length(Susceptible_Stratification_H5N1),1);
    
    Parameter_Name=cell(22,3);
    Parameter_Name{1,1}='Susceptibility';
    Parameter_Name{1,2}='Intercept';
    Parameter_Name{2,2}='Inventory';
    Parameter_Name{3,2}='Poultry Operations';
    Parameter_Name{4,2}='Connectivity';
    Parameter_Name{5,2}='Stratified operations';
    Parameter_Name{5,3}='Cattle Inventory: Less than 50';
    Parameter_Name{6,3}='Cattle Inventory: 50+';
    Parameter_Name{7,2}='Stratified operations';
    Parameter_Name{7,3}='Cattle Inventory: Less than 100';
    Parameter_Name{8,3}='Cattle Inventory: 100+';
    Parameter_Name{9,2}='Stratified operations';
    Parameter_Name{9,3}='Cattle Inventory: Less than 200';
    Parameter_Name{10,3}='Cattle Inventory: 200+';
    Parameter_Name{11,2}='Stratified operations';
    Parameter_Name{11,3}='Cattle Inventory: Less than 500';
    Parameter_Name{12,3}='Cattle Inventory: 500+';
    Parameter_Name{13,2}='Stratified operations';
    Parameter_Name{13,3}='Cattle Inventory: Less than 10';
    Parameter_Name{14,3}='Cattle Inventory: 10 to 19';
    Parameter_Name{15,3}='Cattle Inventory: 20 to 49';
    Parameter_Name{16,3}='Cattle Inventory: 50 to 99';
    Parameter_Name{17,3}='Cattle Inventory: 100 to 199';
    Parameter_Name{18,3}='Cattle Inventory: 200 to 499';
    Parameter_Name{19,3}='Cattle Inventory: 500+';
    Parameter_Name{20,1}='Exposure';
    Parameter_Name{20,2}='Intercept';
    Parameter_Name{21,2}='H5N1 incidence among wild birds';
    Parameter_Name{22,2}='Migratory bird density';
    Parameter_Name{23,2}='Detection of H5N1 among birds';
    Parameter_Name{24,2}='Detection of H5N1 among birds and bird density';



    for mm=1:height(Dairy_Model)
        [X_County,~,~,~,~,~,~,~,~,~,logic_par] = Dairy_Covariates(Dairy_Model.Model_H5N1{mm},Dairy_Model.Model_Farm{mm},Dairy_Model.Model_Stratified_Operations{mm});
        Indicator_Full(:,mm)=double(logic_par(1:end-1));
        x=par_est{mm};
        beta_x=x(1:(1+size(X_County,1)));
        if(length(beta_x)>1)
            beta_x(2:end)=10.^beta_x(2:end);
        end
        beta_y=x((1+length(beta_x)):(end-1));
        if(length(beta_y)>1)
            beta_y(2:end)=10.^beta_y(2:end);
        end
        x_temp=[beta_x(:); beta_y(:)];
        
        Parameter_Full(logic_par(1:end-1),mm)=x_temp;
    
    end
elseif(strcmp(Var_Table,'Poultry'))
    load('Poultry_Models_Fit.mat',"par_est","w_AIC","Poultry_Model");
    Table_Summary=Poultry_Model;

    Parameter_Full=zeros(16,length(w_AIC));
    Indicator_Full=zeros(16,length(w_AIC));

    Parameter_Name=cell(16,3);

    Parameter_Name{1,1}='Susceptibility';
    Parameter_Name{1,2}='Intercept';
    Parameter_Name{2,2}='Poultry_Operations';
    Parameter_Name{3,2}='Turkey_Operations';
    Parameter_Name{4,2}='Total Inventory:';
    Parameter_Name{5,2}='Stratified inventory';
    Parameter_Name{5,3}='Layer invenotry';
    Parameter_Name{6,3}='Broiler and Other inventory';
    Parameter_Name{7,2}='Stratified inventory';
    Parameter_Name{7,3}='Broiler invenotry';
    Parameter_Name{8,3}='Layer and Other inventory';
    Parameter_Name{9,2}='Stratified inventory';
    Parameter_Name{9,3}='Other invenotry';
    Parameter_Name{10,3}='Layer inventory';
    Parameter_Name{11,3}='Broiler inventory';
    Parameter_Name{12,1}='Exposure';
    Parameter_Name{12,2}='Intercept';
    Parameter_Name{13,2}='H5N1 surviellance in migratory birds';
    Parameter_Name{14,2}='Migratory bird intensity';
    Parameter_Name{15,2}='Binary indicator of H5N1 in migratory birds';
    Parameter_Name{16,2}='Combination of binary indicator of H5N1 in migratory birds and their intensity';

   
    for mm=1:height(Poultry_Model)
        [X_County,~,~,~,~,logic_par] = Poultry_Covariates(Poultry_Model.Model_H5N1{mm},Poultry_Model.Model_Farm{mm},Poultry_Model.Model_Stratified_Chicken_Inventory{mm});
        Indicator_Full(:,mm)=double(logic_par(1:end-1));
        
        x=par_est{mm};
        beta_x=x(1:(1+size(X_County,1)));
        if(length(beta_x)>1)
            beta_x(2:end)=10.^beta_x(2:end);
        end
        beta_y=x((1+length(beta_x)):(end-1));
        if(length(beta_y)>1)
            beta_y(2:end)=10.^beta_y(2:end);
        end
        x_temp=[beta_x(:); beta_y(:)];
        
        Parameter_Full(logic_par(1:end-1),mm)=x_temp;
    end
end

Avg_Parameter=Parameter_Full*w_AIC;

Model_Inclusion=Indicator_Full*w_AIC;


r=rand(10^6,1);
wc=cumsum(w_AIC);
indx_m=zeros(10^6,1);
for ii=1:10^6
    indx_m(ii)=find(r(ii)<=wc,1);
end

Dist_Parameter=Parameter_Full(:,indx_m);

Lower_Bound=prctile(Dist_Parameter,2.5,2);
Upper_Bound=prctile(Dist_Parameter,97.5,2);
T=table(Parameter_Name,Model_Inclusion,Avg_Parameter,Lower_Bound,Upper_Bound);

writetable(T,[Var_Table '_Risk_Model_Summary.csv']);



end


