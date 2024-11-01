clear;

load([pwd '/Data/Data_US_County.mat'],'US_County');
load('Dairy_Models_Fit.mat',"par_est","w_AIC");

% Define the variables to loop through
H5N1_Variable_v={'Migratory_Bird_Density','Migratory_Birds','Mammals'};
Farm_Variables_v={'Inventory','Poultry_Operations','Connectivity'};
Stratified_Operations_Variables_v={'Cattle_Inventory_50','Cattle_Inventory_100','Cattle_Inventory_200','Cattle_Inventory_500','All'};

bin_farm=dec2bin([0:2^11-1]',11)=='1';
bin_farm=bin_farm(sum(bin_farm(:,7:11),2)<=1,:);

overall_risk_dairy_farm_County=zeros(height(US_County),length(par_est));
exposure_risk_dairy_farm_County=zeros(height(US_County),length(par_est));
susceptible_risk_dairy_farm_County=zeros(height(US_County),length(par_est));
spillover_risk_dairy_farm_County=zeros(height(US_County),length(par_est));

for mm=1:size(bin_farm,1)
    if(sum(bin_farm(mm,1:3))>0)
        H5N1_Variable=H5N1_Variable_v(bin_farm(mm,1:3));
    else
        H5N1_Variable={};
    end
    if(sum(bin_farm(mm,4:6))>0)
        Farm_Variables=Farm_Variables_v(bin_farm(mm,4:6));
    else
        Farm_Variables={};
    end
    if(sum(bin_farm(mm,7:11))>0)
        Stratified_Operations_Variables=Stratified_Operations_Variables_v(bin_farm(mm,7:11));
    else
        Stratified_Operations_Variables={};
    end

    [X_County,Y_County,County_Farms,Affected_County_Farms,~,~,~,~,~,Dairy_Network,~] = Dairy_Covariates(H5N1_Variable,Farm_Variables,Stratified_Operations_Variables);
    x=par_est{mm};
    beta_x=x(1:(1+size(X_County,1)));
    if(length(beta_x)>1)
        beta_x(2:end)=10.^beta_x(2:end);
    end
    beta_y=x((1+length(beta_x)):(end-1));
    if(length(beta_y)>1)
        beta_y(2:end)=10.^beta_y(2:end);
    end
    if(bin_farm(mm,6))
        beta_x(1+sum(bin_farm(mm,4:6)))=0; % set the connectivity parameter to zero for the susceptbility map
    end
    
    t_connect=sum(isinf(X_County),2)==size(X_County,2);
    if(~isempty(Dairy_Network))
        beta_x_temp=beta_x([1 1+find(~t_connect)']);
        r_farm_temp = Risk_Assesment_Farms(beta_y,Y_County).*Risk_Assesment_Farms(beta_x_temp,X_County(~t_connect,:));
        r_farm_temp(County_Farms==0)=0;
        X_County(t_connect,:)=(r_farm_temp(:)'*Dairy_Network);
    end

    r_nbin=10.^x(end);
    
    overall_risk_dairy_farm_County(:,mm) = Risk_Assesment_Farms(beta_y,Y_County).*Risk_Assesment_Farms(beta_x,X_County);                
    overall_risk_dairy_farm_County(County_Farms==0,mm)=0;

    exposure_risk_dairy_farm_County(:,mm) = Risk_Assesment_Farms(beta_y,Y_County);
    exposure_risk_dairy_farm_County(County_Farms==0,mm)=0;

    susceptible_risk_dairy_farm_County(:,mm) = Risk_Assesment_Farms(beta_x,X_County);                
    susceptible_risk_dairy_farm_County(County_Farms==0,mm)=0;

    spillover_risk_dairy_farm_County(:,mm)=1-nbincdf(0,r_nbin,1-overall_risk_dairy_farm_County(:,mm));
    spillover_risk_dairy_farm_County(County_Farms==0,mm)=0;
end

avg_overall_risk_dairy_farm_County=overall_risk_dairy_farm_County*w_AIC;
avg_exposure_risk_dairy_farm_County=exposure_risk_dairy_farm_County*w_AIC;
avg_susceptible_risk_dairy_farm_County=susceptible_risk_dairy_farm_County*w_AIC;
avg_spillover_risk_dairy_farm_County=spillover_risk_dairy_farm_County*w_AIC;

avg_overall_risk_dairy_farm_County(County_Farms==0)=NaN;
avg_exposure_risk_dairy_farm_County(County_Farms==0)=NaN;
avg_susceptible_risk_dairy_farm_County(County_Farms==0)=NaN;
avg_spillover_risk_dairy_farm_County(County_Farms==0)=NaN;

save('Average_Risk_Dairy_No_Connectivity.mat','avg_overall_risk_dairy_farm_County','avg_exposure_risk_dairy_farm_County','avg_susceptible_risk_dairy_farm_County','overall_risk_dairy_farm_County','exposure_risk_dairy_farm_County','susceptible_risk_dairy_farm_County','w_AIC','avg_spillover_risk_dairy_farm_County','spillover_risk_dairy_farm_County');