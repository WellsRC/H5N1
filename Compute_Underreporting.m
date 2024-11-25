clear;
clc;

load('Average_Risk_Dairy.mat','avg_overall_risk_dairy_farm_County','avg_spillover_risk_dairy_farm_County');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd '/Data/Data_US_County.mat'],'US_County');

Cattle_Out_All=US_County.CONNECT_DAIRY;
Cattle_Out=zeros(height(US_County),1);
County_Farms_All=US_County.TOTAL_DAIRY_OPERATIONS;

for ss=1:height(US_County)
    t_out_state=~strcmp(US_County.STATE_NAME{ss},US_County.STATE_NAME);
    if(County_Farms_All(ss)>0)
        Cattle_Out(ss)=sum(Cattle_Out_All(ss,t_out_state))./US_County.CATTLE_INVENTORY(ss);
    end
end
avg_overall_risk_dairy_farm_County(County_Farms_All==0)=0;

Affected_County_Farms_All = US_County.DAIRY_HPAI_OUTBREAK_KNOWN;
ur_level=0.05:0.05:0.95;


Affected_County_Farms = US_County.DAIRY_HPAI_OUTBREAK_KNOWN+US_County.DAIRY_HPAI_OUTBREAK_UNKNOWN;


County_Farms=US_County.TOTAL_DAIRY_OPERATIONS;

Affected_County_Farms(County_Farms==0)=0;

Tot_Affected_Farms=sum(Affected_County_Farms);


MC_Samp=10^3;
Prob_UR=zeros(size(avg_overall_risk_dairy_farm_County,1),length(ur_level));
p_risk_all=zeros(size(avg_overall_risk_dairy_farm_County,1),length(ur_level));
new_spillover_risk_all=zeros(size(avg_overall_risk_dairy_farm_County,1),length(ur_level));
for uu=1:length(ur_level)
        x0=linspace(-2,2,50001);
        dx=x0(2)-x0(1);
        fval=zeros(50001,1);
        parfor xx=1:50001
            fval(xx)=Objective_Function_Expontential_Risk(x0(xx),avg_overall_risk_dairy_farm_County,Cattle_Out,Tot_Affected_Farms,County_Farms,ur_level(uu));
        end
        lambda_r=min(10.^x0(fval==min(fval)));
        n_risk_county=(avg_overall_risk_dairy_farm_County(:)-min(avg_overall_risk_dairy_farm_County(:)))./(max(avg_overall_risk_dairy_farm_County(:))-min(avg_overall_risk_dairy_farm_County(:)));
        n_Cattle_Out=(Cattle_Out(:)-min(Cattle_Out(~isnan(avg_overall_risk_dairy_farm_County))))./(max(Cattle_Out(~isnan(avg_overall_risk_dairy_farm_County)))-min(Cattle_Out(~isnan(avg_overall_risk_dairy_farm_County))));

        r=n_risk_county.^2+(1-n_Cattle_Out).^2;

        p_risk=1-exp(-lambda_r.*r);
        p_risk_all(:,uu)=p_risk;
        temp_MC=zeros(size(avg_overall_risk_dairy_farm_County,1),MC_Samp);
        floor_a_farm=floor(Affected_County_Farms);
        add_farm_w=Affected_County_Farms-floor_a_farm;
        parfor mc=1:MC_Samp
            r=add_farm_w-rand(size(add_farm_w));
            r(r>0)=1;
            r(r<=0)=0;
            temp_MC(:,mc)=1-binocdf(floor_a_farm+r,County_Farms,p_risk);
        end   
        Prob_UR(:,uu)=mean(temp_MC,2);
end    

Prob_UR(County_Farms==0,:)=NaN;
p_risk_all(County_Farms==0,:)=NaN;

save('Underreporting_Dairy_Farms.mat',"Prob_UR",'ur_level','County_Farms','p_risk_all')
