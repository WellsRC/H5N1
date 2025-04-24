clear;
clc;

load('Uncertainty_AIC_Poultry_Model.mat');
[F_County,X_County,P_County,County_Farms,Affected_County_Farms,state_weight_matrix,State_Spillover_Events,logic_par,logic_temperature] = Poultry_Covariates(Poultry_Model.Model_H5N1{1},Poultry_Model.Model_Farm{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NED TO ADJUST FOR SELECTED MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H5N1_Variable_v={'Temperature'};
Other_Variables_v={'Turkey_Operations','Turkey_Inventory'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Parameter=cell(sum(logic_par)+9,1);

Parameter(1)={'Zero Inflated Intercept: Atlantic Flyway'};
Parameter(2)={'Zero Inflated Intercept: Missippi Flyway'};
Parameter(3)={'Zero Inflated Intercept: Pacific Flyway'};
Parameter(4)={'Zero Inflated Intercept: Central Flyway'};

Parameter(5+size(P_County,1))={'Outbreak Intercept: Atlantic Flyway'};
Parameter(6+size(P_County,1))={'Outbreak Intercept: Missippi Flyway'};
Parameter(7+size(P_County,1))={'Outbreak Intercept: Pacific Flyway'};
Parameter(8+size(P_County,1))={'Outbreak Intercept: Central Flyway'};
for ii=1:length(H5N1_Variable_v)
    Parameter(ii+4)={['Zero Inflated: ' H5N1_Variable_v{ii}]};
end

for ff=1:length(Other_Variables_v)
    Parameter(8+size(P_County,1)+ff)={['Outbreak: ' Other_Variables_v{ff}]};
end

Parameter(sum(logic_par)+9)={'Spillover per outbreak'};

Point_Estimate=zeros(sum(logic_par)+9,1);

x=par_mle;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
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

Point_Estimate(1:length(beta_p))=beta_p;
Point_Estimate((1+length(beta_p)):(length(beta_p)+length(beta_x)))=beta_x;

Point_Estimate(end)=10.^x(end);


Point_Samp=zeros(5000,sum(logic_par)+9);

for jj=1:5000
    x=par_samp(jj,:);  

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

    

    Point_Samp(jj,1:length(beta_p))=beta_p;
    Point_Samp(jj,(1+length(beta_p)):(length(beta_p)+length(beta_x)))=beta_x;
    
    Point_Samp(jj,end)=10.^x(end);

end

Point_Samp=Point_Samp';

Par_95UR=prctile(Point_Samp,[2.5 97.5],2);


MLE=cell(size(Point_Estimate));
Uncertanty_Range=cell(size(Point_Estimate));

for jj=1:length(MLE)
    if(abs(Point_Estimate(jj))<10^(-4))
        MLE{jj}=[num2str(Point_Estimate(jj),'%3.2e')];
        Uncertanty_Range{jj}=['(' num2str(Par_95UR(jj,1),'%3.2e') char(8211) (num2str(Par_95UR(jj,2),'%3.2e')) ')'];
    else
        MLE{jj}=[num2str(Point_Estimate(jj),'%5.4f')];
        Uncertanty_Range{jj}=['(' num2str(Par_95UR(jj,1),'%5.4f') char(8211) (num2str(Par_95UR(jj,2),'%5.4f')) ')'];
    end
end

T=table(Parameter,MLE,Uncertanty_Range);

writetable(T,'H5N1_Poultry_Risk_Model.xlsx','Sheet','Paramter_Uncertainty');