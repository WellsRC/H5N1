clear;
clc;
parpool(24);
% Define the variables to loop through
Exposure_Variables_v={'Connectivity','Population','Urban'};
Suceptibility_Variables_v={'Education','Gini','Income','Density','Eldery','Under_15'};
bin_suc=dec2bin([0:2^9-1]',9)=='1';
[~,srt_indx]=sort(sum(bin_suc,2),'ascend');
bin_suc=bin_suc(srt_indx,:);
par_est_H1N1=cell(size(bin_suc,1),1);
L_H1N1=zeros(size(bin_suc,1),1);
AIC_H1N1=zeros(size(bin_suc,1),1);
par_est_COVID=cell(size(bin_suc,1),1);
L_COVID=zeros(size(bin_suc,1),1);
AIC_COVID=zeros(size(bin_suc,1),1);

Model_H1N1_Exposure=cell(size(bin_suc,1),1);
Model_H1N1_Suceptibility=cell(size(bin_suc,1),1);
Model_COVID_Exposure=cell(size(bin_suc,1),1);
Model_COVID_Suceptibility=cell(size(bin_suc,1),1);


Sub_Model_Fit_H1N1=NaN.*zeros(size(bin_suc,1),12);
Sub_Model_Fit_COVID=NaN.*zeros(size(bin_suc,1),13);
logic_temp_H1N1=cell(size(bin_suc,1),1);
logic_temp_COVID=cell(size(bin_suc,1),1);

for ss=0:9
    m_indx=find(sum(bin_suc,2)==ss);
    m_start=min(m_indx);
    m_end=max(m_indx);
    
    if(ss>0)
        x0_pot_H1N1=Sub_Model_Fit_H1N1(sum(bin_suc,2)<ss,:);
        x0_pot_COVID=Sub_Model_Fit_COVID(sum(bin_suc,2)<ss,:);
    else
        x0_pot_H1N1=[];
        x0_pot_COVID=[];
    end
    
    parfor mm=m_start:m_end
        if(sum(bin_suc(mm,1:3))>0)
            Exposure_Variables=Exposure_Variables_v(bin_suc(mm,1:3));
        else
            Exposure_Variables={};
        end
        if(sum(bin_suc(mm,4:9))>0)
            Suceptibility_Variables=Suceptibility_Variables_v(bin_suc(mm,4:9));
        else
            Suceptibility_Variables={};
        end
        % Obtain the covariates and data needed for the model
        [X_County_H1N1,Y_County_H1N1,Cases_County_H1N1,Population_H1N1,logic_par_H1N1,X_County_COVID,Y_County_COVID,Deaths_County_COVID,Population_COVID,logic_par_COVID] = Population_Covariates(Exposure_Variables,Suceptibility_Variables);
        
        logic_temp_H1N1{mm}=logic_par_H1N1;
        if(~isempty(x0_pot_H1N1))
            xt=max(abs(x0_pot_H1N1(:,~logic_par_H1N1)),[],2);
            x0=x0_pot_H1N1(isnan(xt),:);
            x0(isnan(x0))=0;
            x0=x0(:,logic_par_H1N1);
        else
            x0=[];
        end
        [par_est_H1N1{mm},L_H1N1(mm),AIC_H1N1(mm)]=Optimize_Population_Risk(X_County_H1N1,Y_County_H1N1,Cases_County_H1N1,Population_H1N1,x0);
        Model_H1N1_Exposure{mm}=Exposure_Variables;
        Model_H1N1_Suceptibility{mm}=Suceptibility_Variables;
    
        logic_temp_COVID{mm}=logic_par_COVID;
        if(~isempty(x0_pot_COVID))
            xt=max(abs(x0_pot_COVID(:,~logic_par_COVID)),[],2);
            x0=x0_pot_COVID(isnan(xt),:);
            x0(isnan(x0))=0;
            x0=x0(:,logic_par_COVID);
        else
            x0=[];
        end
        [par_est_COVID{mm},L_COVID(mm),AIC_COVID(mm)]=Optimize_Population_Risk(X_County_COVID,Y_County_COVID,Deaths_County_COVID,Population_COVID,x0);
        Model_COVID_Exposure{mm}=Exposure_Variables;
        Model_COVID_Suceptibility{mm}=Suceptibility_Variables;
    end

    for mm= m_start:m_end
        Sub_Model_Fit_H1N1(mm,logic_temp_H1N1{mm})=par_est_H1N1{mm};
        Sub_Model_Fit_COVID(mm,logic_temp_COVID{mm})=par_est_COVID{mm};
    end
end

dAIC=AIC_H1N1-min(AIC_H1N1);
w_AIC_H1N1=exp(-dAIC./2)./sum(exp(-dAIC./2));
Population_Model_H1N1=table(Model_H1N1_Exposure,Model_H1N1_Suceptibility,L_H1N1,AIC_H1N1,w_AIC_H1N1);

dAIC=AIC_COVID-min(AIC_COVID);
w_AIC_COVID=exp(-dAIC./2)./sum(exp(-dAIC./2));
Population_Model_COVID=table(Model_COVID_Exposure,Model_COVID_Suceptibility,L_COVID,AIC_COVID,w_AIC_COVID);

Model_Combined_Suceptibility=Model_COVID_Suceptibility;
Model_Combined_Exposure=Model_COVID_Exposure;
L_Combined=L_COVID+L_H1N1;
AIC_Combined=AIC_H1N1+AIC_COVID;
dAIC=AIC_Combined-min(AIC_Combined);
w_AIC_Combined=exp(-dAIC./2)./sum(exp(-dAIC./2));
Population_Model_Combined=table(Model_Combined_Exposure,Model_Combined_Suceptibility,L_Combined,AIC_Combined,w_AIC_Combined);

save('Population_Models_Fit_H1N1.mat',"Population_Model_H1N1","par_est_H1N1","L_H1N1","AIC_H1N1","w_AIC_H1N1");
save('Population_Models_Fit_COVID.mat',"Population_Model_COVID","par_est_COVID","L_COVID","AIC_COVID","w_AIC_COVID");
save('Population_Models_Fit_Combined.mat',"Population_Model_Combined","par_est_COVID","par_est_H1N1","L_Combined","AIC_Combined","w_AIC_Combined");
