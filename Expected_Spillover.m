function Expected_Spillover
load('Average_Risk_Dairy.mat','avg_spillover_risk_scenario_dairy_farm_State','State_Name');

C=[0:size(avg_spillover_risk_scenario_dairy_farm_State,1)-1];


Spillover_Data=readtable([pwd '\Data\H5N1_Outbreaks\H5N1_State_Human.csv']);
Spillover_Data=Spillover_Data(Spillover_Data.Dairy>0,:);

P=zeros(height(Spillover_Data),1);

for ii=1:length(P)
    t_state=strcmp(Spillover_Data.StateName{ii},State_Name);
    t_case=C==Spillover_Data.Dairy(ii);
    P(ii)=avg_spillover_risk_scenario_dairy_farm_State(t_case,t_state);
end

Trim_State=State_Name(~ismember(State_Name,Spillover_Data.StateName) & ~isnan(sum(avg_spillover_risk_scenario_dairy_farm_State,1)'));

Est_Spillover=zeros(length(Trim_State),1);

    for ii=1:length(Trim_State)
        t_state=strcmp(Trim_State{ii},State_Name);
        p_temp=avg_spillover_risk_scenario_dairy_farm_State(:,t_state);
        L=zeros(length(C),1);
        for jj=1:length(C)
            L(jj)=sum((log10(p_temp(jj))-log10(P)).^2);
        end
        Est_Spillover(ii)=C(L==min(L));
    end


bar(Trim_State,Est_Spillover);

end