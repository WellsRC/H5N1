clear;
clc;

County_Strat=readtable('County_Operations_with_Inventory_Cattle_Milk_2022.csv');
County_Strat=County_Strat(~strcmp(County_Strat.State,'ALASKA')& ~strcmp(County_Strat.State,'HAWAII'),:);

State_Inv=readtable('State_Milk_Cow_Inventory_2022.csv');
State_Inv=State_Inv(~strcmp(State_Inv.State,'ALASKA')& ~strcmp(State_Inv.State,'HAWAII'),:);

County_Inv=readtable('County_Milk_Cow_Inventory_2022.csv');
County_Inv=County_Inv(~strcmp(County_Inv.State,'ALASKA')& ~strcmp(County_Inv.State,'HAWAII'),:);

Strat=unique(County_Strat.DomainCategory);
Strat=Strat([1 2 4 6 3 5 7]);

lb_v=[1 10 20 50 100 200 500];
ub_v=[9 19 49 99 199 499 50000];

for ss=1:height(State_Inv)
    f_state=strcmp(County_Inv.State,State_Inv.State(ss));
    Tot_state=State_Inv.Value(ss);
    temp_v=County_Inv.Value(f_state);
    indx_need=find(strcmp(County_Inv.State,State_Inv.State(ss)) & isnan(County_Inv.Value) );
    County_Needed=County_Inv.County(indx_need);
    
    county_reported=sum(temp_v(~isnan(temp_v)));
    
    Rem_state=Tot_state-county_reported;
    % Go through the stratifications to count those that are needed
    count_strat=zeros(length(lb_v),length(County_Needed));
    for cc=1:length(County_Needed)
        for ii=1:length(Strat)
            f_count=strcmp(County_Strat.State,State_Inv.State{ss}) & strcmp(County_Strat.County,County_Needed{cc}) & strcmp(County_Strat.DomainCategory,Strat{ii});
            if(sum(f_count)>0)
                count_strat(ii,cc)=County_Strat.Value(f_count);
            end
        end
    end

    count_strat_tot=sum(count_strat,2);
    lb=lb_v*count_strat_tot;
    ub=ub_v*count_strat_tot;
    x=(Rem_state-lb)./(ub-lb);

    for ii=1:length(indx_need)
        County_Inv.Value(indx_need(ii))=((1-x).*lb_v+x.*ub_v)*count_strat(:,ii);
    end
end


writetable(County_Inv,'Inferred_County_Milk_Cow_Inventory_2022.csv');