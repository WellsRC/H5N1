clear;
clc;

Data_Strat=readtable('County_Operations_with_Inventory_Layer_Chickens_2022.csv');
Strat=unique(Data_Strat.DomainCategory);
Strat=Strat([1 8 3 7 6 2 5 9 4]);

Data_Strat=Data_Strat(~strcmp(Data_Strat.State,'ALASKA')& ~strcmp(Data_Strat.State,'HAWAII'),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Layer_State=readtable('Layer_Inventory_State_2022.csv');
Layer_State=Layer_State(~strcmp(Layer_State.State,'ALASKA')& ~strcmp(Layer_State.State,'HAWAII'),:);

Broiler_State=readtable('Broiler_Inventory_State_2022.csv');
Broiler_State=Broiler_State(~strcmp(Broiler_State.State,'ALASKA')& ~strcmp(Broiler_State.State,'HAWAII'),:);

Pullet_State=readtable('Pullet_Inventory_State_2022.csv');
Pullet_State=Pullet_State(~strcmp(Pullet_State.State,'ALASKA')& ~strcmp(Pullet_State.State,'HAWAII'),:);

Turkey_State=readtable('Turkey_Inventory_State_2022.csv');
Turkey_State=Turkey_State(~strcmp(Turkey_State.State,'ALASKA')& ~strcmp(Turkey_State.State,'HAWAII'),:);

Num_Nan=[sum(isnan(Layer_State.Value)) sum(isnan(Broiler_State.Value)) sum(isnan(Pullet_State.Value)) sum(isnan(Turkey_State.Value))];
[Num_Nan,S_indx]=sort(Num_Nan);

All_Data=[Layer_State.Value Broiler_State.Value Pullet_State.Value Turkey_State.Value];

All_Data_National=All_Data(1,:);
All_Data=All_Data(2:end,:);
for nn=1:length(S_indx)
    if(Num_Nan(nn)>0)
        poultry_data=S_indx(nn);
        Y=All_Data(:,poultry_data);
        X=All_Data(:,~ismember(1:4,poultry_data));
        t_nan= isnan(Y+sum(X,2));
        Y_fit=Y(~t_nan);
        X_fit=X(~t_nan,:);
        mdl_est=fitlm(X_fit,log(Y_fit));
        y_pred=exp(predict(mdl_est,X(isnan(Y),:)));
        if(sum(isnan(y_pred))==0)
            y_state=Y;
            YR=All_Data_National(poultry_data)-sum(y_state(~isnan(y_state)));
            Y(isnan(Y))=YR.*y_pred./sum(y_pred);
        elseif(sum(isnan(y_pred))==1)
            Y(isnan(Y))=y_pred;
            y_state=Y;
            Y(isnan(Y))=All_Data_National(poultry_data)-sum(y_state(~isnan(y_state)));
        else
            Y=All_Data(:,poultry_data);
            X=All_Data(:,~ismember(1:4,[poultry_data S_indx(end)]));
            t_nan= isnan(Y+sum(X,2));
            Y_fit=Y(~t_nan);
            X_fit=X(~t_nan,:);
            mdl_est=fitlm(X_fit,log(Y_fit));
            y_pred=exp(predict(mdl_est,X(isnan(Y),:)));
            if(sum(isnan(y_pred))==0)
                y_state=Y;
                YR=All_Data_National(poultry_data)-sum(y_state(~isnan(y_state)));
                Y(isnan(Y))=YR.*y_pred./sum(y_pred);
            elseif(sum(isnan(y_pred))==1)
                Y(isnan(Y))=y_pred;
                y_state=Y;
                Y(isnan(Y))=All_Data_National(poultry_data)-sum(y_state(~isnan(y_state)));
            end
        end

        All_Data(:,poultry_data)=Y;
    end
end

Layer_State.Value=[All_Data_National(1); All_Data(:,1)];
Broiler_State.Value=[All_Data_National(2); All_Data(:,2)];
Pullet_State.Value=[All_Data_National(3) ; All_Data(:,3)];
Turkey_State.Value=[All_Data_National(4); All_Data(:,4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pullet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pullet=readtable('County_Chicken_Pullets_Inventory_2022.csv');
Pullet=Pullet(~strcmp(Pullet.State,'ALASKA')& ~strcmp(Pullet.State,'HAWAII'),:);

Pullet_Opp=readtable('County_Level_Pullet_Operations_with_Inventory_2022.csv');
Pullet_Opp=Pullet_Opp(~strcmp(Pullet_Opp.State,'ALASKA')& ~strcmp(Pullet_Opp.State,'HAWAII'),:);


X_Pullet=zeros(height(Pullet),9+1);
X_State=categorical(Pullet.State);
Y_Pullet=Pullet.Value;

for cc=1:height(Pullet)
    tf=strcmpi(Pullet.State{cc},Data_Strat.State) & strcmpi(Pullet.County{cc},Data_Strat.County);
    for ss=1:length(Strat)
        ts=strcmp(Strat{ss},Data_Strat.DomainCategory);
        if(sum(ts & tf)>0)
            X_Pullet(cc,ss)=Data_Strat.Value(ts & tf);
        end
    end
    tf=strcmpi(Pullet.State{cc},Pullet_Opp.State) & strcmpi(Pullet.County{cc},Pullet_Opp.County);
    X_Pullet(cc,10)=Pullet_Opp.Value(tf);

end

X_Table=[array2table(X_Pullet) table(X_State)];

X_fit=X_Table(~isnan(Y_Pullet) & sum(X_Pullet,2)>0,:);
Y_fit=Y_Pullet(~isnan(Y_Pullet) & sum(X_Pullet,2)>0);

mdl_Pullet=fitlm(X_fit,log(Y_fit));

y_pred=predict(mdl_Pullet,X_Table(isnan(Y_Pullet),:));

Y_temp=Y_Pullet;
Y_temp(isnan(Y_Pullet))=exp(y_pred);

for ss=1:height(Pullet_State)
    tf=strcmpi(Pullet_State.State(ss),Pullet.State);
    t_nan=isnan(Y_Pullet);
    if(~isempty(Y_temp(tf & t_nan)))
        Y_State=Pullet_State.Value(ss);
        Y_Remain=Y_State-sum(Y_Pullet(tf & ~t_nan));
        Y_Pullet(tf & t_nan)=Y_Remain.*Y_temp(tf & t_nan)./sum(Y_temp(tf & t_nan));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Broilers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Broiler=readtable('County_Chicken_Broilers_Inventory_2022.csv');
Broiler=Broiler(~strcmp(Broiler.State,'ALASKA')& ~strcmp(Broiler.State,'HAWAII'),:);

Broiler_Opp=readtable('County_Level_Broiler_Operations_with_Inventory_2022.csv');
Broiler_Opp=Broiler_Opp(~strcmp(Broiler_Opp.State,'ALASKA')& ~strcmp(Broiler_Opp.State,'HAWAII'),:);

X_Broiler=zeros(height(Broiler),9+1);
X_State=categorical(Broiler.State);
Y_Broiler=Broiler.Value;

for cc=1:height(Broiler)
    tf=strcmpi(Broiler.State{cc},Data_Strat.State) & strcmpi(Broiler.County{cc},Data_Strat.County);
    for ss=1:length(Strat)
        ts=strcmp(Strat{ss},Data_Strat.DomainCategory);
        if(sum(ts & tf)>0)
            X_Broiler(cc,ss)=Data_Strat.Value(ts & tf);
        end
    end
    tf=strcmpi(Broiler.State{cc},Broiler_Opp.State) & strcmpi(Broiler.County{cc},Broiler_Opp.County);
    X_Broiler(cc,10)=Broiler_Opp.Value(tf);
end

X_Table=[array2table(X_Broiler) table(X_State)];

X_fit=X_Table(~isnan(Y_Broiler) & sum(X_Broiler,2)>0,:);
Y_fit=Y_Broiler(~isnan(Y_Broiler) & sum(X_Broiler,2)>0);

mdl_Broiler=fitlm(X_fit,log(Y_fit));

y_pred=predict(mdl_Broiler,X_Table(isnan(Y_Broiler),:));

Y_temp=Y_Broiler;
Y_temp(isnan(Y_Broiler))=exp(y_pred);

for ss=1:height(Broiler_State)
    tf=strcmpi(Broiler_State.State(ss),Broiler.State);
    t_nan=isnan(Y_Broiler);
    if(~isempty(Y_temp(tf & t_nan)))
        Y_State=Broiler_State.Value(ss);
        Y_Remain=Y_State-sum(Y_Broiler(tf & ~t_nan));
        Y_Broiler(tf & t_nan)=Y_Remain.*Y_temp(tf & t_nan)./sum(Y_temp(tf & t_nan));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turkey
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Turkey=readtable('Turkey_Inventory_County_2022.csv');
Turkey=Turkey(~strcmp(Turkey.State,'ALASKA')& ~strcmp(Turkey.State,'HAWAII'),:);

Turkey_Opp=readtable('County_Turkey_Farms_with_Inventory_2022.csv');
Turkey_Opp=Turkey_Opp(~strcmp(Turkey_Opp.State,'ALASKA')& ~strcmp(Turkey_Opp.State,'HAWAII'),:);

X_Turkey=zeros(height(Turkey),9+1);
X_State=categorical(Turkey.State);
Y_Turkey=Turkey.Value;

for cc=1:height(Turkey)
    tf=strcmpi(Turkey.State{cc},Data_Strat.State) & strcmpi(Turkey.County{cc},Data_Strat.County);
    for ss=1:length(Strat)
        ts=strcmp(Strat{ss},Data_Strat.DomainCategory);
        if(sum(ts & tf)>0)
            X_Turkey(cc,ss)=Data_Strat.Value(ts & tf);
        end
    end
    tf=strcmpi(Turkey.State{cc},Turkey_Opp.State) & strcmpi(Turkey.County{cc},Turkey_Opp.County);
    X_Turkey(cc,10)=Turkey_Opp.Value(tf);
end

X_Table=[array2table(X_Turkey) table(X_State)];

X_fit=X_Table(~isnan(Y_Turkey) & sum(X_Turkey,2)>0,:);
Y_fit=Y_Turkey(~isnan(Y_Turkey) & sum(X_Turkey,2)>0);

mdl_Turkey=fitlm(X_fit,log(Y_fit));

y_pred=predict(mdl_Turkey,X_Table(isnan(Y_Turkey),:));

Y_temp=Y_Turkey;
Y_temp(isnan(Y_Turkey))=exp(y_pred);

for ss=1:height(Turkey_State)
    tf=strcmpi(Turkey_State.State(ss),Turkey.State);
    t_nan=isnan(Y_Turkey);
    if(~isempty(Y_temp(tf & t_nan)))
        Y_State=Turkey_State.Value(ss);
        Y_Remain=Y_State-sum(Y_Turkey(tf & ~t_nan));
        Y_Turkey(tf & t_nan)=Y_Remain.*Y_temp(tf & t_nan)./sum(Y_temp(tf & t_nan));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb_v=[1 50 100 400 3200 10000 20000 50000 100000];
ub_v=[lb_v(2:end)-1 1000000];

Layers=readtable('County_Chicken_Layers_Inventory_2022.csv');
Layers=Layers(~strcmp(Layers.State,'ALASKA')& ~strcmp(Layers.State,'HAWAII'),:);


Y_Layers=Layers.Value;

for ss=1:height(Layer_State)
    f_state=strcmp(Layers.State,Layer_State.State(ss));
    Tot_state=Layer_State.Value(ss);
    temp_v=Layers.Value(f_state);
    indx_need=find(strcmp(Layers.State,Layer_State.State(ss)) & isnan(Layers.Value) );
    County_Needed=Layers.County(indx_need);
    
    county_reported=sum(temp_v(~isnan(temp_v)));
    
    Rem_state=Tot_state-county_reported;
    % Go through the stratifications to count those that are needed
    count_strat=zeros(length(lb_v),length(County_Needed));
    for cc=1:length(County_Needed)
        for ii=1:length(Strat)
            f_count=strcmp(Data_Strat.State,Layer_State.State{ss}) & strcmp(Data_Strat.County,County_Needed{cc}) & strcmp(Data_Strat.DomainCategory,Strat{ii});
            if(sum(f_count)>0)
                count_strat(ii,cc)=Data_Strat.Value(f_count);
            end
        end
    end

    count_strat_tot=sum(count_strat,2);
    lb=lb_v*count_strat_tot;
    ub=ub_v*count_strat_tot;
    x=(Rem_state-lb)./(ub-lb);

    for ii=1:length(indx_need)
        Y_Layers(indx_need(ii))=((1-x).*lb_v+x.*ub_v)*count_strat(:,ii);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Write to files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

Turkey.Value=Y_Turkey;
writetable(Turkey,'Inferred_Turkey_Inventory_County_2022.csv');

Pullet.Value=Y_Pullet;
writetable(Pullet,'Inferred_Pullet_Inventory_County_2022.csv');

Layers.Value=Y_Layers;
writetable(Layers,'Inferred_Layers_Inventory_County_2022.csv');

Broiler.Value=Y_Broiler;
writetable(Broiler,'Inferred_Broiler_Inventory_County_2022.csv');


