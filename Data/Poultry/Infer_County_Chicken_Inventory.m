clear;
clc;

Data_Strat=readtable('County_Operations_with_Inventory_Layer_Chickens_2022.csv');
Strat=unique(Data_Strat.DomainCategory);
Strat=Strat([1 8 3 7 6 2 5 9 4]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rooster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rooster=readtable('County_Chicken_Roosters_Inventory_2022.csv');

X_Rooster=zeros(height(Rooster),9);
Y_Rooster=Rooster.Value;

for cc=1:height(Rooster)
    tf=strcmpi(Rooster.State{cc},Data_Strat.State) & strcmpi(Rooster.County{cc},Data_Strat.County);
    for ss=1:length(Strat)
        ts=strcmp(Strat{ss},Data_Strat.DomainCategory);
        if(sum(ts & tf)>0)
            X_Rooster(cc,ss)=Data_Strat.Value(ts & tf);
        end
    end
end

X_fit=X_Rooster(~isnan(Y_Rooster) & sum(X_Rooster,2)>0,:);
Y_fit=Y_Rooster(~isnan(Y_Rooster) & sum(X_Rooster,2)>0);

mdl_Rooster=fitlm(X_fit,log(Y_fit));

y_pred=predict(mdl_Rooster,X_Rooster(isnan(Y_Rooster),:));

Y_Rooster(isnan(Y_Rooster))=exp(y_pred);
f_Rooster=sum(X_Rooster,2)==0;
Y_Rooster(f_Rooster)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pullet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pullet=readtable('County_Chicken_Pullets_Inventory_2022.csv');

X_Pullet=zeros(height(Pullet),9);
Y_Pullet=Pullet.Value;

for cc=1:height(Pullet)
    tf=strcmpi(Pullet.State{cc},Data_Strat.State) & strcmpi(Pullet.County{cc},Data_Strat.County);
    for ss=1:length(Strat)
        ts=strcmp(Strat{ss},Data_Strat.DomainCategory);
        if(sum(ts & tf)>0)
            X_Pullet(cc,ss)=Data_Strat.Value(ts & tf);
        end
    end
end

X_fit=X_Pullet(~isnan(Y_Pullet) & sum(X_Pullet,2)>0,:);
Y_fit=Y_Pullet(~isnan(Y_Pullet) & sum(X_Pullet,2)>0);

mdl_Pullet=fitlm(X_fit,log(Y_fit));

y_pred=predict(mdl_Pullet,X_Pullet(isnan(Y_Pullet),:));

Y_Pullet(isnan(Y_Pullet))=exp(y_pred);
f_Pullet=sum(X_Pullet,2)==0;
Y_Pullet(f_Pullet)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb=[1 50 100 400 3200 10000 20000 50000 100000];
ub=[lb(3:end)-1 1000000];
lb=log10(lb);
ub=log10(ub);
Layers=readtable('County_Chicken_Layers_Inventory_2022.csv');

X_Layers=zeros(height(Layers),9);
Y_Layers=str2double(Layers.Value);

for cc=1:height(Layers)
    tf=strcmpi(Layers.State{cc},Data_Strat.State) & strcmpi(Layers.County{cc},Data_Strat.County);
    for ss=1:length(Strat)
        ts=strcmp(Strat{ss},Data_Strat.DomainCategory);
        if(sum(ts & tf)>0)
            X_Layers(cc,ss)=Data_Strat.Value(ts & tf);
        end
    end
end

X_fit=X_Layers(~isnan(Y_Layers),:);
Y_fit=Y_Layers(~isnan(Y_Layers));

ga_options=optimoptions('ga','FunctionTolerance',10^(-8),'MaxGenerations',10^4,'PlotFcn','gaplotbestf','UseParallel',true);

[par_est]=ga(@(x)mean(((log10(X_fit*(10.^x')))-log10(Y_fit)).^2),length(lb),[],[],[],[],lb,ub,[],[],ga_options);

par_Layers=10.^par_est;

Y_Layers(isnan(Y_Layers))=X_Layers(isnan(Y_Layers),:)*(par_Layers');
Y_Layers(sum(X_Layers,2)==0)=NaN;
f_Layers=sum(X_Layers,2)==0;
Y_Layers(f_Layers)=NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Broilers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb=-12.*ones(1,length(Strat));
ub=6.*ones(1,length(Strat));

Broiler=readtable('County_Chicken_Broilers_Inventory_2022.csv');

X_Broiler=zeros(height(Broiler),9);
Y_Broiler=str2double(Broiler.Value);

for cc=1:height(Broiler)
    tf=strcmpi(Broiler.State{cc},Data_Strat.State) & strcmpi(Broiler.County{cc},Data_Strat.County);
    for ss=1:length(Strat)
        ts=strcmp(Strat{ss},Data_Strat.DomainCategory);
        if(sum(ts & tf)>0)
            X_Broiler(cc,ss)=Data_Strat.Value(ts & tf);
        end
    end
end


X_fit=X_Broiler(~isnan(Y_Broiler) & sum(X_Broiler,2)>0,:);
Y_fit=Y_Broiler(~isnan(Y_Broiler) & sum(X_Broiler,2)>0);

mdl_Broiler=fitlm(X_fit,log(Y_fit));

y_pred=predict(mdl_Broiler,X_Broiler(isnan(Y_Broiler),:));

Y_Broiler(isnan(Y_Broiler))=exp(y_pred);

f_Broiler=sum(X_Broiler,2)==0;
Y_Broiler(f_Broiler)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill in Gaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(sum(f_Rooster)>0)
    f_indx=find(f_Rooster);

    for ff=1:length(f_indx)
        tf=strcmpi(Rooster.State{cc},Rooster.State) & ~isnan(Y_Rooster);
        Y_Rooster(f_indx(ff))=median(Y_Rooster(tf));
    end
end

Rooster.Value=Y_Rooster;
writetable(Rooster,'Inferred_Rooster_Inventory_County_2022.csv');

if(sum(f_Pullet)>0)
    f_indx=find(f_Pullet);

    for ff=1:length(f_indx)
        tf=strcmpi(Pullet.State{cc},Pullet.State) & ~isnan(Y_Pullet);
        Y_Pullet(f_indx(ff))=median(Y_Pullet(tf));
    end
end

Pullet.Value=Y_Pullet;
writetable(Pullet,'Inferred_Pullet_Inventory_County_2022.csv');

if(sum(f_Layers)>0)
    f_indx=find(f_Layers);

    for ff=1:length(f_indx)
        tf=strcmpi(Layers.State{cc},Layers.State) & ~isnan(Y_Layers);
        Y_Layers(f_indx(ff))=median(Y_Layers(tf));
    end
end

Layers.Value=Y_Layers;
writetable(Layers,'Inferred_Layers_Inventory_County_2022.csv');

if(sum(f_Broiler)>0)
    f_indx=find(f_Broiler);

    for ff=1:length(f_indx)
        tf=strcmpi(Broiler.State{cc},Broiler.State) & ~isnan(Y_Broiler);
        Y_Broiler(f_indx(ff))=median(Y_Broiler(tf));
    end
end

Broiler.Value=Y_Broiler;
writetable(Broiler,'Inferred_Broiler_Inventory_County_2022.csv');


