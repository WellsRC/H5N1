clear;
clc;

Data_Strat=readtable('County_Level_Operations_Inventory_2022.csv');
Strat=unique(Data_Strat.DomainCategory);
Strat=Strat([1 5 6 3 4 7 2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hogs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb=[1 25 50 100 200 500 1000];
ub=[lb(3:end)-1 10000];
lb=log10(lb);
ub=log10(ub);
Hogs=readtable('County_Level_Hog_Inventory_2022.csv');

X_Hogs=zeros(height(Hogs),7);
Y_Hogs=Hogs.Value;

for cc=1:height(Hogs)
    tf=strcmpi(Hogs.State{cc},Data_Strat.State) & strcmpi(Hogs.County{cc},Data_Strat.County);
    for ss=1:length(Strat)
        ts=strcmp(Strat{ss},Data_Strat.DomainCategory);
        if(sum(ts & tf)>0)
            X_Hogs(cc,ss)=Data_Strat.Value(ts & tf);
        end
    end
end

X_fit=X_Hogs(~isnan(Y_Hogs),:);
Y_fit=Y_Hogs(~isnan(Y_Hogs));

ga_options=optimoptions('ga','FunctionTolerance',10^(-8),'MaxGenerations',10^4,'PlotFcn','gaplotbestf','UseParallel',true);

[par_est,fval]=ga(@(x)mean(((log10(X_fit*(10.^x')))-log10(Y_fit)).^2),length(lb),[],[],[],[],lb,ub,[],[],ga_options);

par_Hogs=10.^par_est;

Y_Hogs(isnan(Y_Hogs))=X_Hogs(isnan(Y_Hogs),:)*(par_Hogs');
Y_Hogs(sum(X_Hogs,2)==0)=NaN;
f_Hogs=sum(X_Hogs,2)==0;
Y_Hogs(f_Hogs)=NaN;

if(sum(f_Hogs)>0)
    f_indx=find(f_Hogs);

    for ff=1:length(f_indx)
        tf=strcmpi(Hogs.State{cc},Hogs.State) & ~isnan(Y_Hogs);
        Y_Hogs(f_indx(ff))=median(Y_Hogs(tf));
    end
end

Hogs.Value=Y_Hogs;
writetable(Hogs,'Inferred_Hogs_Inventory_County_2022.csv');
