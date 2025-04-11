clear;
clc;
rng(5120410);
load('Average_Risk_Poultry.mat','w_AIC','post_spillover_poultry_farm_State');
w_AICp=cumsum(w_AIC);
load('Average_Risk_Dairy.mat','w_AIC','post_spillover_dairy_farm_State');
w_AICd=cumsum(w_AIC);

r_p=rand(10^4,1);
r_d=rand(10^4,1);

f_indxp=zeros(10^4,1);
f_indxd=zeros(10^4,1);

for ii=1:10^4
    f_indxp(ii)=find(r_p(ii)<=w_AICp,1,"first");
    f_indxd(ii)=find(r_d(ii)<=w_AICd,1,"first");
end

spill_dairy=post_spillover_dairy_farm_State(:,:,f_indxd);
spill_poultry=post_spillover_poultry_farm_State(:,:,f_indxp);

k_onward_transmission=2.69;
R0=0.05;
p_no_onward_transmission=nbinpdf(0,k_onward_transmission,k_onward_transmission./(k_onward_transmission+R0));

