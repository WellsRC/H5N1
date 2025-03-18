function F = Objective_High_Risk_Regions(rx,Average_County_Estimate,D_County)

x=double(D_County<=10.^rx);

c_z=sum(Average_County_Estimate(x==1));
n_z=sum(x==1);

c_out=sum(Average_County_Estimate(x==0));
n_out=sum(x==0);

c_tot=sum(Average_County_Estimate);
n_tot=length(x);

if(c_z/n_z>c_out/n_out)
    % L=(c_z./n_z).^c_z.*((n_z-c_z)./n_z).^(n_z-c_z).*(c_out./n_out).^c_out.*((n_out-c_out)./n_out).^(n_out-c_out);
    LL=c_z.*(log(c_z)-log(n_z))+(n_z-c_z).*(log((n_z-c_z))-log(n_z))+c_out.*(log(c_out)-log(n_out))+(n_out-c_out).*(log((n_out-c_out))-log(n_out));
else
    % L=c_tot.^c_tot.*(n_tot-c_tot).^(n_tot-c_tot)./(n_tot.^n_tot);
    LL=c_tot.*log(c_tot)+(n_tot-c_tot).*log(n_tot-c_tot)-n_tot.*log(n_tot);
end
F=-LL;
end

