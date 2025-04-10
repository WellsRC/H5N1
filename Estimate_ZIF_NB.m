function [z_nb,k_nb,p_nb,L_Nan]=Estimate_ZIF_NB(mu_County_NB,var_County_NB,p_zero_county)

r=10.^linspace(-4,3,501);
K_Const=(var_County_NB-mu_County_NB)./mu_County_NB+mu_County_NB;

p_nb=(r+1)./(r+1+K_Const);
m_nb=r.*(1-p_nb)./p_nb;
z_nb=1-mu_County_NB./m_nb;

rt=r(z_nb>=0);
p_nb=p_nb(z_nb>=0);
z_nb=z_nb(z_nb>=0);
p_temp_state=log(z_nb+(1-z_nb).*nbinpdf(0,rt,p_nb));

[p_temp_state,ia]=unique(p_temp_state);
rt=rt(ia);
rt=rt(~isinf(p_temp_state) & ~isnan(p_temp_state));
p_temp_state=p_temp_state(~isinf(p_temp_state) & ~isnan(p_temp_state));
L_Nan=false;
if(length(rt)>=2)
    k_nb=interp1(p_temp_state,rt,log(p_zero_county),"pchip");
    if(k_nb<=0)
        k_nb=10.^(interp1(p_temp_state,log10(rt),log(p_zero_county),"pchip"));
    end
    p_nb=(k_nb+1)./(k_nb+1+K_Const);
    m_nb=k_nb.*(1-p_nb)./p_nb;
    z_nb=max(1-mu_County_NB./m_nb,0);
else
    % Arbitrary since we will be returning a NaN value for the lieklihood
    k_nb=1;
    p_nb=(k_nb+1)./(k_nb+1+K_Const);
    z_nb=0;
    L_Nan=true;
end

end
