function [k_nb,p_nb,L_nan]=Approx_NB(mu_County_NB,p_zero_county)

r=10.^linspace(-4,3,501);
p_temp_state=log(nbinpdf(0,r,r./(r+mu_County_NB)));
[p_temp_state,ia]=unique(p_temp_state);
rt=r(ia);
rt=rt(~isinf(p_temp_state) & ~isnan(p_temp_state));

if(length(rt)>=2)
    p_temp_state=p_temp_state(~isinf(p_temp_state) & ~isnan(p_temp_state));

    k_nb=interp1(p_temp_state,rt,log(p_zero_county),"pchip");
    p_nb=k_nb./(k_nb+mu_County_NB);
    L_nan=false;
else
    L_nan=true;
    k_nb=1;
    p_nb=0;
end

end