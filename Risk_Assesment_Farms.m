function mu_farm = Risk_Assesment_Farms(beta_x,Xt)
% A regression model to quantify the relative level of risk among
% a spatial regions for a specified type of farm

% Compute the logit transformed value 
z_farm = beta_x*Xt;
% Transform to value between 0 and 1
mu_farm = exp(z_farm);
end

