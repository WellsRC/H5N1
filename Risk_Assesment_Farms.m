function mu_farm = Risk_Assesment_Farms(beta_x,X)
% A logisitic regression model to quantify the relative level of risk among
% a spatial regions for a specified type of farm

% Append ones to the covariate matrix for the cotanstan
Xt=[ones(1,size(X,2));X];
% Compute the logit transformed value 
z_farm = beta_x*Xt;
% Transform to value between 0 and 1
mu_farm = exp(z_farm);
end

