function p_inf=Zero_Inflation(beta_p,Xt)

% Compute the logit transformed value 
z = beta_p*Xt;
% Transform to value between 0 and 1
p_inf = 1./(1+exp(-z(:)));

end