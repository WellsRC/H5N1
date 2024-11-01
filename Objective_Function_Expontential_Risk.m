function F = Objective_Function_Expontential_Risk(x,risk_county,Cattle_Out,Tot_Affected_Farms,County_Farms,ur_level)


n_risk_county=(risk_county(:)-min(risk_county(:)))./(max(risk_county(:))-min(risk_county(:)));
n_Cattle_Out=(Cattle_Out(:)-min(Cattle_Out(~isnan(risk_county))))./(max(Cattle_Out(~isnan(risk_county)))-min(Cattle_Out(~isnan(risk_county))));


r=n_risk_county.^2+(1-n_Cattle_Out).^2;

lambda_risk=10.^x;

p_risk_county=1-exp(-lambda_risk.*r);

p_risk_US=p_risk_county.*County_Farms./sum(County_Farms);

p_risk_US=sum(p_risk_US(~isnan(p_risk_US)));

Mode_Affected_US=ceil((sum(County_Farms)+1).*p_risk_US)-1;

est_ur=1-Tot_Affected_Farms./Mode_Affected_US;

F=10.^6.*((est_ur-ur_level).^2);

end