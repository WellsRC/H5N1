clear;
clc;
close all;

load([pwd '/Data/Data_US_County.mat'],'US_County');
state_remove=strcmp(US_County.STATE_NAME,"Alaska");
US_County=US_County(~state_remove,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Poultry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

State=US_County.STATE_NAME;
County=US_County.NAME;
Total_H5N1_Poultry=US_County.POULTRY_HPAI_OUTBREAK+US_County.PULLET_HPAI_OUTBREAK+US_County.TURKEY_HPAI_OUTBREAK+US_County.LAYER_HPAI_OUTBREAK+US_County.BROILER_HPAI_OUTBREAK;
Pullet_H5N1=US_County.PULLET_HPAI_OUTBREAK;
Turkey_H5N1=US_County.TURKEY_HPAI_OUTBREAK;
Layer_H5N1=US_County.LAYER_HPAI_OUTBREAK;
Broiler_H5N1=US_County.BROILER_HPAI_OUTBREAK;
Poultry_H5N1=US_County.POULTRY_HPAI_OUTBREAK;

State=State(Total_H5N1_Poultry>0);
County=County(Total_H5N1_Poultry>0);
Pullet_H5N1=Pullet_H5N1(Total_H5N1_Poultry>0);
Turkey_H5N1=Turkey_H5N1(Total_H5N1_Poultry>0);
Layer_H5N1=Layer_H5N1(Total_H5N1_Poultry>0);
Broiler_H5N1=Broiler_H5N1(Total_H5N1_Poultry>0);
Poultry_H5N1=Poultry_H5N1(Total_H5N1_Poultry>0);
Total_H5N1_Poultry=Total_H5N1_Poultry(Total_H5N1_Poultry>0);

T=table(State,County,Total_H5N1_Poultry,Poultry_H5N1,Pullet_H5N1,Turkey_H5N1,Layer_H5N1,Broiler_H5N1);

writetable(T,'Data_Summary.xlsx','Sheet','H5N1_Poultry_Farms')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Dairy Cattle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

State=US_County.STATE_NAME;
County=US_County.NAME;
Total_H5N1_Dairy=US_County.DAIRY_HPAI_OUTBREAK_KNOWN;


State=State(Total_H5N1_Dairy>0);
County=County(Total_H5N1_Dairy>0);
Total_H5N1_Dairy=Total_H5N1_Dairy(Total_H5N1_Dairy>0);

T=table(State,County,Total_H5N1_Dairy);


State=unique(US_County.STATE_NAME);

Total_H5N1_Dairy=zeros(size(State));

for ss=1:length(State)
    tf=strcmp(State{ss},US_County.STATE_NAME);
    Total_H5N1_Dairy(ss)=sum(US_County.DAIRY_HPAI_OUTBREAK_UNKNOWN(tf));
end

State=State(Total_H5N1_Dairy>0);
County=cell(size(State));
for cc=1:length(County)
    County{cc}='Unknown';
end
Total_H5N1_Dairy=Total_H5N1_Dairy(Total_H5N1_Dairy>0);

T=[T; table(State,County,Total_H5N1_Dairy)];

writetable(T,'Data_Summary.xlsx','Sheet','H5N1_Dairy_Farms')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Spillover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55



State=unique(US_County.STATE_NAME);

Spillover_Poultry=zeros(size(State)); 
Spillover_Dairy=zeros(size(State)); 

for ss=1:length(State)
    tf=strcmp(US_County.STATE_NAME,State{ss});
    Spillover_Poultry(ss)=sum(US_County.SPILLOVER_POULTRY(tf));
    Spillover_Dairy(ss)=sum(US_County.SPILLOVER_DAIRY(tf));
end
tz=Spillover_Dairy>0 | Spillover_Poultry>0;

State=State(tz);
County=cell(size(State));
for cc=1:length(County)
    County{cc}='Unknown';
end
Spillover_Poultry=Spillover_Poultry(tz);
Spillover_Dairy=Spillover_Dairy(tz);

T=[table(State,County,Spillover_Poultry,Spillover_Dairy)];

writetable(T,'Data_Summary.xlsx','Sheet','H5N1_Spillover_to_Humans')
