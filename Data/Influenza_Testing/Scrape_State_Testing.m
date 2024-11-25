clear;
clc;



temp_dir=pwd;
temp_dir=temp_dir(1:(length(temp_dir)-length('Data\Influenza_Testing')));

load([temp_dir '/Data/Data_US_County.mat'],'US_County');
US_County=US_County(:,[1:8 57]);
PHL=readtable("WHO_NREVSS_Public_Health_Labs.csv");
CL=readtable("WHO_NREVSS_Clinical_Labs.csv");

PHL_Pos=sum(table2array(PHL(:,5:end)),2);
CL_Pos=sum(table2array(CL(:,[6 7])),2);

States=unique(US_County.STATE_NAME);

TEST_PER_CAPITA=zeros(height(US_County),1);
PERCENT_POS_INFV=zeros(height(US_County),1);

for ss=1:length(States)
    t_us=strcmp(States{ss},US_County.STATE_NAME);

    t_phl=strcmp(States{ss},PHL.REGION);
    t_cl=strcmp(States{ss},CL.REGION);

    test_cl=CL.TOTALSPECIMENS(t_cl);
    test_cl=test_cl(~isnan(test_cl));

    test_phl=PHL.TOTALSPECIMENS(t_phl);
    test_phl=test_phl(~isnan(test_phl));
    Tests_State=sum(test_phl)+sum(test_cl);
    
    test_cl=CL_Pos(t_cl);
    test_cl=test_cl(~isnan(test_cl));

    test_phl=PHL_Pos(t_phl);
    test_phl=test_phl(~isnan(test_phl));

    Pos_test=sum(test_phl)+sum(test_cl);

    Population_State=sum(US_County.POPULATION_SIZE_2022(t_us));
    TEST_PER_CAPITA(t_us)=Tests_State/Population_State;
    PERCENT_POS_INFV(t_us)=Pos_test./Tests_State;
end

US_County_Influenza_Test=[US_County array2table(TEST_PER_CAPITA) array2table(PERCENT_POS_INFV)];

save('State_Level_Influenza_Testing.mat','US_County_Influenza_Test');