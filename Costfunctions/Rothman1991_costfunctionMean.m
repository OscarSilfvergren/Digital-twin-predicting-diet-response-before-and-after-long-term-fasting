function Cost_Model = Rothman1991_costfunctionMean(Rothman1991,time,params,modelName,body_information, meal_information,DatapointsGly,DatapointsGNG)

model = str2func(modelName);
params=exp(params);

params = AssignParameter(params, body_information, meal_information);

try
    sim = model(time,[], params);
    
catch error
    
    sim = ones(size(time))*10^5;
    Cost_Model = 1e60;
    
    return
end


%% Data
simGly  = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
simGly  = [simGly(7245,1); simGly(7416,1); simGly(7594,1); simGly(7776,1); simGly(7941,1); simGly(8315,1); simGly(8666,1); simGly(9398,1); simGly(9747,1); simGly(10102,1); simGly(10835,1)];

GlycogenSEM = Rothman1991.glycogenpSEM(1:11);
Value = Rothman1991.glycogenpALL(1:11);

CostGly = ((Value(1:DatapointsGly) - simGly(1:DatapointsGly)).^2)./(GlycogenSEM(1:DatapointsGly).^2);

Cost_Model = nansum(CostGly,"all");

Gluconeogenesis1    = sum(sim.reactionvalues(7245:8565,ismember(sim.reactions,'Gluconeogenesis'))) + sum(sim.reactionvalues(7245:8565,ismember(sim.reactions,'EGP_Kidneys')));
Gluconeogenesis2    = sum(sim.reactionvalues(8565:9405,ismember(sim.reactions,'Gluconeogenesis'))) + sum(sim.reactionvalues(8565:9405,ismember(sim.reactions,'EGP_Kidneys')));
Gluconeogenesis3    = sum(sim.reactionvalues(9405:10485,ismember(sim.reactions,'Gluconeogenesis'))) + sum(sim.reactionvalues(9405:10485,ismember(sim.reactions,'EGP_Kidneys')));
EGP1                = sum(sim.reactionvalues(7245:8565,ismember(sim.reactions,'EGP')));
EGP2                = sum(sim.reactionvalues(8565:9405,ismember(sim.reactions,'EGP')));
EGP3                = sum(sim.reactionvalues(9405:10485,ismember(sim.reactions,'EGP')));
Procentage1         = Gluconeogenesis1/EGP1;
Procentage2         = Gluconeogenesis2/EGP2;
Procentage3         = Gluconeogenesis3/EGP3;

data       = [0.64, 0.82,  0.96];
SEM        = [0.05, 0.05, 0.01];
Simulation = [Procentage1, Procentage2, Procentage3];

CostEGP = ((data(1:DatapointsGNG) - Simulation(1:DatapointsGNG)).^2)./(SEM(1:DatapointsGNG).^2);

if Procentage3 > 1
    Cost_Model = 500;
end

Cost_Model = Cost_Model + nansum(CostEGP,"all");

end

