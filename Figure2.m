%% Simulate old Model and compare to Lerche data
load('Lerche2009_data');
load('Lerche2009_Param');

body_information = [1 , 0, 180, 85 ];  % female, male, height, weight
meal_information = [10, 0, 6.7, 73];

modelName = 'DallamanLerche';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

 Dallaman_simulation;

%% Simulate old Model and compare to Dalla Man data
load('DallaMan2007_data');

body_information  = [0 , 1, 170, 78 ];  % female, male, height, weight
meal_information  = [5 , 0, 78 , 0];

modelName = 'ModelSimpleMetabolismSteadyState';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

 DallaMan2007_Trained