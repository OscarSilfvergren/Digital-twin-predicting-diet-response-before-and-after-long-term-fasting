load('Lerche2009_data');
load('Krssak2004_dataDiabetes.mat');
load('Krssak2004_dataHealthy.mat');
load('Magnusson1992_data');
load('ModelValidationHealthy_params');
load('Silfvergren2021Everyday_data');

modelName = 'ModelSimpleMetabolismSteadyState';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);
  %%
  Trained_Final
  