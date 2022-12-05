LoadEverything

%% Mex and index

modelName = 'SimpleHepaticMetabolism';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

[pNames, params]               = IQMparameters(model);
[sNames, ODEs, initialvalues]  = IQMstates(model);

%% ParameterEstimation

% [minCostPS] = ParameterEstimation(i,Krssak2004Healthy_data,Rothman1991_data,Silfvergren2021_data,params,pNames,initialvalues,sNames,model,options);

%% Plot

 PlotEstimationStudies

%% Kill all paralell workers and clusters

poolobj = gcp('nocreate');
delete(poolobj)
myCluster = parcluster('local');
delete(myCluster.Jobs)