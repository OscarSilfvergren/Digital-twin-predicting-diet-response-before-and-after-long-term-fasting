clear
close all
clc

LoadEverything

modelName = 'ModelSimpleMetabolismSteadyState';
model = IQMmodel([modelName '.txt']);
IQMmakeMEXmodel(model);
model = str2func(modelName);

% We advice users to plot one figure at a time due to memory requirements.
% Remove "%" sign before plot script and press "run" above to run code.
%% Plot figures from article

%  Figure2   % Plot takes about 40 sec
%  Figure3   % Plot takes about 70 sec
%  Figure4   % Plot takes about 80 sec
%  Figure5   % Plot takes about 25 sec
%  Figure6   % plot takes about 50 sec
%  Figure7   % plot takes about 30 sec







%%
%% Model estimation optimization

% ModelEstimationHealthy_Optimization
% ModelEstimationT2DM_Optimization

%% Model prediction calibration

% Rothman1991_Calibration
% Firth1986_Calibration
% Taylor1996_Calibration
% Silfvergren2021_Calibration
