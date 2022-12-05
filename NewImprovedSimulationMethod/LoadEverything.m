clear
close all
clc
format longG

addpath('Data')
addpath('IQMtools')

load('Rothman1991_data.mat')
load('Silfvergren2021_data');
load('Krssak2004Healthy_data.mat');

options = optimoptions(@particleswarm,'PlotFcn','pswplotbestf','MaxStallIterations',300,'FunctionTolerance',1.0e-7,'UseParallel', true,'UseVectorized', false);
i = 1;

if ~exist('IQMsimulate','file')
    fprintf('\n\nFor these scripts to work, the IQM tools toolbox have to be compiled.\n')
    disp('To do this, a valid C-compiler is necessary')
    disp('If the compilation of the toolbox does not, it is likely that a valid C-compiler is missing.')
    choice = input('press enter to continue and compile the toolbox.');
    run('./IQMtools/installIQMtoolsInitial.m')
end
