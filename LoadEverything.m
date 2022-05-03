addpath('Costfunctions')
addpath('Data')
addpath('Models')
addpath('OptimizationScripts')
addpath('Parameters')
addpath('PlotScripts')
addpath('IQMtools')
addpath('GeneralFunctions')

load('ParameterBounds');
load('ModelStartGuess');
load('ModelValidationHealthyA');
load('Lerche2009_data');
load('Lerche2009_Param');
load('DallaMan2007_data');
load('Silfvergren2021Everyday_data')
load('Rothman1991_data.mat')
load('Taylor1996_data');
load('Firth1986_data');
load('Silfvergren2021_data');
load('Lerche2009_data');
load('Krssak2004_dataDiabetes.mat');
load('Krssak2004_dataHealthy.mat');
load('Magnusson1992_data');
load('ModelValidationHealthy_params');
load('Silfvergren2021Everyday_data');
load('OptimizeDallaManOld_params');

options2 = optimoptions(@particleswarm,'PlotFcn','pswplotbestf');
time = linspace(0,12500,12501);
pic_count = 1;

colorOrange = [0.8500, 0.3250, 0.0980];	% Orange
colorGrey   = [0.5 0.5 0.5];
LineWidthValue = 15;
MarkerSizeValue = 160;
CapSize = 30;
set(0, 'DefaultFigureRenderer', 'painters');

if ~exist('IQMsimulate','file')
    fprintf('\n\nFor these scripts to work, the IQM tools toolbox have to be compiled.\n')
    disp('To do this, a valid C-compiler is necessary')
    disp('If the compilation of the toolbox does not, it is likely that a valid C-compiler is missing.')
    choice = input('press enter to continue and compile the toolbox.');
    run('./IQMtools/installIQMtoolsInitial.m')
end
