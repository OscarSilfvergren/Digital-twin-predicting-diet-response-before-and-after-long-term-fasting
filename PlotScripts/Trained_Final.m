%% pre-calculations
modelName = 'ModelSimpleMetabolismSteadyState';
time2 = [time, fliplr(time)];

%% Sort all unrejected parameters for healthy estimation data based on cost

load('ModelValidationHealthy_params');
[row column] = size(ModelValidationHealthy);
ModelValidationHealthy = sortrows(ModelValidationHealthy,column);

disp(' ')
if ModelValidationHealthy(1,column) < chi2inv(0.95,190)
    disp('----- Best fit to healthy estimation data is below treshold -----')
end

fprintf('Best fit to healthy data: %.2f, Statistical Limit: %.2f (dgf = %i)', ModelValidationHealthy(1,column), chi2inv(0.95,190), 190)

%% Healthy Krsak

for i = 1:row
    optimizedParamTemp = ModelValidationHealthy(i,1:(column-1));
    optimizedParamTemp(99:101)    = [1,0,0];                         % Steady State Meals are solid and there is no meal in t0
    optimizedParamTemp(102:105)   = [0 , 1, 175, 75];                % female, male, height, weight
    optimizedParamTemp(106:109)   = [1,0,0,0];                       % This is study A
    optimizedParamTemp(110:117)   = [10, 0, 87, 23,1,0,0,0];         % Meal Information
    try
        sim = model(time,[],optimizedParamTemp);
        
    catch error
    end
    
    if i == 1
        simBest = sim;
    end
    
    %Glucose
    if i == 1
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18; % /18 to change unit to mM from mg/dl
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
    else
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = max(maxG2,maxG1);
        minG2 = min(minG2,minG1);
    end
    
    % Insulin
    if i == 1
        maxI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        minI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        maxI2 = sim.variablevalues(:,ismember(sim.variables,'I'));
        minI2 = sim.variablevalues(:,ismember(sim.variables,'I'));
    else
        maxI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        minI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        maxI2 = max(maxI2,maxI1);
        minI2 = min(minI2,minI1);
    end
    
    % Glycogen
    if i == 1
        maxGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        maxGly2 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly2 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
    else
        maxGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        maxGly2 = max(maxGly2,maxGly1);
        minGly2 = min(minGly2,minGly1);
    end
    
    % EGP
    if i == 1
        maxEGP1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
        minEGP1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
        maxEGP2 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
        minEGP2 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
    else
        maxEGP1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
        minEGP1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
        maxEGP2 = max(maxEGP2,maxEGP1);
        minEGP2 = min(minEGP2,minEGP1);
    end
end
%%
% Glycogen
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,10,20],'ytick',[100,200,300],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-127,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004_Healthy.timeGlycogen(1:22)/60-127,Krssak2004_Healthy.glycogen(1:22),Krssak2004_Healthy.glycogenSEM(1:22),' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([1 1.2], [150 150],'Color','k','LineWidth',8);
xlim([0 20])
ylim([150 360])
hold off
%%
% EGP
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,5,10,15,20],'ytick',[0,0.5,1,1.5,2,2.5],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxEGP2', fliplr(minEGP2')];
fill(time2/60-127,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127,simBest.reactionvalues(:,ismember(simBest.reactions,'EGP')),'b-.','LineWidth',LineWidthValue);
ylabel({'EGP' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004_Healthy.timeEGP(1:20)/60-127,Krssak2004_Healthy.EGP(1:20),Krssak2004_Healthy.EGPSEM(1:20),' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([1 1.2], [0.4 0.4],'Color','k','LineWidth',8);
xlim([0 10])
ylim([0.4 1.85])
hold off

% G
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,5,10,15,20],'ytick',[0,5,10,15,20],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxG2', fliplr(minG2')];
fill(time2/60-127,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue)
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004_Healthy.timeGlucose(1:25)/60-127,Krssak2004_Healthy.glucose(1:25)/18,Krssak2004_Healthy.glucoseSEM(1:25)/18,' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([1 1.2], [2 2],'Color','k','LineWidth',8);
xlim([0 10])
ylim([2 11])
hold off

% I
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,5,10,15,20],'ytick',[0,200,400,600,800,1000,1200],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxI2', fliplr(minI2')];
fill(time2/60-127,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBest.variablevalues(:,ismember(simBest.variables,'I')),'b-.','LineWidth',LineWidthValue)
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004_Healthy.timeInsulin(1:10)/60-127,Krssak2004_Healthy.insulin(1:10),Krssak2004_Healthy.insulinSEM(1:10),' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([1 1.2], [0 0],'Color','k','LineWidth',8);
xlim([0 10])
ylim([0 850])
hold off

%% Lerche Healthy

for i = 1:row
    optimizedParamTemp = ModelValidationHealthy(i,1:(column-1));
    optimizedParamTemp(99:101)    = [1,0,0];                         % Steady State Meals are solid and there is no meal in t0
    optimizedParamTemp(102:105)   = [0 , 1, 175, 75];                % female, male, height, weight
    optimizedParamTemp(106:109)   = [0,1,0,0];                       % This is study B
    optimizedParamTemp(110:117)   = [10, 2420, 75, 0, 1, 0,0,0];     % Meal Information
    
    try
        sim = model(time,[],optimizedParamTemp);
    catch error
    end
    
    if i == 1
        simBest = sim;
    end
    
    %Glucose
    if i == 1
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18; % /18 to change unit to mM from mg/dl
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
    else
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = max(maxG2,maxG1);
        minG2 = min(minG2,minG1);
    end
    
    % Insulin
    if i == 1
        maxI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        minI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        maxI2 = sim.variablevalues(:,ismember(sim.variables,'I'));
        minI2 = sim.variablevalues(:,ismember(sim.variables,'I'));
    else
        maxI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        minI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        maxI2 = max(maxI2,maxI1);
        minI2 = min(minI2,minI1);
    end
end

% Glucose
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[0,3,4,5,6,8],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxG2', fliplr(minG2')];
fill(time2/1440-5,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/1440-5,simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue)
errorbar(Lerche2009_data.time_glucose_healthy/1440-5,Lerche2009_data.glucose_healthy/18,Lerche2009_data.glucoseSEM_healthy/18,' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
%xlim([46 52])
xlim([0 2])
ylim([3 5.2])
hold off

% Insulin
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,1,2],'ytick',[0,20,40],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxI2', fliplr(minI2')];
fill(time2/1440-5,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/1440-5,simBest.variablevalues(:,ismember(simBest.variables,'I')),'b-.','LineWidth',LineWidthValue)
errorbar(Lerche2009_data.time_insulin_healthy/1440-5,Lerche2009_data.insulin_healthy,Lerche2009_data.insulinSEMFIXED_healthy,' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (days)",'FontSmoothing','on','fontname','Arial');
xlim([0 2])
ylim([0 40])
hold off

% Insulin meal
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,5,10],'ytick',[0,300,600],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxI2', fliplr(minI2')];
fill(time2/60-164.5,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-164.5,simBest.variablevalues(:,ismember(simBest.variables,'I')),'b-.','LineWidth',LineWidthValue)
errorbar(Lerche2009_data.time_insulin_healthy/60-164.5,Lerche2009_data.insulin_healthy,Lerche2009_data.insulinSEMFIXED_healthy,' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
line([4 4.2], [2 2],'Color','k','LineWidth',6);
xlim([0 10])
ylim([0 600])
hold off

% Glucose meal
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,5,10],'ytick',[0,2,6,10],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxG2', fliplr(minG2')];
fill(time2/60-164.5,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-164.5,simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue)
errorbar(Lerche2009_data.time_glucose_healthy/60-164.5,Lerche2009_data.glucose_healthy/18,Lerche2009_data.glucoseSEM_healthy/18,' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
line([4 4.2], [2 2],'Color','k','LineWidth',6);
xlim([0 10])
ylim([2 12])
hold off

%% Healthy magnusson
clear ratio_Gluconeogenesis

for i = 1:row
    optimizedParamTemp = ModelValidationHealthy(i,1:(column-1));
    optimizedParamTemp(99:101)     = [1,0,0];                         % Steady State Meals are solid and there is no meal in t0
    optimizedParamTemp(102:105)    = [0 , 1, 175, 75 ];               % female, male, height, weight
    optimizedParamTemp(106:109)    = [0,0,1,0];                       % This is study C
    optimizedParamTemp(110:117)    = [5, 0, 98.3125, 26, 1, 0,0,0];   % Meal Information
    try
        sim = model(time,[],optimizedParamTemp);
        
    catch error
    end
    
    if i == 1
        simBest = sim;
    end
    
    if i == 1
        maxGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        maxGly2 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly2 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
    else
        maxGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        maxGly2 = max(maxGly2,maxGly1);
        minGly2 = min(minGly2,minGly1);
    end
    
    Gluconeogenesis          = sum(sim.reactionvalues(7914:9294,ismember(sim.reactions,'Gluconeogenesis'))) + sum(simBest.reactionvalues(7914:9294,ismember(simBest.reactions,'EGP_Kidneys')));
    EGP                      = sum(sim.reactionvalues(7914:9294,ismember(sim.reactions,'EGP')));
    ratio_Gluconeogenesis(i) = Gluconeogenesis/EGP;
    
end

% Glycogen
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,12,24],'ytick',[100,250,400],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-130,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-130, sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Magnusson1992_data.time_glycogen_healthy/60-130,Magnusson1992_data.glycogen_healthy, Magnusson1992_data.glycogenSEM_healthy,' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
xlim([0 24])
ylim([80 400])
hold off

X = categorical({'Simulation', 'Data'});
X = reordercats(X,{'Simulation', 'Data'});

ratio_Gluconeogenesis = BestMinMax(ratio_Gluconeogenesis,70/100);

% All
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'ytick',[0,50,100],'FontSize', 65,'fontname','Arial')%,'FontSmoothing','on')
bar(X,[ratio_Gluconeogenesis(1);70] ,'b','FaceAlpha',0.4,'EdgeAlpha',0.6);
errorbar(X,[ratio_Gluconeogenesis(1);70],[ratio_Gluconeogenesis(2);6],' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
ylabel({'Gluconeogenesis' ; 'contribution to EGP (%)'},'FontSmoothing','on','fontname','Arial');
ylim([0 100])
hold off

%% Silfvergren

for i = 1:row
    optimizedParamTemp = ModelValidationHealthy(i,1:(column-1));
    optimizedParamTemp(99:101)    = [1,0,0];                           % Steady State Meals are solid and there is no meal in t0
    optimizedParamTemp(102:105)   = [0 , 1, 175, 75 ];                 % female, male, height, weight
    optimizedParamTemp(106:109)   = [0,0,0,1];                         % This is study D
    optimizedParamTemp(110:117)   = [10, 0, 155.35, 20.7, 1, 0,0,0];   % Meal Information
    
    try
        sim = model(time,[],optimizedParamTemp);
    catch error
    end
    
    if i == 1
        simBest = sim;
    end
    
    %Glucose
    if i == 1
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18; % /18 to change unit to mM from mg/dl
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
    else
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = max(maxG2,maxG1);
        minG2 = min(minG2,minG1);
    end
    
end

% Glucose
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,2,4,6,8],'ytick',[2,6,10],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxG2', fliplr(minG2')];
fill(time2/60-127,y,'b','LineStyle','none','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127,simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue) %155
plot(Silfvergren2021Everyday_data.time/60-127,Silfvergren2021Everyday_data.value,'kx','MarkerSize',20,'LineWidth',6);
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
xlim([0 7])
line([1 1.2], [2 2],'Color','k','LineWidth',8);
ylim([2 10])
hold off

%% pre-calculations T2D
modelName = 'ModelSimpleMetabolismSteadyState';

load('ModelValidationT2D')
[row column] = size(ModelValidationT2D);
ModelValidationT2D = sortrows(ModelValidationT2D,column);

disp(' ')
disp(' ')
if ModelValidationT2D(1,column) < chi2inv(0.95,83)
    disp('----- Best fit to diabetic estimation data is below treshold -----')
end

fprintf('Best fit to diabetic data: %.2f, Statistical Limit: %.2f (dgf = %i)', ModelValidationT2D(1,column), chi2inv(0.95,83), 83)

%% T2DM Krsak population

for i = 1:row
    optimizedParamTemp = ModelValidationT2D(i,1:(column-1));
    
    optimizedParamTemp(99:101)  = [1,0,0];                           % Steady State Meals are solid and there is no meal in t0
    optimizedParamTemp(102:105)   = [0 , 1, 175, 75];                % female, male, height, weight
    optimizedParamTemp(106:109)   = [1,0,0,0];                       % This is study A
    optimizedParamTemp(110:117)   = [10, 0, 87, 23,1,0,0,0];         % Meal Information
    
    try
        sim = model(time,[],optimizedParamTemp);
    catch error
    end
    
    if i == 1
        simBest = sim;
    end
    
    %Glucose
    if i == 1
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18; % /18 to change unit to mM from mg/dl
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG2 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
    else
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        minG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
        maxG2 = max(maxG2,maxG1);
        minG2 = min(minG2,minG1);
    end
    
    % Insulin
    if i == 1
        maxI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        minI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        maxI2 = sim.variablevalues(:,ismember(sim.variables,'I'));
        minI2 = sim.variablevalues(:,ismember(sim.variables,'I'));
    else
        maxI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        minI1 = sim.variablevalues(:,ismember(sim.variables,'I'));
        maxI2 = max(maxI2,maxI1);
        minI2 = min(minI2,minI1);
    end
    
    % Glycogen
    if i == 1
        maxGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        maxGly2 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly2 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
    else
        maxGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        maxGly2 = max(maxGly2,maxGly1);
        minGly2 = min(minGly2,minGly1);
    end
    
    % EGP
    if i == 1
        maxEGP1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
        minEGP1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
        maxEGP2 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
        minEGP2 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
    else
        maxEGP1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
        minEGP1 = sim.reactionvalues(:,ismember(sim.reactions,'EGP'));
        maxEGP2 = max(maxEGP2,maxEGP1);
        minEGP2 = min(minEGP2,minEGP1);
    end
end

% Glycogen
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,10,20],'ytick',[0,150,225,300],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-127,y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'r-.','LineWidth',LineWidthValue);
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004_Diabetes.timeGlycogen(1:22)/60-127,Krssak2004_Diabetes.glycogen(1:22),Krssak2004_Diabetes.glycogenSEM(1:22),' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
line([1 1.2], [100 100],'Color','k','LineWidth',8);
xlim([0 20])
ylim([150 300])
hold off

% EGP
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,10,20],'ytick',[0,1,2],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxEGP2', fliplr(minEGP2')];
fill(time2/60-127,y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127,simBest.reactionvalues(:,ismember(simBest.reactions,'EGP')),'r-.','LineWidth',LineWidthValue)
ylabel({'EGP' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004_Diabetes.timeEGP(1:20)/60-127,Krssak2004_Diabetes.EGP(1:20),Krssak2004_Diabetes.EGPSEM(1:20),' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
line([1 1.2], [0 0],'Color','k','LineWidth',8);
xlim([0 10])
ylim([0 2.2])

% G
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,5,10,15,20,25],'ytick',[5,10,15,20],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxG2', fliplr(minG2')];
fill(time2/60-127,y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'r-.','LineWidth',LineWidthValue)
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004_Diabetes.timeGlucose(1:25)/60-127,Krssak2004_Diabetes.glucose(1:25)/18,Krssak2004_Diabetes.glucoseSEM(1:25)/18,' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
line([1 1.2], [5 5],'Color','k','LineWidth',8);
xlim([0 10])
ylim([5 20])
hold off

% I
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,5,10,15,20,25],'ytick',[0,200,400],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxI2', fliplr(minI2')];
fill(time2/60-127,y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBest.variablevalues(:,ismember(simBest.variables,'I')),'r-.','LineWidth',LineWidthValue)
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Krssak2004_Diabetes.timeInsulin(1:10)/60-127,Krssak2004_Diabetes.insulin(1:10),Krssak2004_Diabetes.insulinSEM(1:10),' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
line([1 1.2], [0 0],'Color','k','LineWidth',8);
xlim([0 10])
ylim([0 400])
hold off


%% T2DM magnusson
clear ratio_Gluconeogenesis

for i = 1:row
    optimizedParamTemp = ModelValidationT2D((i),1:(column-1));
    optimizedParamTemp(99:101)     = [1,0,0];                         % Steady State Meals are solid and there is no meal in t0
    optimizedParamTemp(102:105)    = [0 , 1, 175, 75 ];               % female, male, height, weight
    optimizedParamTemp(106:109)    = [0,0,1,0];                       % This is study C
    optimizedParamTemp(110:117)    = [5, 0, 98.3125, 26, 1, 0,0,0];   % Meal Information
    
    try
        sim = model(time,[],optimizedParamTemp);
        
    catch error
    end
    
    if i == 1
        simBest = sim;
    end
    
    if i == 1
        maxGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        maxGly2 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly2 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
    else
        maxGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        minGly1 = sim.variablevalues(:,ismember(sim.variables,'Glycogen_liver'));
        maxGly2 = max(maxGly2,maxGly1);
        minGly2 = min(minGly2,minGly1);
    end
    
    Gluconeogenesis          = sum(sim.reactionvalues(7914:9294,ismember(sim.reactions,'Gluconeogenesis'))) + sum(simBest.reactionvalues(7914:9294,ismember(simBest.reactions,'EGP_Kidneys')));
    EGP                      = sum(sim.reactionvalues(7914:9294,ismember(sim.reactions,'EGP')));
    ratio_Gluconeogenesis(i) = Gluconeogenesis/EGP;
end

% Glycogen
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,12,24],'ytick',[0,100,200,300,400],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-130,y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-130, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'r-.','LineWidth',LineWidthValue);
errorbar(Magnusson1992_data.time_glycogen_diabetes/60-130,Magnusson1992_data.glycogen_diabetes, Magnusson1992_data.glycogenSEM_diabetes,' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
xlim([0 24])
ylim([0 200])
hold off

% GNG
X = categorical({'Simulation', 'Data'});
X = reordercats(X,{'Simulation', 'Data'});

ratio_Gluconeogenesis = BestMinMax(ratio_Gluconeogenesis,88/100);
%%
% All
figure('Name', " ", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'ytick',[0,50,100],'FontSize', 65)%,'FontSmoothing','on')
ratio_EGP = 1;
bar(X,[ratio_Gluconeogenesis(1);88],'r');
errorbar(X,[ratio_Gluconeogenesis(1);88],[ratio_Gluconeogenesis(2);2],' k .','MarkerSize',1,'LineWidth',10,'CapSize',20');
ylabel({'Gluconeogenesis' ; ' contribution to EGP %'},'FontSmoothing','on','fontname','Arial');
ylim([0 100])
hold off