%% Pre-Calculations
modelName = 'ModelSimpleMetabolismSteadyState';
time2 = [time, fliplr(time)];

body_information = [0 , 1, 180, 75 ];  % female, male, height, weight
meal_information = [1, 0, 0, 0];

%% P1
 load('Rothman1991P1CalibratedDone');

[row column] = size(Rothman1991_paramP1Calibrated);
clear optimizedParamTemp2

for i = 1:row
    optimizedParamTemp  = Rothman1991_paramP1Calibrated((i),1:(column-1));
    optimizedParamTemp  = log(optimizedParamTemp);
    cost                = Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp1(1:9),Rothman1991_data.timep1(1:9),time,optimizedParamTemp,modelName,body_information, meal_information);
    optimizedParamTemp  = exp(optimizedParamTemp);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Rothman1991_paramP1Calibrated = sortrows(optimizedParamTemp2,column);


for i = 1:row
    optimizedParamTemp  = Rothman1991_paramP1Calibrated(i,1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
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
    
end


% Gly
figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,30,60],'ytick',[0,250,500],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-120,y ,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-120, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot(Rothman1991_data.timep1(1:1)/60-120,Rothman1991_data.glycogenp1(1:1),'kx','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep1(1:10)/60-120,Rothman1991_data.glycogenp1(1:10),'k.','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
xlim([0 62])
ylim([0 500])
hold off

%% P2
load('Rothman1991P2CalibratedDone');

[row column] = size(Rothman1991_paramP2Calibrated);
clear optimizedParamTemp2

for i = 1:row
    optimizedParamTemp  = Rothman1991_paramP2Calibrated(i,1:(column-1));
    optimizedParamTemp  = log(optimizedParamTemp);
    cost                = Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp2(1:11),Rothman1991_data.timep2(1:11),time,optimizedParamTemp,modelName,body_information, meal_information);
    optimizedParamTemp  = exp(optimizedParamTemp);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Rothman1991_paramP2Calibrated = sortrows(optimizedParamTemp2,column);


for i = 1:row
    
    optimizedParamTemp = Rothman1991_paramP2Calibrated(i,1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
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
end

figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,30,60],'ytick',[0,250,500],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-120,y ,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-120, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot(Rothman1991_data.timep2(1:1)/60-120,Rothman1991_data.glycogenp2(1:1),'kx','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep2(1:11)/60-120,Rothman1991_data.glycogenp2(1:11),'k.','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
xlim([0 62])
ylim([0 500])
hold off

%% P3
load('Rothman1991P3CalibratedDone');

[row column] = size(Rothman1991_paramP3Calibrated);
clear optimizedParamTemp2

for i = 1:row
    optimizedParamTemp  = Rothman1991_paramP3Calibrated(i,1:(column-1));
    optimizedParamTemp  = log(optimizedParamTemp);
    cost                = Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp3(1:10),Rothman1991_data.timep3(1:10),time,optimizedParamTemp,modelName,body_information, meal_information);
    optimizedParamTemp  = exp(optimizedParamTemp);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Rothman1991_paramP3Calibrated = sortrows(optimizedParamTemp2,column);


for i = 1:row
    
    optimizedParamTemp = Rothman1991_paramP3Calibrated(i,1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
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
end

figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,30,60],'ytick',[0,250,500],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-120,y ,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-120, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot(Rothman1991_data.timep3(1)/60-120,Rothman1991_data.glycogenp3(1),'kx','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep3(1:10)/60-120,Rothman1991_data.glycogenp3(1:10),'k.','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
xlim([0 62])
ylim([0 500])
hold off

%% P4
load('Rothman1991P4CalibratedDone');


[row column] = size(Rothman1991_paramP4Calibrated);
clear optimizedParamTemp2

for i = 1:row
    optimizedParamTemp  = Rothman1991_paramP4Calibrated(i,1:(column-1));
    optimizedParamTemp  = log(optimizedParamTemp);
    cost                = Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp4(1:10),Rothman1991_data.timep4(1:10),time,optimizedParamTemp,modelName,body_information, meal_information);
    optimizedParamTemp  = exp(optimizedParamTemp);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Rothman1991_paramP4Calibrated = sortrows(optimizedParamTemp2,column);


for i = 1:row
    
    optimizedParamTemp = Rothman1991_paramP4Calibrated(i,1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
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
end

figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,30,60],'ytick',[0,250,500],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-120,y ,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-120, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot(Rothman1991_data.timep4(1:1)/60-120,Rothman1991_data.glycogenp4(1:1),'kx','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep4(1:11)/60-120,Rothman1991_data.glycogenp4(1:11),'k.','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
xlim([0 62])
ylim([0 500])
hold off

%% P5
load('Rothman1991P5CalibratedDone');

[row column] = size(Rothman1991_paramP5Calibrated);
clear optimizedParamTemp2

for i = 1:row
    optimizedParamTemp  = Rothman1991_paramP5Calibrated((i),1:(column-1));
    optimizedParamTemp  = log(optimizedParamTemp);
    cost                = Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp5(1:8),Rothman1991_data.timep5(1:8),time,optimizedParamTemp,modelName,body_information, meal_information);
    optimizedParamTemp  = exp(optimizedParamTemp);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Rothman1991_paramP5Calibrated = sortrows(optimizedParamTemp2,column);


for i = 1:row
    
    optimizedParamTemp = Rothman1991_paramP5Calibrated((i),1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
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
end

figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,30,60],'ytick',[0,250,500],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-120,y ,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-120, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot(Rothman1991_data.timep5(1:1)/60-120,Rothman1991_data.glycogenp5(1:1),'kx','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep5(1:8)/60-120,Rothman1991_data.glycogenp5(1:8),'k.','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
xlim([0 62])
ylim([0 500])
hold off

%% P6
load('Rothman1991P6CalibratedDone');

[row column] = size(Rothman1991_paramP6Calibrated);
clear optimizedParamTemp2

for i = 1:row
    optimizedParamTemp  = Rothman1991_paramP6Calibrated((i),1:(column-1));
    optimizedParamTemp  = log(optimizedParamTemp);
    cost                = Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp6(1:11),Rothman1991_data.timep6(1:11),time,optimizedParamTemp,modelName,body_information, meal_information);
    optimizedParamTemp  = exp(optimizedParamTemp);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Rothman1991_paramP6Calibrated = sortrows(optimizedParamTemp2,column);

for i = 1:row
    
    optimizedParamTemp = Rothman1991_paramP6Calibrated((i),1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
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
end

figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,30,60],'ytick',[0,250,500],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-120,y ,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-120, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot(Rothman1991_data.timep6(1:1)/60-120,Rothman1991_data.glycogenp6(1:1),'kx','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep6(1:11)/60-120,Rothman1991_data.glycogenp6(1:11),'k.','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
xlim([0 62])
ylim([0 500])
hold off

%% P7
load('Rothman1991P7CalibratedDone');

[row column] = size(Rothman1991_paramP7Calibrated);
clear optimizedParamTemp2

for i = 1:row
    optimizedParamTemp  = Rothman1991_paramP7Calibrated((i),1:(column-1));
    optimizedParamTemp  = log(optimizedParamTemp);
    cost                = Rothman1991_costfunctionIndividual(Rothman1991_data.glycogenp7(1:9),Rothman1991_data.timep7(1:9),time,optimizedParamTemp,modelName,body_information, meal_information);
    optimizedParamTemp  = exp(optimizedParamTemp);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Rothman1991_paramP7Calibrated = sortrows(optimizedParamTemp2,column);


for i = 1:row
    
    optimizedParamTemp = Rothman1991_paramP7Calibrated((i),1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
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
end

figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,30,60],'ytick',[0,250,500],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-120,y ,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-120, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
plot(Rothman1991_data.timep7(1)/60-120,Rothman1991_data.glycogenp7(1),'kx','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
plot(Rothman1991_data.timep7(1:9)/60-120,Rothman1991_data.glycogenp7(1:9),'k.','MarkerSize',MarkerSizeValue,'LineWidth',LineWidthValue)
xlim([0 62])
ylim([0 500])
hold off

%% Pmean
load('Rothman1991MeanCalibrated');

[row column] = size(Rothman1991_paramPMeanCalibrated);
clear optimizedParamTemp2

for i = 1:row
    optimizedParamTemp  = Rothman1991_paramPMeanCalibrated((i),1:(column-1));
    optimizedParamTemp  = log(optimizedParamTemp);
    cost                = Rothman1991_costfunctionMean(Rothman1991_data,time,optimizedParamTemp,modelName,body_information, meal_information,11,3);
    optimizedParamTemp  = exp(optimizedParamTemp);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Rothman1991_paramPMeanCalibrated = sortrows(optimizedParamTemp2,column);

disp(' ')
if Rothman1991_paramPMeanCalibrated(1,column) < chi2inv(0.95,14)
    disp('----- Best prediction of Rothman data is below treshold -----')
end

fprintf('Best prediction of Rothman data: %.2f, Statistical Limit: %.2f (dgf = %i)', Rothman1991_paramPMeanCalibrated(1,column), chi2inv(0.95,14), 14)
disp(' ')

for i = 1:row
    optimizedParamTemp = Rothman1991_paramPMeanCalibrated(i,1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
    if i == 1
        simBest = sim;
    end
    
    Gluconeogenesis1    = sum(sim.reactionvalues(7245:8565,ismember(sim.reactions,'Gluconeogenesis'))) + sum(sim.reactionvalues(7245:8565,ismember(sim.reactions,'EGP_Kidneys')));
    Gluconeogenesis2    = sum(sim.reactionvalues(8565:9405,ismember(sim.reactions,'Gluconeogenesis'))) + sum(sim.reactionvalues(8565:9405,ismember(sim.reactions,'EGP_Kidneys')));
    Gluconeogenesis3    = sum(sim.reactionvalues(9405:10485,ismember(sim.reactions,'Gluconeogenesis'))) + sum(sim.reactionvalues(9405:10485,ismember(sim.reactions,'EGP_Kidneys')));
    EGP1                = sum(sim.reactionvalues(7245:8565,ismember(sim.reactions,'EGP')));
    EGP2                = sum(sim.reactionvalues(8565:9405,ismember(sim.reactions,'EGP')));
    EGP3                = sum(sim.reactionvalues(9405:10485,ismember(sim.reactions,'EGP')));
    Procentage1         = Gluconeogenesis1/EGP1;
    Procentage2         = Gluconeogenesis2/EGP2;
    Procentage3         = Gluconeogenesis3/EGP3;
    Procentage          = [Procentage1, Procentage2, Procentage3];
    
    if i == 1
        maxsim1 = Procentage;
        minsim1 = Procentage;
        maxsim2 = Procentage;
        minsim2 = Procentage;
    else
        maxsim1 = Procentage;
        minsim1 = Procentage;
        maxsim2 = max(maxsim2,maxsim1);
        minsim2 = min(minsim2,minsim1);
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
    
end


%Glycogen
figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,30,60],'ytick',[0,250,500],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-120,y ,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-120, simBest.variablevalues(:,ismember(simBest.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
GlycogenSEM = Rothman1991_data.glycogenpSEM;
TimeSEM = Rothman1991_data.timepSEM/60;
errorbar(Rothman1991_data.timepALL(1)/60-120,Rothman1991_data.glycogenpALL(1),GlycogenSEM(1),' k x','MarkerSize',60,'LineWidth',10,'CapSize',40');
errorbar(Rothman1991_data.timepALL(2:10)/60-120,Rothman1991_data.glycogenpALL(2:10),GlycogenSEM(2:10),' k .','MarkerSize',1,'LineWidth',10,'CapSize',40');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlim([0 62])
ylim([0 500])
hold off

%Gluconeogenesis
temptime   = [8565, 9405, 10485];
temptime   = [temptime, fliplr(temptime)];

figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[20,40,60],'ytick',[0,20,40,60,80,100],'FontSize', 70,'FontSmoothing','on','fontname','Arial')
y        = [maxsim2, fliplr(minsim2)];
fill(temptime/60-120,y*100 ,'b','FaceAlpha',0.2,'EdgeAlpha',0);

Gluconeogenesis1    = sum(simBest.reactionvalues(7245:8565,ismember(simBest.reactions,'Gluconeogenesis'))) + sum(simBest.reactionvalues(7245:8565,ismember(simBest.reactions,'EGP_Kidneys')));
Gluconeogenesis2    = sum(simBest.reactionvalues(8565:9405,ismember(simBest.reactions,'Gluconeogenesis'))) + sum(simBest.reactionvalues(8565:9405,ismember(simBest.reactions,'EGP_Kidneys')));
Gluconeogenesis3    = sum(simBest.reactionvalues(9405:10485,ismember(simBest.reactions,'Gluconeogenesis'))) + sum(simBest.reactionvalues(9405:10485,ismember(simBest.reactions,'EGP_Kidneys')));
EGP1                = sum(simBest.reactionvalues(7245:8565,ismember(simBest.reactions,'EGP')));
EGP2                = sum(simBest.reactionvalues(8565:9405,ismember(simBest.reactions,'EGP')));
EGP3                = sum(simBest.reactionvalues(9405:10485,ismember(simBest.reactions,'EGP')));
Procentage1         = Gluconeogenesis1/EGP1;
Procentage2         = Gluconeogenesis2/EGP2;
Procentage3         = Gluconeogenesis3/EGP3;
Procentage          = [Procentage1, Procentage2, Procentage3];
temptime            = [8565, 9405, 10485];
plot(temptime/60-120, Procentage*100,'b-.','LineWidth',LineWidthValue)

errorbar([8565]/60-120,[0.64]*100,[0.05]*100,' k x','MarkerSize',60,'LineWidth',10,'CapSize',40');
errorbar([9405, 10485]/60-120,[0.82, 0.96]*100,[0.05, 0.01]*100,' k .','MarkerSize',1,'LineWidth',10,'CapSize',40');
xlim([20 60])
ylim([40 100])
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Gluconeogenesis' ; 'contribution to EGP(%)'},'FontSmoothing','on','fontname','Arial');
hold off