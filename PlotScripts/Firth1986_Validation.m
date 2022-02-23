%% Pre-Calculations
modelName = 'ModelSimpleMetabolismSteadyState';
time2 = [time, fliplr(time)];

body_information = [1 , 0, 175, 90 ];  % female, male, height, weight
meal_information = [10, 0, 50, 0];

%% Healthy
 load('Firth1986_ParamHealthyCalibratedDone');

[row column] = size(Firth1986_ParamHealthyCalibrated);
clear optimizedParamTemp2

for i = 1:row
    optimizedParamTemp  = Firth1986_ParamHealthyCalibrated(i,1:(column-1));
    cost                = Firth1986_costfunctionHealthy(Firth1986_data,time,log(optimizedParamTemp),modelName,meal_information,body_information,14,15,15);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Firth1986_ParamHealthyCalibrated = sortrows(optimizedParamTemp2,column);

disp(' ')
if Firth1986_ParamHealthyCalibrated(1,column) < chi2inv(0.95,44)
    disp('----- Best prediction of Firth data is below treshold -----')
end

fprintf('Best fit Firth data: %.2f, Statistical Limit: %.2f (dgf = %i)', Firth1986_ParamHealthyCalibrated(1,column), chi2inv(0.95,44), 44)
disp(' ')

%%
for i = 1:row
    optimizedParamTemp = Firth1986_ParamHealthyCalibrated(i,1:(column-1));
    sim = model(time,[],optimizedParamTemp);
    
    if i==1
        simBest = model(time,[],optimizedParamTemp);
    end
    
    % Glucose
    if i == 1
        maxG1 = sim.variablevalues(:,ismember(sim.variables,'G'))/18;
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

% EGP
figure('Name', "EGP", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,4,8],'ytick',[0,1,2,],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y   = [maxEGP2', fliplr(minEGP2')];
fill(time2/60-127,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127,simBest.reactionvalues(:,ismember(simBest.reactions,'EGP')),'b-.','LineWidth',LineWidthValue)
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'EGP' ; '(mg/kg/min)'},'FontSmoothing','on','fontname','Arial');
errorbar(Firth1986_data.time_EGP_healthy(1:3)/60-127,Firth1986_data.EGP_healthy((1:3)),Firth1986_data.EGPSEM_healthy((1:3)),' k x','MarkerSize',70,'LineWidth',10,'CapSize',CapSize');
errorbar(Firth1986_data.time_EGP_healthy(1:14)/60-127,Firth1986_data.EGP_healthy(1:14),Firth1986_data.EGPSEM_healthy(1:14),' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([1 1.16], [0 0],'Color','k','LineWidth',12);
xlim([0 8.5])
ylim([0 2.5])
hold off
%%
% Glucose
figure('Name', "Glucose", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,4,8],'ytick',[0,4,8],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxG2', fliplr(minG2')];
fill(time2/60-127,y ,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBest.variablevalues(:,ismember(simBest.variables,'G'))/18,'b-.','LineWidth',LineWidthValue)
errorbar(Firth1986_data.time_glucose_healthy((1:5))/60-127,Firth1986_data.glucose_healthy((1:5))/18,Firth1986_data.glucoseSEM_healthy((1:5))/18,' k x','MarkerSize',70,'LineWidth',10,'CapSize',CapSize');
errorbar(Firth1986_data.time_glucose_healthy(1:15)/60-127,Firth1986_data.glucose_healthy(1:15)/18,Firth1986_data.glucoseSEM_healthy(1:15)/18,' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
line([1 1.16], [2 2],'Color','k','LineWidth',12);
xlim([0 8.5])
hold off
%%
% Insulin
figure('Name', "Insulin", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,4,8],'ytick',[0,200,400],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxI2', fliplr(minI2')];
fill(time2/60-127,y ,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBest.variablevalues(:,ismember(simBest.variables,'I')),'b-.','LineWidth',LineWidthValue)
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on','fontname','Arial');
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
errorbar(Firth1986_data.time_insulin_healthy((1:5))/60-127,Firth1986_data.insulin_healthy((1:5)),Firth1986_data.insulinSEM_healthy((1:5)),' k x','MarkerSize',70,'LineWidth',10,'CapSize',CapSize');
errorbar(Firth1986_data.time_insulin_healthy(1:15)/60-127,Firth1986_data.insulin_healthy(1:15),Firth1986_data.insulinSEM_healthy(1:15),' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([1 1.16], [0 0],'Color','k','LineWidth',12);
xlim([0 8.5])
hold off
