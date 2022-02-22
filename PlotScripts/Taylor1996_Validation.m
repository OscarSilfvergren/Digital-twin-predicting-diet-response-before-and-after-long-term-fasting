%% Pre-Calculations
load('Taylor1996_ParamHealthyCalibrated');

body_information = [0 , 1, 180, 75.6 ];  % female, male, height, weight
meal_information = [10, 0, 138.64, 29.25];

[row column] = size(Taylor1996_ParamHealthyCalibrated);

modelName = 'ModelSimpleMetabolismSteadyState';
time2 = [time, fliplr(time)];

%% Taylor

clear optimizedParamTemp2

for i = 1:row
    optimizedParamTemp  = Taylor1996_ParamHealthyCalibrated((i),1:(column-1));
    optimizedParamTemp(16) = optimizedParamTemp(15); % Identical meals, but consumed at different times. Params are set to not allow differences in insulin response etc.
    optimizedParamTemp(34) = optimizedParamTemp(33);
    optimizedParamTemp(56) = optimizedParamTemp(55);
    optimizedParamTemp(60) = optimizedParamTemp(59);
    optimizedParamTemp(64) = optimizedParamTemp(63);
    cost                = Taylor1996_costfunction(Taylor1996_data,time,log(optimizedParamTemp),modelName,body_information, meal_information,9,13,13);
    optimizedParamTemp(column) = cost;
    optimizedParamTemp2(i,:)   = optimizedParamTemp;
end
Taylor1996_ParamHealthyCalibrated = sortrows(optimizedParamTemp2,column);

disp(' ')
if Taylor1996_ParamHealthyCalibrated(1,column) < chi2inv(0.95,35)
    disp('----- Best prediction of Taylor data is below treshold -----')
end

fprintf('Best prediction of Taylor data: %.2f, Statistical Limit: %.2f (dgf = %i)', Taylor1996_ParamHealthyCalibrated(1,column), chi2inv(0.95,35), 35)
disp(' ')
%%
for i = 1:row
    
    try
        optimizedParamTemp = Taylor1996_ParamHealthyCalibrated(i,1:(column-1));
        
        optimizedParamTemp(87:101)    = [0,0,1,0,0,0,0,0,1,0,0,0,1,0,0];   % Steady State Meals are solid and there is no meal in t0
        optimizedParamTemp(102:105)   = [0 , 1, 180, 75.6];                % female, male, height, weight
        optimizedParamTemp(106:109)   = [1,0,0,0];                         % This is meal A; Identical meals
        optimizedParamTemp(110:117)   = [10, 0, 138.64, 29.25,1,0,0,0];    % Meal Information
        
        simA = model(time,[],optimizedParamTemp);
        
        optimizedParamTemp(106:109)   = [0,1,0,0];                         % This is meal B; Identical meals
        simB = model(time,[],optimizedParamTemp);
        
        
        
        if i == 1
            simBestA = simA;
            simBestB = simB;
        end
        
        % Glucose
        if i == 1
            maxG1 = simB.variablevalues(:,ismember(simB.variables,'G'))/18;
            minG1 = simB.variablevalues(:,ismember(simB.variables,'G'))/18;
            maxG2 = simB.variablevalues(:,ismember(simB.variables,'G'))/18;
            minG2 = simB.variablevalues(:,ismember(simB.variables,'G'))/18;
        else
            maxG1 = simB.variablevalues(:,ismember(simB.variables,'G'))/18;
            minG1 = simB.variablevalues(:,ismember(simB.variables,'G'))/18;
            maxG2 = max(maxG2,maxG1);
            minG2 = min(minG2,minG1);
        end
        
        % Insulin
        if i == 1
            maxI1 = simB.variablevalues(:,ismember(simB.variables,'I'));
            minI1 = simB.variablevalues(:,ismember(simB.variables,'I'));
            maxI2 = simB.variablevalues(:,ismember(simB.variables,'I'));
            minI2 = simB.variablevalues(:,ismember(simB.variables,'I'));
        else
            maxI1 = simB.variablevalues(:,ismember(simB.variables,'I'));
            minI1 = simB.variablevalues(:,ismember(simB.variables,'I'));
            maxI2 = max(maxI2,maxI1);
            minI2 = min(minI2,minI1);
        end
        
        % Glycogen
        if i == 1
            maxGly1 = simA.variablevalues(:,ismember(simA.variables,'Glycogen_liver'));
            minGly1 = simA.variablevalues(:,ismember(simA.variables,'Glycogen_liver'));
            maxGly2 = simA.variablevalues(:,ismember(simA.variables,'Glycogen_liver'));
            minGly2 = simA.variablevalues(:,ismember(simA.variables,'Glycogen_liver'));
        else
            maxGly1 = simA.variablevalues(:,ismember(simA.variables,'Glycogen_liver'));
            minGly1 = simA.variablevalues(:,ismember(simA.variables,'Glycogen_liver'));
            maxGly2 = max(maxGly2,maxGly1);
            minGly2 = min(minGly2,minGly1);
        end
    end
end

%%
% Glycogen
figure('Name', "Glycogen", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
set(gca,'xtick',[0,6,12],'ytick',[150,250,350],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxGly2', fliplr(minGly2')];
fill(time2/60-127,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBestA.variablevalues(:,ismember(simBestA.variables,'Glycogen_liver')),'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Hepatic glycogen' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
errorbar(Taylor1996_data.time_glycogen_healthy(1:2)/60-127,Taylor1996_data.glycogen_healthy(1:2), Taylor1996_data.glycogenSEM_healthy(1:2),' k x','MarkerSize',70,'LineWidth',10,'CapSize',CapSize');
errorbar(Taylor1996_data.time_glycogen_healthy/60-127,Taylor1996_data.glycogen_healthy, Taylor1996_data.glycogenSEM_healthy,' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([1 1.16], [150 150],'Color','k','LineWidth',12);
xlim([0 12])
ylim([150 350])
hold off

%%
% Glucose
figure('Name', "Glucose", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,6,12],'ytick',[0,2,4,6,8,10],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [minG2', fliplr(maxG2')];
fill(time2/60-127,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBestB.variablevalues(:,ismember(simBestB.variables,'G'))/18,'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Plasma glucose' ; '(mM)'},'FontSmoothing','on','fontname','Arial');
errorbar(Taylor1996_data.time_glucose_healthy(1:3)/60-127,Taylor1996_data.glucose_healthy(1:3)/18, Taylor1996_data.glucoseSEM_healthy(1:3)/18,' k x','MarkerSize',70,'LineWidth',10,'CapSize',CapSize');
errorbar(Taylor1996_data.time_glucose_healthy(1:13)/60-127,Taylor1996_data.glucose_healthy(1:13)/18, Taylor1996_data.glucoseSEM_healthy(1:13)/18,' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([1 1.16], [4 4],'Color','k','LineWidth',12);
xlim([0 12])
ylim([4 11])
hold off

%%
% Insulin
figure('Name', "Insulin", 'units', 'normalized', 'outerposition', [0 0 1 1])
hold on
a = gca;
set(a,'xtick',[0,6,12],'ytick',[0,400,800],'FontSize', 70,'fontname','Arial')%,'FontSmoothing','on')
y        = [maxI2', fliplr(minI2')];
fill(time2/60-127,y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
plot(time/60-127, simBestB.variablevalues(:,ismember(simBestB.variables,'I')),'b-.','LineWidth',LineWidthValue);
xlabel("Time (h)",'FontSmoothing','on','fontname','Arial');
ylabel({'Plasma insulin' ; '(pM)'},'FontSmoothing','on','fontname','Arial');
errorbar(Taylor1996_data.time_Insulin_healthy(1:3)/60-127,Taylor1996_data.Insulin_healthy(1:3), Taylor1996_data.InsulinSEM_healthy(1:3),' kx','MarkerSize',70,'LineWidth',10,'CapSize',CapSize');
errorbar(Taylor1996_data.time_Insulin_healthy(1:13)/60-127,Taylor1996_data.Insulin_healthy(1:13), Taylor1996_data.InsulinSEM_healthy(1:13),' k .','MarkerSize',1,'LineWidth',10,'CapSize',CapSize');
line([1 1.16], [0 0],'Color','k','LineWidth',12);
xlim([0 12])
hold off