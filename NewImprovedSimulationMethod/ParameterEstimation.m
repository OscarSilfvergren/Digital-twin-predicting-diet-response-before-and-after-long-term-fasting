function [minCostPS] = ParameterEstimation(seed, Krssak2004Healthy_data, Rothman1991_data, Silfvergren2021_data,paramsModelValues, pNames, initialvalues, sNames, model, options)

rng('shuffle') % set the stochastic seed
s = rng;
rng(s.Seed+seed);
AmountParametersOptimized = 58;

%% Pre-Optimization - create bounds

params = 1:1:AmountParametersOptimized;
lb     = log(paramsModelValues(1:AmountParametersOptimized)/3);
ub     = log(paramsModelValues(1:AmountParametersOptimized)*3);

% Meals
lb(ismember(pNames,'Meal_CarbohydratesFlow'))    = log(100*1000/1440);  % 50g/day
ub(ismember(pNames,'Meal_CarbohydratesFlow'))    = log(300*1000/1440); % 300g/day etc.
lb(ismember(pNames,'Meal_ProteinFlow'))          = log(40*1000/1440);
ub(ismember(pNames,'Meal_ProteinFlow'))          = log(200*1000/1440);

% Exp parameters fixed to not be too big
ub(ismember(pNames,'GlyDepEXP_Meal'))             = log(5);
lb(ismember(pNames,'GlyDepEXP_Meal'))             = log(2);
ub(ismember(pNames,'GlyDepEXP_TCA'))              = log(1.1);
lb(ismember(pNames,'GlyDepEXP_TCA'))              = log(0.1);
ub(ismember(pNames,'GlyDepEXP_Glugoneogenesis'))  = log(1.1);
lb(ismember(pNames,'GlyDepEXP_Glugoneogenesis'))  = log(0.1);
ub(ismember(pNames,'GlyDepEXP_KidneysEGP'))       = log(1.1);
lb(ismember(pNames,'GlyDepEXP_KidneysEGP'))       = log(0.1);

% Fix f - loss of glucose in convertion from carbohydrates
ub(ismember(pNames,'f'))             = log(0.99);
lb(ismember(pNames,'f'))             = log(0.8);

% Fix Aminoprofile_k - Procentage of destination of aminoacids to liver
ub(ismember(pNames,'Aminoprofile_k'))             = log(0.95);
lb(ismember(pNames,'Aminoprofile_k'))             = log(0.05);

% Fix blood volume uncertainties from Nadler Eq.
ub(ismember(pNames,'BloodVolumeUncertainty'))            = log(1/0.7);
lb(ismember(pNames,'BloodVolumeUncertainty'))            = log(0.7);
ub(ismember(pNames,'BloodLiverUncerteinty'))             = log(1/0.7);
lb(ismember(pNames,'BloodLiverUncerteinty'))             = log(0.7);

%% Pre-Optimization - Add parameters for multiple studies

AmountExtraParametersOptimized = 9;
AmountExtraStudies             = 2;


for i = 2:AmountExtraStudies+1
    params(end+1:end+AmountExtraParametersOptimized)  = [i];
    
    % Meals
    lb(end+1)     = lb(ismember(pNames,'Meal_CarbohydratesFlow'));    % Steady state carbs
    ub(end+1)     = ub(ismember(pNames,'Meal_CarbohydratesFlow'));
    lb(end+1)     = lb(ismember(pNames,'Meal_ProteinFlow'));          % Steady state protein
    ub(end+1)     = ub(ismember(pNames,'Meal_ProteinFlow'));
    
    % Basal Values: Personal Scaling -  Param 3 & 4
    lb(end+1)     = log(0.7);               % Basal Insulin
    ub(end+1)     = log(1/0.7);
    lb(end+1)     = log(0.7);               % Basal Glucose
    ub(end+1)     = log(1/0.7);
    
    % InsulinProduction: Personal Scaling - Param 5
    lb(end+1)     = log(0.7); % InsulinProduction - all parameters scaled together
    ub(end+1)     = log(1/0.7);
    
    % InsulinClearance: Personal Scaling -  Param 6
    lb(end+1)     = log(0.7); % InsulinClearance - all parameters scaled together
    ub(end+1)     = log(1/0.7);
    
    % InsulinResponse: Personal Scaling -  Param 7
    lb(end+1)     = log(0.7); % InsulinResponse - all parameters scaled together
    ub(end+1)     = log(1/0.7);
    
    % BloodVolumeUncertainty: Personal Scaling -  Param 8 & 9
    lb(end+1)     = log(0.7);         % BloodVolumeUncertainty
    ub(end+1)     = log(1/0.7);
    lb(end+1)     = log(0.7);         % BloodLiverUncerteinty
    ub(end+1)     = log(1/0.7);
    
end

func =@(params)EstimationData_costfunction(Krssak2004Healthy_data,Rothman1991_data,Silfvergren2021_data,params,model,initialvalues, pNames,sNames,AmountParametersOptimized);

%% Optimization

[EstimationData_params, minCostPS] = particleswarm(func, length(lb), lb, ub, options);

% Save
EstimationData_params = exp(EstimationData_params);
save(['EstimationData_params' datestr(now, 'yymmdd-HHMMSS')],'EstimationData_params');

end

