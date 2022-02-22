function [lb,ub] = DifferentPopulationCalibration(lb,ub,ParameterBounds)

% Basal Glucose
    lb(15)   = log(ParameterBounds.LowerBoundHealthy(15));
    ub(15)   = log(ParameterBounds.UpperBoundHealthy(15));
    
%Basal Insulin    
    lb(33)   = log(ParameterBounds.LowerBoundHealthy(33));
    ub(33)   = log(ParameterBounds.UpperBoundHealthy(33));
    
% Different insulin resistance    
    lb(55)   = log(ParameterBounds.LowerBoundHealthy(55));
    ub(55)   = log(ParameterBounds.UpperBoundHealthy(55));
    
%Different insulin clearance    
    lb(59)   = log(ParameterBounds.LowerBoundHealthy(59));
    ub(59)   = log(ParameterBounds.UpperBoundHealthy(59));
    
%Different insulin production    
    lb(63)   = log(ParameterBounds.LowerBoundHealthy(63));
    ub(63)   = log(ParameterBounds.UpperBoundHealthy(63));
    
end

