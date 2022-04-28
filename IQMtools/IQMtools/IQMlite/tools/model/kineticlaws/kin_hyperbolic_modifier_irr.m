function [R] = kin_hyperbolic_modifier_irr(V,substrate,b,modifier,a,Kd,Km)
% kin_hyperbolic_modifier_irr: Hyperbolic modifier (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_hyperbolic_modifier_irr(V,substrate,b,modifier,a,Kd,Km)
%
% Output Arguments:
% =================
% R = V*substrate*(1+b*modifier/(a*Kd)) / ( Km*(1+modifier/Kd) + substrate*(1+modifier/(a*Kd)) )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = V*substrate*(1+b*modifier/(a*Kd)) / ( Km*(1+modifier/Kd) + substrate*(1+modifier/(a*Kd)) );

