function [R] = kin_hill_cooperativity_irr(V,substrate,h,Shalve)
% kin_hill_cooperativity_irr: Hill type (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_hill_cooperativity_irr(V,substrate,h,Shalve)
%
% Output Arguments:
% =================
% R = V*substrate^h / ( Shalve^h + substrate^h )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = V*substrate^h / ( Shalve^h + substrate^h );

