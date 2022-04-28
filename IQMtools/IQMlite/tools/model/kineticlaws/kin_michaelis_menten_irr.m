function [R] = kin_michaelis_menten_irr(V,substrate,Km)
% kin_michaelis_menten_irr: Michaelis Menten (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_michaelis_menten_irr(V,substrate,Km)
%
% Output Arguments:
% =================
% R = V*substrate / ( Km + substrate )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = V*substrate / ( Km + substrate );

