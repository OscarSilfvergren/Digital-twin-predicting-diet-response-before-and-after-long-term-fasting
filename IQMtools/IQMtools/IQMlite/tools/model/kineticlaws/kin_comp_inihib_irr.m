function [R] = kin_comp_inihib_irr(V,substrate,Km,inhibitor,Ki)
% kin_comp_inihib_irr: Competitive inhibition (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_comp_inihib_irr(V,substrate,Km,inhibitor,Ki)
%
% Output Arguments:
% =================
% R = V*substrate / ( Km*(1+inhibitor/Ki) + substrate )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = V*substrate / ( Km*(1+inhibitor/Ki) + substrate );

