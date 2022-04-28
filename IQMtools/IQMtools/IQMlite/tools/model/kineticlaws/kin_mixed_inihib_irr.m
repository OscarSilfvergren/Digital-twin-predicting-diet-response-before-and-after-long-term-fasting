function [R] = kin_mixed_inihib_irr(V,substrate,Km,inhibitor,Kis,Kic)
% kin_mixed_inihib_irr: Mixed inhibition (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_mixed_inihib_irr(V,substrate,Km,inhibitor,Kis,Kic)
%
% Output Arguments:
% =================
% R = V*substrate / ( (Km*(1+inhibitor/Kis) + substrate*(1+inhibitor/Kic) )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = V*substrate / ( Km*(1+inhibitor/Kis) + substrate*(1+inhibitor/Kic) );

