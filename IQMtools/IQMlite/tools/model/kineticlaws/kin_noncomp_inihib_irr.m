function [R] = kin_noncomp_inihib_irr(V,substrate,Km,inhibitor,Ki)
% kin_noncomp_inihib_irr: Noncompetitive inhibition (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_noncomp_inihib_irr(V,substrate,Km,inhibitor,Ki)
%
% Output Arguments:
% =================
% R = V*substrate / ( (Km+substrate)*(1+inhibitor/Ki) )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = V*substrate / ( (Km+substrate)*(1+inhibitor/Ki) );

