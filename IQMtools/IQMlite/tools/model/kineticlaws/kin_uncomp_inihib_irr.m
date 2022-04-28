function [R] = kin_uncomp_inihib_irr(V,substrate,Km,inhibitor,Ki)
% kin_uncomp_inihib_irr: Uncompetitive inhibition (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_uncomp_inihib_irr(V,substrate,Km,inhibitor,Ki)
%
% Output Arguments:
% =================
% R = V*substrate / ( Km + substrate*(1+inhibitor/Ki) )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = V*substrate / ( Km + substrate*(1+inhibitor/Ki) );

