function [R] = kin_specific_activation_irr(V,substrate,activator,Kms,Ka)
% kin_specific_activation_irr: Specific activation (irreversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_specific_activation_irr(V,substrate,activator,Kms,Ka)
%
% Output Arguments:
% =================
% R = V*substrate*activator/( Kms*Ka+(Kms+substrate)*activator )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = V*substrate*activator/( Kms*Ka+(Kms+substrate)*activator );

