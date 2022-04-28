function [R] = kin_mixed_activation_irr(V,substrate,activator,Kms,Kas,Kac)
% kin_mixed_activation_irr: Mixed activation irreversible
%
% Application restrictions: 
% =========================
%   - Only irreversible
%   - exactly 1 substrate
%
% USAGE:
% ======
% R = kin_mixed_activation_irr(V,substrate,activator,Kms,Kas,Kac)
%
% Output Arguments:
% =================
% R = V*substrate*activator / ( Kms*(Kas+activator) + substrate*(Kac+activator) )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = V*substrate*activator / ( Kms*(Kas+activator) + substrate*(Kac+activator) );

