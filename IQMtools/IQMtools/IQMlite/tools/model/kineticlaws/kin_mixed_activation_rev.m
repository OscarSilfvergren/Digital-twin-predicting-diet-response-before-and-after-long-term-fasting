function [R] = kin_mixed_activation_rev(Vf,substrate,Kms,Vr,product,Kmp,activator,Kas,Kac)
% kin_mixed_activation_rev: Mixed activation reversible
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_mixed_activation_rev(Vf,substrate,Kms,Vr,product,Kmp,activator,Kas,Kac)
%
% Output Arguments:
% =================
% R = (Vf*substrate/Kms - Vr*product/Kmp)*activator / ( Kas+activator+(substrate/Kms+product/Kmp)*(Kac+activator) )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = (Vf*substrate/Kms - Vr*product/Kmp)*activator / ( Kas+activator+(substrate/Kms+product/Kmp)*(Kac+activator) );

