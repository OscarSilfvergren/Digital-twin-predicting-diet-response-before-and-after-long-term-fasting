function [R] = kin_catalytic_activation_rev(Vf,substrate,Kms,Vr,product,Kmp,activator,Ka)
% kin_catalytic_activation_rev: Catalytic activation (reversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_catalytic_activation_rev(Vf,substrate,Kms,Vr,product,Kmp,activator,Ka)
%
% Output Arguments:
% =================
% R = ( Vf*substrate/Kms - Vr*product/Kmp )*activator / ( (1+substrate/Kms+product/Kmp)*(Ka+activator) )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = ( Vf*substrate/Kms - Vr*product/Kmp )*activator / ( (1+substrate/Kms+product/Kmp)*(Ka+activator) );

