function [R] = kin_specific_activation_rev(Vf,substrate,Kms,Vr,product,Kmp,activator,Ka)
% kin_specific_activation_rev: Specific activation (reversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_specific_activation_rev(Vf,substrate,Kms,Vr,product,Kmp,activator,Ka)
%
% Output Arguments:
% =================
% R = ( Vf*substrate/Kms-Vr*product/Kmp )*activator / ( Ka+(1+substrate/Kms+product/Kmp)*activator )  

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = ( Vf*substrate/Kms-Vr*product/Kmp )*activator / ( Ka+(1+substrate/Kms+product/Kmp)*activator );

