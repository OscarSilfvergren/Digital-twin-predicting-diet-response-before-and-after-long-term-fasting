function [R] = kin_uni_uni_rev(Vf,substrate,product,Keq,Kms,Kmp)
% kin_uni_uni_rev: uni uni reversible kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_uni_uni_rev(Vf,substrate,product,Keq,Kms,Kmp)
%
% Output Arguments:
% =================
% R = Vf*( substrate-product/Keq ) / ( substrate+Kms*(1+product/Kmp) )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = Vf*( substrate-product/Keq ) / ( substrate+Kms*(1+product/Kmp) );

