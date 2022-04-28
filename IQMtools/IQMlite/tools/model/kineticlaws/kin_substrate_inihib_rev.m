function [R] = kin_substrate_inihib_rev(Vf,substrate,Kms,Vr,product,Kmp,Ki)
% kin_substrate_inihib_rev: Substrate inhibition (reversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_substrate_inihib_rev(Vf,substrate,Kms,Vr,product,Kmp,Ki)
%
% Output Arguments:
% =================
% R = (Vf*substrate/Kms-Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp+(substrate/Ki)^2 ) 

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = (Vf*substrate/Kms-Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp+(substrate/Ki)^2 );

