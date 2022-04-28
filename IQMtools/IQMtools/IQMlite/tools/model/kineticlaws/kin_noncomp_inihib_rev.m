function [R] = kin_noncomp_inihib_rev(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Ki)
% kin_noncomp_inihib_rev: Noncompetitive inhibition (reversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_noncomp_inihib_rev(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Ki)
%
% Output Arguments:
% =================
% R = (Vf*substrate/Kms-Vr*product/Kmp) / ( (1+substrate/Kms+product/Kmp)*(1+inhibitor/Ki) )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = (Vf*substrate/Kms-Vr*product/Kmp) / ( (1+substrate/Kms+product/Kmp)*(1+inhibitor/Ki) );

