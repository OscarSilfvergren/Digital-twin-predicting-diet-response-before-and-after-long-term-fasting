function [R] = kin_comp_inihib_rev(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Ki)
% kin_comp_inihib_rev: Competitive inhibition (reversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_comp_inihib_rev(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Ki)
%
% Output Arguments:
% =================
% R = (Vf*substrate/Kms - Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp+inhibitor/Ki )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = (Vf*substrate/Kms - Vr*product/Kmp) / ( 1+substrate/Kms+product/Kmp+inhibitor/Ki );

