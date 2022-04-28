function [R] = kin_mixed_inihib_rev(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Kis,Kic)
% kin_mixed_inihib_rev: Mixed inhibition (reversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_mixed_inihib_rev(Vf,substrate,Kms,Vr,product,Kmp,inhibitor,Kis,Kic)
%
% Output Arguments:
% =================
% R = (Vf*substrate/Kms - Vr*product/Kmp) / ( 1+inhibitor/Kis+(substrate/Kms+product/Kmp)*(1+inhibitor/Kic) )

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = (Vf*substrate/Kms - Vr*product/Kmp) / ( 1+inhibitor/Kis+(substrate/Kms+product/Kmp)*(1+inhibitor/Kic) );

