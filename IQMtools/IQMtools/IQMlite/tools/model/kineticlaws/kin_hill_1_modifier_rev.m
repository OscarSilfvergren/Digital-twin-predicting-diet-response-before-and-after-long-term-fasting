function [R] = kin_hill_1_modifier_rev(Vf,substrate,Shalve,product,Keq,Phalve,h,modifier,Mhalve,alpha)
% kin_hill_1_modifier_rev: Reversible Hill type kinetics with one modifier
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_hill_1_modifier_rev(Vf,substrate,Shalve,product,Keq,Phalve,h,modifier,Mhalve,alpha)
%
% Output Arguments:
% =================
% R = Vf*substrate/Shalve*(1-product/(substrate*Keq))*(substrate/Shalve+product/Phalve)^(h-1) / ( (1+(modifier/Mhalve)^h)/(1+alpha*(modifier/Mhalve)^h)+(substrate/Shalve + product/Phalve)^h ) 

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

R = Vf*substrate/Shalve*(1-product/(substrate*Keq))*(substrate/Shalve+product/Phalve)^(h-1) / ( (1+(modifier/Mhalve)^h)/(1+alpha*(modifier/Mhalve)^h)+(substrate/Shalve + product/Phalve)^h );
