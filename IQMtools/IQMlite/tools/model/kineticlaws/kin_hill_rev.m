function [R] = kin_hill_rev(Vf,substrate,Shalve,product,Keq,Phalve,h)
% kin_hill_rev: Hill type (reversible) kinetics
%
% Application restrictions: 
% =========================
%   - Only reversible
%   - exactly 1 substrate
%   - exactly 1 product
%
% USAGE:
% ======
% R = kin_hill_rev(Vf,substrate,Shalve,product,Keq,Phalve,h)
%
% Output Arguments:
% =================
% R = Vf*substrate/Shalve*(1-product/(substrate*Keq))*(substrate/Shalve+product/Phalve)^(h-1) / ( 1+(substrate/Shalve + product/Phalve)^h ) 

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>> 

R = Vf*substrate/Shalve*(1-product/(substrate*Keq))*(substrate/Shalve+product/Phalve)^(h-1) / ( 1+(substrate/Shalve + product/Phalve)^h );

