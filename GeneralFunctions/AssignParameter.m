function [newparam] = AssignParameter(params, body_information, meal_information)

newparam         = params;
newparam(81:98)  = 0;
newparam(83)     = 1;
newparam(89)     = 1;
newparam(95)     = 1;
newparam(99:101)  = [1,0,0];
newparam(102:105)  = body_information;
newparam(106:109)  = [1,0,0,0];
newparam(110:113)  = meal_information;
newparam(114:117)  = [1,0, 0,0];
return
end

