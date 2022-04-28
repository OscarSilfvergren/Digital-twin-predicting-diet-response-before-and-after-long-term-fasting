function [ ispresent ] = isIQMproPresent()
% isIQMproPresent: Checks if the IQMpro tools are installed. 

ispresent = 1;
if exist('IQMPsimulate') ~= 2,
    ispresent = 0;
end  

