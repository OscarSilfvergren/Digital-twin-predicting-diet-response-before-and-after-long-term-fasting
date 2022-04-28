function [boolean] = isSOSign(character)
% isSOSign
% checks, wether character is sign of "Second Order" (* or /)
%
% USAGE:
% ======
%
% [boolean] = isSOSign(character)
%
% character: test character
%
% boolean: true if character is mathematical operator of second order
%          (* or /) otherwise false

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>


boolean = false;

if (strcmp(character, '*') || strcmp(character, '/')),
    boolean = true;
end

return