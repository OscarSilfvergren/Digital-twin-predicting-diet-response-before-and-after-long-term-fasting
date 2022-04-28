function [username] = usernameIQM()
% usernameIQM: Get the name of the current user

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

if ispc,
    username = getenv('UserName');
else
    username = getenv('USER');
end
return
