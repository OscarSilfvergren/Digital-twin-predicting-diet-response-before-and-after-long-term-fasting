function [output] = delayIQM(input,tau,time,queuename)
% delayIQM: realizes a time delay of "tau" time units
%
% input:     the input that is to be delayed
% tau:       the delay
% time:      the time of the input
% queuename: unique name for the variable storing the queue data

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

% get the queuedata for the current delay
eval(sprintf('global %s',queuename));
eval(sprintf('queue = %s;',queuename));

if isempty(queue),
% initialize the queue if first call
    queue(1,1) = -1e10;
    queue(1,2) = input;
    queue(2,1) = time+tau;  % add delay to time in order to store the output time
    queue(2,2) = input;
else
% add new time point to the queue
    % if last time+tau < max stored time ... delete all time points larger
    % than time + tau (can happen during event handling)
    queue(find(queue(:,1)>=time+tau),:) = [];
    % add the new point
    queue(end+1,1) = time+tau;   % add delay to time in order to store the output time
    queue(end,2) = input;
end

% % make unique time values (not necessary due to line 23 above???)
% [dummy,indexunique] = unique(queue(:,1),'first');
% queue = queue(indexunique,:);

% interpolate to find the correct value for the current time. 
% (linear interpolation)
output = interp1IQM(queue(:,1),queue(:,2),time);

% store the queue under the right name again
eval(sprintf('%s = queue;',queuename));


